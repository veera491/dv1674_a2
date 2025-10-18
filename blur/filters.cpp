/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <vector>

#ifndef FILTER_SMALL_R_UNROLL
#define FILTER_SMALL_R_UNROLL 8
#endif

// Tile size for blocked transpose (tune 16..64 depending on cache)
#ifndef FILTER_TRANSPOSE_TILE
#define FILTER_TRANSPOSE_TILE 32
#endif

namespace Filter {

namespace Gauss {
// Implemented here to satisfy the linker and keep everything self-contained.
// weights[i] is the coefficient for offset ±i (0..n) before normalization.
void get_weights(int n, double* weights_out)
{
    // same constants as your original
    constexpr float max_x{1.33f};
    constexpr float pi{3.14159f};

    for (int i = 0; i <= n; ++i) {
        double x = static_cast<double>(i) * max_x / static_cast<double>(n);
        weights_out[i] = std::exp(-x * x * pi);
    }
}
} // namespace Gauss

// ------------------------ utilities ------------------------

static inline unsigned char f2u8(float v) {
    if (v < 0.0f)   return 0;
    if (v > 255.0f) return 255;
    return static_cast<unsigned char>(v + 0.5f);
}

static inline void matrix_to_planar(const Matrix& m,
                                    float* __restrict r,
                                    float* __restrict g,
                                    float* __restrict b)
{
    const unsigned W = m.get_x_size();
    const unsigned H = m.get_y_size();
    for (unsigned y = 0; y < H; ++y) {
        const unsigned row = y * W;
        for (unsigned x = 0; x < W; ++x) {
            const unsigned idx = row + x;
            r[idx] = static_cast<float>(m.r(x, y));
            g[idx] = static_cast<float>(m.g(x, y));
            b[idx] = static_cast<float>(m.b(x, y));
        }
    }
}

static inline void planar_to_matrix(const float* __restrict r,
                                    const float* __restrict g,
                                    const float* __restrict b,
                                    Matrix& out)
{
    const unsigned W = out.get_x_size();
    const unsigned H = out.get_y_size();
    for (unsigned y = 0; y < H; ++y) {
        const unsigned row = y * W;
        for (unsigned x = 0; x < W; ++x) {
            const unsigned idx = row + x;
            out.r(x, y) = f2u8(r[idx]);
            out.g(x, y) = f2u8(g[idx]);
            out.b(x, y) = f2u8(b[idx]);
        }
    }
}

// Blocked transpose: src is (H x W) in row-major, dst becomes (W x H) row-major.
static inline void transpose_plane_blocked(const float* __restrict src,
                                           float* __restrict dst,
                                           unsigned W, unsigned H,
                                           unsigned B = FILTER_TRANSPOSE_TILE)
{
    for (unsigned y0 = 0; y0 < H; y0 += B) {
        const unsigned y1 = std::min(H, y0 + B);
        for (unsigned x0 = 0; x0 < W; x0 += B) {
            const unsigned x1 = std::min(W, x0 + B);
            for (unsigned y = y0; y < y1; ++y) {
                const unsigned srcRow = y * W;
                for (unsigned x = x0; x < x1; ++x) {
                    // (x, y) in src -> (y, x) in dst of size WxH -> index y*H + x
                    dst[y * H + x] = src[srcRow + x];
                }
            }
        }
    }
}

// Pre-pad with reflected borders so inner loops have no bounds checks.
static inline void reflect_pad_plane(const float* __restrict in,
                                     float* __restrict out,
                                     unsigned W, unsigned H, int R)
{
    const unsigned PW = W + 2u * static_cast<unsigned>(R);

    // 1) center copy
    for (unsigned y = 0; y < H; ++y) {
        const unsigned srcRow = y * W;
        const unsigned dstRow = (y + static_cast<unsigned>(R)) * PW + static_cast<unsigned>(R);
        std::copy(in + srcRow, in + srcRow + W, out + dstRow);
    }
    // 2) left/right reflect per row
    for (unsigned y = 0; y < H; ++y) {
        const unsigned dstRow = (y + static_cast<unsigned>(R)) * PW;
        for (int k = 0; k < R; ++k) {
            out[dstRow + (static_cast<unsigned>(R) - 1u - static_cast<unsigned>(k))] =
                out[dstRow + (static_cast<unsigned>(R) + static_cast<unsigned>(k))];
            out[dstRow + (static_cast<unsigned>(R) + W + static_cast<unsigned>(k))] =
                out[dstRow + (static_cast<unsigned>(R) + W - 1u - static_cast<unsigned>(k))];
        }
    }
    // 3) top/bottom reflect full padded rows
    for (int k = 0; k < R; ++k) {
        const unsigned PWstride = PW;
        const unsigned srcTop    = (static_cast<unsigned>(R) + static_cast<unsigned>(k)) * PWstride;
        const unsigned srcBottom = (static_cast<unsigned>(R) + H - 1u - static_cast<unsigned>(k)) * PWstride;
        const unsigned dstTop    = (static_cast<unsigned>(R) - 1u - static_cast<unsigned>(k)) * PWstride;
        const unsigned dstBottom = (static_cast<unsigned>(R) + H + static_cast<unsigned>(k)) * PWstride;
        std::copy(out + srcTop,    out + srcTop    + PWstride, out + dstTop);
        std::copy(out + srcBottom, out + srcBottom + PWstride, out + dstBottom);
    }
}

// Convolve one padded row (branchless) and write **transposed** into dstT.
//  - row_pad points to the *start of the padded row* (left padding at +0).
//  - We write result for output row y_out into transposed buffer at (x,y)->dstT[x*H + y_out].
static inline void gauss_row_padded_to_transposed(const float* __restrict row_pad,
                                                  float* __restrict dstT,
                                                  unsigned W, unsigned H, unsigned y_out,
                                                  const std::vector<float>& gw,
                                                  int radius)
{
    const unsigned base = static_cast<unsigned>(radius);   // first valid sample
    for (unsigned x = 0; x < W; ++x) {
        const unsigned cx = base + x;
        float acc = gw[0] * row_pad[cx];
        #pragma omp simd
        for (int k = 1; k <= radius; ++k) {
            const float w = gw[k];
            acc += w * ( row_pad[cx - static_cast<unsigned>(k)] +
                         row_pad[cx + static_cast<unsigned>(k)] );
        }
        // transpose write (x,y_out) -> [x*H + y_out]
        dstT[x * H + y_out] = acc;
    }
}

// Convolve one padded row of the **transposed** image and write to another
// transposed buffer (i.e., both src and dst are HxW row-major).
static inline void gauss_row_padded_transposed_to_transposed(const float* __restrict row_pad,
                                                             float* __restrict dstT,
                                                             unsigned TW, unsigned TH, unsigned y_tr,
                                                             const std::vector<float>& gw,
                                                             int radius)
{
    // y_tr is the row index in the transposed grid (range [0..TH-1])
    // TW is the width of the transposed grid (= original H)
    const unsigned base = static_cast<unsigned>(radius);
    const unsigned row_start = y_tr * TW;
    for (unsigned x = 0; x < TW; ++x) {
        const unsigned cx = base + x;
        float acc = gw[0] * row_pad[cx];
        #pragma omp simd
        for (int k = 1; k <= radius; ++k) {
            const float w = gw[k];
            acc += w * ( row_pad[cx - static_cast<unsigned>(k)] +
                         row_pad[cx + static_cast<unsigned>(k)] );
        }
        dstT[row_start + x] = acc; // row-wise contiguous write
    }
}

// ------------------------ main blur -----------------------------------------

Matrix blur(Matrix m, const int radius)
{
    const unsigned W = m.get_x_size();
    const unsigned H = m.get_y_size();
    if (W == 0 || H == 0 || radius <= 0) {
        return m;
    }

    // 1) Precompute Gaussian weights (0..R) and **normalize to sum==1**
    //    With reflected padding, every pixel has a full neighborhood,
    //    so we can normalize once and avoid per-pixel division.
    std::vector<double> gw_d(static_cast<std::size_t>(radius) + 1);
    Gauss::get_weights(radius, gw_d.data());

    // sum_full = w0 + 2*sum_{k=1..R} wk
    double sum_full = gw_d[0];
    for (int k = 1; k <= radius; ++k) sum_full += 2.0 * gw_d[k];
    const float inv_sum = static_cast<float>(1.0 / sum_full);

    // convert to normalized float weights
    std::vector<float> gw(gw_d.size());
    gw[0] = static_cast<float>(gw_d[0]) * inv_sum;
    for (int k = 1; k <= radius; ++k) {
        gw[k] = static_cast<float>(gw_d[k]) * inv_sum;
    }

    // 2) Convert to planar
    const std::size_t N = static_cast<std::size_t>(W) * static_cast<std::size_t>(H);
    std::vector<float> r_in(N), g_in(N), b_in(N);
    matrix_to_planar(m, r_in.data(), g_in.data(), b_in.data());

    // 3) Reflected padding for the original (W x H)
    const unsigned PW = W + 2u * static_cast<unsigned>(radius);
    const unsigned PH = H + 2u * static_cast<unsigned>(radius);
    const std::size_t NP = static_cast<std::size_t>(PW) * static_cast<std::size_t>(PH);
    std::vector<float> rp(NP), gp(NP), bp(NP);
    reflect_pad_plane(r_in.data(), rp.data(), W, H, radius);
    reflect_pad_plane(g_in.data(), gp.data(), W, H, radius);
    reflect_pad_plane(b_in.data(), bp.data(), W, H, radius);

    // 4) Pass 1: horizontal on padded rows, **write transposed** directly (HxW buffers)
    std::vector<float> r_tr(N), g_tr(N), b_tr(N); // dims H x W
    for (unsigned y = 0; y < H; ++y) {
        const unsigned prow = (y + static_cast<unsigned>(radius)) * PW; // start of that padded row
        gauss_row_padded_to_transposed(rp.data() + prow, r_tr.data(), W, H, y, gw, radius);
        gauss_row_padded_to_transposed(gp.data() + prow, g_tr.data(), W, H, y, gw, radius);
        gauss_row_padded_to_transposed(bp.data() + prow, b_tr.data(), W, H, y, gw, radius);
    }

    // 5) Reflected padding for the transposed image (TW x TH) = (H x W)
    const unsigned TW = H, TH = W; // transposed dims
    const unsigned TPW = TW + 2u * static_cast<unsigned>(radius);
    const unsigned TPH = TH + 2u * static_cast<unsigned>(radius);
    const std::size_t TNP = static_cast<std::size_t>(TPW) * static_cast<std::size_t>(TPH);
    std::vector<float> rtp(TNP), gtp(TNP), btp(TNP);
    reflect_pad_plane(r_tr.data(), rtp.data(), TW, TH, radius);
    reflect_pad_plane(g_tr.data(), gtp.data(), TW, TH, radius);
    reflect_pad_plane(b_tr.data(), btp.data(), TW, TH, radius);

    // 6) Pass 2: horizontal on padded transposed rows, **write to transposed output** (HxW)
    std::vector<float> r_tr2(N), g_tr2(N), b_tr2(N);
    for (unsigned y_tr = 0; y_tr < TH; ++y_tr) {
        const unsigned prow = (y_tr + static_cast<unsigned>(radius)) * TPW;
        gauss_row_padded_transposed_to_transposed(rtp.data() + prow, r_tr2.data(), TW, TH, y_tr, gw, radius);
        gauss_row_padded_transposed_to_transposed(gtp.data() + prow, g_tr2.data(), TW, TH, y_tr, gw, radius);
        gauss_row_padded_transposed_to_transposed(btp.data() + prow, b_tr2.data(), TW, TH, y_tr, gw, radius);
    }

    // 7) Blocked transpose back to WxH (fast) & pack to Matrix
    std::vector<float> r_out(N), g_out(N), b_out(N);
    transpose_plane_blocked(r_tr2.data(), r_out.data(), TW, TH); // (H x W) -> (W x H)
    transpose_plane_blocked(g_tr2.data(), g_out.data(), TW, TH);
    transpose_plane_blocked(b_tr2.data(), b_out.data(), TW, TH);

    Matrix dst{m};
    planar_to_matrix(r_out.data(), g_out.data(), b_out.data(), dst);
    return dst;
}

} // namespace Filter
