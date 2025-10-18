/*
Author: David Holmqvist <daae19@student.bth.se>
Optimized by: (your names)
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <vector>

// Optional backends (compile-time). Default keeps exact Gaussian FIR.
//// #define FILTER_USE_BOX3
//// #define FILTER_USE_IIR
//// #define FILTER_USE_FFT

#ifndef FILTER_SMALL_R_UNROLL
#define FILTER_SMALL_R_UNROLL 8
#endif

namespace Filter {

namespace Gauss {

// IMPLEMENTATION INCLUDED to satisfy the linker (same math as your original).
// Builds weights[0..n] where weights[i] corresponds to offset ±i.
void get_weights(int n, double* weights_out)
{
    // constants from your original code
    constexpr float max_x{1.33f};
    constexpr float pi{3.14159f};

    for (int i = 0; i <= n; ++i)
    {
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

static inline void transpose_plane(const float* __restrict src,
                                   float* __restrict dst,
                                   unsigned W, unsigned H)
{
    for (unsigned y = 0; y < H; ++y) {
        const unsigned row = y * W;
        for (unsigned x = 0; x < W; ++x) {
            dst[x * H + y] = src[row + x];
        }
    }
}

// Pre-pad a plane with reflected borders so the blur inner loop has no bounds checks.
static inline void reflect_pad_plane(const float* __restrict in,
                                     float* __restrict out,
                                     unsigned W, unsigned H, int R)
{
    const unsigned PW = W + 2u * static_cast<unsigned>(R);
    // center copy
    for (unsigned y = 0; y < H; ++y) {
        const unsigned srcRow = y * W;
        const unsigned dstRow = (y + static_cast<unsigned>(R)) * PW + static_cast<unsigned>(R);
        std::copy(in + srcRow, in + srcRow + W, out + dstRow);
    }
    // left/right reflect inside each padded row
    for (unsigned y = 0; y < H; ++y) {
        const unsigned dstRow = (y + static_cast<unsigned>(R)) * PW;
        for (int k = 0; k < R; ++k) {
            out[dstRow + (static_cast<unsigned>(R) - 1u - static_cast<unsigned>(k))] =
                out[dstRow + (static_cast<unsigned>(R) + static_cast<unsigned>(k))];
            out[dstRow + (static_cast<unsigned>(R) + W + static_cast<unsigned>(k))] =
                out[dstRow + (static_cast<unsigned>(R) + W - 1u - static_cast<unsigned>(k))];
        }
    }
    // top/bottom reflect full padded rows
    const unsigned PWbytes = PW;
    for (int k = 0; k < R; ++k) {
        const unsigned srcTop    = (static_cast<unsigned>(R) + static_cast<unsigned>(k)) * PWbytes;
        const unsigned srcBottom = (static_cast<unsigned>(R) + H - 1u - static_cast<unsigned>(k)) * PWbytes;
        const unsigned dstTop    = (static_cast<unsigned>(R) - 1u - static_cast<unsigned>(k)) * PWbytes;
        const unsigned dstBottom = (static_cast<unsigned>(R) + H + static_cast<unsigned>(k)) * PWbytes;
        std::copy(out + srcTop,    out + srcTop    + PWbytes, out + dstTop);
        std::copy(out + srcBottom, out + srcBottom + PWbytes, out + dstBottom);
    }
}

// Convolve a single *padded* row branchlessly with symmetric Gaussian weights.
static inline void horizontal_gauss_padded_row(const float* __restrict row_padded,
                                               float* __restrict row_out,
                                               unsigned W,
                                               const std::vector<float>& gw,
                                               int radius,
                                               float norm_full)
{
    const unsigned base = static_cast<unsigned>(radius); // first valid index into padded row
    for (unsigned x = 0; x < W; ++x) {
        const unsigned cx = base + x;
        float acc = gw[0] * row_padded[cx];

        // Pair symmetric taps; interior is branch-free due to padding.
        #pragma omp simd
        for (int k = 1; k <= radius; ++k) {
            const float w = gw[k];
            acc += w * ( row_padded[cx - static_cast<unsigned>(k)]
                       + row_padded[cx + static_cast<unsigned>(k)] );
        }
        row_out[x] = acc / norm_full;
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

    // Precompute Gaussian weights (exact FIR, your formula), convert to float.
    std::vector<double> gw_d(static_cast<std::size_t>(radius) + 1);
    Gauss::get_weights(radius, gw_d.data());
    std::vector<float> gw(gw_d.size());
    for (std::size_t i = 0; i < gw.size(); ++i) gw[i] = static_cast<float>(gw_d[i]);

    // Full normalization for interior (symmetric kernel).
    float norm_full = gw[0];
    for (int i = 1; i <= radius; ++i) norm_full += 2.0f * gw[i];

    // Convert to planar (SoA)
    const std::size_t N = static_cast<std::size_t>(W) * static_cast<std::size_t>(H);
    std::vector<float> r_in(N), g_in(N), b_in(N);
    matrix_to_planar(m, r_in.data(), g_in.data(), b_in.data());

    // ---- Pass 1: horizontal on rows with reflected padding, write transposed ----
    const unsigned PW = W + 2u * static_cast<unsigned>(radius);
    const std::size_t NP = static_cast<std::size_t>(PW) * static_cast<std::size_t>(H + 2u * static_cast<unsigned>(radius));
    std::vector<float> rp(NP), gp(NP), bp(NP);
    reflect_pad_plane(r_in.data(), rp.data(), W, H, radius);
    reflect_pad_plane(g_in.data(), gp.data(), W, H, radius);
    reflect_pad_plane(b_in.data(), bp.data(), W, H, radius);

    // temp row buffers (untransposed outputs of pass 1)
    std::vector<float> rowR(W), rowG(W), rowB(W);
    // transposed after pass 1: dims H x W
    std::vector<float> r_tr(N), g_tr(N), b_tr(N);

    for (unsigned y = 0; y < H; ++y) {
        // padded row start (skip the left padding by +radius)
        const unsigned prow = (y + static_cast<unsigned>(radius)) * PW;

        horizontal_gauss_padded_row(rp.data() + prow, rowR.data(), W, gw, radius, norm_full);
        horizontal_gauss_padded_row(gp.data() + prow, rowG.data(), W, gw, radius, norm_full);
        horizontal_gauss_padded_row(bp.data() + prow, rowB.data(), W, gw, radius, norm_full);

        // transpose write: (x,y) -> (y,x) : dst[x*H + y]
        for (unsigned x = 0; x < W; ++x) {
            r_tr[x * H + y] = rowR[x];
            g_tr[x * H + y] = rowG[x];
            b_tr[x * H + y] = rowB[x];
        }
    }

    // ---- Pass 2: run the same horizontal kernel over the transposed image ----
    const unsigned TW = H, TH = W;

    const unsigned TPW = TW + 2u * static_cast<unsigned>(radius);
    const std::size_t TNP = static_cast<std::size_t>(TPW) * static_cast<std::size_t>(TH + 2u * static_cast<unsigned>(radius));
    std::vector<float> rtp(TNP), gtp(TNP), btp(TNP);
    reflect_pad_plane(r_tr.data(), rtp.data(), TW, TH, radius);
    reflect_pad_plane(g_tr.data(), gtp.data(), TW, TH, radius);
    reflect_pad_plane(b_tr.data(), btp.data(), TW, TH, radius);

    std::vector<float> rowRT(TW), rowGT(TW), rowBT(TW);
    std::vector<float> r_tr2(N), g_tr2(N), b_tr2(N);

    for (unsigned y = 0; y < TH; ++y) {
        const unsigned prow = (y + static_cast<unsigned>(radius)) * TPW;

        horizontal_gauss_padded_row(rtp.data() + prow, rowRT.data(), TW, gw, radius, norm_full);
        horizontal_gauss_padded_row(gtp.data() + prow, rowGT.data(), TW, gw, radius, norm_full);
        horizontal_gauss_padded_row(btp.data() + prow, rowBT.data(), TW, gw, radius, norm_full);

        // write back “untransposed on transposed grid”
        const unsigned row = y * TW;
        std::copy(rowRT.begin(), rowRT.end(), r_tr2.begin() + row);
        std::copy(rowGT.begin(), rowGT.end(), g_tr2.begin() + row);
        std::copy(rowBT.begin(), rowBT.end(), b_tr2.begin() + row);
    }

    // ---- Transpose back to WxH and pack to Matrix ----
    std::vector<float> r_out(N), g_out(N), b_out(N);
    transpose_plane(r_tr2.data(), r_out.data(), TW, TH);
    transpose_plane(g_tr2.data(), g_out.data(), TW, TH);
    transpose_plane(b_tr2.data(), b_out.data(), TW, TH);

    Matrix dst{m};
    planar_to_matrix(r_out.data(), g_out.data(), b_out.data(), dst);
    return dst;
}

} // namespace Filter
