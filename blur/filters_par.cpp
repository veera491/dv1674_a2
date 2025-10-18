#include "filters_par.hpp"
#include "matrix.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <pthread.h>
#include <vector>

namespace Filter {

// ---- Original Gauss::get_weights (kept here to avoid linker issues) ----
namespace Gauss {
static inline void get_weights(int n, double* weights_out) {
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

static inline void reflect_pad_plane(const float* __restrict in,
                                     float* __restrict out,
                                     unsigned W, unsigned H, int R)
{
    const unsigned PW = W + 2u * static_cast<unsigned>(R);
    // center
    for (unsigned y = 0; y < H; ++y) {
        const unsigned srcRow = y * W;
        const unsigned dstRow = (y + static_cast<unsigned>(R)) * PW + static_cast<unsigned>(R);
        std::copy(in + srcRow, in + srcRow + W, out + dstRow);
    }
    // left/right reflect
    for (unsigned y = 0; y < H; ++y) {
        const unsigned dstRow = (y + static_cast<unsigned>(R)) * PW;
        for (int k = 0; k < R; ++k) {
            out[dstRow + (static_cast<unsigned>(R) - 1u - static_cast<unsigned>(k))] =
                out[dstRow + (static_cast<unsigned>(R) + static_cast<unsigned>(k))];
            out[dstRow + (static_cast<unsigned>(R) + W + static_cast<unsigned>(k))] =
                out[dstRow + (static_cast<unsigned>(R) + W - 1u - static_cast<unsigned>(k))];
        }
    }
    // top/bottom reflect
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
            acc += w * ( row_pad[cx - static_cast<unsigned>(k)]
                       + row_pad[cx + static_cast<unsigned>(k)] );
        }
        dstT[x * H + y_out] = acc; // transpose write
    }
}

static inline void gauss_row_padded_transposed_to_final(const float* __restrict row_pad,
                                                        float* __restrict dst_final,
                                                        unsigned TW, unsigned TH, unsigned y_tr,
                                                        const std::vector<float>& gw,
                                                        int radius,
                                                        unsigned W /*final width*/)
{
    const unsigned base = static_cast<unsigned>(radius);
    for (unsigned x = 0; x < TW; ++x) {            // x == original y
        const unsigned cx = base + x;
        float acc = gw[0] * row_pad[cx];
        #pragma omp simd
        for (int k = 1; k <= radius; ++k) {
            const float w = gw[k];
            acc += w * ( row_pad[cx - static_cast<unsigned>(k)]
                       + row_pad[cx + static_cast<unsigned>(k)] );
        }
        // Correct mapping: (x,y_tr) in transposed => (y=x, x=y_tr) in final
        const unsigned out_index = x * W + y_tr;
        dst_final[out_index] = acc;
    }
}

// ------------------------ pthreads plumbing ------------------------

struct Pass1Args {
    const float* rp; const float* gp; const float* bp;
    float* r_tr; float* g_tr; float* b_tr;
    unsigned W, H, PW;
    const std::vector<float>* gw;
    int R;
    unsigned y0, y1; // [y0, y1)
};

static void* pass1_worker(void* argp) {
    auto* a = static_cast<Pass1Args*>(argp);
    for (unsigned y = a->y0; y < a->y1; ++y) {
        const unsigned prow = (y + static_cast<unsigned>(a->R)) * a->PW;
        gauss_row_padded_to_transposed(a->rp + prow, a->r_tr, a->W, a->H, y, *a->gw, a->R);
        gauss_row_padded_to_transposed(a->gp + prow, a->g_tr, a->W, a->H, y, *a->gw, a->R);
        gauss_row_padded_to_transposed(a->bp + prow, a->b_tr, a->W, a->H, y, *a->gw, a->R);
    }
    return nullptr;
}

struct Pass2Args {
    const float* rtp; const float* gtp; const float* btp;
    float* r_out; float* g_out; float* b_out;
    unsigned TW, TH, TPW, W;
    const std::vector<float>* gw;
    int R;
    unsigned y0, y1; // [y0, y1)
};

static void* pass2_worker(void* argp) {
    auto* a = static_cast<Pass2Args*>(argp);
    for (unsigned y_tr = a->y0; y_tr < a->y1; ++y_tr) {
        const unsigned prow = (y_tr + static_cast<unsigned>(a->R)) * a->TPW;
        gauss_row_padded_transposed_to_final(a->rtp + prow, a->r_out, a->TW, a->TH, y_tr, *a->gw, a->R, a->W);
        gauss_row_padded_transposed_to_final(a->gtp + prow, a->g_out, a->TW, a->TH, y_tr, *a->gw, a->R, a->W);
        gauss_row_padded_transposed_to_final(a->btp + prow, a->b_out, a->TW, a->TH, y_tr, *a->gw, a->R, a->W);
    }
    return nullptr;
}

// ------------------------ main (parallel) -----------------------------------

Matrix blur_par(Matrix m, int radius, int num_threads)
{
    const unsigned W = m.get_x_size();
    const unsigned H = m.get_y_size();
    if (W == 0 || H == 0 || radius <= 0) {
        return m;
    }
    if (num_threads < 1) num_threads = 1;
    if ((unsigned)num_threads > H) num_threads = (int)H; // cap

    // 1) Precompute Gaussian weights and normalize to 1.0 (avoid per-pixel division).
    std::vector<double> gw_d((std::size_t)radius + 1);
    Gauss::get_weights(radius, gw_d.data());
    double sum_full = gw_d[0];
    for (int k = 1; k <= radius; ++k) sum_full += 2.0 * gw_d[k];
    const float inv_sum = static_cast<float>(1.0 / sum_full);
    std::vector<float> gw(gw_d.size());
    for (int k = 0; k <= radius; ++k) gw[k] = static_cast<float>(gw_d[k]) * inv_sum;

    // 2) Planar
    const std::size_t N = (std::size_t)W * (std::size_t)H;
    std::vector<float> r_in(N), g_in(N), b_in(N);
    matrix_to_planar(m, r_in.data(), g_in.data(), b_in.data());

    // 3) Padding for original (W x H)
    const unsigned PW = W + 2u * (unsigned)radius;
    const unsigned PH = H + 2u * (unsigned)radius;
    const std::size_t NP = (std::size_t)PW * (std::size_t)PH;
    std::vector<float> rp(NP), gp(NP), bp(NP);
    reflect_pad_plane(r_in.data(), rp.data(), W, H, radius);
    reflect_pad_plane(g_in.data(), gp.data(), W, H, radius);
    reflect_pad_plane(b_in.data(), bp.data(), W, H, radius);

    // 4) Pass 1: horizontal rows → write transposed (HxW)
    std::vector<float> r_tr(N), g_tr(N), b_tr(N);

    {
        std::vector<pthread_t> th(num_threads);
        std::vector<Pass1Args> args(num_threads);

        unsigned rows_per = H / (unsigned)num_threads;
        unsigned rem = H % (unsigned)num_threads;
        unsigned y = 0;
        for (int t = 0; t < num_threads; ++t) {
            unsigned take = rows_per + (rem ? 1u : 0u);
            if (rem) --rem;
            args[t] = Pass1Args{
                rp.data(), gp.data(), bp.data(),
                r_tr.data(), g_tr.data(), b_tr.data(),
                W, H, PW, &gw, radius,
                y, y + take
            };
            pthread_create(&th[t], nullptr, pass1_worker, &args[t]);
            y += take;
        }
        for (int t = 0; t < num_threads; ++t) pthread_join(th[t], nullptr);
    }

    // 5) Padding for transposed (TW x TH) = (H x W)
    const unsigned TW = H, TH = W;
    const unsigned TPW = TW + 2u * (unsigned)radius;
    const unsigned TPH = TH + 2u * (unsigned)radius;
    const std::size_t TNP = (std::size_t)TPW * (std::size_t)TPH;
    std::vector<float> rtp(TNP), gtp(TNP), btp(TNP);
    reflect_pad_plane(r_tr.data(), rtp.data(), TW, TH, radius);
    reflect_pad_plane(g_tr.data(), gtp.data(), TW, TH, radius);
    reflect_pad_plane(b_tr.data(), btp.data(), TW, TH, radius);

    // 6) Pass 2: horizontal rows on transposed → write final (WxH)
    std::vector<float> r_out(N), g_out(N), b_out(N);

    {
        std::vector<pthread_t> th(num_threads);
        std::vector<Pass2Args> args(num_threads);

        unsigned rows_per = TH / (unsigned)num_threads;
        unsigned rem = TH % (unsigned)num_threads;
        unsigned y = 0;
        for (int t = 0; t < num_threads; ++t) {
            unsigned take = rows_per + (rem ? 1u : 0u);
            if (rem) --rem;
            args[t] = Pass2Args{
                rtp.data(), gtp.data(), btp.data(),
                r_out.data(), g_out.data(), b_out.data(),
                TW, TH, TPW, W, &gw, radius,
                y, y + take
            };
            pthread_create(&th[t], nullptr, pass2_worker, &args[t]);
            y += take;
        }
        for (int t = 0; t < num_threads; ++t) pthread_join(th[t], nullptr);
    }

    Matrix dst{m};
    planar_to_matrix(r_out.data(), g_out.data(), b_out.data(), dst);
    return dst;
}

} // namespace Filter
