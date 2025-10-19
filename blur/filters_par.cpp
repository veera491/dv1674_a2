/*
Author: David Holmqvist <daae19@student.bth.se>
Parallelized by: Nithin, VEERA
*/

#include "filters_par.hpp"
#include "matrix.hpp"
#include <pthread.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <vector>

namespace FilterPar {

namespace Gauss {

void get_weights(int n, double* weights_out)
{
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
    // 1) center copy
    for (unsigned y = 0; y < H; ++y) {
        const unsigned srcRow = y * W;
        const unsigned dstRow = (y + static_cast<unsigned>(R)) * PW + static_cast<unsigned>(R);
        std::copy(in + srcRow, in + srcRow + W, out + dstRow);
    }
    // 2) left/right reflect in each row
    for (unsigned y = 0; y < H; ++y) {
        const unsigned dstRow = (y + static_cast<unsigned>(R)) * PW;
        for (int k = 0; k < R; ++k) {
            out[dstRow + (static_cast<unsigned>(R) - 1u - static_cast<unsigned>(k))] =
                out[dstRow + (static_cast<unsigned>(R) + static_cast<unsigned>(k))];
            out[dstRow + (static_cast<unsigned>(R) + W + static_cast<unsigned>(k))] =
                out[dstRow + (static_cast<unsigned>(R) + W - 1u - static_cast<unsigned>(k))];
        }
    }
    // 3) top/bottom reflect entire padded rows
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
    const unsigned base = static_cast<unsigned>(radius);
    for (unsigned x = 0; x < W; ++x) {
        const unsigned cx = base + x;
        float acc = gw[0] * row_pad[cx];
        for (int k = 1; k <= radius; ++k) {
            const float w = gw[k];
            acc += w * ( row_pad[cx - static_cast<unsigned>(k)] +
                         row_pad[cx + static_cast<unsigned>(k)] );
        }
        dstT[x * H + y_out] = acc;
    }
}

static inline void gauss_row_padded_transposed_to_final(const float* __restrict row_pad,
                                                        float* __restrict dst_final,
                                                        unsigned TW, unsigned TH, unsigned y_tr,
                                                        const std::vector<float>& gw,
                                                        int radius, unsigned W)
{
    const unsigned base = static_cast<unsigned>(radius);
    for (unsigned x = 0; x < TW; ++x) {
        const unsigned cx = base + x;
        float acc = gw[0] * row_pad[cx];
        for (int k = 1; k <= radius; ++k) {
            const float w = gw[k];
            acc += w * ( row_pad[cx - static_cast<unsigned>(k)] +
                         row_pad[cx + static_cast<unsigned>(k)] );
        }
        const unsigned out_index = x * W + y_tr;
        dst_final[out_index] = acc;
    }
}

// ------------------------ Thread Data Structures ------------------------

struct Pass1ThreadData {
    const float* rp;
    const float* gp;
    const float* bp;
    float* r_tr;
    float* g_tr;
    float* b_tr;
    unsigned W;
    unsigned H;
    unsigned PW;
    const std::vector<float>* gw;
    int radius;
    unsigned start_row;
    unsigned end_row;
};

struct Pass2ThreadData {
    const float* rtp;
    const float* gtp;
    const float* btp;
    float* r_out;
    float* g_out;
    float* b_out;
    unsigned TW;
    unsigned TH;
    unsigned TPW;
    unsigned W;  // final width
    const std::vector<float>* gw;
    int radius;
    unsigned start_row;
    unsigned end_row;
};

// ------------------------ Thread Workers ------------------------

void* pass1_worker(void* arg) {
    Pass1ThreadData* data = static_cast<Pass1ThreadData*>(arg);

    for (unsigned y = data->start_row; y < data->end_row; ++y) {
        const unsigned prow = (y + static_cast<unsigned>(data->radius)) * data->PW;

        gauss_row_padded_to_transposed(
            data->rp + prow, data->r_tr,
            data->W, data->H, y,
            *data->gw, data->radius
        );

        gauss_row_padded_to_transposed(
            data->gp + prow, data->g_tr,
            data->W, data->H, y,
            *data->gw, data->radius
        );

        gauss_row_padded_to_transposed(
            data->bp + prow, data->b_tr,
            data->W, data->H, y,
            *data->gw, data->radius
        );
    }

    return nullptr;
}

void* pass2_worker(void* arg) {
    Pass2ThreadData* data = static_cast<Pass2ThreadData*>(arg);

    for (unsigned y_tr = data->start_row; y_tr < data->end_row; ++y_tr) {
        const unsigned prow = (y_tr + static_cast<unsigned>(data->radius)) * data->TPW;

        gauss_row_padded_transposed_to_final(
            data->rtp + prow, data->r_out,
            data->TW, data->TH, y_tr,
            *data->gw, data->radius, data->W
        );

        gauss_row_padded_transposed_to_final(
            data->gtp + prow, data->g_out,
            data->TW, data->TH, y_tr,
            *data->gw, data->radius, data->W
        );

        gauss_row_padded_transposed_to_final(
            data->btp + prow, data->b_out,
            data->TW, data->TH, y_tr,
            *data->gw, data->radius, data->W
        );
    }

    return nullptr;
}

// ------------------------ main blur -----------------------------------------

Matrix blur(Matrix m, const int radius, int num_threads)
{
    const unsigned W = m.get_x_size();
    const unsigned H = m.get_y_size();
    if (W == 0 || H == 0 || radius <= 0) {
        return m;
    }

    // Clamp num_threads to reasonable values
    if (num_threads < 1) num_threads = 1;
    if (num_threads > static_cast<int>(H)) num_threads = H;

    // 1) Precompute Gaussian weights and normalize
    std::vector<double> gw_d(static_cast<std::size_t>(radius) + 1);
    Gauss::get_weights(radius, gw_d.data());

    double sum_full = gw_d[0];
    for (int k = 1; k <= radius; ++k) sum_full += 2.0 * gw_d[k];
    const float inv_sum = static_cast<float>(1.0 / sum_full);

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

    // 4) Pass 1: horizontal on padded rows (PARALLEL)
    std::vector<float> r_tr(N), g_tr(N), b_tr(N);

    std::vector<pthread_t> threads(num_threads);
    std::vector<Pass1ThreadData> thread_data(num_threads);

    unsigned rows_per_thread = H / num_threads;
    unsigned remaining_rows = H % num_threads;

    unsigned current_row = 0;
    for (int t = 0; t < num_threads; ++t) {
        thread_data[t].rp = rp.data();
        thread_data[t].gp = gp.data();
        thread_data[t].bp = bp.data();
        thread_data[t].r_tr = r_tr.data();
        thread_data[t].g_tr = g_tr.data();
        thread_data[t].b_tr = b_tr.data();
        thread_data[t].W = W;
        thread_data[t].H = H;
        thread_data[t].PW = PW;
        thread_data[t].gw = &gw;
        thread_data[t].radius = radius;
        thread_data[t].start_row = current_row;

        unsigned rows_for_this_thread = rows_per_thread + (t < static_cast<int>(remaining_rows) ? 1 : 0);
        thread_data[t].end_row = current_row + rows_for_this_thread;
        current_row += rows_for_this_thread;

        pthread_create(&threads[t], nullptr, pass1_worker, &thread_data[t]);
    }

    for (int t = 0; t < num_threads; ++t) {
        pthread_join(threads[t], nullptr);
    }

    // 5) Reflected padding for the transposed image (TW x TH) = (H x W)
    const unsigned TW = H, TH = W;
    const unsigned TPW = TW + 2u * static_cast<unsigned>(radius);
    const unsigned TPH = TH + 2u * static_cast<unsigned>(radius);
    const std::size_t TNP = static_cast<std::size_t>(TPW) * static_cast<std::size_t>(TPH);
    std::vector<float> rtp(TNP), gtp(TNP), btp(TNP);
    reflect_pad_plane(r_tr.data(), rtp.data(), TW, TH, radius);
    reflect_pad_plane(g_tr.data(), gtp.data(), TW, TH, radius);
    reflect_pad_plane(b_tr.data(), btp.data(), TW, TH, radius);

    // 6) Pass 2: horizontal on padded transposed rows (PARALLEL)
    std::vector<float> r_out(N), g_out(N), b_out(N);

    std::vector<Pass2ThreadData> thread_data2(num_threads);

    // Clamp threads for pass 2
    int num_threads_pass2 = num_threads;
    if (num_threads_pass2 > static_cast<int>(TH)) num_threads_pass2 = TH;

    rows_per_thread = TH / num_threads_pass2;
    remaining_rows = TH % num_threads_pass2;

    current_row = 0;
    for (int t = 0; t < num_threads_pass2; ++t) {
        thread_data2[t].rtp = rtp.data();
        thread_data2[t].gtp = gtp.data();
        thread_data2[t].btp = btp.data();
        thread_data2[t].r_out = r_out.data();
        thread_data2[t].g_out = g_out.data();
        thread_data2[t].b_out = b_out.data();
        thread_data2[t].TW = TW;
        thread_data2[t].TH = TH;
        thread_data2[t].TPW = TPW;
        thread_data2[t].W = W;
        thread_data2[t].gw = &gw;
        thread_data2[t].radius = radius;
        thread_data2[t].start_row = current_row;

        unsigned rows_for_this_thread = rows_per_thread + (t < static_cast<int>(remaining_rows) ? 1 : 0);
        thread_data2[t].end_row = current_row + rows_for_this_thread;
        current_row += rows_for_this_thread;

        pthread_create(&threads[t], nullptr, pass2_worker, &thread_data2[t]);
    }

    for (int t = 0; t < num_threads_pass2; ++t) {
        pthread_join(threads[t], nullptr);
    }

    // 7) Pack back to Matrix
    Matrix dst{m};
    planar_to_matrix(r_out.data(), g_out.data(), b_out.data(), dst);
    return dst;
}

} // namespace FilterPar