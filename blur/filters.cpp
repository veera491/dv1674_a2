/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cmath>
#include <vector>
#include <algorithm> // std::min, std::max

namespace Filter
{

    namespace Gauss
    {
        // As in your original: builds weights[0..radius] where weights[i] is the i-th offset weight
        void get_weights(int n, double *weights_out)
        {
            for (auto i{0}; i <= n; i++)
            {
                double x{static_cast<double>(i) * max_x / n};
                weights_out[i] = std::exp(-x * x * pi);
            }
        }
    }

    // --- helpers: planar conversion and transpose (float planes for speed) ---

    struct PlanesF {
        float* __restrict r;
        float* __restrict g;
        float* __restrict b;
        std::size_t size; // W*H
    };

    static inline void matrix_to_planar(const Matrix& m, float* __restrict r,
                                        float* __restrict g, float* __restrict b)
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
                out.r(x, y) = static_cast<unsigned char>(std::clamp(r[idx], 0.0f, 255.0f));
                out.g(x, y) = static_cast<unsigned char>(std::clamp(g[idx], 0.0f, 255.0f));
                out.b(x, y) = static_cast<unsigned char>(std::clamp(b[idx], 0.0f, 255.0f));
            }
        }
    }

    // transpose: (W x H) -> (H x W) buffers
    static inline void transpose_plane(const float* __restrict src, float* __restrict dst,
                                       unsigned W, unsigned H)
    {
        for (unsigned y = 0; y < H; ++y) {
            const unsigned row = y * W;
            for (unsigned x = 0; x < W; ++x) {
                dst[x * H + y] = src[row + x];
            }
        }
    }

    // Horizontal 1D Gaussian pass on a single plane (row-wise).
    // Edges handled separately so the interior has no bounds checks.
    static inline void horizontal_gauss_plane(const float* __restrict src,
                                              float* __restrict dst,
                                              unsigned W, unsigned H,
                                              const std::vector<double>& gw,
                                              int radius,
                                              double norm_full)
    {
        // Left edge region width (clamped access)
        const unsigned edge = std::min<unsigned>(W, static_cast<unsigned>(radius));

        for (unsigned y = 0; y < H; ++y) {
            const unsigned row = y * W;

            // ---- Left edge: x in [0, edge-1]
            for (unsigned x = 0; x < edge; ++x) {
                double acc = gw[0] * src[row + x];
                double nrm = gw[0];

                // guarded neighbors
                for (int k = 1; k <= radius; ++k) {
                    const unsigned xl = (x >= static_cast<unsigned>(k)) ? (x - static_cast<unsigned>(k)) : 0u;
                    const unsigned xr = std::min<unsigned>(x + static_cast<unsigned>(k), W - 1);
                    const double w = gw[k];
                    if (xl != x) { acc += w * src[row + xl]; nrm += w; }
                    if (xr != x) { acc += w * src[row + xr]; nrm += w; }
                }
                dst[row + x] = static_cast<float>(acc / nrm);
            }

            // ---- Interior: x in [radius, W-1-radius] (no bounds checks)
            if (W > static_cast<unsigned>(2 * radius)) {
                const unsigned x0 = static_cast<unsigned>(radius);
                const unsigned x1 = W - static_cast<unsigned>(radius);
                for (unsigned x = x0; x < x1; ++x) {
                    double acc = gw[0] * src[row + x];
                    // symmetric taps: pair left/right
                    #pragma omp simd
                    for (int k = 1; k <= radius; ++k) {
                        const double w = gw[k];
                        acc += w * ( src[row + (x - static_cast<unsigned>(k))] +
                                     src[row + (x + static_cast<unsigned>(k))] );
                    }
                    dst[row + x] = static_cast<float>(acc / norm_full);
                }
            }

            // ---- Right edge: x in [max(radius, W-radius) .. W-1]
            const unsigned xr_start = (W > static_cast<unsigned>(radius)) ? (W - static_cast<unsigned>(radius)) : 0u;
            for (unsigned x = xr_start; x < W; ++x) {
                double acc = gw[0] * src[row + x];
                double nrm = gw[0];

                for (int k = 1; k <= radius; ++k) {
                    const unsigned xl = (x >= static_cast<unsigned>(k)) ? (x - static_cast<unsigned>(k)) : 0u;
                    const unsigned xr = std::min<unsigned>(x + static_cast<unsigned>(k), W - 1);
                    const double w = gw[k];
                    if (xl != x) { acc += w * src[row + xl]; nrm += w; }
                    if (xr != x) { acc += w * src[row + xr]; nrm += w; }
                }
                dst[row + x] = static_cast<float>(acc / nrm);
            }
        }
    }

    Matrix blur(Matrix m, const int radius)
    {
        // === (already done earlier) Precompute Gaussian weights once per call ===
        std::vector<double> gw(static_cast<std::size_t>(radius) + 1);
        Gauss::get_weights(radius, gw.data());

        const unsigned W = m.get_x_size();
        const unsigned H = m.get_y_size();

        // Quick outs
        if (W == 0 || H == 0 || radius <= 0) {
            return m;
        }

        // Precompute full-window normalization for interior (symmetric kernel)
        double norm_full = gw[0];
        for (int i = 1; i <= radius; ++i) norm_full += 2.0 * gw[i];

        // === Planar (SoA) buffers ===
        const std::size_t N = static_cast<std::size_t>(W) * static_cast<std::size_t>(H);
        std::vector<float> r_in(N), g_in(N), b_in(N);
        std::vector<float> r_tmp(N), g_tmp(N), b_tmp(N);               // after horizontal
        std::vector<float> r_tr(N),  g_tr(N),  b_tr(N);                // transposed
        std::vector<float> r_out_tr(N), g_out_tr(N), b_out_tr(N);      // second pass (on transposed)
        std::vector<float> r_out(N), g_out(N), b_out(N);               // final after transpose back

        matrix_to_planar(m, r_in.data(), g_in.data(), b_in.data());

        // === Pass 1: Horizontal on rows ===
        horizontal_gauss_plane(r_in.data(), r_tmp.data(), W, H, gw, radius, norm_full);
        horizontal_gauss_plane(g_in.data(), g_tmp.data(), W, H, gw, radius, norm_full);
        horizontal_gauss_plane(b_in.data(), b_tmp.data(), W, H, gw, radius, norm_full);

        // === Transpose so that columns become rows ===
        transpose_plane(r_tmp.data(), r_tr.data(), W, H);
        transpose_plane(g_tmp.data(), g_tr.data(), W, H);
        transpose_plane(b_tmp.data(), b_tr.data(), W, H);

        // Now dimensions are (H x W)
        const unsigned TW = H;
        const unsigned TH = W;

        // === Pass 2: "Vertical" as another horizontal over transposed ===
        horizontal_gauss_plane(r_tr.data(), r_out_tr.data(), TW, TH, gw, radius, norm_full);
        horizontal_gauss_plane(g_tr.data(), g_out_tr.data(), TW, TH, gw, radius, norm_full);
        horizontal_gauss_plane(b_tr.data(), b_out_tr.data(), TW, TH, gw, radius, norm_full);

        // === Transpose back to original layout ===
        transpose_plane(r_out_tr.data(), r_out.data(), TW, TH);
        transpose_plane(g_out_tr.data(), g_out.data(), TW, TH);
        transpose_plane(b_out_tr.data(), b_out.data(), TW, TH);

        // === Pack back to Matrix ===
        Matrix dst{m}; // keep same metadata (sizes, color_max)
        planar_to_matrix(r_out.data(), g_out.data(), b_out.data(), dst);
        return dst;
    }

}
