/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cmath>
#include <vector>

namespace Filter
{

    namespace Gauss
    {
        void get_weights(int n, double *weights_out)
        {
            for (auto i{0}; i <= n; i++)
            {
                double x{static_cast<double>(i) * max_x / n};
                weights_out[i] = std::exp(-x * x * pi);
            }
        }
    }

    // Helper: make a full symmetric kernel of length (2R+1) from half-weights [0..R]
    static inline void make_full_kernel(const std::vector<double>& halfW, int radius,
                                        std::vector<double>& fullW)
    {
        fullW.resize(static_cast<std::size_t>(2 * radius + 1));
        // halfW is already normalized for sum over [-R..R] to be 1:
        // fullW[k+R] = halfW[abs(k)]
        for (int k = -radius; k <= radius; ++k)
            fullW[static_cast<std::size_t>(k + radius)] = halfW[static_cast<std::size_t>(std::abs(k))];
    }

    Matrix blur(Matrix m, const int radius)
    {
        // === Step 1 (kept): precompute Gaussian half-kernel once per call ===
        std::vector<double> halfW(static_cast<std::size_t>(radius) + 1);
        Gauss::get_weights(radius, halfW.data());

        // Normalize so sum over k=-R..R equals 1.0
        if (radius > 0) {
            double total = halfW[0];
            for (int i = 1; i <= radius; ++i) total += 2.0 * halfW[i];
            const double inv = 1.0 / total;
            for (double& w : halfW) w *= inv;
        } else {
            halfW[0] = 1.0;
        }

        // Build a full symmetric kernel of length (2R+1) for branch-free inner loops
        std::vector<double> fullW;
        make_full_kernel(halfW, radius, fullW);

        // Keep original construction style (Matrix has no (w,h) ctor)
        Matrix scratch{m};
        Matrix dst{m};

        const unsigned W = m.get_x_size();
        const unsigned H = m.get_y_size();

        if (W == 0u || H == 0u) return m;

        // ---- Horizontal pass: per-row padded line buffers, then 1-D conv ----
        {
            const std::size_t pad = static_cast<std::size_t>(radius);
            const std::size_t lineLen = static_cast<std::size_t>(W) + 2 * pad;

            std::vector<double> rline(lineLen), gline(lineLen), bline(lineLen);

            for (unsigned y = 0; y < H; ++y)
            {
                // Fill middle from source row
                for (unsigned x = 0; x < W; ++x) {
                    rline[pad + x] = static_cast<double>(m.r(x, y));
                    gline[pad + x] = static_cast<double>(m.g(x, y));
                    bline[pad + x] = static_cast<double>(m.b(x, y));
                }

                // Edge replicate into the left padding
                const double rLeft = rline[pad + 0];
                const double gLeft = gline[pad + 0];
                const double bLeft = bline[pad + 0];
                for (std::size_t i = 0; i < pad; ++i) {
                    rline[i] = rLeft;
                    gline[i] = gLeft;
                    bline[i] = bLeft;
                }

                // Edge replicate into the right padding
                const double rRight = rline[pad + (W - 1)];
                const double gRight = gline[pad + (W - 1)];
                const double bRight = bline[pad + (W - 1)];
                for (std::size_t i = 0; i < pad; ++i) {
                    rline[pad + W + i] = rRight;
                    gline[pad + W + i] = gRight;
                    bline[pad + W + i] = bRight;
                }

                // Convolve into scratch
                // out(x) = sum_{k=-R..R} fullW[k+R] * line[(x + k) + pad]
                for (unsigned x = 0; x < W; ++x) {
                    const std::size_t base = pad + x;
                    double accR = 0.0, accG = 0.0, accB = 0.0;

                    // Tight inner loop: no clamp/abs, contiguous loads
                    for (int k = -radius; k <= radius; ++k) {
                        const double w = fullW[static_cast<std::size_t>(k + radius)];
                        const std::size_t idx = base + static_cast<std::size_t>(k);
                        accR += w * rline[idx];
                        accG += w * gline[idx];
                        accB += w * bline[idx];
                    }

                    scratch.r(x, y) = static_cast<unsigned char>(accR + 0.5); // rounded
                    scratch.g(x, y) = static_cast<unsigned char>(accG + 0.5);
                    scratch.b(x, y) = static_cast<unsigned char>(accB + 0.5);
                }
            }
        }

        // ---- Vertical pass: per-column padded line buffers, then 1-D conv ----
        {
            const std::size_t pad = static_cast<std::size_t>(radius);
            const std::size_t lineLen = static_cast<std::size_t>(H) + 2 * pad;

            std::vector<double> rcol(lineLen), gcol(lineLen), bcol(lineLen);

            for (unsigned x = 0; x < W; ++x)
            {
                // Fill middle from scratch column
                for (unsigned y = 0; y < H; ++y) {
                    rcol[pad + y] = static_cast<double>(scratch.r(x, y));
                    gcol[pad + y] = static_cast<double>(scratch.g(x, y));
                    bcol[pad + y] = static_cast<double>(scratch.b(x, y));
                }

                // Edge replicate into the top padding
                const double rTop = rcol[pad + 0];
                const double gTop = gcol[pad + 0];
                const double bTop = bcol[pad + 0];
                for (std::size_t i = 0; i < pad; ++i) {
                    rcol[i] = rTop;
                    gcol[i] = gTop;
                    bcol[i] = bTop;
                }

                // Edge replicate into the bottom padding
                const double rBot = rcol[pad + (H - 1)];
                const double gBot = gcol[pad + (H - 1)];
                const double bBot = bcol[pad + (H - 1)];
                for (std::size_t i = 0; i < pad; ++i) {
                    rcol[pad + H + i] = rBot;
                    gcol[pad + H + i] = gBot;
                    bcol[pad + H + i] = bBot;
                }

                // Convolve into dst
                for (unsigned y = 0; y < H; ++y) {
                    const std::size_t base = pad + y;
                    double accR = 0.0, accG = 0.0, accB = 0.0;

                    for (int k = -radius; k <= radius; ++k) {
                        const double w = fullW[static_cast<std::size_t>(k + radius)];
                        const std::size_t idx = base + static_cast<std::size_t>(k);
                        accR += w * rcol[idx];
                        accG += w * gcol[idx];
                        accB += w * bcol[idx];
                    }

                    dst.r(x, y) = static_cast<unsigned char>(accR + 0.5);
                    dst.g(x, y) = static_cast<unsigned char>(accG + 0.5);
                    dst.b(x, y) = static_cast<unsigned char>(accB + 0.5);
                }
            }
        }

        return dst;
    }

}
