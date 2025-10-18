/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cmath>
#include <vector>
#include <algorithm> // std::clamp

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

    Matrix blur(Matrix m, const int radius)
    {
        // === Step 1 (kept): precompute Gaussian weights once per call ===
        std::vector<double> weights(static_cast<std::size_t>(radius) + 1);
        Gauss::get_weights(radius, weights.data());

        // Normalize the symmetric 1D kernel once so sum(-r..r)=1
        // total = w[0] + 2*(w[1]+...+w[r])
        if (radius > 0) {
            double total = weights[0];
            for (int i = 1; i <= radius; ++i) total += 2.0 * weights[i];
            const double inv = 1.0 / total;
            for (double& w : weights) w *= inv;
        } else {
            weights[0] = 1.0;
        }

        // Keep original construction style (Matrix has no (w,h) ctor)
        Matrix scratch{m};
        Matrix dst{m};

        // Cache sizes as signed ints for math, cast back on access
        const int W = static_cast<int>(m.get_x_size());
        const int H = static_cast<int>(m.get_y_size());

        // -----------------------
        // Horizontal 1-D pass
        // -----------------------
        for (int y = 0; y < H; ++y)
        {
            const unsigned uy = static_cast<unsigned>(y);
            for (int x = 0; x < W; ++x)
            {
                double accR = 0.0, accG = 0.0, accB = 0.0;

                // full 1-D kernel with clamped borders; index by abs(k)
                for (int k = -radius; k <= radius; ++k)
                {
                    const int xx  = std::clamp(x + k, 0, W - 1);
                    const unsigned uxx = static_cast<unsigned>(xx);
                    const double w = weights[std::abs(k)];

                    accR += w * static_cast<double>(m.r(uxx, uy));
                    accG += w * static_cast<double>(m.g(uxx, uy));
                    accB += w * static_cast<double>(m.b(uxx, uy));
                }

                const unsigned ux = static_cast<unsigned>(x);
                scratch.r(ux, uy) = static_cast<unsigned char>(std::clamp(accR, 0.0, 255.0));
                scratch.g(ux, uy) = static_cast<unsigned char>(std::clamp(accG, 0.0, 255.0));
                scratch.b(ux, uy) = static_cast<unsigned char>(std::clamp(accB, 0.0, 255.0));
            }
        }

        // -----------------------
        // Vertical 1-D pass
        // -----------------------
        for (int x = 0; x < W; ++x)
        {
            const unsigned ux = static_cast<unsigned>(x);
            for (int y = 0; y < H; ++y)
            {
                double accR = 0.0, accG = 0.0, accB = 0.0;

                for (int k = -radius; k <= radius; ++k)
                {
                    const int yy  = std::clamp(y + k, 0, H - 1);
                    const unsigned uyy = static_cast<unsigned>(yy);
                    const double w = weights[std::abs(k)];

                    accR += w * static_cast<double>(scratch.r(ux, uyy));
                    accG += w * static_cast<double>(scratch.g(ux, uyy));
                    accB += w * static_cast<double>(scratch.b(ux, uyy));
                }

                const unsigned uy = static_cast<unsigned>(y);
                dst.r(ux, uy) = static_cast<unsigned char>(std::clamp(accR, 0.0, 255.0));
                dst.g(ux, uy) = static_cast<unsigned char>(std::clamp(accG, 0.0, 255.0));
                dst.b(ux, uy) = static_cast<unsigned char>(std::clamp(accB, 0.0, 255.0));
            }
        }

        return dst;
    }

}
