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
                weights_out[i] = exp(-x * x * pi);
            }
        }
    }

    Matrix blur(Matrix m, const int radius)
    {
        // Precompute Gaussian weights once per call (avoid per-pixel recomputation)
        std::vector<double> __weights(radius + 1);
        Filter::Gauss::get_weights(radius, __weights.data());

        Matrix scratch(m.get_x_size(), m.get_y_size());
        Matrix dst(m.get_x_size(), m.get_y_size());

        // Horizontal pass
        for (auto y{0u}; y < m.get_y_size(); y++)
        {
            for (auto x{0u}; x < m.get_x_size(); x++)
            {
                // use precomputed __weights instead of recomputing per pixel
                auto r{__weights[0] * m.r(x, y)}, g{__weights[0] * m.g(x, y)}, b{__weights[0] * m.b(x, y)}, n{__weights[0]};

                for (auto wi{1}; wi <= radius; wi++)
                {
                    auto wc{__weights[wi]};
                    // left
                    if (x >= wi)
                    {
                        auto x2{x - wi};
                        r += wc * m.r(x2, y);
                        g += wc * m.g(x2, y);
                        b += wc * m.b(x2, y);
                        n += wc;
                    }
                    // right
                    auto x3{x + wi};
                    if (x3 < m.get_x_size())
                    {
                        r += wc * m.r(x3, y);
                        g += wc * m.g(x3, y);
                        b += wc * m.b(x3, y);
                        n += wc;
                    }
                }
                scratch.r(x, y) = r / n;
                scratch.g(x, y) = g / n;
                scratch.b(x, y) = b / n;
            }
        }

        // Vertical pass
        for (auto x{0u}; x < m.get_x_size(); x++)
        {
            for (auto y{0u}; y < m.get_y_size(); y++)
            {
                // use precomputed __weights instead of recomputing per pixel
                auto r{__weights[0] * scratch.r(x, y)}, g{__weights[0] * scratch.g(x, y)}, b{__weights[0] * scratch.b(x, y)}, n{__weights[0]};

                for (auto wi{1}; wi <= radius; wi++)
                {
                    auto wc{__weights[wi]};
                    // up
                    if (y >= wi)
                    {
                        auto y2{y - wi};
                        r += wc * scratch.r(x, y2);
                        g += wc * scratch.g(x, y2);
                        b += wc * scratch.b(x, y2);
                        n += wc;
                    }
                    // down
                    auto y3{y + wi};
                    if (y3 < m.get_y_size())
                    {
                        r += wc * scratch.r(x, y3);
                        g += wc * scratch.g(x, y3);
                        b += wc * scratch.b(x, y3);
                        n += wc;
                    }
                }
                dst.r(x, y) = r / n;
                dst.g(x, y) = g / n;
                dst.b(x, y) = b / n;
            }
        }

        return dst;
    }

}
