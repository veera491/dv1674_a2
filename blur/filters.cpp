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
        // unchanged: builds weights[0..radius] once per radius
        void get_weights(int n, double *weights_out)
        {
            for (auto i{0}; i <= n; i++)
            {
                double x{static_cast<double>(i) * max_x / n};
                weights_out[i] = std::exp(-x * x * pi);
            }
        }
    }

    // --- Optional: O(1) per-pixel box-blur (rolling sum) fast path -----------
    // OFF by default to keep output identical to your Gaussian implementation.
    // Enable by compiling with:  -DFILTER_BOX_BLUR
#ifdef FILTER_BOX_BLUR
    static void horizontal_box(const Matrix& in, Matrix& tmp, int radius)
    {
        const unsigned W = in.get_x_size();
        const unsigned H = in.get_y_size();
        if (radius == 0 || W == 0 || H == 0) { tmp = in; return; }

        const int win = 2*radius + 1;
        for (unsigned y = 0; y < H; ++y)
        {
            // rolling sums for each channel
            double sr = 0.0, sg = 0.0, sb = 0.0;

            // initialize window [0 .. radius] with clamping to edges
            for (int dx = -radius; dx <= radius; ++dx)
            {
                unsigned xx = (dx < 0) ? 0u : (dx >= (int)W ? W-1 : (unsigned)dx);
                sr += in.r(xx, y);
                sg += in.g(xx, y);
                sb += in.b(xx, y);
            }
            tmp.r(0, y) = sr / win;
            tmp.g(0, y) = sg / win;
            tmp.b(0, y) = sb / win;

            for (unsigned x = 1; x < W; ++x)
            {
                int xl = (int)x - radius - 1;
                int xr = (int)x + radius;
                unsigned x_out = (xr >= (int)W) ? (W - 1) : (unsigned)xr;
                unsigned x_in  = (xl < 0)       ? 0u       : (unsigned)xl;

                sr += in.r(x_out, y) - in.r(x_in, y);
                sg += in.g(x_out, y) - in.g(x_in, y);
                sb += in.b(x_out, y) - in.b(x_in, y);

                tmp.r(x, y) = sr / win;
                tmp.g(x, y) = sg / win;
                tmp.b(x, y) = sb / win;
            }
        }
    }

    static void vertical_box(const Matrix& tmp, Matrix& out, int radius)
    {
        const unsigned W = tmp.get_x_size();
        const unsigned H = tmp.get_y_size();
        if (radius == 0 || W == 0 || H == 0) { out = tmp; return; }

        const int win = 2*radius + 1;
        for (unsigned x = 0; x < W; ++x)
        {
            double sr = 0.0, sg = 0.0, sb = 0.0;

            for (int dy = -radius; dy <= radius; ++dy)
            {
                unsigned yy = (dy < 0) ? 0u : (dy >= (int)H ? H-1 : (unsigned)dy);
                sr += tmp.r(x, yy);
                sg += tmp.g(x, yy);
                sb += tmp.b(x, yy);
            }
            out.r(x, 0) = sr / win;
            out.g(x, 0) = sg / win;
            out.b(x, 0) = sb / win;

            for (unsigned y = 1; y < H; ++y)
            {
                int yu = (int)y - radius - 1;
                int yd = (int)y + radius;
                unsigned y_out = (yd >= (int)H) ? (H - 1) : (unsigned)yd;
                unsigned y_in  = (yu < 0)       ? 0u       : (unsigned)yu;

                sr += tmp.r(x, y_out) - tmp.r(x, y_in);
                sg += tmp.g(x, y_out) - tmp.g(x, y_in);
                sb += tmp.b(x, y_out) - tmp.b(x, y_in);

                out.r(x, y) = sr / win;
                out.g(x, y) = sg / win;
                out.b(x, y) = sb / win;
            }
        }
    }
#endif
    // ------------------------------------------------------------------------

    Matrix blur(Matrix m, const int radius)
    {
        // === from step #1: precompute Gaussian weights once per call ===
        std::vector<double> gw(static_cast<std::size_t>(radius) + 1);
        Gauss::get_weights(radius, gw.data());

        // cache image geometry (step #3: kill accessor overhead for sizes)
        const unsigned W = m.get_x_size();
        const unsigned H = m.get_y_size();

        Matrix scratch{m}; // will hold horizontal pass
        Matrix dst{m};     // final

#ifdef FILTER_BOX_BLUR
        // --- Step #4: optional O(1) box blur path (kept OFF by default) ---
        horizontal_box(m, scratch, radius);
        vertical_box(scratch, dst, radius);
        return dst;
#endif

        // === Steps #2, #5 for Gaussian (separable + branchless interior) ===
        // Precompute full-window normalization for interior (no border loss).
        // For symmetric kernel: norm_full = w0 + 2 * sum(wi, i=1..radius)
        double norm_full = gw[0];
        for (int i = 1; i <= radius; ++i) norm_full += 2.0 * gw[i];

        // ----------------- Horizontal pass -----------------
        for (unsigned y = 0; y < H; ++y)
        {
            // Left edge: x in [0, radius-1] → needs boundary checks
            for (unsigned x = 0; x < std::min<unsigned>(W, (unsigned)radius); ++x)
            {
                double r = gw[0] * m.r(x, y);
                double g = gw[0] * m.g(x, y);
                double b = gw[0] * m.b(x, y);
                double n = gw[0];

                for (int k = 1; k <= radius; ++k)
                {
                    double w = gw[k];

                    // left neighbor (clamped)
                    if (x >= (unsigned)k) {
                        unsigned xx = x - (unsigned)k;
                        r += w * m.r(xx, y);
                        g += w * m.g(xx, y);
                        b += w * m.b(xx, y);
                        n += w;
                    }

                    // right neighbor (clamped)
                    unsigned xr = x + (unsigned)k;
                    if (xr < W) {
                        r += w * m.r(xr, y);
                        g += w * m.g(xr, y);
                        b += w * m.b(xr, y);
                        n += w;
                    }
                }

                scratch.r(x, y) = r / n;
                scratch.g(x, y) = g / n;
                scratch.b(x, y) = b / n;
            }

            // Interior: x in [radius, W-1-radius] → no bounds checks
            if (W > (unsigned)(2*radius))
            {
                for (unsigned x = (unsigned)radius; x < W - (unsigned)radius; ++x)
                {
                    double r = gw[0] * m.r(x, y);
                    double g = gw[0] * m.g(x, y);
                    double b = gw[0] * m.b(x, y);

                    for (int k = 1; k <= radius; ++k)
                    {
                        double w = gw[k];
                        unsigned xl = x - (unsigned)k;
                        unsigned xr = x + (unsigned)k;

                        r += w * ( m.r(xl, y) + m.r(xr, y) );
                        g += w * ( m.g(xl, y) + m.g(xr, y) );
                        b += w * ( m.b(xl, y) + m.b(xr, y) );
                    }

                    // no need to compute n per pixel, use constant full-window sum
                    scratch.r(x, y) = r / norm_full;
                    scratch.g(x, y) = g / norm_full;
                    scratch.b(x, y) = b / norm_full;
                }
            }

            // Right edge: x in [max(radius, W-radius) .. W-1] → boundary checks
            unsigned x_start = (W > (unsigned)radius) ? (W - (unsigned)radius) : 0u;
            for (unsigned x = x_start; x < W; ++x)
            {
                double r = gw[0] * m.r(x, y);
                double g = gw[0] * m.g(x, y);
                double b = gw[0] * m.b(x, y);
                double n = gw[0];

                for (int k = 1; k <= radius; ++k)
                {
                    double w = gw[k];

                    // left
                    if (x >= (unsigned)k) {
                        unsigned xl = x - (unsigned)k;
                        r += w * m.r(xl, y);
                        g += w * m.g(xl, y);
                        b += w * m.b(xl, y);
                        n += w;
                    }

                    // right
                    unsigned xr = x + (unsigned)k;
                    if (xr < W) {
                        r += w * m.r(xr, y);
                        g += w * m.g(xr, y);
                        b += w * m.b(xr, y);
                        n += w;
                    }
                }

                scratch.r(x, y) = r / n;
                scratch.g(x, y) = g / n;
                scratch.b(x, y) = b / n;
            }
        }

        // ----------------- Vertical pass -----------------
        for (unsigned x = 0; x < W; ++x)
        {
            // Top edge: y in [0, radius-1]
            for (unsigned y = 0; y < std::min<unsigned>(H, (unsigned)radius); ++y)
            {
                double r = gw[0] * scratch.r(x, y);
                double g = gw[0] * scratch.g(x, y);
                double b = gw[0] * scratch.b(x, y);
                double n = gw[0];

                for (int k = 1; k <= radius; ++k)
                {
                    double w = gw[k];

                    // up
                    if (y >= (unsigned)k) {
                        unsigned yu = y - (unsigned)k;
                        r += w * scratch.r(x, yu);
                        g += w * scratch.g(x, yu);
                        b += w * scratch.b(x, yu);
                        n += w;
                    }

                    // down
                    unsigned yd = y + (unsigned)k;
                    if (yd < H) {
                        r += w * scratch.r(x, yd);
                        g += w * scratch.g(x, yd);
                        b += w * scratch.b(x, yd);
                        n += w;
                    }
                }

                dst.r(x, y) = r / n;
                dst.g(x, y) = g / n;
                dst.b(x, y) = b / n;
            }

            // Interior: y in [radius, H-1-radius] → no bounds checks
            if (H > (unsigned)(2*radius))
            {
                for (unsigned y = (unsigned)radius; y < H - (unsigned)radius; ++y)
                {
                    double r = gw[0] * scratch.r(x, y);
                    double g = gw[0] * scratch.g(x, y);
                    double b = gw[0] * scratch.b(x, y);

                    for (int k = 1; k <= radius; ++k)
                    {
                        double w = gw[k];
                        unsigned yu = y - (unsigned)k;
                        unsigned yd = y + (unsigned)k;

                        r += w * ( scratch.r(x, yu) + scratch.r(x, yd) );
                        g += w * ( scratch.g(x, yu) + scratch.g(x, yd) );
                        b += w * ( scratch.b(x, yu) + scratch.b(x, yd) );
                    }

                    dst.r(x, y) = r / norm_full;
                    dst.g(x, y) = g / norm_full;
                    dst.b(x, y) = b / norm_full;
                }
            }

            // Bottom edge: y in [max(radius, H-radius) .. H-1]
            unsigned y_start = (H > (unsigned)radius) ? (H - (unsigned)radius) : 0u;
            for (unsigned y = y_start; y < H; ++y)
            {
                double r = gw[0] * scratch.r(x, y);
                double g = gw[0] * scratch.g(x, y);
                double b = gw[0] * scratch.b(x, y);
                double n = gw[0];

                for (int k = 1; k <= radius; ++k)
                {
                    double w = gw[k];

                    // up
                    if (y >= (unsigned)k) {
                        unsigned yu = y - (unsigned)k;
                        r += w * scratch.r(x, yu);
                        g += w * scratch.g(x, yu);
                        b += w * scratch.b(x, yu);
                        n += w;
                    }

                    // down
                    unsigned yd = y + (unsigned)k;
                    if (yd < H) {
                        r += w * scratch.r(x, yd);
                        g += w * scratch.g(x, yd);
                        b += w * scratch.b(x, yd);
                        n += w;
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
