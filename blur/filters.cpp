#include "filters.hpp"
#include "matrix.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <vector>

namespace Filter {
namespace Gauss {

// --------- Kernel cache (from Optimization #1) ----------
static inline float default_sigma(int radius) {
    return (radius <= 0) ? 0.0f : (radius * 0.5f + 0.5f);
}

const std::vector<float>& get_weights(int radius) {
    assert(radius >= 0);

    static thread_local int    cached_radius = -1;
    static thread_local float  cached_sigma  = -1.0f;
    static thread_local std::vector<float> kernel;

    if (radius == 0) {
        if (cached_radius != 0 || kernel.size() != 1) {
            kernel.assign(1, 1.0f);
            cached_radius = 0;
            cached_sigma  = 0.0f;
        }
        return kernel;
    }

    const float sigma = default_sigma(radius);
    const int needed_size = 2 * radius + 1;

    if (cached_radius != radius || cached_sigma != sigma ||
        static_cast<int>(kernel.size()) != needed_size) {

        kernel.resize(needed_size);
        const float twoSigma2 = 2.0f * sigma * sigma;
        float sum = 0.0f;

        for (int i = -radius; i <= radius; ++i) {
            const float w = std::exp(-(static_cast<float>(i * i)) / twoSigma2);
            kernel[i + radius] = w;
            sum += w;
        }
        const float inv = 1.0f / sum;
        for (float& w : kernel) w *= inv;

        cached_radius = radius;
        cached_sigma  = sigma;
    }

    return kernel;
}

} // namespace Gauss

// ----------------- Separable blur implementation -----------------

// Clamp helper to keep border handling branch-light
static inline int clampi(int v, int lo, int hi) {
    return (v < lo) ? lo : (v > hi ? hi : v);
}

// Convert float to uint8 [0,255], with rounding
static inline uint8_t f2u8(float x) {
    x = std::min(std::max(x, 0.0f), 255.0f);
    return static_cast<uint8_t>(x + 0.5f);
}

// Horizontal pass (into a temp Matrix of floats) -> Vertical pass (into output Matrix of u8).
// To minimize accessor overhead *without* refactoring Matrix, we store the intermediate pass in
// three float buffers (planar), then write back to Matrix once per pixel.

struct PlanarFloat {
    std::vector<float> r, g, b;
    void resize(size_t n) { r.resize(n); g.resize(n); b.resize(n); }
};

static void horizontal_pass(const Matrix& src, PlanarFloat& tmp, int radius, const std::vector<float>& k) {
    const int W = static_cast<int>(src.get_x_size());
    const int H = static_cast<int>(src.get_y_size());
    const int K = static_cast<int>(k.size());
    const int R = radius;

    // For each row, convolve horizontally with clamped borders
    for (int y = 0; y < H; ++y) {
        const int rowOff = y * W;
        for (int x = 0; x < W; ++x) {
            float accR = 0.f, accG = 0.f, accB = 0.f;

            // 1D kernel centered at (x)
            // Note: if performance requires, you can split into "interior" vs "border" loops
            // to remove clamps in the interior. This version keeps it simple & correct.
            for (int t = -R; t <= R; ++t) {
                const int xx = clampi(x + t, 0, W - 1);
                const float w = k[t + R];

                // Accessors: Matrix::r/g/b(x,y)
                accR += w * static_cast<float>(src.r(xx, y));
                accG += w * static_cast<float>(src.g(xx, y));
                accB += w * static_cast<float>(src.b(xx, y));
            }

            const int i = rowOff + x;
            tmp.r[i] = accR;
            tmp.g[i] = accG;
            tmp.b[i] = accB;
        }
    }
}

static void vertical_pass(const PlanarFloat& tmp, Matrix& dst, int radius, const std::vector<float>& k) {
    const int W = static_cast<int>(dst.get_x_size());
    const int H = static_cast<int>(dst.get_y_size());
    const int R = radius;

    // For each column, convolve vertically with clamped borders
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            float accR = 0.f, accG = 0.f, accB = 0.f;

            for (int t = -R; t <= R; ++t) {
                const int yy = clampi(y + t, 0, H - 1);
                const float w = k[t + R];

                const int idx = yy * W + x;
                accR += w * tmp.r[idx];
                accG += w * tmp.g[idx];
                accB += w * tmp.b[idx];
            }

            // Write back once per pixel
            dst.r(x, y) = f2u8(accR);
            dst.g(x, y) = f2u8(accG);
            dst.b(x, y) = f2u8(accB);
        }
    }
}

Matrix blur(const Matrix& src, int radius) {
    assert(radius >= 0);

    // Trivial case: radius 0 returns a copy
    if (radius == 0) {
        Matrix out = src;
        return out;
    }

    const int W = static_cast<int>(src.get_x_size());
    const int H = static_cast<int>(src.get_y_size());
    Matrix out(W, H);                   // assumes Matrix(w,h) ctor exists
    PlanarFloat tmp; tmp.resize(static_cast<size_t>(W) * static_cast<size_t>(H));

    // Reuse cached Gaussian kernel
    const std::vector<float>& k = Gauss::get_weights(radius);

    // Pass 1: horizontal into tmp
    horizontal_pass(src, tmp, radius, k);

    // Pass 2: vertical into out
    vertical_pass(tmp, out, radius, k);

    return out;
}

} // namespace Filter
