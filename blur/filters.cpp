#include "filters.hpp"
#include <cmath>
#include <cassert>
#include <vector>

namespace Filter {
namespace Gauss {

// A simple rule-of-thumb mapping from radius -> sigma.
static inline float default_sigma(int radius) {
    return (radius <= 0) ? 0.0f : (radius * 0.5f + 0.5f);
}

const std::vector<float>& get_weights(int radius) {
    assert(radius >= 0);

    static thread_local int    cached_radius = -1;
    static thread_local float  cached_sigma  = -1.0f;
    static thread_local std::vector<float> kernel;

    if (radius == 0) {
        // Degenerate case: identity kernel
        if (cached_radius != 0 || kernel.size() != 1) {
            kernel.assign(1, 1.0f);
            cached_radius = 0;
            cached_sigma  = 0.0f;
        }
        return kernel;
    }

    const float sigma = default_sigma(radius);

    // Only (re)build when parameters change (or the buffer doesn't match size).
    const int needed_size = 2 * radius + 1;
    const bool need_rebuild = (cached_radius != radius) ||
                              (cached_sigma  != sigma)  ||
                              (static_cast<int>(kernel.size()) != needed_size);

    if (need_rebuild) {
        kernel.resize(needed_size);

        const float twoSigma2 = 2.0f * sigma * sigma;
        float sum = 0.0f;

        // Symmetric weights centered at 0
        for (int i = -radius; i <= radius; ++i) {
            // Using float math is faster and fully adequate for 8-bit images
            const float w = std::exp(-(static_cast<float>(i * i)) / twoSigma2);
            kernel[i + radius] = w;
            sum += w;
        }

        // Normalize so sum(weights) == 1.0
        const float inv = 1.0f / sum;
        for (float& w : kernel) w *= inv;

        cached_radius = radius;
        cached_sigma  = sigma;
    }

    return kernel;
}

} // namespace Gauss
} // namespace Filter
