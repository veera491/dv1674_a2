#include "filters.hpp"
#include <cmath>
#include <cassert>
#include <vector>

namespace Filter {
    namespace Gauss {

        // Compute a default sigma based on the blur radius.
        // You can adjust this formula if your course defines another relationship.
        static inline float default_sigma(int radius) {
            return (radius <= 0) ? 0.0f : (radius * 0.5f + 0.5f);
        }

        const std::vector<float>& get_weights(int radius) {
            assert(radius >= 0);

            // Thread-local cache to avoid recomputation and keep things thread-safe.
            static thread_local int cached_radius = -1;
            static thread_local float cached_sigma = -1.0f;
            static thread_local std::vector<float> kernel;

            if (radius == 0) {
                // Degenerate kernel for radius 0
                if (cached_radius != 0 || kernel.size() != 1) {
                    kernel.assign(1, 1.0f);
                    cached_radius = 0;
                    cached_sigma = 0.0f;
                }
                return kernel;
            }

            const float sigma = default_sigma(radius);

            // Only rebuild when radius or sigma changes
            const int size = 2 * radius + 1;
            if (cached_radius != radius || cached_sigma != sigma ||
                static_cast<int>(kernel.size()) != size) {

                kernel.resize(size);
                const float twoSigma2 = 2.0f * sigma * sigma;
                float sum = 0.0f;

                // Build Gaussian weights symmetrically
                for (int i = -radius; i <= radius; ++i) {
                    float w = std::exp(-(i * i) / twoSigma2);
                    kernel[i + radius] = w;
                    sum += w;
                }

                // Normalize so the weights sum to 1
                const float inv_sum = 1.0f / sum;
                for (float& w : kernel) w *= inv_sum;

                cached_radius = radius;
                cached_sigma = sigma;
                }

            return kernel;
        }

    } // namespace Gauss
} // namespace Filter
