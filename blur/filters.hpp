#pragma once
#include <vector>

//
// Drop-in header to provide a cached Gaussian kernel for the blur.
// This matches the symbol seen in your profiles: Filter::Gauss::get_weights(...)
//
namespace Filter {
    namespace Gauss {

        /// Return a normalized 1D Gaussian kernel for a given radius.
        /// IMPORTANT:
        /// - This function now caches the kernel per-thread and per-radius.
        /// - If the radius hasn't changed since the last call in the same thread,
        ///   it returns the previously computed kernel without recomputation.
        /// - The returned reference remains valid until the next call in the same thread
        ///   with a different radius.
        ///
        /// @param radius  Non-negative blur radius (kernel width = 2*radius + 1).
        /// @return        Reference to normalized weights (sum == 1.0f).
        const std::vector<float>& get_weights(int radius);

    } // namespace Gauss
} // namespace Filter
