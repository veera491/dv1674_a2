#pragma once
#include <vector>

namespace Filter {
    namespace Gauss {

        /// Returns a normalized 1D Gaussian kernel for a given blur radius.
        /// The kernel is cached per thread, so it is computed only once per radius.
        /// Thread-safe and efficient for both sequential and parallel execution.
        ///
        /// @param radius  Non-negative blur radius (kernel width = 2*radius + 1)
        /// @return        Reference to the cached normalized weights (sum == 1.0f)
        const std::vector<float>& get_weights(int radius);

    } // namespace Gauss
} // namespace Filter
