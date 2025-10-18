#pragma once
#include <vector>

// Public surface remains minimal & stable.
// - Filter::Gauss::get_weights(radius) is the same symbol your code already calls.
// - Filter::blur(...) is now implemented as a *separable* (horizontal + vertical) pass.

namespace Filter {

    // Forward-declare Matrix to avoid dragging headers here; include in .cpp.
    class Matrix;

    /// Apply a Gaussian blur with the given radius (separable: horizontal then vertical).
    /// Returns a new Matrix with the same dimensions as the input.
    Matrix blur(const Matrix& src, int radius);

    namespace Gauss {
        /// Return a normalized 1D Gaussian kernel for a given radius.
        /// This function caches the kernel per-thread and only recomputes when the radius changes.
        const std::vector<float>& get_weights(int radius);
    }

} // namespace Filter
