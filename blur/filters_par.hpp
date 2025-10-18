#pragma once
#include "matrix.hpp"

namespace Filter {
    // Parallel Gaussian blur with pthreads
    Matrix blur_par(Matrix m, int radius, int num_threads);
}
