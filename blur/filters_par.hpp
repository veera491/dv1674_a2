/*
Author: David Holmqvist <daae19@student.bth.se>
Parallelized by: Nithin
*/

#ifndef FILTERS_PAR_HPP
#define FILTERS_PAR_HPP

#include "matrix.hpp"

namespace FilterPar
{
    namespace Gauss
    {
        constexpr unsigned max_radius{1000};
        constexpr float max_x{1.33};
        constexpr float pi{3.14159};

        void get_weights(int n, double *weights_out);
    }

    Matrix blur(Matrix m, const int radius, int num_threads);
}

#endif