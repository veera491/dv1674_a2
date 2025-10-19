/*
Parallel Pearson (pthreads) — packed upper triangle + fused 4x-unrolled kernel (index-based).
CLI: ./pearson_par [dataset] [outfile] [num_threads]
*/

#include <pthread.h>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include "vector.hpp"
#include "dataset.hpp"

struct ThreadCtx {
    unsigned tid;
    unsigned T;
    unsigned N;
    const std::vector<Vector>* datasets;
    const std::vector<double>* means;
    const std::vector<double>* invstds;
    const std::vector<unsigned>* lens;
    const std::vector<size_t>* row_start;
    std::vector<double>* R_ut; // packed upper triangle
};

// 4x-unrolled fused kernel using operator[] (no raw pointers needed)
static inline double pearson_pair_kernel_unrolled_idx(
    const Vector& a, const Vector& b,
    unsigned n,
    const double mean_a, const double invstd_a,
    const double mean_b, const double invstd_b,
    const double invN)
{
    double acc = 0.0;

    unsigned k = 0;
    const unsigned n4 = n & ~3u;
    for (; k < n4; k += 4) {
        const double a0 = (a[k+0] - mean_a) * invstd_a; const double b0 = (b[k+0] - mean_b) * invstd_b;
        const double a1 = (a[k+1] - mean_a) * invstd_a; const double b1 = (b[k+1] - mean_b) * invstd_b;
        const double a2 = (a[k+2] - mean_a) * invstd_a; const double b2 = (b[k+2] - mean_b) * invstd_b;
        const double a3 = (a[k+3] - mean_a) * invstd_a; const double b3 = (b[k+3] - mean_b) * invstd_b;
        acc += a0*b0 + a1*b1 + a2*b2 + a3*b3;
    }
    for (; k < n; ++k) {
        const double da = (a[k] - mean_a) * invstd_a;
        const double db = (b[k] - mean_b) * invstd_b;
        acc += da * db;
    }
    return acc * invN;
}

static void* worker(void* arg) {
    ThreadCtx* C = reinterpret_cast<ThreadCtx*>(arg);
    const unsigned N   = C->N;
    const unsigned T   = C->T;
    const unsigned tid = C->tid;

    // Static block partition over rows i of the upper triangle
    const unsigned rows_per = (N + T - 1) / T;
    const unsigned i_begin  = tid * rows_per;
    const unsigned i_end    = std::min(N, i_begin + rows_per);

    for (unsigned i = i_begin; i < i_end; ++i) {
        const unsigned ni   = (*(C->lens))[i];
        const double   mi   = (*(C->means))[i];
        const double   invsi= (*(C->invstds))[i];
        const bool     ok_i = (ni > 0 && invsi != 0.0);

        const Vector& Ai = (*(C->datasets))[i];

        for (unsigned j = i + 1; j < N; ++j) {
            double rij = 0.0;

            const unsigned nj   = (*(C->lens))[j];
            const double   mj   = (*(C->means))[j];
            const double   invsj= (*(C->invstds))[j];
            const bool     ok_j = (nj > 0 && invsj != 0.0);

            if (ok_i && ok_j && ni == nj) {
                const double invN = 1.0 / static_cast<double>(ni);
                const Vector& Bj = (*(C->datasets))[j];
                rij = pearson_pair_kernel_unrolled_idx(Ai, Bj, ni, mi, invsi, mj, invsj, invN);
            } else {
                rij = 0.0; // mismatched length or zero variance → define 0
            }

            const size_t p = (*(C->row_start))[i] + static_cast<size_t>(j - i - 1);
            (*(C->R_ut))[p] = rij;
        }
    }
    return nullptr;
}

// Mean (population)
static inline double compute_mean_vec(const Vector& v) {
    const unsigned n = v.get_size();
    if (n == 0) return 0.0;
    double s = 0.0;
    for (unsigned k = 0; k < n; ++k) s += v[k];
    return s / static_cast<double>(n);
}

// Inverse stddev from mean (population)
static inline double compute_invstd_from_mean_vec(const Vector& v, double mean) {
    const unsigned n = v.get_size();
    if (n == 0) return 0.0;
    double sumsq = 0.0;
    for (unsigned k = 0; k < n; ++k) {
        const double d = v[k] - mean;
        sumsq += d * d;
    }
    const double var = sumsq / static_cast<double>(n);
    if (var <= 0.0) return 0.0;
    return 1.0 / std::sqrt(var);
}

int main(int argc, char const* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " [dataset] [outfile] [num_threads]\n";
        std::exit(1);
    }
    const char* infile   = argv[1];
    const char* outfile  = argv[2];
    const int   T_arg    = std::atoi(argv[3]);
    if (T_arg <= 0) {
        std::cerr << "num_threads must be positive\n";
        std::exit(1);
    }
    const unsigned T = static_cast<unsigned>(T_arg);

    // Read datasets
    std::vector<Vector> datasets = Dataset::read(infile);
    const unsigned N = static_cast<unsigned>(datasets.size());

    // Precompute stats
    std::vector<double>   means(N, 0.0), invstds(N, 0.0);
    std::vector<unsigned> lens(N, 0);
    for (unsigned i = 0; i < N; ++i) {
        const unsigned n = datasets[i].get_size();
        lens[i] = n;
        const double m = compute_mean_vec(datasets[i]);
        means[i]   = m;
        invstds[i] = compute_invstd_from_mean_vec(datasets[i], m);
    }

    // Packed upper-triangle storage
    const size_t P = static_cast<size_t>(N) * static_cast<size_t>(N - 1) / 2;
    std::vector<double> R_ut(P, 0.0);

    // Precompute row starts for O(1) packed indexing: p(i,j) = row_start[i] + (j - i - 1)
    std::vector<size_t> row_start(N, 0);
    for (unsigned i = 0; i < N; ++i) {
        row_start[i] = static_cast<size_t>(i) * static_cast<size_t>(N)
                     - static_cast<size_t>(i) * static_cast<size_t>(i + 1) / 2;
    }

    // Launch threads
    std::vector<pthread_t> threads(T);
    std::vector<ThreadCtx> ctx(T);
    for (unsigned t = 0; t < T; ++t) {
        ctx[t] = ThreadCtx{t, T, N, &datasets, &means, &invstds, &lens, &row_start, &R_ut};
        const int rc = pthread_create(&threads[t], nullptr, worker, &ctx[t]);
        if (rc != 0) {
            std::cerr << "pthread_create failed: " << rc << "\n";
            std::exit(1);
        }
    }
    for (unsigned t = 0; t < T; ++t) {
        pthread_join(threads[t], nullptr);
    }

    // Expand to full N×N flattened matrix and set diagonals
    std::vector<double> R_full(static_cast<size_t>(N) * static_cast<size_t>(N), 0.0);
    for (unsigned i = 0; i < N; ++i) {
        R_full[static_cast<size_t>(i) * N + i] = 1.0;
        for (unsigned j = i + 1; j < N; ++j) {
            const size_t p = row_start[i] + static_cast<size_t>(j - i - 1);
            const double v = R_ut[p];
            R_full[static_cast<size_t>(i) * N + j] = v;
            R_full[static_cast<size_t>(j) * N + i] = v;
        }
    }

    // Write result
    Dataset::write(R_full, outfile);
    return 0;
}
