#include <cmath>
#include <algorithm>
#include <vector>
#include "analysis.hpp"
#include "vector.hpp"

// ---- Vector helpers bound to your API ----
static inline unsigned VLEN(const Vector& v) { return v.get_size(); }

// Mean (population)
static inline double compute_mean(const Vector& v) {
    const unsigned n = VLEN(v);
    if (n == 0) return 0.0;
    double s = 0.0;
    for (unsigned k = 0; k < n; ++k) s += v[k];
    return s / static_cast<double>(n);
}

// Inverse stddev from mean (population variance ÷ n)
// If your verifier expects sample variance, switch to (n-1) with guard.
static inline double compute_invstd_from_mean(const Vector& v, double mean) {
    const unsigned n = VLEN(v);
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

// ---- 4x-unrolled fused kernel using operator[] (no raw pointers) ----
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

// ---------------- Pearson for two vectors (API-compat) ----------------
double Analysis::pearson(Vector vec1, Vector vec2)
{
    const unsigned n1 = VLEN(vec1);
    const unsigned n2 = VLEN(vec2);
    if (n1 == 0 || n2 == 0 || n1 != n2) return 0.0;

    const double m1    = compute_mean(vec1);
    const double invs1 = compute_invstd_from_mean(vec1, m1);
    const double m2    = compute_mean(vec2);
    const double invs2 = compute_invstd_from_mean(vec2, m2);

    if (invs1 == 0.0 || invs2 == 0.0) return 0.0;

    const double invN = 1.0 / static_cast<double>(n1);
    return pearson_pair_kernel_unrolled_idx(vec1, vec2, n1, m1, invs1, m2, invs2, invN);
}

// ---------------- Correlation matrix (Hybrid Opt3) ----------------
//
// - Precompute means/invstd once.
// - Compute only upper triangle into packed array (size P = N*(N-1)/2) using fused kernel.
// - Expand to full N×N flattened matrix and set diagonals.
//
std::vector<double> Analysis::correlation_coefficients(std::vector<Vector> datasets)
{
    const unsigned N = static_cast<unsigned>(datasets.size());
    std::vector<double> R_full;
    if (N == 0) return R_full;

    // Stats per vector
    std::vector<double> means(N, 0.0);
    std::vector<double> invstds(N, 0.0);
    std::vector<unsigned> lens(N, 0);

    for (unsigned i = 0; i < N; ++i) {
        const unsigned n = VLEN(datasets[i]);
        lens[i] = n;
        const double m = compute_mean(datasets[i]);
        means[i]   = m;
        invstds[i] = compute_invstd_from_mean(datasets[i], m);
    }

    // Packed upper triangle
    const size_t P = static_cast<size_t>(N) * static_cast<size_t>(N - 1) / 2;
    std::vector<double> R_ut(P, 0.0);

    // Precompute row starts for O(1) packed indexing: p(i,j) = row_start[i] + (j - i - 1)
    std::vector<size_t> row_start(N, 0);
    for (unsigned i = 0; i < N; ++i) {
        row_start[i] = static_cast<size_t>(i) * static_cast<size_t>(N)
                     - static_cast<size_t>(i) * static_cast<size_t>(i + 1) / 2;
    }

    // Compute upper triangle
    for (unsigned i = 0; i < N; ++i) {
        const unsigned ni = lens[i];
        const double mi = means[i];
        const double invsi = invstds[i];
        const bool ok_i = (ni > 0 && invsi != 0.0);

        for (unsigned j = i + 1; j < N; ++j) {
            double rij = 0.0;

            const unsigned nj = lens[j];
            const double mj = means[j];
            const double invsj = invstds[j];
            const bool ok_j = (nj > 0 && invsj != 0.0);

            if (ok_i && ok_j && ni == nj) {
                const double invN = 1.0 / static_cast<double>(ni);
                rij = pearson_pair_kernel_unrolled_idx(
                    datasets[i], datasets[j], ni,
                    mi, invsi, mj, invsj, invN
                );
            } else {
                rij = 0.0; // mismatched length or zero variance → define 0
            }

            const size_t p = row_start[i] + static_cast<size_t>(j - i - 1);
            R_ut[p] = rij;
        }
    }

    // Expand to full N×N flattened matrix and set diagonals
    R_full.assign(static_cast<size_t>(N) * static_cast<size_t>(N), 0.0);

    for (unsigned i = 0; i < N; ++i) {
        R_full[static_cast<size_t>(i) * N + i] = 1.0;
        for (unsigned j = i + 1; j < N; ++j) {
            const size_t p = row_start[i] + static_cast<size_t>(j - i - 1);
            const double v = R_ut[p];
            R_full[static_cast<size_t>(i) * N + j] = v;
            R_full[static_cast<size_t>(j) * N + i] = v;
        }
    }

    return R_full;
}
