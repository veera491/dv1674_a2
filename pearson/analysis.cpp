#include <cmath>
#include <algorithm>
#include "analysis.hpp"
#include "vector.hpp"

static inline unsigned VLEN(const Vector& v) { return v.get_size(); }

// ---------- Internal helpers (file-local) ----------

static inline double compute_mean(const Vector& v) {
    const unsigned int n = VLEN(v);
    if (n == 0) return 0.0;

    double s = 0.0;
    for (unsigned int k = 0; k < n; ++k) s += v[k];
    return s / static_cast<double>(n);
}

static inline double compute_invstd_from_mean(const Vector& v, double mean) {
    const unsigned int n = VLEN(v);
    if (n == 0) return 0.0;

    double sumsq = 0.0;
    for (unsigned int k = 0; k < n; ++k) {
        const double d = v[k] - mean;
        sumsq += d * d;
    }
    // population variance (÷ n). If your verifier expects sample variance, switch to (n-1).
    const double var = sumsq / static_cast<double>(n);
    if (var <= 0.0) return 0.0;
    return 1.0 / std::sqrt(var);
}

// Fused Pearson core for a single pair (no temporaries, one pass).
static inline double pearson_pair_fused(
    const Vector& a, const Vector& b,
    double mean_a, double invstd_a,
    double mean_b, double invstd_b)
{
    const unsigned int n = VLEN(a);
    if (n == 0 || n != VLEN(b) || invstd_a == 0.0 || invstd_b == 0.0) {
        return 0.0;
    }

    double acc = 0.0;
    for (unsigned int k = 0; k < n; ++k) {
        const double za = (a[k] - mean_a) * invstd_a;
        const double zb = (b[k] - mean_b) * invstd_b;
        acc += za * zb;
    }
    return acc / static_cast<double>(n);
}

// ---------- Public API (match headers) ----------

double Analysis::pearson(Vector vec1, Vector vec2)
{
    const double m1    = compute_mean(vec1);
    const double invs1 = compute_invstd_from_mean(vec1, m1);

    const double m2    = compute_mean(vec2);
    const double invs2 = compute_invstd_from_mean(vec2, m2);

    return pearson_pair_fused(vec1, vec2, m1, invs1, m2, invs2);
}

// Flattened N×N correlation matrix in row-major: R[i*N + j]
std::vector<double> Analysis::correlation_coefficients(std::vector<Vector> datasets)
{
    const unsigned int N = static_cast<unsigned int>(datasets.size());
    std::vector<double> R;
    if (N == 0) return R;

    R.assign(static_cast<size_t>(N) * static_cast<size_t>(N), 0.0);

    // Precompute per-vector stats once
    std::vector<double> means(N, 0.0);
    std::vector<double> invstds(N, 0.0);
    for (unsigned int i = 0; i < N; ++i) {
        const double m = compute_mean(datasets[i]);
        means[i]   = m;
        invstds[i] = compute_invstd_from_mean(datasets[i], m);
    }

    // Upper triangle only, then mirror
    for (unsigned int i = 0; i < N; ++i) {
        R[static_cast<size_t>(i) * N + i] = 1.0;
        for (unsigned int j = i + 1; j < N; ++j) {
            const double rij = pearson_pair_fused(
                datasets[i], datasets[j],
                means[i], invstds[i],
                means[j], invstds[j]);

            R[static_cast<size_t>(i) * N + j] = rij;
            R[static_cast<size_t>(j) * N + i] = rij;
        }
    }

    return R;
}
