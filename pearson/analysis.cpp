#include <cmath>
#include <algorithm>
#include <vector>
#include "analysis.hpp"
#include "vector.hpp"

// --- small helpers bound to your Vector API ---
static inline unsigned VLEN(const Vector& v) { return v.get_size(); }

// Compute mean (population) of v
static inline double compute_mean(const Vector& v) {
    const unsigned n = VLEN(v);
    if (n == 0) return 0.0;
    double s = 0.0;
    for (unsigned k = 0; k < n; ++k) s += v[k];
    return s / static_cast<double>(n);
}

// From mean, compute inverse stddev (population variance ÷ n)
// If your verifier expects sample variance, switch to n-1 with guard.
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

// ---------------- Pearson for two vectors (kept simple & compatible) ----------------
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

    // 8-way unrolled fused pass (no temporaries)
    double acc = 0.0;
    unsigned k = 0;
    const unsigned n8 = n1 & ~7u;
    for (; k < n8; k += 8) {
        const double a0 = (vec1[k+0] - m1) * invs1; const double b0 = (vec2[k+0] - m2) * invs2;
        const double a1 = (vec1[k+1] - m1) * invs1; const double b1 = (vec2[k+1] - m2) * invs2;
        const double a2 = (vec1[k+2] - m1) * invs1; const double b2 = (vec2[k+2] - m2) * invs2;
        const double a3 = (vec1[k+3] - m1) * invs1; const double b3 = (vec2[k+3] - m2) * invs2;
        const double a4 = (vec1[k+4] - m1) * invs1; const double b4 = (vec2[k+4] - m2) * invs2;
        const double a5 = (vec1[k+5] - m1) * invs1; const double b5 = (vec2[k+5] - m2) * invs2;
        const double a6 = (vec1[k+6] - m1) * invs1; const double b6 = (vec2[k+6] - m2) * invs2;
        const double a7 = (vec1[k+7] - m1) * invs1; const double b7 = (vec2[k+7] - m2) * invs2;
        acc += a0*b0 + a1*b1 + a2*b2 + a3*b3 + a4*b4 + a5*b5 + a6*b6 + a7*b7;
    }
    for (; k < n1; ++k) {
        acc += ((vec1[k] - m1) * invs1) * ((vec2[k] - m2) * invs2);
    }
    return acc * invN;
}

// ---------------- Correlation matrix (Opt3) ----------------
//
// Strategy:
// 1) Precompute means & invstd for each vector.
// 2) Pre-normalize once: Z[i][k] = (X[i][k] - mean_i) * invstd_i.
//    If invstd_i == 0 or length 0, fill Z[i] with zeros (pairs with i -> 0 corr).
// 3) Compute only the upper triangle into packed array R_ut of size P=N*(N-1)/2.
//    Use precomputed row_start[i] = i*N - (i*(i+1))/2; index p(i,j)=row_start[i] + (j - i - 1).
// 4) Expand packed to full N×N flattened row-major vector and set diagonals to 1.
//
std::vector<double> Analysis::correlation_coefficients(std::vector<Vector> datasets)
{
    const unsigned N = static_cast<unsigned>(datasets.size());
    std::vector<double> R_full;
    if (N == 0) return R_full;

    // --- Step 1: stats per vector ---
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

    // --- Step 2: pre-normalize to Z (z-scores) ---
    std::vector<std::vector<double>> Z;
    Z.resize(N);
    for (unsigned i = 0; i < N; ++i) {
        const unsigned n = lens[i];
        Z[i].assign(n, 0.0);
        if (n == 0 || invstds[i] == 0.0) continue; // zero-variance or empty → leave zeros
        const double m  = means[i];
        const double is = invstds[i];
        unsigned k = 0;
        const unsigned n8 = n & ~7u;
        for (; k < n8; k += 8) {
            Z[i][k+0] = (datasets[i][k+0] - m) * is;
            Z[i][k+1] = (datasets[i][k+1] - m) * is;
            Z[i][k+2] = (datasets[i][k+2] - m) * is;
            Z[i][k+3] = (datasets[i][k+3] - m) * is;
            Z[i][k+4] = (datasets[i][k+4] - m) * is;
            Z[i][k+5] = (datasets[i][k+5] - m) * is;
            Z[i][k+6] = (datasets[i][k+6] - m) * is;
            Z[i][k+7] = (datasets[i][k+7] - m) * is;
        }
        for (; k < n; ++k) {
            Z[i][k] = (datasets[i][k] - m) * is;
        }
    }

    // --- Step 3: compute upper triangle into packed array ---
    const size_t P = static_cast<size_t>(N) * static_cast<size_t>(N - 1) / 2;
    std::vector<double> R_ut(P, 0.0);

    // Precompute row starts so p(i,j) is O(1)
    std::vector<size_t> row_start(N, 0);
    for (unsigned i = 0; i < N; ++i) {
        row_start[i] = static_cast<size_t>(i) * static_cast<size_t>(N)
                     - static_cast<size_t>(i) * static_cast<size_t>(i + 1) / 2;
    }

    for (unsigned i = 0; i < N; ++i) {
        const unsigned ni = lens[i];
        const bool zi_valid = (ni > 0 && invstds[i] != 0.0);
        const std::vector<double>& Zi = Z[i];

        for (unsigned j = i + 1; j < N; ++j) {
            double rij = 0.0;
            const unsigned nj = lens[j];
            const bool zj_valid = (nj > 0 && invstds[j] != 0.0);

            if (zi_valid && zj_valid && ni == nj) {
                const double invN = 1.0 / static_cast<double>(ni);
                const std::vector<double>& Zj = Z[j];

                double acc = 0.0;
                unsigned k = 0;
                const unsigned n8 = ni & ~7u;
                for (; k < n8; k += 8) {
                    acc += Zi[k+0]*Zj[k+0] + Zi[k+1]*Zj[k+1]
                         + Zi[k+2]*Zj[k+2] + Zi[k+3]*Zj[k+3]
                         + Zi[k+4]*Zj[k+4] + Zi[k+5]*Zj[k+5]
                         + Zi[k+6]*Zj[k+6] + Zi[k+7]*Zj[k+7];
                }
                for (; k < ni; ++k) {
                    acc += Zi[k] * Zj[k];
                }
                rij = acc * invN;
            } else {
                rij = 0.0; // mismatched length or zero-variance → define as 0
            }

            const size_t p = row_start[i] + static_cast<size_t>(j - i - 1);
            R_ut[p] = rij;
        }
    }

    // --- Step 4: expand to full N×N flattened matrix and set diagonals ---
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
