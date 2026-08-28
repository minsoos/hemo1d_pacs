#include "hemo1d/dg/quadrature.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace hemo1d::dg {

namespace {

constexpr Real kPi = 3.14159265358979323846;
constexpr Real kMaxNewtonIterations = 100;
constexpr Real kNewtonTolerance = 1e-15;

std::vector<Real> legendreValues(unsigned n, Real x) {
    std::vector<Real> P(n + 1);
    P[0] = 1.0;
    if (n >= 1) P[1] = x;
    for (unsigned k = 2; k <= n; ++k) {
        P[k] = ((2.0 * k - 1.0) * x * P[k - 1] - (k - 1.0) * P[k - 2]) / k;
    }

    return P;   
}

} // namespace


QuadratureRule gaussLegendre(unsigned numPoints) {
    if (numPoints < 1) {
        throw std::invalid_argument("gaussLegendre: Number of points must be at least 1.");
    }

    const unsigned n = numPoints;

    std::vector<std::pair<Real, Real>> pointsAndWeights(n);
    for (unsigned i = 0; i < n; ++i) {
        Real x = std::cos(kPi * (i + 0.75) / (n + 0.5));
        Real Pn = 0.0, dPn = 0.0;
        for (int it = 0; it < kMaxNewtonIterations; ++it) {
            const std::vector<Real> P = legendreValues(n, x);
            Pn = P[n];
            const Real Pnm1 = P[n - 1];
            dPn = n * (x * Pn - Pnm1) / (x * x - 1.0);
            const Real dx = -Pn / dPn;
            x += dx;
            if (std::abs(dx) < kNewtonTolerance) break;
        }
        pointsAndWeights[i] = {x, 2.0 / ((1.0 - x * x) * dPn * dPn)};
    }

    std::sort(pointsAndWeights.begin(), pointsAndWeights.end());

    QuadratureRule rule;
    rule.points.resize(n);
    rule.weights.resize(n);
    for (unsigned i = 0; i < n; ++i) {
        rule.points[i] = pointsAndWeights[i].first;
        rule.weights[i] = pointsAndWeights[i].second;
    }

    return rule;
}


QuadratureRule gaussLobattoLegendre(unsigned numPoints) {
    if (numPoints < 2) {
        throw std::invalid_argument("gaussLobattoLegendre: numPoints must be >= 2");
    }

    const unsigned N = numPoints - 1;
    const unsigned N1 = numPoints;

    std::vector<Real> x(N1);
    for (unsigned j = 0; j < N1; ++j) {
        x[j] = -std::cos(kPi * j / N);
    }

    for (int it = 0; it < kMaxNewtonIterations; ++it) {
        Real maxDelta = 0.0;
        for (unsigned j = 0; j < N1; ++j) {
            const std::vector<Real> P = legendreValues(N, x[j]);
            const Real delta = (x[j] * P[N] - P[N - 1]) / (N1 * P[N]);
            x[j] -= delta;
            maxDelta = std::max(maxDelta, std::abs(delta));
        }
        if (maxDelta < kNewtonTolerance) break;
    }

    x.front() = -1.0;
    x.back() = 1.0;

    QuadratureRule rule;
    rule.points = x;
    rule.weights.resize(N1);
    for (unsigned j = 0; j < N1; ++j) {
        const std::vector<Real> P = legendreValues(N, x[j]);
        rule.weights[j] = 2.0 / (N * N1 * P[N] * P[N]);
    }
    return rule;
}

} // namespace hemo1d::dg