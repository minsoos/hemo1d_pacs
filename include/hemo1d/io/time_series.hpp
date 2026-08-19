#pragma once

#include <vector>
#include <filesystem>

#include "hemo1d/core/types.hpp"

namespace hemo1d::io{

class TimeSeries{
private:
    std::vector<Real> times_;
    std::vector<Real> values_;

public:
    TimeSeries(std::vector<Real> times, std::vector<Real> values);
    static TimeSeries fromCsv(const std::filesystem::path& path);

    // value: Gives the value of the time series at time t
    // It can take values between times().front() and times().back()
    // It uses interpolation for non sampled values.
    Real value(Real t) const;

    const std::vector<Real>& times() const noexcept { return times_;}
    const std::vector<Real>& values() const noexcept { return values_;}

};
} // namespace hemo1d::io