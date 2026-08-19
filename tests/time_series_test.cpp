#include <filesystem>
 
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "hemo1d/core/types.hpp"
#include "hemo1d/io/time_series.hpp"

using namespace hemo1d;
using hemo1d::io::TimeSeries;
using Catch::Approx;

TEST_CASE("TimeSeries: Creates successfully", "[time_series]") {
    std::vector<Real> t({0.0, 1.0, 2.0});
    std::vector<Real> v({0.0, 10.0, 10.0});
    TimeSeries ts(t, v);

    CHECK(ts.times() == t);
    CHECK(ts.values() == v);
    CHECK(ts.value(1.0) == Approx(10.0));
}

TEST_CASE("TimeSeries: Interpolates linearly between samples", "[time_series]") {
    TimeSeries ts({0.0, 1.0, 2.0}, {0.0, 10.0, 10.0});
 
    CHECK(ts.value(0.0) == Approx(0.0));
    CHECK(ts.value(0.5) == Approx(5.0));
    CHECK(ts.value(1.0) == Approx(10.0));
    CHECK(ts.value(1.5) == Approx(10.0));
    CHECK(ts.value(2.0) == Approx(10.0));
}

TEST_CASE("TimeSeries: Non non-decreasing time values", "[time_series]") {
    std::vector<Real> t({0.0, 1, 0.5});
    std::vector<Real> v({0.0, 2.0, 3.0});
    CHECK_THROWS(TimeSeries(t, v));
}

TEST_CASE("TimeSeries: Non strictly increasing time values", "[time_series]") {
    std::vector<Real> t({0.0, 0.5, 0.5});
    std::vector<Real> v({0.0, 2.0, 3.0});
    CHECK_THROWS(TimeSeries(t, v));
}

TEST_CASE("TimeSeries: Different lengths of times and values", "[time_series]") {
    std::vector<Real> t({0.0, 0.5, 0.5});
    std::vector<Real> v({0.0, 2.0});
    CHECK_THROWS(TimeSeries(t, v));
}

TEST_CASE("TimeSeries: 0 length vectors", "[time_series]") {
    std::vector<Real> t({});
    std::vector<Real> v({});
    CHECK_THROWS(TimeSeries(t, v));
}


TEST_CASE("TimeSeries::fromCsv search at a missing file", "[time_series]") {
    CHECK_THROWS(TimeSeries::fromCsv("does/not/exist.csv"));
}

TEST_CASE("TimeSeries evaluate outside the range", "[time_series]") {
    TimeSeries ts({0.0, 1.0, 2.0}, {0.0, 10.0, 10.0});

    CHECK(ts.value(0.0)  == Approx(0.0));     // exactly at front
    CHECK(ts.value(2.0)  == Approx(10.0));    // exactly at back
    CHECK(ts.value(2.0 + 1e-15) == Approx(10.0));  // within tolerance -> clamps
    CHECK_THROWS(ts.value(2.1));                    // clearly outside -> throws
    CHECK_THROWS(ts.value(-0.5));
}
