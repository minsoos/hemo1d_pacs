#include "hemo1d/io/time_series.hpp"
 
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace hemo1d::io {

namespace{

std::string trim(const std::string& s){
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last-first+1);
}

} // namespace

TimeSeries::TimeSeries(std::vector<Real> times, std::vector<Real> values):
        times_(std::move(times)), values_(std::move(values)){
            if (times_.size() != values_.size()){
                throw std::runtime_error("TimeSeries: times and values must have the same size");
            }
            if (times_.empty()) {
                throw std::runtime_error("TimeSeries: At least one sample must been given");
            }
            if (std::adjacent_find(times_.begin(), times_.end(),
            [](Real a, Real b) { return a>=b; }) !=times_.end()){
                throw std::runtime_error("TimeSeries: times must be strictly increasing");
            }
        }

TimeSeries TimeSeries::fromCsv(const std::filesystem::path& path){
    std::ifstream file(path);
    if (!file){
        throw std::runtime_error("TimeSeries: fromCsv: Cannot open file" + path.string());
    }

    std::vector<Real> times;
    std::vector<Real> values;

    std::string line;

    std::size_t lineNo = 0;
    bool headerAllowed = true;
    while (std::getline(file, line)) {
        ++lineNo;
        const std::string trimmed = trim(line);

        // Empty line or comment
        if (trimmed.empty() || trimmed[0] == '#') continue;

        std::stringstream ss(trimmed);
        std::string timeToken, valueToken;

        // Malformed line
        if (!std::getline(ss, timeToken, ',') || !std::getline(ss, valueToken, ',')) {
            throw std::invalid_argument("TimeSeries: fromCsv: Expected two argument splitted by a comma");
        };

        timeToken = trim(timeToken);
        valueToken = trim(valueToken);

        try {
            const Real t = std::stod(timeToken);
            const Real v = std::stod(valueToken);
            times.push_back(t);
            values.push_back(v);
            headerAllowed = false;
        } catch (const std::exception&){
            // Header: Non numeric row.
            if (headerAllowed) { headerAllowed = false; continue; }
            // Malformed data row.
            throw std::runtime_error("TimeSeries: fromCsv: Malformed row in path: " + path.string() + 
                ": line " + std::to_string(lineNo) + ": with content: \"" + trimmed + "\"");
        }
    }

    if (times.empty()) {
        throw std::runtime_error("TimeSeries: fromCsv: No data rows found in " + path.string());
    }
    return TimeSeries(std::move(times), std::move(values));
}


Real TimeSeries::value(Real t) const {
    const Real tFront = times_.front();
    const Real tBack = times_.back();

    if (t <= tFront) return values_.front();
    if (t >= tBack)  return values_.back();
    const auto it = std::upper_bound(times_.begin(), times_.end(), t);
    const Index i = static_cast<Index>(it - times_.begin());
    const Real t0 = times_[i-1];
    const Real t1 = times_[i];
    const Real v0 = values_[i-1];
    const Real v1 = values_[i];
    const Real frac = (t-t0) / (t1-t0);
    return v0 + frac * (v1-v0);
}

} // namespace hemo1d::io