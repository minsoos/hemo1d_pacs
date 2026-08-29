#pragma once

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include "hemo1d/dg/mesh.hpp"
#include "hemo1d/physics/state.hpp"
#include "hemo1d/physics/tube_law.hpp"

namespace hemo1d::output {

void writeVtk(const std::filesystem::path& path, const dg::Mesh& mesh, 
            const physics::State& state, const physics::TubeLaw& tubeLaw,
            Real rho, Real time);

class VtkSeriesWriter{

private:
    void writeCollectionFile() const;
    std::filesystem::path directory_;
    std::string baseName_;
    int nextIndex_ = 0;
    std::vector<std::pair<int, Real>> writtenSnapshots_;

public:
    VtkSeriesWriter(std::filesystem::path directory, std::string basName);

    void write(const dg::Mesh& mesh, const physics::State& state, 
        const physics::TubeLaw& tubeLaw, Real rho, Real time);
};

} // namespace hemo1d::output