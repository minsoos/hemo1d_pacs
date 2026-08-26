#pragma once
 
#include <filesystem>
 
#include "hemo1d/core/network.hpp"

namespace hemo1d::io {

// Parses a JSON network file into a Network object
Network loadNetwork(const std::filesystem::path& jsonfile);

} // namespace hemo1d::io