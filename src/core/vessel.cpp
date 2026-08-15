#include "hemo1d/core/vessel.hpp"


namespace hemo1d{

Vessel::Vessel(Id id, std::string name, VesselParameters params):
    id_(id), name_(std::move(name)), params_(params){}
} //namespace hemo1d