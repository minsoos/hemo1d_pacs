#include "hemo1d/couplings/register.hpp"

#include "hemo1d/couplings/windkessel_coupling.hpp"
#include "hemo1d/physics/terminal_coupling.hpp"

namespace hemo1d::couplings {

void registerBuiltinCouplings() {
    static bool done = false;
    if (done) return;
    done = true;

    physics::TerminalCouplingRegistry::instance().add(
        "windkessel", [](const std::string& paramsJson) {
            return makeWindkesselCoupling(parseWindkesselParams(paramsJson));
        }
    );
}

} // namespace hemo1d::couplings