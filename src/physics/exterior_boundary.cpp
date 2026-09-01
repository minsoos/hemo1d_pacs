#include "hemo1d/physics/exterior_boundary.hpp"

#include <utility>
#include <Eigen/Dense>

#include "hemo1d/physics/blood_flow_model.hpp"

namespace hemo1d::physics{

SectionState closeTerminal(
	VesselEnd end, const TerminalRow& physical, SectionState trace,
    SectionGradient grad, const VesselParameters& p,
    const BloodFlowModel& model, Real dt
) {
	const LeftEigenvectors eig = model.leftEigenvectors(trace, p);
	const std::pair<Real, Real> lOut = (end == VesselEnd::Proximal) ? eig.lMinus : eig.lPlus;

	Eigen::Matrix2d M = Eigen::Matrix2d::Zero();
	Eigen::Vector2d rhs = Eigen::Vector2d::Zero();

	// Row 0: the caller-supplied physical / coupling condition
	M(0, 0) = physical.cA;
	M(0, 1) = physical.cQ;
	rhs(0) = physical.rhs;

	// Row 1: the compatibility condition on the outgoing characteristic.
	const SectionState cc = model.compatibilityPrediction(trace, grad, p, dt);
	M(1, 0) = lOut.first;
	M(1, 1) = lOut.second;
	rhs(1) = lOut.first * cc.A + lOut.second * cc.Q;

	const Eigen::Vector2d solution = M.partialPivLu().solve(rhs);
	return {solution(0), solution(1)};
}

SectionState solveExteriorBoundary(
	VesselEnd end, const BoundaryConditionSpec& bc, Real prescribedValue,
    SectionState trace, SectionGradient grad, const VesselParameters& p,
    const BloodFlowModel& model, Real dt
){
    TerminalRow row;

    if (bc.type == BoundaryConditionType::NonReflecting) {
        const LeftEigenvectors eig = model.leftEigenvectors(trace, p);
        const std::pair<Real, Real> lIn = (end == VesselEnd::Proximal) ? eig.lPlus : eig.lMinus;
        const SectionState nr = model.compatibilityPrediction(trace, SectionGradient{}, p, dt);
        row.cA = lIn.first;
        row.cQ = lIn.second;
        row.rhs = lIn.first * nr.A + lIn.second * nr.Q;
    } else {
        switch (bc.quantity) {
            case PrescribedQuantity::FlowRate:
                row.cQ = 1.0;
                row.rhs = prescribedValue;
                break;
            case PrescribedQuantity::Area:
                row.cA = 1.0;
                row.rhs = prescribedValue;
                break;
            case PrescribedQuantity::Velocity:
                row.cA = -prescribedValue;
                row.cQ = 1.0;
                row.rhs = 0.0;
                break;
            case PrescribedQuantity::Pressure:
                row.cA = 1.0;
                row.rhs = model.tubeLaw().areaFromPressure(prescribedValue, p);
                break;
        }
    }

    return closeTerminal(end, row, trace, grad, p, model, dt);
}

} // namespace hemo1d::physics