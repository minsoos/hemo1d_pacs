#include "hemo1d/physics/exterior_boundary.hpp"

#include <utility>
#include <Eigen/Dense>

#include "hemo1d/physics/blood_flow_model.hpp"

namespace hemo1d::physics{

SectionState solveExteriorBoundary(
	VesselEnd end, const BoundaryConditionSpec& bc, Real prescribedValue,
    SectionState trace, SectionGradient grad, const VesselParameters& p,
    const BloodFlowModel& model, Real dt
){
	const LeftEigenvectors eig = model.leftEigenvectors(trace, p);
	const std::pair<Real, Real> lIn = (end == VesselEnd::Proximal) ? eig.lPlus : eig.lMinus;
	const std::pair<Real, Real> lOut = (end == VesselEnd::Proximal) ? eig.lMinus : eig.lPlus;

	Eigen::Matrix2d M = Eigen::Matrix2d::Zero();
	Eigen::Vector2d rhs = Eigen::Vector2d::Zero();
	if (bc.type == BoundaryConditionType::NonReflecting){
		const SectionState nr = model.compatibilityPrediction(trace, SectionGradient{}, p, dt);

		M(0, 0) = lIn.first;
		M(0, 1) = lIn.second;
		rhs(0) = lIn.first * nr.A + lIn.second * nr.Q;
	} else{
		switch (bc.quantity){
			case PrescribedQuantity::FlowRate:
				M(0,1) = 1.0;
				rhs(0) = prescribedValue;
				break;		
			case PrescribedQuantity::Area:
				M(0,0) = 1.0;
				rhs(0) = prescribedValue;
				break;	
			case PrescribedQuantity::Velocity:
				M(0,0) = -prescribedValue;
				M(0,1) = 1.0;
				rhs(0) = 0.0;
				break;
			case PrescribedQuantity::Pressure:
				M(0,0) = 1.0;
				rhs(0) = model.tubeLaw().areaFromPressure(prescribedValue, p);
				break;	
		}
	}
	
	const SectionState cc = model.compatibilityPrediction(trace, grad, p, dt);
	
	M(1,0) = lOut.first;
	M(1,1) = lOut.second;
	rhs(1) = lOut.first * cc.A + lOut.second * cc.Q;

	const Eigen::Vector2d solution = M.partialPivLu().solve(rhs);
	return {solution(0), solution(1)};
	}

} // namespace hemo1d::physics