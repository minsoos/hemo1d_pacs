#include "hemo1d/physics/exterior_boundary.hpp"

#include <stdexcept>
#include <vector>

#include "hemo1d/core/linear_algebra.hpp"
#include "hemo1d/physics/compatibility.hpp"

namespace hemo1d::physics
{
std::pair<Real, Real> solveExteriorBoundary(VesselEnd end, const BoundaryConditionSpec &bc,
		Real prescribedValue, Real A, Real Q, Real dAdz, Real dQdz, const VesselParameters &params,
		Real rho, const TubeLaw &tubeLaw, Real dt){
			const LeftEigenvectors eig = computeLeftEigenvectors(A, Q, params.A0, params.beta, params.alpha,
																	rho, tubeLaw);
			const std::pair<Real, Real> lIn = (end == VesselEnd::Proximal) ? eig.lPlus : eig.lMinus;
			const std::pair<Real, Real> lOut = (end == VesselEnd::Proximal) ? eig.lMinus : eig.lPlus;

			DenseMatrix M(2, 2, 0.0);
			std::vector<Real> rhs(2, 0.0);
			if (bc.type == BoundaryConditionType::NonReflecting){
				const auto [nrA, nrQ] = compatibilityPrediction(A, Q, 0.0, 0.0, params.A0,
									params.beta, params.alpha, rho, params.frictionKr, tubeLaw, dt);

				M(0, 0) = lIn.first;
				M(0, 1) = lIn.second;
				rhs[0] = lIn.first * nrA + lIn.second * nrQ;
			} else{
				switch (bc.quantity){
					case PrescribedQuantity::FlowRate:
						M(0,1) = 1.0;
						rhs[0] = prescribedValue;
						break;		
					case PrescribedQuantity::Area:
						M(0,0) = 1.0;
						rhs[0] = prescribedValue;
						break;	
					case PrescribedQuantity::Velocity:
						M(0,0) = -prescribedValue;
						M(0,1) = 1.0;
						rhs[0] = 0.0;
						break;
					case PrescribedQuantity::Pressure:
						M(0,0) = 1.0;
						rhs[0] = tubeLaw.areaFromPressure(prescribedValue, params.A0, params.beta);
						break;	
				}
			}
	
	const auto [ccA, ccQ] = compatibilityPrediction(A, Q, dAdz, dQdz, params.A0,
						params.beta, params.alpha, rho, params.frictionKr, tubeLaw, dt);
	
	M(1,0) = lOut.first;
	M(1,1) = lOut.second;
	rhs[1] = lOut.first * ccA + lOut.second * ccQ;

	const std::vector<Real> solution = solveLinearSystem(M, rhs);
	return {solution[0], solution[1]};
	}

} // namespace hemo1d::physics