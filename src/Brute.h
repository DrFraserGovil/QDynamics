#pragma once
#include "Integrator.h"

namespace QDynamics
{
	
	/*!
	*  The dumbest possible way to integrate. Overrides Integrator::UpdatePosition() with the brute force method: \verbatim embed:rst
	  .. math::
		\mathsf{q}  &\to \mathsf{q} + \delta t \dot{\mathsf{q}} 
		~\\
		\mathsf{p} & \to \mathsf{p} + \delta t \dot{\mathsf{p}} 
		~\\
		\dot{\mathsf{q}} & = \frac{1}{2} \mathsf{q} \otimes J^{-1} \overline{\mathsf{q}} \otimes \mathsf{p}
		~\\
		\dot{\mathsf{p}} & = \frac{1}{2} \mathsf{p} \otimes J^{-1} \overline{\mathsf{p}} \otimes \mathsf{q} + \nabla_\mathsf{q} U(\mathsf{q},t)
	  \endverbatim 
	  There is also a normalisation step on q, to prevent changes in the magnitude.
	*/
	class BruteInt : public Integrator
	{
		public:
			//! Identical to Integrator::Integrator(double,double), but also assigns Integrator::Name to "Brute" \param T The total duration of the integration \param deltaT the width of timesteps used
			BruteInt(double T,double deltaT) : Integrator(T,deltaT)
			{
				Name = "Brute";
			}
			//! Identical to Integrator::Integrator(double,double,int), but also assigns Integrator::Name to "Brute" \param T The total duration of the integration \param deltaT the width of timesteps used \param skipper The number of epochs which pass in between updates to the progess Integrator::Buffer
			BruteInt(double T, double deltaT, int skipper) : Integrator(T,deltaT, skipper)
			{
				Name = "Brute";
			}
		
		private:
		
		//! Executes a brute-force linear update to the position and momentum terms. A dumb, first order Euler integrator with a normalistion step, nothing more. \param t The current time (used for the potential U())
		void virtual UpdatePosition(double t)
		{
			w =  InvMult(J,L);
			w.Scalar() = 0;
			Quaternion qDot = 0.5 * q * w;
			Quaternion pDot = 0.5 * p * w - GradU(t);
			
	
			p = p + TimeStep * pDot;
			q = q + TimeStep * qDot;
			
			q/=q.Norm();
			L = 0.5 * q.Conjugate() * p;
		}
	};
}
