#pragma once
#include "Integrator.h"


namespace QDynamics
{
	
	/*! A simple subclass of QDynamics::Magi which overloads the Magi::Magnus() function with an exact analytical expressions, removing the need for numerical estimation of the terms. We therefore set Magi::Resolution to 0, as this is no longer needed. 
	 * 
	 * When working with Symi, note that you have to be explicit about member variables because of template based inheritance: Dependent names vs Independent names. It's inbuilt into the specification for reasons I don't quite know. Therefore, you must refer to member variables as ``Symi::x``, rather than just ``x``.
	
	*/
	template<unsigned int order, UpdateType step>
	class Symi : public Magi<order, step>
	{
		public:
		
			//!Constructor, as per Integrator::Integrator(double,double). \param T The total duration of the integration \param deltaT the width of timesteps used 
			Symi(double T,double deltaT) : Magi<order,step>(T,deltaT,0)
			{
				Symi::Name = "Symi" + std::to_string(Symi::Order);
			}
			//!Constructor, as per Integrator::Integrator(double,double,int). \param T The total duration of the integration \param deltaT the width of timesteps used \param skipper The number of epochs which pass in between updates to the progess Integrator::Buffer
			Symi(double T, double deltaT, int skipper) : Magi<order,step>(T,deltaT,skipper,0)
			{
				Symi::Name = "Symi" + std::to_string(Symi::Order);
			}
		
		protected:
				
		//! The primary purpose of this subclass. This heavily mathematical function computes exactly the terms of the Magnus expansion exactly (up to numerical precision), in the case where the moment of inertia matrix is aligned such that J = (I,I,J_z). \param t The duration of the torque-free precession
		Quaternion virtual Magnus(double duration)
		{
			if (Symi::Order == 0 || abs(Symi::L.Vector(2)) < 1e-8)
			{
				return duration * InvMult(Symi::J,Symi::L);
			}
			
			double lParallel = sqrt(Symi::L.Vector(0)*Symi::L.Vector(0) + Symi::L.Vector(1) * Symi::L.Vector(1));
			double varphi = atan2(Symi::L.Vector(1),Symi::L.Vector(0));
			double zeta =  Symi::L.Vector(2) * (Symi::J[1] - Symi::J[3])/(Symi::J[1] * Symi::J[3]);
			
			double s = sin(zeta * duration - varphi);
			double c = cos(zeta * duration - varphi);
			
			//Update momentum
			Symi::L.Vector(0) = lParallel * c;
			Symi::L.Vector(1) = -lParallel * s;
			
			//first order
			double pre =lParallel/(zeta * Symi::J[1]); 
			double m1x = pre * (s + sin(varphi));
			double m1y= pre * (c - cos(varphi));
			double m1z = Symi::L.Vector(2)/Symi::J[3] * duration;
			
			Quaternion M(0,m1x,m1y,m1z);
			
			if (Symi::Order > 1)
			{
				double sHalf = sin(zeta * duration/2);
				double cHalf = cos(zeta * duration/2);
				double fac =  ( zeta*duration*cHalf - 2*sHalf);
				double pre2 = 2 * lParallel *Symi::L.Vector(2)/(zeta*zeta * Symi::J[1]*Symi::J[3]);
				double m2x = pre2 * cos(zeta * duration/2 - varphi) *fac;
				double m2y = - pre2 * sin(zeta * duration/2 - varphi) * fac;
				double m2z = lParallel*lParallel/(zeta*zeta * Symi::J[1]*Symi::J[1]) * (zeta * duration - sin(zeta*duration));
				
				M.Vector(0) += m2x;
				M.Vector(1) += m2y;
				M.Vector(2) += m2z;
				
			}
			return M;
		}
		
		//! Generates a print-and-filename-friendly name depending on the properties of the integrator. 
		virtual void MakeName()
		{
			Symi::Name = "Sym";
			switch (Symi::StepMode)
			{
				case (Brute):
					Symi::Name = "Brute";
					break;
				case (Euler):
					Symi::Name += "_Euler";
					break;
				case (Leapfrog):
					Symi::Name += "_Frog";
					break;
			}
			Symi::Name +=  "_" + std::to_string(Symi::Order);
		}
	};
	
	
}
