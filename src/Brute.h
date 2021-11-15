#pragma once
#include "Integrator.h"

namespace QDynamics
{
	
	/*
	  The dumbest possible way to integrate 
	*/
	class Brute : public Integrator
	{
		public:
			Brute(double T,double deltaT) : Integrator(T,deltaT)
			{
				Name = "Brute";
			}
			Brute(double T, double deltaT, int skipper) : Integrator(T,deltaT, skipper)
			{
				Name = "Brute";
			}
		
		private:
		
		void virtual UpdatePosition(double t)
		{
			w = 0.5 * InvMult(J,q.Conjugate() * p);
			w.Scalar() = 0;
			Quaternion qDot = 0.5 * q * w;
			Quaternion pDot = 0.5 * p * w + GradU(t);
			
	
			p = p + TimeStep * pDot;
			q = q + TimeStep * qDot;
			
			q/=q.Norm();
		
		}
	};
}
