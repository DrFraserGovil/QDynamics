#pragma once
#include "Integrator.h"


namespace QDynamics
{
	template<unsigned int order, StepType step>
	class Symi : public Magi<order, step>
	{
		public:
		
			Symi(double T,double deltaT) : Magi<order,step>(T,deltaT,0)
			{
				Symi::Name = "Symi" + std::to_string(Symi::Order);
				//~ MagnusBuffers = std::vector<std::vector<Quaternion>>(order,std::vector<Quaternion>(resolution,Quaternion::Zero()));
			}
			Symi(double T, double deltaT, int skipper) : Magi<order,step>(T,deltaT,skipper,0)
			{
				Symi::Name = "Symi" + std::to_string(Symi::Order);
				//~ MagnusBuffers = std::vector<std::vector<Quaternion>>(order,std::vector<Quaternion>(resolution,Quaternion::Zero()));
			}
		
		protected:
				
		
		Quaternion virtual Magnus_Compute(double duration)
		{
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
	};
}
