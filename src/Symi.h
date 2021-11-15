#pragma once
#include "Integrator.h"


namespace QDynamics
{
	class Symi : public Magi
	{
		public:
			Symi(int order, double T,double deltaT) : Magi(order,0,T,deltaT)
			{
				Name = "Symi" + std::to_string(order);
				//~ MagnusBuffers = std::vector<std::vector<Quaternion>>(order,std::vector<Quaternion>(resolution,Quaternion::Zero()));
			}
			Symi(int order, double T, double deltaT, int skipper) : Magi(order,0,T,deltaT)
			{
				Name = "Symi" + std::to_string(order);
				//~ MagnusBuffers = std::vector<std::vector<Quaternion>>(order,std::vector<Quaternion>(resolution,Quaternion::Zero()));
			}
		
		protected:
				
		
		Quaternion virtual Magnus_Compute(double duration)
		{
			double lParallel = sqrt(L.Vector(0)*L.Vector(0) + L.Vector(1) * L.Vector(1));
			double varphi = atan2(L.Vector(1),L.Vector(0));
			double zeta =  L.Vector(2) * (J[1] - J[3])/(J[1] * J[3]);
			
			double s = sin(zeta * duration - varphi);
			double c = cos(zeta * duration - varphi);
			
			//Update momentum
			L.Vector(0) = lParallel * c;
			L.Vector(1) = -lParallel * s;
			
			//first order
			double pre =lParallel/(zeta * J[1]); 
			double m1x = pre * (s + sin(varphi));
			double m1y= pre * (c - cos(varphi));
			double m1z = L.Vector(2)/J[3] * duration;
			
			Quaternion M(0,m1x,m1y,m1z);
			
			if (Order > 1)
			{
				double sHalf = sin(zeta * duration/2);
				double cHalf = cos(zeta * duration/2);
				double fac =  ( zeta*duration*cHalf - 2*sHalf);
				double pre2 = 2 * lParallel *L.Vector(2)/(zeta*zeta * J[1]*J[3]);
				double m2x = pre2 * cos(zeta * duration/2 - varphi) *fac;
				double m2y = - pre2 * sin(zeta * duration/2 - varphi) * fac;
				double m2z = lParallel*lParallel/(zeta*zeta * J[1]*J[1]) * (zeta * duration - sin(zeta*duration));
				
				M.Vector(0) += m2x;
				M.Vector(1) += m2y;
				M.Vector(2) += m2z;
				
				//~ std::cout << m2x << "  " << m2y << "  " << m2z << std::endl;
			}
			//~ std::cout << M[0] << std::endl;
			return M;
		}
	};
	
	
	class SymiL : public Symi
	{
		public:
			SymiL(int order, double T,double deltaT) : Symi(order,T,deltaT)
			{
				Name = "SymiL" + std::to_string(order);
			}
			SymiL(int order,double T, double deltaT, int skipper) : Symi(order,T,deltaT,skipper)
			{
				Name = "SymiL" + std::to_string(order);
			}
		
		private:
		
		void virtual UpdatePosition(double t)
		{
			//drift
			L = 0.5 * q.Conjugate() * p;
			Quaternion wMean = Magnus(TimeStep/2); 
			Quaternion mod = exp(0.5 * wMean);
			
			q = q * mod;
			p = 2 * q * L;
			
			//kick
			Quaternion gU = GradU(t + TimeStep/2);
			p = p + TimeStep * gU;
			L = 0.5 * q.Conjugate() * p;
			
			wMean = Magnus(TimeStep/2);
			
			mod = exp(0.5 * wMean);
			q = q * mod;
			p = 2 * q * L;
			
		}
	};
}
