#pragma once
#include "Integrator.h"


namespace QDynamics
{
	class Magi : public Integrator
	{
		public:
			double unlock = 2600.07;
				bool unlocked = false;
			
			Magi(int order, int resolution, double T,double deltaT) : Integrator(T,deltaT), Order(order), Resolution(resolution)
			{
				Name = "MagM" + std::to_string(order);
				//~ MagnusBuffers = std::vector<std::vector<Quaternion>>(order,std::vector<Quaternion>(resolution,Quaternion::Zero()));
			}
			Magi(int order, int resolution, double T, double deltaT, int skipper) : Integrator(T,deltaT, skipper), Order(order), Resolution(resolution)
			{
				Name = "MagM" + std::to_string(order);
				//~ MagnusBuffers = std::vector<std::vector<Quaternion>>(order,std::vector<Quaternion>(resolution,Quaternion::Zero()));
			}
		
		protected:
			const int Order;
			//~ std::vector<std::vector<Quaternion>> MagnusBuffers;
			const int Resolution;
			
		void virtual UpdatePosition(double t)
		{
			if (t > unlock)
			{
				unlocked = true;
			}
			if ( t > unlock + 0.1)
			{
				unlocked = false;
			}
			Quaternion tau = 0.5 * q.Conjugate() * GradU(t);
			L = L  + TimeStep * tau ;
	
			
			
			if (unlocked)
			{
				const Quaternion z(0,0,0,1);
				double cosTheta =  (q * z * q.Conjugate()).Dot(z);
				
				cosTheta = std::min(1.0,std::max(-1.0,cosTheta));
				double theta = acos(cosTheta);
			
				Quaternion gradVerticalAngle = -2 * z * q * z;
				
				double sincInv = theta / sin(theta);
				if (abs(theta) < 1e-8)
				{
					sincInv = 1;
				}
				Quaternion gradV = 2 * sincInv * gradVerticalAngle;
				Quaternion gradU = gradV - (gradV.Dot(q))* q;
				std::cout << "\n\n" << t << " q = " << q << "p = " << p << " L = " << L << "tau = " << tau << std::endl;
				std::cout << "cosTheta = " << cosTheta << " theta = " << theta << "  gradCosTheta = " << gradVerticalAngle << "  invSinc = " << sincInv << "  gradV = " << gradV << "  gradU = " << gradU << std::endl;
			}
			
			Quaternion wMean =  Magnus(TimeStep);
	
			Quaternion mod = exp(0.5 * wMean);
			
			
			if (Order == 0)
			{
				p = (p + 2 * q *tau) * mod;
				q = q * mod;
				L = 0.5 * q.Conjugate() * p;
			}
			else
			{
				q = q * mod;
				p = 2 * q * L;
			}
			
			if (unlocked)
			{
				std::cout << "w = " << wMean << " mod = " << mod << " newL = " << L << " new q = " << q  << "  new p = " << p << std::endl;
			}
		}
		
		
	
		
		Quaternion virtual Magnus(double duration)
		{
			//~ double zeta = L.Vector(2) * (J[1] - J[3])/(J[1] * J[3]);
			
			if (Order == 0 || abs(L.Vector(2)) < 1e-8)
			{
				if (unlocked)
				{
					std::cout << "Zerp" << std::endl;
				}
				return duration * InvMult(J,L);
			}
			else
			{
				if (unlocked)
				{
					std::cout << "Marg" << std::endl;
				}
				return Magnus_Compute(duration);
			}
		}
		
		Quaternion virtual Magnus_Compute(double duration)
		{
			Quaternion Lt = L;
			Quaternion wt = InvMult(J,Lt);
			
			
			Quaternion wdot;
			
			double ddt = duration/Resolution;
			
			std::vector<Quaternion> M = std::vector<Quaternion>(Order,Quaternion::Zero());
			double LNorm = 1.0/L.SqNorm();
			for (int i = 0; i < Resolution; ++i)
			{
				//~ JSL::Vector ldot = Lt.Vector().Cross(wt.Vector());
				L = L * exp(0.5 * (wt - LNorm * L.Conjugate() * wt* L)*ddt);
				//~ std::cout << LNorm << L << std::endl;
				//~ L = L + 0.5 * (L * wt - wt* L)*ddt;
				
				wt = InvMult(J,L);
				
				M[0] = M[0] + wt;
				if (Order > 1)
				{
					M[1] = M[1] + Quaternion(0,wt.Vector().Cross(M[0].Vector()));
				}
			}
			Quaternion output = M[0] * ddt;
			for (int i = 1; i < Order; ++i)
			{
				output = output + M[i] * pow(ddt,i+1);
			}
			//~ std::cout << M[0] << std::endl;
			return output;
		}
	};
	
	
	class Leapi : public Magi
	{
		public:
			Leapi(int order, int resolution, double T,double deltaT) : Magi(order,resolution,T,deltaT)
			{
				Name = "MagiML" + std::to_string(order);
			}
			Leapi(int order, int resolution, double T, double deltaT, int skipper) : Magi(order,resolution,T,deltaT,skipper)
			{
				Name = "MagiML" + std::to_string(order);
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
};
