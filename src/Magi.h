#pragma once
#include "Integrator.h"


namespace QDynamics
{
	enum StepType {Brute, Euler, Leapfrog};
	 
	template<unsigned int order, StepType step>
	class Magi : public Integrator
	{
		public:
			const unsigned int Order = order;
			const StepType StepMode = step;
			
			Magi(double T,double deltaT, int resolution) : Integrator(T,deltaT), Resolution(resolution)
			{
				Name = "MagM" + std::to_string(Order);
			}
			Magi(double T, double deltaT, int resolution, int skipper) : Integrator(T,deltaT, skipper),  Resolution(resolution)
			{
				Name = "MagM" + std::to_string(Order);
			}
		
		protected:
			const int Resolution;
			
		void virtual UpdatePosition(double t)
		{		
			switch (StepMode)
			{
				case Euler:
					Euler_Update(t);
					break;
				case Leapfrog:
					Leapfrog_Update(t);
					break;	
				case Brute:
					Brute_Update(t);
					break;
					
			}
			

		}
		
		void Brute_Update(double t)
		{
			w = 0.5 * InvMult(J,q.Conjugate() * p);
			w.Scalar() = 0;
			Quaternion qDot = 0.5 * q * w;
			Quaternion pDot = 0.5 * p * w + GradU(t);
			
	
			p = p + TimeStep * pDot;
			q = q + TimeStep * qDot;
			
			q/=q.Norm();
			
		}
		
		void Euler_Update(double t)
		{
			Quaternion tau = 0.5 * q.Conjugate() * GradU(t);
			L = L  + TimeStep * tau ;
				
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
		}
		
		void Leapfrog_Update(double t)
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
	
		Quaternion virtual Magnus(double duration)
		{
			if (Order == 0 || abs(L.Vector(2)) < 1e-8)
			{
				return duration * InvMult(J,L);
			}
			else
			{
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
				L = L * exp(0.5 * (wt - LNorm * L.Conjugate() * wt* L)*ddt);	
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
	
	
	//~ class Leapi : public Magi
	//~ {
		//~ public:
			//~ Leapi(int order, int resolution, double T,double deltaT) : Magi(order,resolution,T,deltaT)
			//~ {
				//~ Name = "MagiML" + std::to_string(order);
			//~ }
			//~ Leapi(int order, int resolution, double T, double deltaT, int skipper) : Magi(order,resolution,T,deltaT,skipper)
			//~ {
				//~ Name = "MagiML" + std::to_string(order);
			//~ }
		
		//~ private:
		
		//~ void virtual UpdatePosition(double t)
		//~ {
			//~ //drift
			//~ L = 0.5 * q.Conjugate() * p;
			//~ Quaternion wMean = Magnus(TimeStep/2); 
			//~ Quaternion mod = exp(0.5 * wMean);
			
			//~ q = q * mod;
			//~ p = 2 * q * L;
			
			//~ //kick
			//~ Quaternion gU = GradU(t + TimeStep/2);
			//~ p = p + TimeStep * gU;
			//~ L = 0.5 * q.Conjugate() * p;
			
			//~ wMean = Magnus(TimeStep/2);
			
			//~ mod = exp(0.5 * wMean);
			//~ q = q * mod;
			//~ p = 2 * q * L;
			
		//~ }
	//~ };
};
