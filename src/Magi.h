#pragma once
#include "Integrator.h"


namespace QDynamics
{
	
	enum UpdateType {
			
		Brute, //!< The dumbest way to update the formula, directly evaluates linear update formula. Note that Brute integrators do **not** compute the Magnus series.
		Euler, //!< A kick-drift integrator. Updates the momentum with an instantaneous torque-impulse, and then undergoes varying degrees of torque free precession (depending on the Order of the Magnus() series)
		Leapfrog //! A drift-kick-drift integrator, distinct from the Euler in that the torque is evaluated directly in the middle of a timestep, making the integrator time-reversible.
	};
	 
	 
	/*
	 * Computes the mean 
	 * 
	*/
	template<unsigned int order, UpdateType step>
	class Magi : public Integrator
	{
		public:
			
			
			//!Constructor, as per Integrator::Integrator(double,double). All arguments except *resolution* are used to build the parent class. \param T The total duration of the integration \param deltaT the width of timesteps used \param resolution The number of points to used whenever a numerical integral is required. 
			Magi(double T,double deltaT, int resolution) : Integrator(T,deltaT), Resolution(resolution)
			{
				MakeName();
			}
			//!Constructor, as per Integrator::Integrator(double,double,int) \param T The total duration of the integration \param deltaT the width of timesteps used \param  resolution The number of points to used whenever a numerical integral is required. \param skipper The number of epochs which pass in between updates to the progess Integrator::Buffer
			Magi(double T, double deltaT, int resolution, int skipper) : Integrator(T,deltaT, skipper),  Resolution(resolution)
			{
				MakeName();
			}
		
		protected:
			const unsigned int Order = order; //!< The first template argument passed to a ``Magi`` object. Determines the order of the computation in the Magnus() calculation.
			const UpdateType StepMode = step; //!< The second template argument passed to a ``Magi`` object. Determines the update formula used.
			const int Resolution; //!< As per the assignation in Magi(), this is the number of points to used whenever a numerical integral is required. 
			
			//! An override of the generic Integrator::UpdatePosition(). Depending on the value of Magi::StepMode, calls an appropriate function. The intent is that as Magi::StepMode is a template value, logic of which function is called is determined at compile time, without additional overhead. \param t The current time, used by Integrator::GradU()
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
			
			//! The function called by UpdatePosition() if ``StepMode == Brute``. See QDynamics::UpdateType for more information. \param t The current time
			void Brute_Update(double t)
			{
				w = 0.5 * InvMult(J,q.Conjugate() * p);
				w.Scalar() = 0;
				Quaternion qDot = 0.5 * q * w;
				Quaternion pDot = 0.5 * p * w - GradU(t);
				
		
				p = p + TimeStep * pDot;
				q = q + TimeStep * qDot;
				
				q/=q.Norm();
				
			}
			
			//! The function called by UpdatePosition() if ``StepMode == Euler``. See QDynamics::UpdateType for more information. \param t The current time
			void Euler_Update(double t)
			{
				Quaternion gU = GradU(t);
				p = p - TimeStep * gU;
				L = 0.5 * q.Conjugate() * p;
				
				Quaternion wMean =  Magnus(TimeStep);
				Quaternion mod = exp(wMean);

				q = q * mod;
				p = p * mod;	
				
				//~ if (Order == 0)
				//~ {
					//~ p = (p + 2 * TimeStep * q *tau) * mod;
					//~ q = q * mod;
					//~ L = 0.5 * q.Conjugate() * p;
				//~ }
				//~ else
				//~ {
								
				//~ }
			}
			
			//! The function called by UpdatePosition() if ``StepMode == Leapfrog``. See QDynamics::UpdateType for more information. \param t The current time
			void Leapfrog_Update(double t)
			{
				//drift
				//~ L = 0.5 * q.Conjugate() * p;
				Quaternion wMean = Magnus(TimeStep/2); 
				Quaternion mod = exp(wMean);
				
				q = q * mod;
				p = p * mod;
				
				//kick
				Quaternion gU = GradU(t + TimeStep/2);
				//~ p = p - TimeStep * gU;
				L = L - 0.5 * TimeStep * q.Conjugate()* gU;
				
				wMean = Magnus(TimeStep/2);
				
				mod = exp(wMean);
				q = q * mod;
				p = 2 * q * L;
			}
		
			//! Generates an estimate of the Magnus series for torque-free drift over the provided duration, given the current value of the angular momentum, Integrator::L. The duration of this compuation is determined by Magi::Order (which determines how many terms are included) and by Magi::Resolution (which determines the accuracy with which each term is computed). \param duration The length of the torque-free drift to be computed
			Quaternion virtual Magnus(double duration)
			{
				if (Order == 0 || abs(L.Vector(2)) < 1e-8)
				{
					return duration * InvMult(J,L) /2;
				}
				
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
					
					M[0] = M[0] + 0.5*wt;
					if (Order > 1)
					{
						M[1] = M[1] - 0.5*Quaternion(0,wt.Vector().Cross(M[0].Vector()));
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
		
		private:
			//! Generates a print-and-filename-friendly name depending on the properties of the integrator. 
			virtual void MakeName()
			{
				Name = "Mag";
				switch (StepMode)
				{
					case (Brute):
						Name = "Brute";
						break;
					case (Euler):
						Name += "_Euler";
						break;
					case (Leapfrog):
						Name += "_Frog";
						break;
				}
				Name +=  "_" + std::to_string(Order);
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
