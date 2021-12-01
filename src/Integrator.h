#pragma once
#include <iostream>
#include <vector>
#include "JSL.h"
#include "DynamicOperations.h"


namespace QDynamics
{
	const int buffersize = 100;
	const int hash = 32;
	
	/*!
	 * A superclass from which all individual integrators will inherit. Handles the mundane stuff such as the progress tracking and saving. Importantly, it provides virtual overload methods for Gradu() and U() which allows external definition of the dynamics, and of UpdatePosition() which allows specific subclasses of this object choose how they want to move from q_i -> q_i+1.
	*/
	class Integrator
	{
		public:
			
			//! Constructor function \param T The total duration of the integration \param deltaT the width of timesteps used
			Integrator(double T,double deltaT)
			{
				Name = "Unassigned";
				Skips = 1;
				TotalTime = T;
				TimeStep = deltaT;
				BufferSize = buffersize;
				NHashes = hash;
			}
			
			//! Constructor function \param T The total duration of the integration \param deltaT the width of timesteps used \param skipper The number of epochs which pass in between updates to the progess Integrator::Buffer
			Integrator(double T, double deltaT, int skipper)
			{
				Name = "Unassigned";
				Skips = std::max(1,skipper);
				TotalTime = T;
				TimeStep = deltaT;
				NHashes = hash;
				BufferSize = buffersize;
			}
			
			/*! The main computation loop. Calls UpdatePosition() at each timestep until the TotalTime is reached.
				\param q0 The initial position Quaternion
				\param p0 The initial momentum Quaternion
				\param J The diagonals of the moment of inertia vector (we have implicitly assumed the system triaxial is such that J is diagonal)
				\param saveFolder The location into which the Integrator::Buffer output is saved.
			*/
			void Evolve(Quaternion q0, Quaternion p0, JSL::Vector J,std::string saveFolder)
			{
				Initialise(q0,p0,J,saveFolder);
	
				double t = 0;
				while (t < TotalTime)
				{
					UpdatePosition(t);
					t+=TimeStep;
					UpdateBuffer(t);
					UpdateProgressBar(t);
					
					if (q.isnan() || p.isnan())
					{
						std::cout << "Quitting due to erroneous state" << std::endl;
						t = TotalTime;
					}
				}
				
				Buffer.resize(BufferPos);
				FlushBuffer();
				UpdateProgressBar(TotalTime);
			}
			
			//! In the Integrator base class, this is an empty function. Child members will mainly focus on overloading this function to their own update formula \param t The current time of the integrator
			void virtual UpdatePosition(double t)
			{
				
			}
		protected:
			Quaternion q; //!< The current position variable
			Quaternion p; //!< The current momentum variable
			Quaternion w; //!< The current body-fixed angular speed
			Quaternion L;  //!< The current body-fixed Angular momentum
			JSL::Vector J; //!< The saved location of the moment of inertia matrix (see Evolve())
			
			double TimeStep; //!< The timestep taken by the integrator
			double TotalTime; //!< The total length of time that the integrator will search for.
			
			std::string Name; //!< The assigned name of the integrator. This is set to "Unassigned" by default: child classes should modify this name. 
			std::string FileName; //!< The location written to by CreateFullName(). A file location which will be used to write the output buffer to.
			
			//! Returns the quaternionic derivative of the neo-potential V (equal to the derivative of the rotation potential U), given the current value of Integrator::q and the current time \param t The current time \return The value of the quaternionic gradient
			virtual Quaternion GradU(double t);
			
			//! Returns the current value of the rotation potential U, given the current value of Integrator::q and the current time. Unlike GradU(), it is not critical that this function be 100% correct, as it is only used by the Integrator::Hamiltonian function as bookkeeping for the current energy. If you do not wish to track the current energy, you may neglect this function.  \param t The current time \return The current potential energy of the system
			virtual double U(double t);
		private:
			int BufferPos; //!< The current write-index of Integrator::Buffer
			int BufferSize; //!< The length of the Integrator::Buffer
			
			int NHashes; //!< The number of hashes a full progress bar is made up of
			bool FinalHash; //!< A simple lock tracking if the final hash has been written, preventing multiple termination characters being printed
			int CurrentHashes;//!< The current number of progress hashes which have been written to the terminal
			
			std::vector<std::string> Buffer; //!< A string-buffer containing the save-values of several past values of the integrator. Once full, the buffer is written to file.
			int SkipID; //!< The current number of files that the buffer has skipped (as per the skipper variable in Integrator(double,double,int))
			int Skips; //!< The maximum value of SkipID before the buffer saves a value.
			
			
			//! Called during every loop of Evolve(), saves some values of interest to the Integrator::Buffer. \param t The current time, needed as part of the buffer input
			void UpdateBuffer(double t)
			{
				
				int prec = 15;
				++SkipID;
				if (t==0 || SkipID == Skips)
				{
					std::ostringstream out;
					out.precision(prec);
					out << std::fixed << t << "," << q.to_string_precision(prec) << "," << p.to_string_precision(prec) << "," << std::fixed <<Hamiltonian(t) << "," << std::fixed  << (0.5 * q.Conjugate() * p).Norm();
					
					
					std::string b = out.str();
					Buffer[BufferPos] = b;
					
					++BufferPos;
					if (BufferPos >= BufferSize)
					{
						FlushBuffer();
					}
					SkipID = 0;
				}	
			}
			
			//! Called when the Integrator::Buffer is full, writes the contents of the buffer to the file determined by Integrator::FileName, and resets the buffer.
			void FlushBuffer()
			{
				JSL::writeVectorToFile(FileName,Buffer,"\n",true);
				BufferPos = 0;
			}
			
			//! Writes a series of # to the screen to act as a progress bar: [#####             ] \param t The current time, so that the current progress can be computed.
			void UpdateProgressBar(double t)
			{
				int hashesRequired = std::min((int)(t/TotalTime * NHashes),NHashes);
				while (hashesRequired > CurrentHashes)
				{
					std::cout << "#" << std::flush;
					++CurrentHashes;
				}
				if (CurrentHashes == NHashes && !FinalHash)
				{
					FinalHash = true;
					std::cout << "]\n";
				}		
			}
			
			
			/*!
			 * Prepares the Integrator for a new computation loop. Called at the start of an Evolve() process, it wipes the buffer and prepares the output file. The inputs are simply those passed to Evolve()
			 * \param q0 The initial position Quaternion
				\param p0 The initial momentum Quaternion
				\param Jin The diagonals of the moment of inertia vector (we have implicitly assumed the system triaxial is such that J is diagonal)
				\param saveFolder The location into which the Integrator::Buffer output is saved.
			 * 
			*/
			void Initialise(Quaternion q0, Quaternion p0,JSL::Vector Jin,std::string saveFolder)
			{
				CreateFullName(saveFolder);
				JSL::initialiseFile(FileName);
				std::vector<std::string> headerNames = {"t","q0","q1","q2","q3","p0","p1","p2","p3","H","L"};
				JSL::writeVectorToFile(FileName,headerNames,",",true);
				Buffer = std::vector<std::string>(BufferSize);
				BufferPos = 0;
				SkipID = 0;
				CurrentHashes = 0;
				FinalHash = false;
				q = q0;
				p = p0;
				L = 0.5 * q.Conjugate() * p;
				w = InvMult(Jin,L);
				
				if (Jin.Size() == 4)
				{
					J = Jin;
				}
				else
				{
					exit(-1);
				}
				
				std::cout << Name << " Integrator: ["<< std::flush;
				UpdateBuffer(0);
				
			}
			
			//! Computes the Integrator::FileName including the filepath from from the Integrator::Name and the provided directory. \param saveFolder The corresponding value passed to Intitialise() 
			void CreateFullName(std::string saveFolder)
			{
				JSL::mkdir(saveFolder);
				FileName = saveFolder + "/" + Name +"_N" + std::to_string((int)round(10*log10(TotalTime/TimeStep))) + ".dat";
			}
			
			//! Computes the current value of the energy of the system given the kinetic energy and the value of U(). \param t The current time, for time-dependent potentials.
			double Hamiltonian(double t)
			{	
				
				//~ Quaternion Lt = 0.5 * q.Conjugate() * p;
				//~ L = 2 * L;
				double T = 0.5 * L.Dot( InvMult(J,L));
			
				double E = T + U(t);	
				//~ std::cout <<t << "   " << T << "  " << U(t) << " " << E << "   " << L << std::endl;
				//~ if( rand() < RAND_MAX /10)
				//~ {
					//~ exit(10);
				//~ }
				//~ std::cout << t << "  " << T << "  " << U(t) << "   " << E << std::endl;
				return E;	
			}
				
	};
}
