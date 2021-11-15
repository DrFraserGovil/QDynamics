#pragma once
#include <iostream>
#include <vector>
#include "JSL.h"
#include "DynamicOperations.h"


namespace QDynamics
{
	const int buffersize = 100;
	const int hash = 32;
	
	
	class Integrator
	{
		public:
			
			Integrator(double T,double deltaT)
			{
				Name = "Unassigned";
				Skips = 1;
				TotalTime = T;
				TimeStep = deltaT;
				BufferSize = buffersize;
				NHashes = hash;
			}
			
			Integrator(double T, double deltaT, int skipper)
			{
				Name = "Unassigned";
				Skips = skipper;
				TotalTime = T;
				TimeStep = deltaT;
				NHashes = hash;
				BufferSize = buffersize;
			}
			
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
						t = TotalTime;
					}
				}
				
				Buffer.resize(BufferPos);
				FlushBuffer();
				UpdateProgressBar(TotalTime);
			}
			
			
			void virtual UpdatePosition(double t)
			{
				
			}
		protected:
			Quaternion q;
			Quaternion p;
			Quaternion w;
			Quaternion L; 
			JSL::Vector J;
			
			double TimeStep;
			
			std::string Name;
			std::string FileName;
			
			virtual Quaternion GradU(double t);
			virtual double U(double t);
		private:
			int BufferPos;
			int BufferSize;
			
			int NHashes;
			bool FinalHash;
			int CurrentHashes;
			std::vector<std::string> Buffer;
			int SkipID;
			int Skips;
			double TotalTime;
			
			void UpdateBuffer(double t)
			{
				int prec = 10;
				++SkipID;
				if (t==0 || SkipID == Skips)
				{
					std::ostringstream out;
					out.precision(10);
					out << std::fixed << t;
					
					std::string b = out.str() + "," + q.to_string_precision(prec) + "," + p.to_string_precision(prec) + ", " + std::to_string(Hamiltonian(t)) + ", " + std::to_string(LabAngularMomentum(q,p).Norm());
					Buffer[BufferPos] = b;
					
					++BufferPos;
					if (BufferPos >= BufferSize)
					{
						FlushBuffer();
					}
					SkipID = 0;
				}	
			}
			
			
			void FlushBuffer()
			{
				JSL::writeVectorToFile(FileName,Buffer,"\n",true);
				BufferPos = 0;
			}
			
			
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
				
				std::cout << "New " << Name << " Integrator Initialised" << "\n ["<< std::flush;
				UpdateBuffer(0);
			}
			
			void CreateFullName(std::string saveFolder)
			{
				JSL::mkdir(saveFolder);
				FileName = saveFolder + "/" + Name +"_N" + std::to_string((int)log10(TotalTime/TimeStep)) + ".dat";
			}
			
			double Hamiltonian(double t)
			{	
				Quaternion Lt = 0.5 * q.Conjugate() * p;
				double T = 0.5 * Lt.Dot( InvMult(J,Lt));
				double E = T + U(t);	
				return E;	
			}
				
	};
}
