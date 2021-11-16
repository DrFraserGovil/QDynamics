#include <iostream>
#include <vector>
#include "JSL.h"
#include "QDynamics.h"


using QDynamics::Quaternion;

const Quaternion z(0,0,0,1);
double U0 = 3;
Quaternion QDynamics::Integrator::GradU(double t)
{
	double cosTheta =  (q * z * q.Conjugate()).Dot(z);
	cosTheta = std::min(1.0,std::max(-1.0,cosTheta));
	double theta = acos(cosTheta);

	Quaternion gradVerticalAngle = -2 * U0* z * q * z;
	
	double s = sin(theta);
	double sincInv; 
	if (abs(s) < 1e-8)
	{
		sincInv = 1;
	}
	else
	{
		sincInv = theta / sin(theta);
	}
	Quaternion gradV = 2 * sincInv * gradVerticalAngle;
	
	
	return gradV - (gradV.Dot(q)) * q;
}


double QDynamics::Integrator::U(double t)
{
	//~ return U0 * (q * z * q.Conjugate()).Dot(z);
	double cosTheta =  (q * z * q.Conjugate()).Dot(z);
	cosTheta = std::min(1.0,std::max(-1.0,cosTheta));
	double theta = acos(cosTheta);
	return U0 * theta * theta;
}

void TrivialTest(double T, double dT,int skipResolution,bool bruteOnly)
{
	double omega = 1;
	double phi = 0;
	U0 = 0;
	JSL::Vector J({2,2,2,1});
	
	Quaternion wInit(0,3*omega,0,omega);
	Quaternion r = Quaternion::Random();
	
	Quaternion qInit(cos(phi/2),sin(phi/2),0,0);
	Quaternion pInit= 2 * qInit * Mult(J,wInit);
	
	std::string folder = "Output/Trivial/";
	QDynamics::BruteInt B(T,dT/10,skipResolution);
	B.Evolve(qInit,pInit,J,folder);
	
	if (!bruteOnly)
	{
		
	
		//~ QDynamics::Magi M1(1,100,T,dT,skipResolution);
		//~ M1.Evolve(qInit,pInit,J,folder);
	
		//~ QDynamics::Magi M2(2,100,T,dT,skipResolution);
		//~ M2.Evolve(qInit,pInit,J,folder);
		
		//~ QDynamics::Leapi L1(1,50,T,dT,skipResolution);
		//~ L1.Evolve(qInit,pInit,J,folder);
		
		//~ QDynamics::Leapi L2(2,50,T,dT,skipResolution);
		//~ L2.Evolve(qInit,pInit,J,folder);
		
		//~ QDynamics::Symi S1(1,T,dT,skipResolution);
		//~ S1.Evolve(qInit,pInit,J,folder);
		
		//~ QDynamics::SymiL SL1(1,T,dT,skipResolution);
		//~ SL1.Evolve(qInit,pInit,J,folder);
	}
}

void HarmonicTest(double T, double dT,bool bruteOnly,int skipResolution)
{
	double omega = 25;
	double phi = 2.2;
	U0 = 3;
	JSL::Vector J({2,2,2,1});
	
	Quaternion wInit(0,omega,0,omega);
	Quaternion r = Quaternion::Random();
	
	Quaternion qInit(cos(phi/2),sin(phi/2),0,0);
	Quaternion pInit= 2 * qInit * Mult(J,wInit);
	

	std::string folder = "Output/Harmonic/";
	
	
	QDynamics::BruteInt B(T,dT,skipResolution);
	B.Evolve(qInit,pInit,J,folder);
	if (!bruteOnly)
	{
	
		QDynamics::Magi<1, QDynamics::Euler> M1(T,dT,50,skipResolution);
		M1.Evolve(qInit,pInit,J,folder);

		//~ QDynamics::Leapi L1(1,100,T,dT,skipResolution);
		//~ L1.Evolve(qInit,pInit,J,folder);
		
		//~ Leapi L2(2,50,T,dT,skipResolution);
		//~ L2.Evolve(qInit,pInit,J,folder);
		
		//~ QDynamics::Symi S1(1,T,dT,skipResolution);
		//~ S1.Evolve(qInit,pInit,J,folder);
		
		//~ Symi S2(2,T,dT,skipResolution);
		//~ S2.Evolve(qInit,pInit,J,folder);
		
		//~ SymiL SL1(1,T,dT,skipResolution);
		//~ SL1.Evolve(qInit,pInit,J,folder);
		
		//~ SymiL SL2(2,T,dT,skipResolution);
		//~ SL2.Evolve(qInit,pInit,J,folder);
		
		//~ SymiL SL3(2,T,dT/10,skipResolution*2);
		//~ SL3.Evolve(qInit,pInit,J,folder);
	}
}

int main(int argc, char * argv[])
{
	double T = 30;
	double dT = 0.0001;

	
	for (int i = 0; i < 2;++i)
	{
		bool bruteOnly = true;
		if (i  == 0)
		{
			bruteOnly = false;
		}
		//~ TrivialTest(T,dT/pow(10,i),1*pow(10,i),bruteOnly);

		//~ HarmonicTest(T,dT/pow(10,i),bruteOnly,10*pow(10,i));
	}
	
	QDynamics::Symi<3,QDynamics::Brute> B(100,0.1,20);
	std::cout << B.Order << "  " << B.StepMode << std::endl;

}
