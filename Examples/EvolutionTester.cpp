#include <iostream>
#include <vector>
#include "JSL.h"
#include "QDynamics.h"


using QDynamics::Quaternion;

const Quaternion z(0,0,0,1);
double U0 = 13;
Quaternion QDynamics::Integrator::GradU(double t)
{
	double cosTheta =  (q * z * q.Conjugate()).Dot(z);
	cosTheta = std::min(1.0,std::max(-1.0,cosTheta));
	
	Quaternion pendulum = -2 * U0* z * q * z;
	//~ Quaternion gradV = pendulum;
	double theta = acos(cosTheta);
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
	Quaternion gradV = - 2 * sincInv * pendulum;
	
	
	return gradV - (gradV.Dot(q)) * q;
	//~ return Quaternion::Zero();
}


double QDynamics::Integrator::U(double t)
{
	double cosTheta =  (q * z * q.Conjugate()).Dot(z);
	cosTheta = std::min(1.0,std::max(-1.0,cosTheta));
	//~ double theta = acos(cosTheta);
	//~ std::cout << theta << std::endl;
	//~ return U0 * theta * theta;
	return U0 * cosTheta;
	//~ return 0;
}

void TrivialTest(double logResolution)
{
	double T = 1000;
	double dTInit = 0.1;
	double omega = 0;

	double phi = 2.0943951;
	U0 = 2;
	JSL::Vector J({1100,1,1,2});
	
	Quaternion wInit(0,omega,0,omega);
	Quaternion r = Quaternion::Random();
	
	Quaternion qInit(cos(phi/2),sin(phi/2),0,0);
	Quaternion pInit= 2 * qInit * Mult(J,wInit);
	
	std::string folder = "Output/Trivial/";
	
	
	double dT = dTInit;
	for (int i = 0; i < 6; ++i)
	{
		
		int skipper = 1;
		double nSteps = log10(T/dT);
		if (nSteps > logResolution)
		{
			skipper = pow(10,nSteps - logResolution);
		}
		
		QDynamics::BruteInt B(T,dT,skipper);
		B.Evolve(qInit,pInit,J,folder);
	
		if (i < 2)
		{
			QDynamics::Magi<0, QDynamics::Euler> M0(T,dT,10,skipper);
			M0.Evolve(qInit,pInit,J,folder);
		}
		dT = dT/10;
	}

}

void OrderTest(double logResolution)
{

	double T = 3000;
	double dTInit = 0.1;
	double omega = 10;

	double phi = 2.0943951;
	U0 = 2;
	JSL::Vector J({1100,1,1,2});
	
	Quaternion wInit(0,omega,0,omega);
	Quaternion r = Quaternion::Random();
	
	Quaternion qInit(cos(phi/2),sin(phi/2),0,0);
	Quaternion pInit= 2 * qInit * Mult(J,wInit);
	
	std::string folder = "Output/Order/";
	
	int N = 4;
	double dT = dTInit;
	for (int i = 0; i < N; ++i)
	{
		
		int skipResolution= 1;
		double nSteps = log10(T/dT);
		if (nSteps > logResolution)
		{
			skipResolution = pow(10,nSteps - logResolution);
		}
		
		
	
		if (i < N -1)
		{
			QDynamics::Magi<0, QDynamics::Euler> M0(T,dT,0,skipResolution);
			M0.Evolve(qInit,pInit,J,folder);
			
			QDynamics::Magi<1, QDynamics::Euler> M1(T,dT,150,skipResolution);
			M1.Evolve(qInit,pInit,J,folder);
			
			QDynamics::Magi<2, QDynamics::Euler> M2(T,dT,150,skipResolution);
			M2.Evolve(qInit,pInit,J,folder);

			
			QDynamics::Symi<1, QDynamics::Euler> S1(T,dT,skipResolution);
			S1.Evolve(qInit,pInit,J,folder);
			
			QDynamics::Symi<2, QDynamics::Euler> S2(T,dT,skipResolution);
			S2.Evolve(qInit,pInit,J,folder);
		}
		else
		{
			//~ dT = dT/10;
			QDynamics::Symi<2, QDynamics::Leapfrog> S2L(T,dT,skipResolution);
			S2L.Evolve(qInit,pInit,J,folder);
		}
		dT = dT/10;
	}
	
	
}



void MagTest(double logResolution)
{

	double T = 100;
	double dTInit = 0.1;
	double omega = 0.5;

	double phi = 2.0943951;
	U0 = 5;
	JSL::Vector J({1100,1,2,3});
	
	Quaternion wInit(0,0,omega,omega);
	Quaternion r = Quaternion::Random();
	
	Quaternion qInit(cos(phi/2),sin(phi/2),0,0);
	Quaternion pInit= 2 * qInit * Mult(J,wInit);
	
	std::string folder = "Output/MagPot/";
	
	int N = 3;
	double dT = dTInit;
	for (int i = 0; i < N; ++i)
	{
		
		int skipResolution= 1;
		double nSteps = log10(T/dT);
		if (nSteps > logResolution)
		{
			skipResolution = pow(10,nSteps - logResolution);
		}
		
		
	
		if (i < N -1)
		{
			QDynamics::BruteInt B(T,dT,skipResolution);
			B.Evolve(qInit,pInit,J,folder);
			
			
			QDynamics::Magi<0, QDynamics::Euler> M0(T,dT,0,skipResolution);
			M0.Evolve(qInit,pInit,J,folder);
			
			QDynamics::Magi<1, QDynamics::Euler> M1(T,dT,10,skipResolution);
			M1.Evolve(qInit,pInit,J,folder);
			
			//~ QDynamics::Magi<1, QDynamics::Euler> M10(T,dT,50,skipResolution);
			//~ M10.Evolve(qInit,pInit,J,folder);
			
			QDynamics::Magi<1, QDynamics::Euler> M100(T,dT,50,skipResolution);
			M100.Evolve(qInit,pInit,J,folder);
			

			
			QDynamics::Magi<1, QDynamics::Euler> M1000(T,dT,250,skipResolution);
			M1000.Evolve(qInit,pInit,J,folder);
			//~ QDynamics::Magi<1, QDynamics::Euler> M2(T,dT,100,skipResolution);
			//~ M2.Evolve(qInit,pInit,J,folder);

			
			//~ QDynamics::Symi<1, QDynamics::Euler> S1(T,dT,skipResolution);
			//~ S1.Evolve(qInit,pInit,J,folder);
			
			//~ QDynamics::Symi<2, QDynamics::Euler> S2(T,dT,skipResolution);
			//~ S2.Evolve(qInit,pInit,J,folder);
		}
		else
		{
			QDynamics::Magi<1, QDynamics::Leapfrog> M100(T,dT/10,150,skipResolution*10);
			M100.Evolve(qInit,pInit,J,folder);
		}
		dT = dT/10;
	}
	
	
}



void FrogTest(double logResolution)
{

	double T = 1000;
	double dTInit = 0.1;
	double omega = 2.5;

	double phi = 0.7;
	U0 = 1;
	JSL::Vector J({1100,1,1,2});
	
	Quaternion wInit(0,0,0,omega);
	Quaternion r = Quaternion::Random();
	
	Quaternion qInit(cos(phi/2),sin(phi/2),0,0);
	Quaternion pInit= 2 * qInit * Mult(J,wInit);
	
	std::string folder = "Output/Frog/";
	
	int N = 3;
	double dT = dTInit;
	for (int i = 0; i < N; ++i)
	{
		
		int skipResolution= 1;
		double nSteps = log10(T/dT);
		if (nSteps > logResolution)
		{
			skipResolution = pow(10,nSteps - logResolution);
		}
		
		
	
		if (i < N -1)
		{
			QDynamics::BruteInt B(T,dT,skipResolution);
			B.Evolve(qInit,pInit,J,folder);
			
			
			QDynamics::Symi<2, QDynamics::Euler> S0(T,dT,skipResolution);
			S0.Evolve(qInit,pInit,J,folder);
			
			QDynamics::Symi<2, QDynamics::Leapfrog> S1(T,dT,skipResolution);
			S1.Evolve(qInit,pInit,J,folder);
			
			
		}
		else
		{
			QDynamics::Symi<2, QDynamics::Leapfrog> M100(T,dT,skipResolution);
			M100.Evolve(qInit,pInit,J,folder);
		}
		dT = dT/10;
	}
	
	
}


void SphereTest(double logResolution)
{

	double T = 100;
	double dTInit = 5e-3;
	double omega = 5;

	double phi = 1.1;
	U0 = 2;
	JSL::Vector J({1100,1,1,2});
	
	
	Quaternion r = Quaternion::Random();
	
	
	
	std::string folder = "Output/Frog/";
	
	int N = 1;
	double dT = dTInit;
	int skipResolution= 1;
		double nSteps = log10(T/dT);
		if (nSteps > logResolution)
		{
			skipResolution = pow(10,nSteps - logResolution);
		}
	QDynamics::Symi<2, QDynamics::Leapfrog> S1(T,dT/100,skipResolution*100);
	QDynamics::BruteInt B(T,dT,skipResolution);
	for (int i = 0; i < N; ++i)
	{
		
		
		
		//~ omega = omega /1.3;
		//~ phi = phi * 1.8;
		Quaternion wInit(0,0,0,omega);
		Quaternion qInit(cos(phi/2),sin(phi/2),0,0);
		Quaternion pInit= 2 * qInit * Mult(J,wInit);
		folder = "Output/Frog" + std::to_string(i) + "/";
		S1.Evolve(qInit,pInit,J,folder);
		B.Evolve(qInit,pInit,J,folder);
	}
	
	
}


int main(int argc, char * argv[])
{
	double logResolution = 5;
	
	SphereTest(logResolution);
	//~ FrogTest(logResolution);
	//~ MagTest(logResolution);
	//~ TrivialTest(logResolution);
	//~ OrderTest(logResolution);
}
