#pragma once
#include <iostream>
#include <vector>
#include "JSL.h"


namespace QDynamics
{
QDynamics::Quaternion Mult(const JSL::Vector &J, const QDynamics::Quaternion & q)
{
	QDynamics::Quaternion m = q;
	m.Scalar() *= J[0];
	for (int i = 0; i < 3; ++i)
	{
		m.Vector(i) *= J[i+1];
	}
	return m;
}
QDynamics::Quaternion InvMult(const JSL::Vector & J, const QDynamics::Quaternion & q)
{
	QDynamics::Quaternion m = q;
	m.Scalar() /= J[0];
	for (int i = 0; i < 3; ++i)
	{
		m.Vector(i) /= J[i+1];
	}
	return m;
}

JSL::Vector LabAngularMomentum(const QDynamics::Quaternion &q, const QDynamics::Quaternion &p)
{
	QDynamics::Quaternion L = 0.5 * q.Conjugate() * p;
	QDynamics::Quaternion LabL = q * L * q.Conjugate();
	
	return LabL.Vector();
	
}
}
