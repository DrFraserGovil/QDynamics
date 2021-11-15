#pragma once
#include <vector>
#include <stdexcept>
#include <string.h>
#include <iostream>
#include <math.h>
#include "JSL.h"
//~ #include "vector.h"
//~ #include "matrix.h"



namespace QDynamics
{
	//forward declarations for circularity
	class Quaternion;
	
	/*!
		A computational replica of the mathematical object *quaternion*, which replicates the expected behaviour of such objects, including addition and scalar multiplication (which it inherits from  \verbatim embed:rst:inline
	`JSL::Vector() <https://jackstandardlibrary.readthedocs.io/en/latest/vectors.html>`_ \endverbatim ), as well as the idiosyncratic multiplication and division native ot this field. We use the notation that  \verbatim embed:rst:inline :math:`\mathsf{q} = (q_0, \vec{q})` \endverbatim , and the interface is based on this assumption. It is intentional that treating quaternions as members of  \verbatim embed:rst:inline :math:`\mathbb{R}^4` \endverbatim is awkward.
	*/
	class Quaternion : public JSL::Vector
	{
		public: 
			
			//! The current value of the Scalar component, returned as a *reference*. This allows it to be assigned to: q.Scalar() = x stores the value of x in the scalar component. \return A reference to the scalar component of the quaternion
			double & Scalar()
			{
				return Data[0];
			}
			
			//! The current value of the ith Vector component, returned as a *reference*. This allows it to be assigned to: q.Vector(i) = x stores the value of x in the scalar component. \return A reference to the ith vector component of the quaternion
			double & Vector(int i)
			{
				return Data[i+1];
			}
			
			//! The current value of the Vector component. This cannot be modified in-place or assigned to as Vector(int i) and Scalar() can be. \return A JSL::Vector object containing the Vector portion of the quaternion.
			JSL::Vector Vector() const
			{
				return JSL::Vector({Data[1],Data[2],Data[3]});
			}
			
			//!Default initialiser. Initialises to  \verbatim embed:rst:inline :math:`\mathsf{q} = (0, \vec{0})` \endverbatim
			Quaternion() : JSL::Vector(4) {	};

			//!  Initialises the object to \verbatim embed:rst:inline :math:`\mathsf{q} = (q_0, \vec{q})` \endverbatim  \param q0 The value \verbatim embed:rst:inline :math:`q_0` \endverbatim \param qVec The value \verbatim embed:rst:inline :math:`\vec{q}` \endverbatim. This must be a vector of length 3, or an error is thrown.
			Quaternion(const double & q0, const JSL::Vector & qVec) : JSL::Vector(4)
			{
				if (qVec.Size() != 3)
				{
					throw std::runtime_error("ERROR in Quaternion: Quaternions(scalar,vector) initialisation only works if the vector has dimension 3");
				}
				Data[0] = q0;
				for (int i = 0; i < 3; ++i)
				{
					Data[i+1] = qVec[i];
				}
			}
			
			//! Initialises the object to the value \verbatim embed:rst:inline :math:`\mathsf{q} = \left(a, \begin{pmatrix} b \\ c \\ d \end{pmatrix}\right)` \endverbatim \param a The scalar value of the quaternion \param b The x-component of the vector \param c The y-component of the vector \param d The z-component of the vector
			Quaternion(const double & a,const double & b,const double & c,const double & d) : JSL::Vector(4)
			{
				Data[0] = a;
				Data[1] = b;
				Data[2] = c;
				Data[3] = d;
			} 
			//! Initialises the quaternion as if it were a member of R^4. \param vec4 A JSL::Vector object of length 4 (else an error is thrown)
			Quaternion(const JSL::Vector & vec4) : JSL::Vector(vec4)
			{
				if (vec4.Size() != 4)
				{
					throw std::runtime_error("ERROR in Quaternion: Quaternions(Vector) initialisation only works if the vector has dimension 4");
				}
			}
			
			//! Initialises the quaternion as if it were a member of R^4. \param vec4 A std::vector<double> object of length 4 (else an error is thrown)
			Quaternion(const std::vector<double> & vec4) : JSL::Vector(vec4)
			{
				if (vec4.size() != 4)
				{
					throw std::runtime_error("ERROR in Quaternion: Quaternions(Vector) initialisation only works if the vector has dimension 4");
				}
			}
			
			//! A static constructor which returns the multiplicative identity \verbatim embed:rst:inline :math:`\mathsf{q} = \left(1, \vec{0}\right)` \endverbatim
			static Quaternion One()
			{
				Quaternion q;
				q[0] = 1;
				return q;
			}
			
			//! A static constructor which returns the additive identity \verbatim embed:rst:inline :math:`\mathsf{q} = \left(0, \vec{0}\right)` \endverbatim. Nominally unneeded, as the default constructor also returns 0, but this is itself based on the default constructor of the JSL::Vector object. To avoid potential future errors, the Zero() function is safer!
			static Quaternion Zero()
			{
				Quaternion q;
				return q;
			}
			
			//! Generates a quaternion which is randomly populated with values between 0 and 1.
			static Quaternion Random()
			{
				Quaternion q;
				
				for (int i = 0; i < 4; ++i)
				{
					q[i] = (double)random() / RAND_MAX;
				}
				
				return q;
			}
			
			

			//! Returns the conjugate of the current quaternion \return The value \verbatim embed:rst:inline :math:`\overline{\mathsf{q}} = (q_0, -\vec{q})` \endverbatim 
			Quaternion Conjugate() const
			{
				return Quaternion(Scalar(), -1*Vector());
			}

			//! Constructs a matrix L(q) such that L(q) b = q*b, replacing the usual quaternion multiplication operation operator*(const Quaternion & lhs, const Quaternion & rhs) \return The JSL::Matrix object which performs the multiplication
			JSL::Matrix LeftMultiplicationMatrix() const
			{
				JSL::Matrix M(4,4);
				
				M(0,0) = Data[0];
				for (int i = 1; i < 4; ++i)
				{
					M(0,i) = -Data[i];
					M(i,0) = Data[i];
					M(i,i) = Data[0];
				}
				
				M(1,2) = -Data[3];
				M(2,1) = Data[3];
				M(1,3) = Data[2];
				M(3,1) = -Data[2];
				M(2,3) = -Data[1];
				M(3,2) = Data[1];
				
				return M;
				
			}
			
			//! Constructs a matrix R(q) such that R(q) b = b*q, replacing the usual quaternion multiplication operation operator*(const Quaternion & lhs, const Quaternion & rhs) \return The JSL::Matrix object which performs the multiplication
			JSL::Matrix RightMultiplicationMatrix() const
			{
				JSL::Matrix M(4,4);
				
				M(0,0) = Data[0];
				for (int i = 1; i < 4; ++i)
				{
					M(0,i) = -Data[i];
					M(i,0) = Data[i];
					M(i,i) = Data[0];
				}
				
				M(1,2) = Data[3];
				M(2,1) = -Data[3];
				M(1,3) = -Data[2];
				M(3,1) = Data[2];
				M(2,3) = Data[1];
				M(3,2) = -Data[1];
				
				return M;
				
			}
			
			//!An annoying double-definition of Scalar(), required for when the quaternion object is a const (and hence the references need to be carefully guarded)
			const double & Scalar() const
			{
				return Data[0];
			}
			//!An annoying double-definition of Scalar(int i), required for when the quaternion object is a const (and hence the references need to be carefully guarded)
			const double & Vector(int id) const
			{
				return Data[id+1];
			}
			
		protected:

			//! We put the Vector brace operator as protected such that external members cannot call Data[i] and mess around with the internals, forcing them to go through the Scalar/vector interface. However, internal functions can still call the brace operations and it works as expected. 
			using JSL::Vector::operator[];
	};
	
	//! The crucial quaternion multiplication operation, defined such that: \verbatim embed:rst:inline :math:`\mathsf{a} \otimes \mathsf{b} = \begin{pmatrix} a_0 b_0 - \vec{a} \cdot \vec{b} \\ a_0 \vec{b} + b_0 \vec{a} + \vec{a} \times \vec{b}\end{pmatrix}` \endverbatim. Note that this product is highly non-commutative in general. \param lhs The first argument of the operation (\verbatim embed:rst:inline :math:`\mathsf{a}` \endverbatim) \param rhs The second argument (\verbatim embed:rst:inline :math:`\mathsf{b}` \endverbatim) \return The quaternion product (\verbatim embed:rst:inline :math:`\mathsf{a}\otimes\mathsf{b}` \endverbatim)
	inline Quaternion operator*(const Quaternion & lhs, const Quaternion & rhs)
	{
		double scalar = lhs.Scalar() * rhs.Scalar() - lhs.Vector().Dot(rhs.Vector());
		JSL::Vector vec = lhs.Scalar() * rhs.Vector() + rhs.Scalar() * lhs.Vector() + lhs.Vector().Cross(rhs.Vector());
		return Quaternion(scalar,vec);
	}
	
	//! The (slightly dodgy) matrix-quaternion product. In reality, casts the quaternion to R^4, multiplies, then casts back. \param lhs A JSL::Matrix object to be multiplied \param rhs A Quaternion object to be multiplied \returns A Quaternion object of the resulting product
	inline Quaternion operator*(const JSL::Matrix & lhs, const Quaternion & rhs)
	{
		//~ std::vector<double> vTemp = ;
		JSL::Vector v({rhs.Scalar(),rhs.Vector(0),rhs.Vector(1),rhs.Vector(2)});
		
		return Quaternion(lhs * v);
	}
	
	//! Quaternion division operation, such that: \verbatim embed:rst:inline :math:`\mathsf{a} \oslash \mathsf{b}  = \frac{1}{|\mathsf{b}|^2} \mathsf{a} \otimes \overline{\mathsf{b}}` \endverbatim.  \param lhs The first argument of the operation (\verbatim embed:rst:inline :math:`\mathsf{a}` \endverbatim) \param rhs The second argument (\verbatim embed:rst:inline :math:`\mathsf{b}` \endverbatim) \return The quaternion division (\verbatim embed:rst:inline :math:`\mathsf{a}\oslash\mathsf{b}` \endverbatim)
	inline Quaternion operator/(const Quaternion & lhs, const Quaternion & rhs)
	{
		double N = rhs.Norm();
		if (N <= 0)
		{
			throw std::runtime_error("Quaternion division is not well defined when the norm is 0");
		}
		
		return 1.0/N * lhs * rhs.Conjugate();
	}
	
	//! The quaternion exponential, defined such that: 	 \verbatim embed:rst:inline :math:`\exp(\mathsf{a}) = \sum_{n = 0}^\infty \frac{\mathsf{a}^n}{n!}` \endverbatim \param lhs The argument of the operation (\verbatim embed:rst:inline :math:`\mathsf{a}` \endverbatim) \return The (analytically computed) quaternion exponential 
	inline Quaternion exp(const Quaternion & a)
	{
		double sc = std::exp(a.Scalar());
		double theta = a.Vector().Norm();
		double sinc = sin(theta)/theta;
		if (theta == 0)
		{
			sinc = 1;
		}
		Quaternion base(cos(theta), sinc * a.Vector());
		return sc * base;
	}
}

