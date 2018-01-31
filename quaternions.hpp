#ifndef QUATERNIONS_H
#define QUATERNIONS_H
/*
This file defines common operations on quaternions 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ~~~~~~~~~~~~~~~~~~~~~~~  Quaternion Convention  ~~~~~~~~~~~~~~~~~~~~~~~ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% - spin axis (3x1):     e     = [qx; qy; qz] (must be unit vector)       %
% - rotation angle:      theta = angle of rotation from Inertial to Body  %
%                                Frame.                                   %
%                                                                         %
% - vector first, scalar last:                                            %
%   - q = [v; q4]                                                         %
%   - v = e.*sin(theta/2) = [q1; q2; q3]                                  %
%   - q4 = cos(theta/2)                                                   %
%                                                                         %
% - conjugate quaternion (inverse = conjugate)                            %
%   - q^-1 = [-v; q4]                                                     %
%   - (p x q)^-1 = q^-1 x p^-1                                            %
%                                                                         %
% - standard transformation                                               %
%   - [vp; 0] = q x [v; 0] x q^-1                                         %
%                                                                         %
% - left-handed convention 					                              %
%   - i^2 = j^2 = k^2 = -1                                                %
%   - i = -jk; j = ik; k = -ij                                            %
%   - quaternion operator transforms a vector from Inertial to Body Frame %
%     coordinates                                                         %
%                                                                         %
% - right-handed convention (currently using)                             %
%   - i^2 = j^2 = k^2 = -1                                                % 
%   - i = jk; j = -ik; k = ij                                             %
%   - quaternion operator transforms a vector from Body to Inertial Frame %
%     coordinates                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Author: Michael Wang
*/

#include "matrix.hpp"
#include <string>
#include <iostream>
#include <math.h>
using namespace std;


// Quaternion class declaration
class Quaternion {
	public:
		Quaternion(void); // default constructor: sets entire quaternion to 0
		Quaternion(double th, const Vec& axis); // main constructor 
		Quaternion(double th, const Matrix& axis); 
		Quaternion(const Matrix& q); // alternate constructor
		Quaternion(const Quaternion& q); // copy constructor 
		Quaternion& operator=(const Quaternion& q); // overload assignment operator for chained operations 
		void print(string quatName);
		double norm(void);
		Matrix getQuat(void); 
		Matrix getQuat(void) const; 
		Quaternion operator*(const Quaternion& q);
		Quaternion conjugate(void);
		bool operator==(const Quaternion& q);

	private: 
		Matrix quat; // quaternion representation
};

// external function declarations
Quaternion purify(const Matrix& v); // construct pure quaternion from 3-element column vector
Matrix transform(const Quaternion& q, const Matrix& v); // transform vector



// function definitions 
Quaternion::Quaternion(void) {
	Matrix q_temp(4, 1); // creates 4-element column zero vector
	quat = q_temp;
}

Quaternion::Quaternion(double th, const Vec& axis) {
	// check dimension of axis
	if (axis.size() != 3) {
		cout << "Axis must be a 3-element vector!" << endl;
		exit(0);
	}

	double temp = 0.0, mag = 0.0;
	for (unsigned int i = 0; i < axis.size(); i++) {
		temp += pow(axis[i], 2);
	}
	mag = sqrt(temp);

	// resize quat
	Vec vec4(4);
	Matrix q_temp(4, 1, vec4);

	// normalize spin axis and assemble quaternion
	for (unsigned int i = 0; i < axis.size(); i++) {
		q_temp(i, 0) = (axis[i] / mag) * sin(th/2);
	}

	// last element of quaternion
	q_temp(3, 0) = cos(th/2);
	quat = q_temp;
}

Quaternion::Quaternion(double th, const Matrix& axis) {
	// check dimension of axis
	if (axis.getRows() != 3 || axis.getCols() != 1) {
		cout << "Axis must be a 3-element vector!" << endl;
		exit(0);
	}

	double temp = 0.0, mag = 0.0;
	for (unsigned int i = 0; i < axis.getRows(); i++) {
		temp += pow(axis(i, 0), 2);
	}
	mag = sqrt(temp);

	// resize quat
	Vec vec4(4);
	Matrix q_temp(4, 1, vec4);

	// normalize spin axis and assemble quaternion
	for (unsigned int i = 0; i < axis.getRows(); i++) {
		q_temp(i, 0) = (axis(i, 0) / mag) * sin(th/2);
	}

	// last element of quaternion
	q_temp(3, 0) = cos(th/2);
	quat = q_temp;
}

Quaternion::Quaternion(const Matrix& q) {
	// check dimension of q
	if (q.getRows() != 4 || q.getCols() != 1) {
		cout << "Vector needs to be a 4-dimensional column vector" << endl;
		exit(0);
	}

	// normalize quaternion
	double temp = 0.0, norm = 0.0;
	for (unsigned int i = 0; i < 4; i++) {
		temp += pow(q(i, 0), 2);
	}
	norm = sqrt(temp);
	quat = q;
	quat *= (1/norm);
}

Quaternion::Quaternion(const Quaternion& q) {
	quat = q.quat;
}

Quaternion& Quaternion::operator=(const Quaternion& q) {
	quat = q.quat;

	return *this;
}

void Quaternion::print(string quatName) {
	cout << endl;
	cout << quatName << " = [";
	for (unsigned int i = 0; i < 4; i++) {
		cout << quat(i, 0);
		if (i == 3) {
			cout << "]" << endl;
		} else {
			cout << endl;
			cout << setw(quatName.size() + 5);
		}
	}
    cout << endl;

    return;
}

double Quaternion::norm(void) {
	double temp = 0.0, norm = 0.0;

	for (unsigned int i = 0; i < 4; i++) {
		temp += pow(quat(i, 0), 2);
	}
	norm = sqrt(temp);

	return norm;
}

Matrix Quaternion::getQuat(void) {
	return quat;
}

Matrix Quaternion::getQuat(void) const {
	return quat;
}

// out = p * q, using right-handed convention
Quaternion Quaternion::operator*(const Quaternion& q) {
	Vec_index vec123(3);
	for (unsigned int i = 0; i < vec123.size(); i++) {
		vec123[i] = i;
	}

	Matrix q_m = q.getQuat();
	double p4 = quat(3, 0), q4 = q_m(3, 0), out4;
	Matrix pvec = quat.slice(vec123, 0), qvec = q_m.slice(vec123, 0), outvec;
	Matrix pvec_t = transpose(pvec);

	outvec = (qvec*p4) + (pvec*q4) + (skew(pvec)*qvec); 
	Matrix temp = pvec_t*qvec;
	out4 = (p4*q4) - temp(0, 0);

	Matrix out_m(4, 1);
	out_m.replace(vec123, 0, outvec);
	out_m(3, 0) = out4;
	Quaternion out(out_m);

	return out;
}

Quaternion purify(const Matrix& v) {
	// check to see if vector is three-dimensional column vector
	if (v.getRows() != 3 || v.getCols() != 1) {
		cout << "Needs to be a three-dimensional column vector! (purify)" << endl;
		exit(0);
	}

	Vec_index vec123(3);
	for (unsigned int i = 0; i < vec123.size(); i++) {
		vec123[i] = i;
	}

	Matrix q_m(4, 1);
	q_m.replace(vec123, 0, v);
	Quaternion q(q_m);

	return q;
}

Quaternion Quaternion::conjugate(void) {
	Vec_index vec123(3);
	for (unsigned int i = 0; i < vec123.size(); i++) {
		vec123[i] = i;
	}

	Matrix q_conj(4, 1);
	for (unsigned int i = 0; i < 3; i++) {
		q_conj(i, 0) = -quat(i, 0);
	}
	q_conj(3, 0) = quat(3, 0);
	Quaternion q(q_conj);

	return q;
}

Matrix transform(const Quaternion& q, const Matrix& v) {
	// check dimension of v
	if (v.getRows() != 3 || v.getCols() != 1) {
		cout << "Input must be a three-dimensional column vector" << endl;
		exit(0);
	}

	Vec_index vec123(3);
	for (unsigned int i = 0; i < vec123.size(); i++) {
		vec123[i] = i;
	}

	Quaternion v_pure = purify(v), q_cur = q;
	Quaternion v_pure_t = (q_cur * v_pure) * q_cur.conjugate();
	Matrix v_t = (v_pure_t.getQuat()).slice(vec123, 0);

	return v_t;
}

bool Quaternion::operator==(const Quaternion& q) {
	return (quat == q.getQuat());
}

#endif // #ifndef QUATERNIONS_H