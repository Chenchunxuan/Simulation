#include "main.hpp"


int main( int argc, const char *argv[] ) {
	Vec temp1(3), temp2(3);
	temp1[0] = 1; temp1[1] = 0; temp1[2] = 0;
	temp2[0] = 0; temp2[1] = 0; temp2[2] = 1;
	Matrix v(3, 1, temp1);
	double theta = M_PI/4;
	Quaternion q(theta, temp2);
	Matrix vp = transform(q, v);

	v.print("v");
	vp.print("v_t");
	

    return 0;
}


