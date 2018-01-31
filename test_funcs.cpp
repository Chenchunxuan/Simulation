#include "matrix.hpp"
#include "quaternions.hpp"
#include <iostream>
#include <string>
#include <cmath>
#include <math.h>
using namespace std;


// global variables
unsigned int success = 0, total = 0; 


// function declarations
Vec getVec(const double d[], unsigned int size);
void test_Matrix_Default_Constructor(void);
void test_Matrix_Zeros_Constructor(void);
void test_Matrix_Main_Constructor(void);
void test_Matrix_Copy_Constructor(void);
void test_Matrix_Get_Functions(void);
void test_Matrix_Element_Get_Operator(void);
void test_Matrix_Element_Set_Operator(void);
void test_Matrix_Replace_Function(void);
void test_Matrix_Assignment_Operator(void);
void test_Matrix_Double_Equal_Sign(void);
void test_Matrix_Addition_Operator(void);
void test_Matrix_Addition_Equal_Operator(void);
void test_Matrix_Subtraction_Operator(void);
void test_Matrix_Subtraction_Equal_Operator(void);
void test_Matrix_Multiplication_Operator(void);
void test_Matrix_Multiplication_Equal_Operator(void);
void test_Matrix_Power_Operator(void);
void test_Matrix_Max_Function(void);
void test_Matrix_Min_Function(void);
void test_Matrix_Slice_Function(void);
void test_Matrix_Identity_Matrix(void);
void test_Matrix_Skew_Matrix(void);
void test_Matrix_Transpose(void);
void test_Matrix_Max_Index(void);
void test_Quaternion_Default_Constructor_Test(void);
void test_Quaternion_Main_Constructor_Test(void);
void test_Quaternion_Alternate_Constructor(void);
void test_Quaternion_Copy_Constructor(void);
void test_Quaternion_Equal_Operator(void);
void test_Quaternion_Norm(void);
void test_Quaternion_Get_Function(void);
void test_Quaternion_Multiplication_Operator(void);
void test_Quaternion_Conjugate(void);
void test_Quaternion_Double_Equal_Operator(void);
void test_Quaternion_Purify_Function(void);
void test_Quaternion_Transform_Function(void);



int main( int argc, const char *argv[] ) {
    // execute tests
    string line = "##########################################################";
    cout << endl << line << endl;
    cout << "Running Tests..." << endl << endl;
    test_Matrix_Default_Constructor();
    test_Matrix_Zeros_Constructor();
    test_Matrix_Main_Constructor();
    test_Matrix_Copy_Constructor();
    test_Matrix_Get_Functions();
    test_Matrix_Element_Get_Operator();
    test_Matrix_Element_Set_Operator();
    test_Matrix_Replace_Function();
    test_Matrix_Assignment_Operator();
    test_Matrix_Double_Equal_Sign();
    test_Matrix_Addition_Operator();
    test_Matrix_Addition_Equal_Operator();
    test_Matrix_Subtraction_Operator();
    test_Matrix_Subtraction_Equal_Operator();
    test_Matrix_Multiplication_Operator();
    test_Matrix_Multiplication_Equal_Operator();
    test_Matrix_Power_Operator();
    test_Matrix_Max_Function();
    test_Matrix_Min_Function();
    test_Matrix_Slice_Function();
    test_Matrix_Identity_Matrix();
    test_Matrix_Skew_Matrix();
    test_Matrix_Transpose();
    test_Matrix_Max_Index();
    test_Quaternion_Default_Constructor_Test();
    test_Quaternion_Main_Constructor_Test();
    test_Quaternion_Alternate_Constructor();
    test_Quaternion_Copy_Constructor();
    test_Quaternion_Equal_Operator();
    test_Quaternion_Norm();
    test_Quaternion_Get_Function();
    test_Quaternion_Multiplication_Operator();
    test_Quaternion_Conjugate();
    test_Quaternion_Double_Equal_Operator();
    test_Quaternion_Purify_Function();
    test_Quaternion_Transform_Function();

    // tally up successful tests
    cout << endl;
    cout << "Success rate: " << success << "/" << total << endl << endl;
    cout << line << endl << endl;

    return 0;
}


// function definitions
Vec getVec(const double d[], unsigned int size) {
	Vec temp(size, 0);
	for (unsigned int i = 0; i < size; i++) {
		temp[i] = d[i];
	}

	return temp;
}

void test_Matrix_Default_Constructor(void) {
	total++;
	Matrix M;
	Mat m(1, Vec(1, 0));
	if (M.getRows() == 1 && M.getCols() == 1 && M.getMatrix() == m) {
		success++;
		cout << "Matrix Default Constructor Test: Success!" << endl;
	} else {
		cout << "Matrix Default Constructor Test: Failure!" << endl;
	}
}

void test_Matrix_Zeros_Constructor(void) {
	total++;
	unsigned int rows = 3, cols = 2;
	Matrix M(rows, cols);
	Mat m(rows, Vec(cols, 0));
	if (M.getRows() == rows && M.getCols() == cols && M.getMatrix() == m) {
		success++;
		cout << "Matrix Zeros Constructor Test: Success!" << endl;
	} else {
		cout << "Matrix Zeros Constructor Test: Failure!" << endl;
	}
}

void test_Matrix_Main_Constructor(void) {
	total++;
	double d[] = {1,2,3,4,5,6,7,8,9};
	Vec v = getVec(d, 9);
	Matrix M(3, 3, v);
	Mat m(3, Vec(3, 0));
	for (unsigned int i = 0; i < 3; i ++) {
		for (unsigned int j = 0; j < 3; j++) {
			m[i][j] = 3*i + j + 1;
		}
	}
	if (M.getRows() == 3 && M.getCols() == 3 && M.getMatrix() == m) {
		success++;
		cout << "Matrix Main Constructor Test: Success!" << endl;
	} else {
		cout << "Matrix Main Constructor Test: Failure!" << endl;
	}
}

void test_Matrix_Copy_Constructor(void) {
	total++;
	double d[] = {0, 1, 0, 1, 0, 1, 0, 1, 0};
	Vec v = getVec(d, 9);
	Matrix M(3, 3, v);
	Matrix M_test(M);
	if (M.getRows() == M_test.getRows() && M.getCols() == M_test.getCols() && M.getMatrix() == M_test.getMatrix()) {
		success++;
		cout << "Matrix Copy Constructor Test: Success!" << endl;
	} else {
		cout << "Matrix Copy Constructor Test: Failure!" << endl;
	}
}

void test_Matrix_Get_Functions(void) {
	total++;
	unsigned int rows = 2, cols = 3;
	double d[] = {3,3,3,3,3,3};
	Vec v = getVec(d, rows*cols);
	Mat m(rows, Vec(cols, 3)); // 2x3 matrix of three's
	Matrix M(rows, cols, v);
	if (M.getCols() == cols && M.getRows() == rows && M.getMatrix() == m) {
		success++;
		cout << "Matrix Get Function Tests: Success!" << endl;
	} else {
		cout << "Matrix Get Function Tests: Failure!" << endl;
	}
}

void test_Matrix_Element_Get_Operator(void) {
	total++;
	double d[] = {1,2,3,4,5,6,7,8,9};
	Vec v = getVec(d, 9);
	Matrix M(3, 3, v);
	if (M(0, 0) == 1 && M(2, 0) == 7 && M(2, 1) == 8) {
		success++;
		cout << "Matrix Get Operator Test: Success!" << endl;
	} else {
		cout << "Matrix Get Operator Test: Failure!" << endl;
	}
}

void test_Matrix_Element_Set_Operator(void) {
	total++;
	double d[] = {1,2,3,4,5,6,7,8,9};
	Vec v = getVec(d, 9);
	Matrix M(3, 3, v);
	M(0, 0) = 10;
	M(2, 0) = 11;
	M(2, 2) = 12;
	if (M(0, 0) == 10 && M(2, 0) == 11 && M(2, 2) == 12) {
		success++;
		cout << "Matrix Set Operator Test: Success!" << endl;
	} else {
		cout << "Matrix Set Operator Test: Failure!" << endl;
	}
}

void test_Matrix_Replace_Function(void) {
	total++;
	double d[] = {1,2,3,4,5,6,7,8,9};
	double d_test[] = {0,0,0,0,5,0,0,0,0};
	Vec v = getVec(d, 9), v_test = getVec(d_test, 9);
	Vec_index v_index(2); 
	v_index[0] = 0; v_index[1] = 2;
	Matrix M(3, 3, v), M_test(3, 3, v_test);
	Matrix M1(2, 2), M2(2, 1), M3(1, 2);
	
	// replace corner elements with 0
	M.replace(v_index, v_index, M1);

	// replace odd elements of second row with 0
	M.replace(1, v_index, M3);

	// replace odd elements of second column with 0
	M.replace(v_index, 1, M2);

	if (M.getRows() == M_test.getRows() && M.getCols() == M_test.getCols() && M.getMatrix() == M_test.getMatrix()) {
		success++;
		cout << "Matrix Replace Function Test: Success!" << endl;
	} else {
		cout << "Matrix Replace Function Test: Failure!" << endl;
	}
}

void test_Matrix_Assignment_Operator(void) {
	total++;
	double d[] = {1,2,3,4,5,6,7,8,9};
	Vec v = getVec(d, 9);
	Matrix M(3, 3, v);
	Matrix M_test1, M_test2; 

	// test assignment operator chaining
	M_test2 = (M_test1 = M);

	if (M.getRows() == M_test1.getRows() && M.getCols() == M_test1.getCols() && M.getMatrix() == M_test1.getMatrix() &&
		M.getRows() == M_test2.getRows() && M.getCols() == M_test2.getCols() && M.getMatrix() == M_test2.getMatrix()) {
		success++;
		cout << "Matrix Assignment Operator Test: Success!" << endl;
	} else {
		cout << "Matrix Assignment Operator Test: Failure!" << endl;
	}
}

void test_Matrix_Double_Equal_Sign(void) {
	total++;
	double d[] = {1,2,3,4,5,6,7,8,9};
	Vec v = getVec(d, 9);
	Matrix M(3, 3, v), M_id(3, 3, v); // create two separate matrices of same value
	Matrix M_zeros(1, 2); 

	bool tf1 = (M == M_id);
	bool tf2 = (M == M_zeros);

	if (tf1 && !tf2) {
		success++;
		cout << "Matrix Double Equal Sign Test: Success!" << endl;
	} else {
		cout << "Matrix Double Equal Sign Test: Failure!" << endl;
	}
}

void test_Matrix_Addition_Operator(void) {
	total++;
	double d[] = {1,2,3,4,5,6,7,8,9};
	Vec v = getVec(d, 9);
	Matrix M(3, 3, v);
	Matrix M_out = (M + M) + 2;
	double d_out[] = {4, 6, 8, 10, 12, 14, 16, 18, 20};
	Vec v_out = getVec(d_out, 9);
	Matrix M_out_test(3, 3, v_out);

	if (M_out.getRows() == M_out_test.getRows() && M_out.getCols() == M_out_test.getCols() && M_out.getMatrix() == M_out_test.getMatrix()) {
		success++;
		cout << "Matrix Addition Operator Test: Success!" << endl;
	} else {
		cout << "Matrix Addition Operator Test: Failure!" << endl;
	}
}

void test_Matrix_Addition_Equal_Operator(void) {
	total++;
	double d[] = {1,2,3,4,5,6,7,8,9};
	Vec v = getVec(d, 9);
	Matrix M(3, 3, v);
	(M += M) += 2; // chaining
	double d_out[] = {4, 6, 8, 10, 12, 14, 16, 18, 20};
	Vec v_out = getVec(d_out, 9);
	Matrix M_out_test(3, 3, v_out);

	if (M.getRows() == M_out_test.getRows() && M.getCols() == M_out_test.getCols() && M.getMatrix() == M_out_test.getMatrix()) {
		success++;
		cout << "Matrix Addition Equal Operator Test: Success!" << endl;
	} else {
		cout << "Matrix Addition Equal Operator Test: Failure!" << endl;
	}
}

void test_Matrix_Subtraction_Operator(void) {
	total++;
	double d[] = {1,2,3,4,5,6,7,8,9};
	Vec v = getVec(d, 9);
	Matrix M(3, 3, v);
	Matrix M_out = (M - M) - 2;
	double d_out[] = {-2, -2, -2, -2, -2, -2, -2, -2, -2};
	Vec v_out = getVec(d_out, 9);
	Matrix M_out_test(3, 3, v_out);

	if (M_out.getRows() == M_out_test.getRows() && M_out.getCols() == M_out_test.getCols() && M_out.getMatrix() == M_out_test.getMatrix()) {
		success++;
		cout << "Matrix Subtraction Operator Test: Success!" << endl;
	} else {
		cout << "Matrix Subtraction Operator Test: Failure!" << endl;
	}
}

void test_Matrix_Subtraction_Equal_Operator(void) {
	total++;
	double d[] = {1,2,3,4,5,6,7,8,9};
	Vec v = getVec(d, 9);
	Matrix M(3, 3, v);
	(M -= M) -= 2; // chaining
	double d_out[] = {-2, -2, -2, -2, -2, -2, -2, -2, -2};
	Vec v_out = getVec(d_out, 9);
	Matrix M_out_test(3, 3, v_out);

	if (M.getRows() == M_out_test.getRows() && M.getCols() == M_out_test.getCols() && M.getMatrix() == M_out_test.getMatrix()) {
		success++;
		cout << "Matrix Subtraction Equal Operator Test: Success!" << endl;
	} else {
		cout << "Matrix Subtraction Equal Operator Test: Failure!" << endl;
	}
}

void test_Matrix_Multiplication_Operator(void) {
	total++;
	double d[] = {1,2,3,4,5,6,7,8,9};
	Vec v = getVec(d, 9);
	Matrix M(3, 3, v);
	Matrix M_out = (M * M) * 2;
	double d_out[] = {60, 72, 84, 132, 162, 192, 204, 252, 300};
	Vec v_out = getVec(d_out, 9);
	Matrix M_out_test(3, 3, v_out);

	if (M_out.getRows() == M_out_test.getRows() && M_out.getCols() == M_out_test.getCols() && M_out.getMatrix() == M_out_test.getMatrix()) {
		success++;
		cout << "Matrix Multiplication Operator Test: Success!" << endl;
	} else {
		cout << "Matrix Multiplication Operator Test: Failure!" << endl;
	}
}

void test_Matrix_Multiplication_Equal_Operator(void) {
	total++;
	double d[] = {1,2,3,4,5,6,7,8,9};
	Vec v = getVec(d, 9);
	Matrix M(3, 3, v);
	(M *= M) *= 2; // chaining
	double d_out[] = {60, 72, 84, 132, 162, 192, 204, 252, 300};
	Vec v_out = getVec(d_out, 9);
	Matrix M_out_test(3, 3, v_out);

	if (M.getRows() == M_out_test.getRows() && M.getCols() == M_out_test.getCols() && M.getMatrix() == M_out_test.getMatrix()) {
		success++;
		cout << "Matrix Multiplication Equal Operator Test: Success!" << endl;
	} else {
		cout << "Matrix Multiplication Equal Operator Test: Failure!" << endl;
	}
}

void test_Matrix_Power_Operator(void) {
	total++;
	double d[] = {1,2,3,4,5,6,7,8,9};
	Vec v = getVec(d, 9);
	Matrix M(3, 3, v);
	Matrix M_out = M^3;
	double d_out[] = {468, 576, 684, 1062, 1305, 1548, 1656, 2034, 2412};
	Vec v_out = getVec(d_out, 9);
	Matrix M_out_test(3, 3, v_out);

	if (M_out.getRows() == M_out_test.getRows() && M_out.getCols() == M_out_test.getCols() && M_out.getMatrix() == M_out_test.getMatrix()) {
		success++;
		cout << "Matrix Power Operator Test: Success!" << endl;
	} else {
		cout << "Matrix Power Operator Test: Failure!" << endl;
	}
}

void test_Matrix_Max_Function(void) {
	total++;
	double d[] = {1,2,3,4,5,6,7,8,9};
	Vec v = getVec(d, 9);
	Matrix M(3, 3, v);
	double maximum = M.max();

	if (maximum == 9) {
		success++;
		cout << "Matrix Max Function Test: Success!" << endl;
	} else {
		cout << "Matrix Max Function Test: Failure!" << endl;
	}
}

void test_Matrix_Min_Function(void) {
	total++;
	double d[] = {1,2,3,4,5,6,7,8,9};
	Vec v = getVec(d, 9);
	Matrix M(3, 3, v);
	double minimum = M.min();

	if (minimum == 1) {
		success++;
		cout << "Matrix Min Function Test: Success!" << endl;
	} else {
		cout << "Matrix Min Function Test: Failure!" << endl;
	}
}

void test_Matrix_Slice_Function(void) {
	total++;
	double d[] = {1,2,3,4,5,6,7,8,9}, d1[] = {1,2}, d2[] = {1,4}, d3[] = {1,2,4,5};
	Vec v = getVec(d, 9), v1 = getVec(d1, 2), v2 = getVec(d2, 2), v3 = getVec(d3, 4);
	Matrix M(3, 3, v), subM1(1, 2, v1), subM2(2, 1, v2), subM3(2, 2, v3);
	Vec_index v_ind(2); v_ind[0] = 0; v_ind[1] = 1;

	if (M.slice(0, v_ind) == subM1 && M.slice(v_ind, 0) == subM2 && M.slice(v_ind, v_ind) == subM3) {
		success++;
		cout << "Matrix Slice Function Test: Success!" << endl;
	} else {
		cout << "Matrix Slice Function Test: Failure!" << endl;
	}
}

void test_Matrix_Identity_Matrix(void) {
	total++;
	double d[] = {1,0,0,0,1,0,0,0,1};
	Vec v = getVec(d, 9);
	Matrix M(3, 3, v);

	if (M == eye(3)) {
		success++;
		cout << "Matrix Identity Matrix Test: Success!" << endl;
	} else {
		cout << "Matrix Identity Matrix Test: Failure!" << endl;
	}
}

void test_Matrix_Skew_Matrix(void) {
	total++;
	double v_arr[] = {-1, 2, -4}, s_arr[] = {0, 4, 2, -4, 0, 1, -2, -1, 0};
	Vec v_vec = getVec(v_arr, 3), s_vec = getVec(s_arr, 9);
	Matrix v(3, 1, v_vec), S(3, 3, s_vec);

	if (skew(v) == S) {
		success++;
		cout << "Matrix Skew Matrix Test: Success!" << endl;
	} else {
		cout << "Matrix Skew Matrix Test: Failure!" << endl;
	}
}

void test_Matrix_Transpose(void) {
	total++;
	double v[] = {1, 2, 3, 4, 5, 6}, vt[] = {1, 3, 5, 2, 4, 6};
	Vec v_vec = getVec(v, 6), vt_vec = getVec(vt, 6);
	Matrix v_mat(3, 2, v_vec), vt_mat(2, 3, vt_vec);

	if (transpose(v_mat) == vt_mat) {
		success++;
		cout << "Matrix Transpose Test: Success!" << endl;
	} else {
		cout << "Matrix Transpose Test: Failure!" << endl;
	}
}

void test_Matrix_Max_Index(void) {
	total++;
	Vec_index v_ind(9);
	for (unsigned int i = 0; i < 9; i++) {
		v_ind[i] = i;
	}

	if (max_index(v_ind) == 8) {
		success++;
		cout << "Matrix Max Index Test: Success!" << endl;
	} else {
		cout << "Matrix Max Index Test: Failure!" << endl;
	}
}

void test_Quaternion_Default_Constructor_Test(void) {
	total++;
	Quaternion q;
	Matrix q_m(4, 1);
	if (q.getQuat() == q_m) {
		success++;
		cout << "Quaternion Default Constructor Test: Success!" << endl;
	} else {
		cout << "Quaternion Default Constructor Test: Failure!" << endl;
	}
}

void test_Quaternion_Main_Constructor_Test(void) {
	total++;
	double d[] = {1, -2, 3}, d_test[] = {0.102276449393203, -0.204552898786406, 0.306829348179609,  0.923879532511287};
	Vec v = getVec(d, 3), v_test = getVec(d_test, 4);
	Matrix v_other(3, 1, v), m_test(4, 1, v_test);
	double theta = M_PI/4;
	Quaternion q_vec(theta, v), q_matrix(theta, v_other);
	Matrix q_vec_m = q_vec.getQuat(), q_matrix_m = q_matrix.getQuat();
	double dec = 100000.;

	for (unsigned int i = 0; i < 4; i++) {
		q_vec_m(i, 0) = trunc(q_vec_m(i, 0) * dec);
		q_matrix_m(i, 0) = trunc(q_matrix_m(i, 0) * dec);
		m_test(i, 0) = trunc(m_test(i, 0) * dec);
	}

	if (q_vec_m == m_test && q_matrix_m == m_test) {
		success++;
		cout << "Quaternion Main Constructor Test: Success!" << endl;
	} else {
		cout << "Quaternion Main Constructor Test: Failure!" << endl;
	}
}

void test_Quaternion_Alternate_Constructor(void) {
	total++;
	double d[] = {4, -4, 4, -4}, dtest[] = {0.5, -0.5, 0.5, -0.5};
	Vec v = getVec(d, 4), vtest = getVec(dtest, 4);
	Matrix m(4, 1, v), mtest(4, 1, vtest);
	Quaternion q(m);

	if (q.getQuat() == mtest) {
		success++;
		cout << "Quaternion Alternate Constructor Test: Success!" << endl;
	} else {
		cout << "Quaternion Alternate Constructor Test: Failure!" << endl;
	}
}

void test_Quaternion_Copy_Constructor(void) {
	total++;
	double d[] = {1, -2, 3};
	Vec v = getVec(d, 3);
	Matrix v_other(3, 1, v);
	double theta = M_PI/4;
	Quaternion q_vec(theta, v);
	Quaternion q_copy(q_vec);

	if (q_vec.getQuat() == q_copy.getQuat()) {
		success++;
		cout << "Quaternion Copy Constructor Test: Success!" << endl;
	} else {
		cout << "Quaternion Copy Constructor Test: Failure!" << endl;
	}
}

void test_Quaternion_Equal_Operator(void) {
	total++;
	double d[] = {1, -2, 3};
	Vec v = getVec(d, 3);
	Matrix v_other(3, 1, v);
	double theta = M_PI/4;
	Quaternion q1(theta, v), q2, q3;

	// chained assignment
	q3 = (q2 = q1);

	if (q1.getQuat() == q2.getQuat() && q1.getQuat() == q3.getQuat()) {
		success++;
		cout << "Quaternion Assignment Operator Test: Success!" << endl;
	} else {
		cout << "Quaternion Assignment Operator Test: Failure!" << endl;
	}
}

void test_Quaternion_Norm(void) {
	total++;
	double d[] = {1, -4, 3, -5};
	Vec v = getVec(d, 4);
	Matrix m(4, 1, v);
	Quaternion q(m);
	double n = 1.0;

	if ((q.norm() - n) < numeric_limits<double>::epsilon()) {
		success++;
		cout << "Quaternion Norm Test: Success!" << endl;
	} else {
		cout << "Quaternion Norm Test: Failure!" << endl;
	}
}

void test_Quaternion_Get_Function(void) {
	total++;
	double d[] = {0.5, -0.5, 0.5, -0.5};
	Vec v = getVec(d, 4);
	Matrix m(4, 1, v);
	Quaternion q(m);

	if (q.getQuat() == m) {
		success++;
		cout << "Quaternion Get Function Test: Success!" << endl;
	} else {
		cout << "Quaternion Get Function Test: Failure!" << endl;
	}
}

void test_Quaternion_Multiplication_Operator(void) {
	total++;
	double d[] = {1, -2, 3}, dd[] = {-3, 1, 2}, dtest[] = {-0.377472641562540, -0.204029331006271, 0.444302874178547, 0.786435879749655};
	Vec pv = getVec(d, 3), qv = getVec(dd, 3), vtest = getVec(dtest, 4);
	Matrix pm(3, 1, pv), qm(3, 1, qv), mtest(4, 1, vtest);
	double thetap = M_PI/4, thetaq = M_PI/3;
	Quaternion p(thetap, pv), q(thetaq, qv), p_times_q = p * q;
	Matrix p_times_q_matrix = p_times_q.getQuat();
	double dec = 100000000.;

	for (unsigned int i = 0; i < 4; i++) {
		p_times_q_matrix(i, 0) = trunc(p_times_q_matrix(i, 0) * dec);
		mtest(i, 0) = trunc(mtest(i, 0) * dec);
	}

	if (p_times_q_matrix == mtest) {
		success++;
		cout << "Quaternion Multiplication Operator Test: Success!" << endl;
	} else {
		cout << "Quaternion Multiplication Operator Test: Failure!" << endl;
	}

}

void test_Quaternion_Conjugate(void) {
	total++;
	double d[] = {1, -2, 3}, dtest[] = {-0.102276449393203, 0.204552898786406, -0.306829348179609, 0.923879532511287};
	Vec v = getVec(d, 3), vtest = getVec(dtest, 4);
	Matrix v_other(3, 1, v), mtest(4, 1, vtest);
	double theta = M_PI/4;
	Quaternion q1(theta, v), q1_conj = q1.conjugate();
	Matrix q1_conj_matrix = q1_conj.getQuat();
	double dec = 100000000.;

	for (unsigned int i = 0; i < 4; i++) {
		q1_conj_matrix(i, 0) = trunc(q1_conj_matrix(i, 0) * dec);
		mtest(i, 0) = trunc(mtest(i, 0) * dec);
	}

	if (q1_conj_matrix == mtest) {
		success++;
		cout << "Quaternion Conjugate Test: Success!" << endl;
	} else {
		cout << "Quaternion Conjugate Test: Failure!" << endl;
	}
}

void test_Quaternion_Double_Equal_Operator(void) {
	total++;
	double d[] = {1, -2, 3}, dfake[] = {1, 1, 2};
	Vec v = getVec(d, 3), vfake = getVec(dfake, 3);
	Matrix v_other(3, 1, v), mfake(3, 1, vfake);
	double theta = M_PI/4;
	Quaternion q(theta, v_other), qfake(theta, mfake), q_copy(q);

	if (q == q_copy && !(q == qfake)) {
		success++;
		cout << "Quaternion Double Equal Operator Test: Success!" << endl;
	} else {
		cout << "Quaternion Double Equal Operator Test: Failure!" << endl;
	}
}

void test_Quaternion_Purify_Function(void) {
	total++;
	double d[] = {1, -2, 3}, dtest[] = {0.267261241912424, -0.534522483824849, 0.801783725737273, 0.0};
	Vec v = getVec(d, 3), vtest = getVec(dtest, 4);
	Matrix v_other(3, 1, v), mtest(4, 1, vtest);
	Quaternion q = purify(v_other);
	Matrix q_matrix = q.getQuat();
	double dec = 100000000.;

	for (unsigned int i = 0; i < 4; i++) {
		q_matrix(i, 0) = trunc(q_matrix(i, 0) * dec);
		mtest(i, 0) = trunc(mtest(i, 0) * dec);
	}

	if (q_matrix == mtest) {
		success++;
		cout << "Quaternion Purify Function Test: Success!" << endl;
	} else {
		cout << "Quaternion Purify Function Test: Failure!" << endl;
	}
}

void test_Quaternion_Transform_Function(void) {
	total++;
	double d[] = {1, -2, 3};
	double dbody[] = {1, 0, 0}; // vector in body frame
	double dtest[] = {0.535714285714286, 0.622936503400842, 0.570052907029133}; // transformed vector in inertial frame

	Vec v = getVec(d, 3), vbody = getVec(dbody, 3), vtest = getVec(dtest, 3);
	Matrix v_other(3, 1, v), mbody(3, 1, vbody), mtest(3, 1, vtest);
	double theta = M_PI/3;
	Quaternion q(theta, v_other);
	Matrix vtrans = transform(q, mbody);
	double dec = 100000000.;

	for (unsigned int i = 0; i < 3; i++) {
		mtest(i, 0) = trunc(mtest(i, 0) * dec);
		vtrans(i, 0) = trunc(vtrans(i, 0) * dec);
	}

	if (mtest == vtrans) {
		success++;
		cout << "Quaternion Transform Test: Success!" << endl;
	} else {
		cout << "Quaternion Transform Test: Failure!" << endl;
	}
}
