/*
This file defines a Matrix class for common matrix operations

Author: Michael Wang
*/

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <string>
#include <iomanip>
#include <limits>
using namespace std;


// definitions for vectors and matrices
typedef vector<double> Vec;
typedef vector<unsigned int> Vec_index;
typedef vector<Vec> Mat;


// Matrix class declaration
class Matrix {
	public:
		Matrix(void); // default constructor 
		Matrix(unsigned int nR, unsigned int nC); // constructor: initialize matrix with zeros
		Matrix(unsigned int nR, unsigned int nC, const Vec& arr); // constructor: fills array row-wise
		Matrix(const Matrix& obj); // copy constructor 
		unsigned int getRows(void);
		unsigned int getCols(void); 
		unsigned int getRows(void) const;
		unsigned int getCols(void) const; 
		Mat getMatrix(void); 
		void print(string matrixName); 
		double& operator()(unsigned int r, unsigned int c); 
		double operator()(unsigned int r, unsigned int c) const;
		void replace(unsigned int r, const Vec_index& col, const Matrix& M);
		void replace(const Vec_index& row, unsigned int c, const Matrix& M);
		void replace(const Vec_index& row, const Vec_index& col, const Matrix& M);
		Matrix& operator=(const Matrix& B); // overload assignment operator  
		bool operator==(const Matrix& B); 
		Matrix operator+(const Matrix& B); 
		Matrix operator+(double s); 
		Matrix& operator+=(const Matrix& B);
		Matrix& operator+=(double s); 
		Matrix operator-(const Matrix& B); 
		Matrix operator-(double s); 
		Matrix& operator-=(const Matrix& B);
		Matrix& operator-=(double s); 
		Matrix operator*(const Matrix& B); 
		Matrix operator*(double s); 
		Matrix& operator*=(const Matrix& B);
		Matrix& operator*=(double s); 
		Matrix operator^(unsigned int s); 
		double max(void); 
		double min(void);
		Matrix slice(unsigned int r, const Vec_index& col);
		Matrix slice(const Vec_index& row, unsigned int c);
		Matrix slice(const Vec_index& row, const Vec_index& col); 

	private: 
		unsigned int numRows;
		unsigned int numCols;
		Mat matrix;
};
// additional function declarations
Matrix eye(unsigned int n);
Matrix skew(const Matrix& v);
Matrix transpose(const Matrix& B);
unsigned int max_index(const Vec_index& M);



// function definitions 
Matrix::Matrix(void) {
	numRows = 1;
	numCols = 1;
	Vec temp(1);
	matrix.resize(1, temp);
}

Matrix::Matrix(unsigned int nR, unsigned int nC) {
	numRows = nR;
	numCols = nC;
	Vec temp(nC);
	matrix.resize(nR, temp);
}

Matrix::Matrix(unsigned int nR, unsigned int nC, const Vec& arr) {
	// check if dimension of arr agrees
	if (nR*nC != arr.size()) {
		cout << "Dimension of arguemnt array does not agree!" << endl;
		exit(0);
	}
	numRows = nR;
	numCols = nC;
	Vec temp(nC);
	matrix.resize(nR, temp);
	for (unsigned int i = 0; i < nR; i++) {
		for (unsigned int j = 0; j < nC; j++) {
			matrix[i][j] = arr[i*nC + j];
		}
	}
}

Matrix::Matrix(const Matrix& obj) {
	numRows = obj.numRows;
	numCols = obj.numCols;
	matrix = obj.matrix;
}

unsigned int Matrix::getRows(void) {
	return numRows;
}

unsigned int Matrix::getCols(void) {
	return numCols;
}

unsigned int Matrix::getRows(void) const {
	return numRows;
}

unsigned int Matrix::getCols(void) const {
	return numCols;
}

Mat Matrix::getMatrix(void) {
	return matrix;
}

void Matrix::print(string matrixName) {
	cout << endl;
	cout << matrixName << " = [";
	for (unsigned int i = 0; i < numRows; i++) {
    	for (unsigned int j = 0; j < numCols; j++) {
    		cout << matrix[i][j] << " ";
    	}
    	if (i == (numRows-1)) {
    		cout << "]" << endl;
    	} else {
    		cout << endl;
    		cout << setw(matrixName.size() + 5);
    	}
    }
    cout << endl;

    return;
}

double& Matrix::operator()(unsigned int r, unsigned int c) {
	// check index out of bounds
	if (r > numRows || c > numCols) {
		cout << "Index out of bounds!" << endl;
		exit(0); 
	}

	return matrix[r][c];
}

double Matrix::operator()(unsigned int r, unsigned int c) const {
	// check index out of bounds
	if (r > numRows || c > numCols) {
		cout << "Index out of bounds!" << endl;
		exit(0); 
	}

	return matrix[r][c];
}

void Matrix::replace(unsigned int r, const Vec_index& col, const Matrix& M) {
	// check index dimensions
	unsigned int c_max = max_index(col);
	if (r > (numRows-1) || c_max > (numCols-1)) {
		cout << "Index out of bounds! (replace)" << endl;
		exit(0);
	}

	// check matrix dimensions
	if (M.getRows() != 1 || M.getCols() != col.size()) {
		cout << "Matrix dimensions must match! (replace)" << endl;
		exit(0);
	}

	for (unsigned int i = 0; i < col.size(); i++) {
		matrix[r][col[i]] = M(0, i);
	}

	return;
}

void Matrix::replace(const Vec_index& row, unsigned int c, const Matrix& M) {
	// check index dimensions
	unsigned int r_max = max_index(row);
	if (r_max > (numRows-1) || c > (numCols-1)) {
		cout << "Index out of bounds! (replace)" << endl;
		exit(0);
	}

	// check matrix dimensions
	if (M.getRows() != row.size() || M.getCols() != 1) {
		cout << "Matrix dimensions must match! (replace)" << endl;
		exit(0);
	}

	for (unsigned int i = 0; i < row.size(); i++) {
		matrix[row[i]][c] = M(i, 0);
	}

	return;
}

void Matrix::replace(const Vec_index& row, const Vec_index& col, const Matrix& M) {
	// check index dimensions
	unsigned int r_max = max_index(row), c_max = max_index(col);
	if (r_max > (numRows-1) || c_max > (numCols-1)) {
		cout << "Index out of bounds! (replace)" << endl;
		exit(0);
	}

	// check matrix dimensions
	if (M.getRows() != row.size() || M.getCols() != col.size()) {
		cout << "Matrix dimensions must match! (replace)" << endl;
		exit(0);
	}

	for (unsigned int i = 0; i < row.size(); i++) {
		for (unsigned int j = 0; j < col.size(); j++) {
			matrix[row[i]][col[j]] = M(i, j);
		}
	}

	return;
}

Matrix& Matrix::operator=(const Matrix& B) {
	numRows = B.numRows;
	numCols = B.numCols;
	matrix = B.matrix;

	return *this;
}

bool Matrix::operator==(const Matrix& B) {
	bool tf = (matrix == B.matrix);

	return tf;
}

Matrix Matrix::operator+(const Matrix& B) {
	// check matrix dimensions
	if (numRows != B.getRows() || numCols != B.getCols()) {
		cout << "Matrix dimensions do not match! (+)." << endl;
		exit(0);
	}

	Vec temp(numRows*numCols);
	Matrix C(numRows, numCols, temp);

	for (unsigned int i = 0; i < numRows; i++) {
		for (unsigned int j = 0; j < numCols; j++) {
			C(i, j) = matrix[i][j] + B(i, j);
		}
	}

	return C;
}

Matrix Matrix::operator+(double s) {
	Vec temp(numRows*numCols);
	Matrix C(numRows, numCols, temp);

	for (unsigned int i = 0; i < numRows; i++) {
		for (unsigned int j = 0; j < numCols; j++) {
			C(i, j) = matrix[i][j] + s;
		}
	}

	return C;
}

Matrix& Matrix::operator+=(const Matrix& B) {
	// check matrix dimensions
	if (numRows != B.getRows() || numCols != B.getCols()) {
		cout << "Matrix dimensions do not match! (+=)." << endl;
		exit(0);
	}

	double temp = 0.0;
	for (unsigned int i = 0; i < numRows; i++) {
		for (unsigned int j = 0; j < numCols; j++) {
			temp = matrix[i][j] + B(i, j);
			matrix[i][j] = temp;
		}
	}

	return *this;
}

Matrix& Matrix::operator+=(double s) {
	double temp = 0.0;
	for (unsigned int i = 0; i < numRows; i++) {
		for (unsigned int j = 0; j < numCols; j++) {
			temp = matrix[i][j] + s;
			matrix[i][j] = temp;
		}
	}

	return *this;
}

Matrix Matrix::operator-(const Matrix& B) {
	// check matrix dimensions
	if (numRows != B.getRows() || numCols != B.getCols()) {
		cout << "Matrix dimensions do not match! (-)." << endl;
		exit(0);
	}

	Vec temp(numRows*numCols);
	Matrix C(numRows, numCols, temp);

	for (unsigned int i = 0; i < numRows; i++) {
		for (unsigned int j = 0; j < numCols; j++) {
			C(i, j) = matrix[i][j] - B(i, j);
		}
	}

	return C;
}

Matrix Matrix::operator-(double s) {
	Vec temp(numRows*numCols);
	Matrix C(numRows, numCols, temp);

	for (unsigned int i = 0; i < numRows; i++) {
		for (unsigned int j = 0; j < numCols; j++) {
			C(i, j) = matrix[i][j] - s;
		}
	}

	return C;
}

Matrix& Matrix::operator-=(const Matrix& B) {
	// check matrix dimensions
	if (numRows != B.getRows() || numCols != B.getCols()) {
		cout << "Matrix dimensions do not match! (-=)." << endl;
		exit(0);
	}

	double temp = 0.0;
	for (unsigned int i = 0; i < numRows; i++) {
		for (unsigned int j = 0; j < numCols; j++) {
			temp = matrix[i][j] - B(i, j);
			matrix[i][j] = temp;
		}
	}

	return *this;
}

Matrix& Matrix::operator-=(double s) {
	double temp = 0.0;
	for (unsigned int i = 0; i < numRows; i++) {
		for (unsigned int j = 0; j < numCols; j++) {
			temp = matrix[i][j] - s;
			matrix[i][j] = temp;
		}
	}

	return *this;
}

Matrix Matrix::operator*(const Matrix& B) {
	// check matrix dimensions
	if (numCols != B.getRows()) {
		cout << "Matrix dimensions do not match! (*)." << endl;
		exit(0);
	}

	Vec temp(numRows*(B.getCols()));
	Matrix C(numRows, B.getCols(), temp);

	double tempNum = 0.0;
	for (unsigned int i = 0; i < numRows; i++) {
		for (unsigned int j = 0; j < B.getCols(); j++) {
			tempNum = 0.0; // re-initialize
			for (unsigned int n = 0; n < numCols; n++) {
				tempNum += matrix[i][n] * B(n, j);
			}
			C(i, j) = tempNum;
		}
	}

	return C;
}

Matrix Matrix::operator*(double s) {
	Vec temp(numRows*numCols);
	Matrix C(numRows, numCols, temp);

	for (unsigned int i = 0; i < numRows; i++) {
		for (unsigned int j = 0; j < numCols; j++) {
			C(i, j) = matrix[i][j] * s;
		}
	}

	return C;
}

Matrix& Matrix::operator*=(const Matrix& B) {
	// check matrix dimensions
	if (numRows != B.getRows() || numCols != B.getCols() || numRows != numCols) {
		cout << "Matrix dimensions do not match! (*=)." << endl;
		exit(0);
	}

	Matrix C = *this * B;
	for (unsigned int i = 0; i < numRows; i++) {
		for (unsigned int j = 0; j < numCols; j++) {
			matrix[i][j] = C(i, j);
		}
	}

	return *this;
}

Matrix& Matrix::operator*=(double s) {
	double temp = 0.0;
	for (unsigned int i = 0; i < numRows; i++) {
		for (unsigned int j = 0; j < numCols; j++) {
			temp = matrix[i][j] * s;
			matrix[i][j] = temp;
		}
	}

	return *this;
}

Matrix Matrix::operator^(unsigned int s) {
	// check matrix dimensions
	if (numRows != numCols) {
		cout << "Matrix dimensions do not match! (^)." << endl;
		exit(0);
	}
	// power of 0
	if (s == 0) {
		Vec temp(numRows * numRows);
		Matrix I(numRows, numRows, temp);

		for (unsigned int i = 0; i < numRows; i++) {
			I(i, i) = 1.0;
		}

		return I;
	}

	Matrix C(*this), temp(C);

	for (unsigned int i = 1; i < s; i++) {
		temp *= C;
	}

	return temp;
}

Matrix eye(unsigned int n) {
	// check n
	if (n == 0) {
		cout << "Invalid Matrix Dimension for Identity Matrix." << endl;
		exit(0);
	}

	Vec temp(n*n);
	Matrix C(n, n, temp);

	for (unsigned int i = 0; i < n; i++) {
		C(i, i) = 1;
	}

	return C;
}

Matrix skew(const Matrix& v) {
	// check of vector v is 3x1
	if (v.getRows() != 3) {
		cout << "Input vector must be a column vector and have dimension of 3!" << endl;
		exit(0);
	}

	Vec temp(9);
	Matrix C(3, 3, temp);

	C(0, 0) = 0;
	C(0, 1) = -v(2, 0);
	C(0, 2) = v(1, 0);
	C(1, 0) = v(2, 0);
	C(1, 1) = 0;
	C(1, 2) = -v(0, 0);
	C(2, 0) = -v(1, 0);
	C(2, 1) = v(0, 0);
	C(2, 2) = 0;

	return C;
}

Matrix transpose(const Matrix& B) {
	Vec temp(B.getRows() * B.getCols());
	Matrix C(B.getCols(), B.getRows(), temp);

	for (unsigned int i = 0; i < B.getCols(); i++) {
		for (unsigned int j = 0; j < B.getRows(); j++) {
			C(i, j) = B(j, i);
		}
	}

	return C;
}

double Matrix::max(void) {
	double cur_max = numeric_limits<double>::min(); 
	for (unsigned int i = 0; i < numRows; i++) {
		for (unsigned int j = 0; j < numCols; j++) {
			if (matrix[i][j] > cur_max) {
				cur_max = matrix[i][j];
			}
		}
	}

	return cur_max;
}

double Matrix::min(void) {
	double cur_min = numeric_limits<double>::max(); 
	for (unsigned int i = 0; i < numRows; i++) {
		for (unsigned int j = 0; j < numCols; j++) {
			if (matrix[i][j] < cur_min) {
				cur_min = matrix[i][j];
			}
		}
	}

	return cur_min;
}

unsigned int max_index(const Vec_index& M) {
	double cur_max = numeric_limits<unsigned int>::min(); 
	for (unsigned int i = 0; i < M.size(); i++) {
		if (M[i] > cur_max) {
			cur_max = M[i];
		}
	}

	return cur_max;
}

Matrix Matrix::slice(unsigned int r, const Vec_index& col) {
	// check index limits
	unsigned int c_max = max_index(col);
	if (r > (numRows-1) || c_max > (numCols-1)) {
		cout << "Index out of bounds! (slice)" << endl;
		exit(0);
	}

	Vec temp(col.size());
	Matrix C(1, col.size(), temp);
	for (unsigned int i = 0; i < col.size(); i++) {
		C(0, i) = matrix[r][col[i]];
	}

	return C;
}

Matrix Matrix::slice(const Vec_index& row, unsigned int c) {
	// check index limits
	unsigned int r_max = max_index(row);
	if (r_max > (numRows-1) || c > (numCols-1)) {
		cout << "Index out of bounds! (slice)" << endl;
		exit(0);
	}

	Vec temp(row.size());
	Matrix C(row.size(), 1, temp);
	for (unsigned int i = 0; i < row.size(); i++) {
		C(i, 0) = matrix[row[i]][c];
	}

	return C;
}

Matrix Matrix::slice(const Vec_index& row, const Vec_index& col) {
	// check index limits
	unsigned int r_max = max_index(row), c_max = max_index(col);
	if (r_max > (numRows-1) || c_max > (numCols-1)) {
		cout << "Index out of bounds! (slice)" << endl;
		exit(0);
	}

	Vec temp(row.size() * col.size());
	Matrix C(row.size(), col.size(), temp);
	for (unsigned int i = 0; i < row.size(); i++) {
		for (unsigned int j = 0; j < col.size(); j++) {
			C(i, j) = matrix[row[i]][col[j]];
		}
	}

	return C;
}

#endif // #ifndef MATRIX_H