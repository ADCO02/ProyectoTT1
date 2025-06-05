/**
 * @file test_orbit_determination.cpp
 * 
 * @brief Unit tests for satellite orbit determination functions and classes.
 */


#include "../include/matrix.hpp"
#include "../include/R_x.hpp"
#include "../include/R_y.hpp"
#include "../include/R_z.hpp"
#include "../include/AccelPointMass.hpp"
#include "../include/Cheb3D.hpp"
#include "../include/EccAnom.hpp"
#include "../include/Frac.hpp"
#include "../include/MeanObliquity.hpp"
#include "../include/Mjday.hpp"
#include "../include/Mjday_TDB.hpp"
#include "../include/Position.hpp"
#include "../include/sign_.hpp"
#include "../include/timediff.hpp"
#include "../include/AzElPa.hpp"
#include "../include/IERS.hpp"
#include "../include/Legendre.hpp"
#include "../include/NutAngles.hpp"
#include "../include/TimeUpdate.hpp"
#include "../include/global.hpp"
#include "../include/AccelHarmonic.hpp"
#include "../include/EqnEquinox.hpp"
#include "../include/LTC.hpp"
#include "../include/NutMatrix.hpp"
#include "../include/PoleMatrix.hpp"
#include "../include/PrecMatrix.hpp"
#include "../include/gmst.hpp"
#include "../include/JPL_Eph_DE430.hpp"
#include "../include/gast.hpp"
#include "../include/G_AccelHarmonic.hpp"
#include "../include/GHAMatrix.hpp"
#include "../include/MeasUpdate.hpp"
#include "../include/Accel.hpp"
#include "../include/VarEqn.hpp"
#include "../include/DEInteg.hpp"
#include "../include/Geodetic.hpp"
#include "../include/angl.hpp"
#include "../include/elements.hpp"
#include "../include/unit.hpp"
#include "../include/hgibbs.hpp"
#include "../include/gibbs.hpp"
#include "../include/anglesg.hpp"
#include <cstdio>
#include <cmath>

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

int m_equals(Matrix& A, Matrix& B, double p) {
	if (A.n_row != B.n_row || A.n_column != B.n_column)
		return 0;
	else
		for(int i = 1; i <= A.n_row; i++)
			for(int j = 1; j <= A.n_column; j++)
				if(fabs(A(i,j)-B(i,j)) > p) {
					printf("%2.20lf %2.20lf\n",A(i,j),B(i,j));
					return 0;
				}
	
	return 1;
}

int m_sum_01() {
    int f = 3;
    int c = 4;
	
	Matrix& A = zeros(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix& B = zeros(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix& C = zeros(f, c);
	C(1,1) = 2; C(1,2) =  2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = 8; C(2,2) = -3; C(2,3) = 1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = -2; C(3,3) = 0; C(3,4) = 7;
	
	Matrix& R = A + B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_sub_01() {
    int f = 3;
    int c = 4;
	
	Matrix& A = zeros(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix& B = zeros(f, c);
	B(1,1) = 2; B(1,2) = 0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix& C = zeros(f, c);
	C(1,1) = -2; C(1,2) = 2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = -6; C(2,2) = 1; C(2,3) = -1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = 4; C(3,3) = 0; C(3,4) = 3;
	
	Matrix& R = A - B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_prod_01() {
    int f = 3;
    int c = 4;
	
	Matrix& A = zeros(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix& B = zeros(c, f);
	B(1,1) = 2; B(1,2) = 7; B(1,3) = 0;
	B(2,1) = 0; B(2,2) = -2; B(2,3) = -3;
	B(3,1) = 0; B(3,2) = 1; B(3,3) = 0;
	B(4,1) = 0; B(4,2) = 0; B(4,3) = 2;
	
	Matrix& C = zeros(f, f);
	C(1,1) = 0; C(1,2) = 4; C(1,3) = -6;
	C(2,1) = 2; C(2,2) = 9; C(2,3) = 3;
	C(3,1) = 0; C(3,2) = -2; C(3,3) = 7;
	
	Matrix& R = A * B;

    A = zeros(2, 2);
	A(1,1) = 0.0; A(1,2) = 0.2;
	A(2,1) = 0.1; A(2,2) = -0.1;
	
	B = zeros(2, 2);
	B(1,1) = 2; B(1,2) = 7;
	B(2,1) = 0; B(2,2) = -2;
	
	C = zeros(2,2);
	C(1,1) = 0; C(1,2) = -0.4;
	C(2,1) = 0.2; C(2,2) = 0.9;
	
	R = A * B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_asign_01(){
	Matrix& A = zeros(2, 2);
    A(1,1) = 1; A(1,2) = 2;
    A(2,1) = 3; A(2,2) = 4;

    Matrix& B = zeros(2, 2);

    B = A;

	_assert(m_equals(A, B, 1e-10));

    A(1,1) = 99;

    _assert(B(1,1) == 1 && A(1,1) == 99);

    return 0;
}

int m_div_01(){
    Matrix& A = zeros(2,2);
	A(1,1) = 1; A(1,2) = 2;
	A(2,1) = 3; A(2,2) = 4;

    Matrix& B = zeros(2,2);
	B(1,1) = 5; B(1,2) = 6;
	B(2,1) = 7; B(2,2) = 8;

    Matrix& R = A / B;

    Matrix& C = zeros(2,2);
	C(1,1) = 3; C(1,2) = -2;
	C(2,1) = 2; C(2,2) = -1;

    _assert(m_equals(C, R, 1e-10));
    return 0;
}

int m_sum_02() {
    int f = 3;
    int c = 4;
	
	Matrix& A = zeros(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix& B = zeros(f, c);
	B(1,1) = 3; B(1,2) = 5; B(1,3) = 11; B(1,4) = 3;
	B(2,1) = 4; B(2,2) = 2; B(2,3) = 3; B(2,4) = 3;
	B(3,1) = 3; B(3,2) = 4; B(3,3) = 3; B(3,4) = 8;
	
	Matrix& R = A + 3;
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_sub_02() {
    int f = 3;
    int c = 4;
	
	Matrix& A = zeros(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix& B = zeros(f, c);
	B(1,1) = -2; B(1,2) = 0; B(1,3) = 6; B(1,4) = -2;
	B(2,1) = -1; B(2,2) = -3; B(2,3) = -2; B(2,4) = -2;
	B(3,1) = -2; B(3,2) = -1; B(3,3) = -2; B(3,4) = 3;
	
	Matrix& R = A - 2;
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_prod_02() {
    int f = 3;
    int c = 4;
	
	Matrix& A = zeros(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix& B = zeros(f, c);
	B(1,1) = 0; B(1,2) = 4; B(1,3) = 16; B(1,4) = 0;
	B(2,1) = 2; B(2,2) = -2; B(2,3) = 0; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = 2; B(3,3) = 0; B(3,4) = 10;
	
	Matrix& R = A * 2;
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_div_02() {
    int f = 3;
    int c = 4;
	
	Matrix& A = zeros(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix& B = zeros(f, c);
	B(1,1) = 0; B(1,2) = 4; B(1,3) = 16; B(1,4) = 0;
	B(2,1) = 2; B(2,2) = -2; B(2,3) = 0; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = 2; B(3,3) = 0; B(3,4) = 10;
	
	Matrix& R = B / 2;
    
    _assert(m_equals(A, R, 1e-10));
    
    return 0;
}

int m_transpose_01() {
    int f = 3;
    int c = 4;
    
    Matrix& A = zeros(f, c);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3; A(1,4) = 4;
    A(2,1) = 5; A(2,2) = 6; A(2,3) = 7; A(2,4) = 8;
    A(3,1) = 9; A(3,2) = 10; A(3,3) = 11; A(3,4) = 12;

    Matrix& T = zeros(c, f);
    T(1,1) = 1; T(1,2) = 5; T(1,3) = 9;
    T(2,1) = 2; T(2,2) = 6; T(2,3) = 10;
    T(3,1) = 3; T(3,2) = 7; T(3,3) = 11;
    T(4,1) = 4; T(4,2) = 8; T(4,3) = 12;

    Matrix& R = A.transpose();

    _assert(m_equals(T, R, 1e-10));

    return 0;
}

int m_extract_vector_01() {
    Matrix& A = zeros(5);
    A(1) = 1; A(2) = 2; A(3) = 3; A(4) = 4; A(5) = 5;

    Matrix& sub_vector = A.extract_vector(2,4);

    Matrix& expected = zeros(3);
    expected(1) = 2; expected(2) = 3; expected(3) = 4;

    _assert(m_equals(sub_vector, expected, 1e-10));
    return 0;
}

int m_union_vector_01(){
	Matrix& m1 = zeros(3);
    Matrix& m2 = zeros(4);

    for (int i = 1; i <= 3; i++) {
        m1(i) = i;
    }
    for (int i = 1; i <= 4; i++) {
        m2(i) = i + 3;
    }

    Matrix& result = m1.union_vector(m2);

	Matrix& expected = zeros(7);
	for (int i = 1; i <= 7; i++) {
        expected(i) = i;
    }
	
	_assert(m_equals(result, expected, 1e-10));
    return 0;
}

int m_extract_row_01() {
    Matrix& A = zeros(3, 3);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;
    A(3,1) = 7; A(3,2) = 8; A(3,3) = 9;

    Matrix& row = A.extract_row(2);

    Matrix& expected = zeros(3);
    expected(1) = 4; expected(2) = 5; expected(3) = 6;

    _assert(m_equals(row, expected, 1e-10));
    return 0;
}

int m_extract_column_01() {
    Matrix& A = zeros(3, 3);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;
    A(3,1) = 7; A(3,2) = 8; A(3,3) = 9;

    Matrix& col = A.extract_column(3);

    Matrix& expected = zeros(3);
    expected(1) = 3; expected(2) = 6; expected(3) = 9;

    _assert(m_equals(col, expected, 1e-10));
    return 0;
}

int m_assign_row_01() {
    Matrix& A = zeros(3, 3);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;
    A(3,1) = 7; A(3,2) = 8; A(3,3) = 9;

    Matrix& row = zeros(3);
    row(1) = -1; row(2) = -2; row(3) = -3;

    A.assign_row(2, row);

    Matrix& expected = zeros(3, 3);
    expected(1,1) = 1; expected(1,2) = 2;  expected(1,3) = 3;
    expected(2,1) = -1; expected(2,2) = -2; expected(2,3) = -3;
    expected(3,1) = 7; expected(3,2) = 8;  expected(3,3) = 9;

    _assert(m_equals(A, expected, 1e-10));
    return 0;
}

int m_assign_column_01() {
    Matrix& A = zeros(3, 3);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;
    A(3,1) = 7; A(3,2) = 8; A(3,3) = 9;

    Matrix& col = zeros(3);
    col(1) = -1; col(2) = -2; col(3) = -3;

    A.assign_column(2, col);

    Matrix& expected = zeros(3, 3);
    expected(1,1) = 1; expected(1,2) = -1; expected(1,3) = 3;
    expected(2,1) = 4; expected(2,2) = -2; expected(2,3) = 6;
    expected(3,1) = 7; expected(3,2) = -3; expected(3,3) = 9;

    _assert(m_equals(A, expected, 1e-10));
    return 0;
}

int m_zeros_01() {
    int f = 3;
    int c = 4;
	
	Matrix& A = zeros(f, c);
	A(1,1) = 0; A(1,2) = 0; A(1,3) = 0; A(1,4) = 0;
	A(2,1) = 0; A(2,2) = 0; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 0; A(3,4) = 0;
	
	Matrix& B = zeros(3, 4);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_zeros_02() {	
	Matrix& A = zeros(4);
	A(1,1) = 0; A(2) = 0; A(3) = 0; A(4) = 0;
	
	Matrix& B = zeros(4);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_eye_01() {		
	Matrix& A = zeros(3,3);
	A(1,1) = 1; A(1,2) = 0; A(1,3) = 0;
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 1;
	
	Matrix& B = eye(3);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_norm_01() {		
	Matrix& A = zeros(3);
	A(1) = 12; A(2) = 3; A(3) = 4;
	
	double result = norm(A);

	double expected_norm = 13;
    
    _assert(fabs(result - expected_norm) < 1e-10);
    
    return 0;
}

int m_dot_01(){
    Matrix& A = zeros(3);
	A(1) = 1; A(2) = 2; A(3) = 3;

    Matrix& B = zeros(3);
	B(1) = 4; B(2) = 5; B(3) = 6;
	
	_assert(fabs(dot(A,B)-32.0)<1e-10);
    return 0;
}

int m_cross_01(){
    Matrix& A = zeros(3);
	A(1) = 1; A(2) = 2; A(3) = 3;

    Matrix& B = zeros(3);
	B(1) = 4; B(2) = 5; B(3) = 6;

    Matrix& expected = zeros(3);
	expected(1) = -3; expected(2) = 6; expected(3) = -3;

    Matrix& R = cross(A,B);
	
	_assert(m_equals(expected, R, 1e-10));
    return 0;
}

int m_inv_01(){
    Matrix& A = zeros(3,3);
	A(1,1) = 2; A(1,2) = 1; A(1,3) = 1;
    A(2,1) = 1; A(2,2) = 1; A(2,3) = 1;
    A(3,1) = 1; A(3,2) = 1; A(3,3) = 2;

    Matrix& R = inv(A);

    Matrix& B = zeros(3,3);
	B(1,1) = 1; B(1,2) = -1; B(1,3) = 0;
    B(2,1) = -1; B(2,2) = 3; B(2,3) = -1;
    B(3,1) = 0; B(3,2) = -1; B(3,3) = 1;
    
    _assert(m_equals(B, R, 1e-10));
    return 0;
}


int i1_R_x_01(){
    Matrix& R = R_x(3);
    
    Matrix& A = zeros(3,3);
    A(1,1) = 1; A(1,2) = 0; A(1,3) = 0;
    A(2,1) = 0; A(2,2) = -0.989992496600445; A(2,3) = 0.141120008059867;
    A(3,1) = 0; A(3,2) = -0.141120008059867; A(3,3) = -0.989992496600445;

    _assert(m_equals(A, R, 1e-10));
    return 0;
}

int i1_R_y_01(){
    Matrix& R = R_y(3);
    
    Matrix& A = zeros(3,3);
    A(1,1) = -0.989992496600445; A(1,2) = 0; A(1,3) = -0.141120008059867;
    A(2,1) = 0; A(2,2) = 1; A(2,3) = 0;
    A(3,1) = 0.141120008059867; A(3,2) = 0; A(3,3) = -0.989992496600445;

    _assert(m_equals(A, R, 1e-10));
    return 0;
}

int i1_R_z_01(){
    Matrix& R = R_z(3);
    
    Matrix& A = zeros(3,3);
    A(1,1) = -0.989992496600445; A(1,2) = 0.141120008059867; A(1,3) = 0;
    A(2,1) = -0.141120008059867; A(2,2) = -0.989992496600445; A(2,3) = 0;
    A(3,1) = 0; A(3,2) = 0; A(3,3) = 1;

    _assert(m_equals(A, R, 1e-10));
    return 0;
}

int i1_AccelPointMass_01(){
    Matrix& r = zeros(3);
    r(1)=1; r(2)=2; r(3)=3;

    Matrix& s = zeros(3);
    s(1)=4; s(2)=5; s(3)=6;

    Matrix& R = AccelPointMass(r,s,5);

    Matrix& expected = zeros(3);
    expected(1)=0.0773165667868213; expected(2)=0.0699165293543773; expected(3)=0.0625164919219332;

    _assert(m_equals(expected, R, 1e-10));
    return 0;
}

int i1_Cheb3D_01(){
    int N = 4;
    double Ta = 0;
    double Tb = 1;
    double t = 0.5;

    Matrix& Cx = zeros(N);
    Cx(1) = 1; Cx(2) = -0.5; Cx(3) = 0.3; Cx(4) = 0.1;
    Matrix& Cy = zeros(N);
    Cy(1) = 0; Cy(2) = 0.7; Cy(3) = -0.2; Cy(4) = 0.05;
    Matrix& Cz = zeros(N);
    Cz(1) = 0.5; Cz(2) = -0.1; Cz(3) = 0.2; Cz(4) = -0.05;

    Matrix& R = Cheb3D(t, N, Ta, Tb, Cx, Cy, Cz);

    Matrix& expected = zeros(3);
    expected(1) = 0.7; expected(2) = 0.2; expected(3) = 0.3;

    _assert(m_equals(expected, R, 1e-10));
    return 0;
}

int i1_EccAnom_01(){
    double R = EccAnom(3,4);
    
    double expected = 3.11327109405639;

    _assert(fabs(expected-R) < 1e-10);
    return 0;
}

int i1_Frac_01(){
    double R = Frac(1.2);
    
    double expected = 0.2;

    _assert(fabs(expected-R) < 1e-10);
    return 0;
}

int i1_MeanObliquity_01(){
    double R = MeanObliquity(3);
    
    double expected = 0.409413051544358;

    _assert(fabs(expected-R) < 1e-10);
    return 0;
}

int i1_Mjday_01(){
    double R = Mjday(1,2,3,4,5,6);
    
    double expected = -678556.829791667;

    _assert(fabs(expected-R) < 1e-9);
    return 0;
}

int i1_Mjday_TDB_01(){
    double R = Mjday_TDB(3.3);
    
    double expected = 3.2999999871282;

    _assert(fabs(expected-R) < 1e-10);
    return 0;
}

int i1_Position_01(){

    Matrix& R = Position(1.1, 2.2, 3.3);

    Matrix& expected = zeros(3);
    expected(1) = -1706329.66806146; expected(2) = -3352527.69377365; expected(3) = 5133425.98442718;

    _assert(m_equals(expected, R, 1e-8));
    return 0;
}

int i1_sign__01(){

    double R = sign_(-3.3,4.3);

    double expected = 3.3;

    _assert(fabs(expected - R)< 1e-10);
    return 0;
}

int i1_timediff_01(){

    double UT1_UTC = 1;
    double TAI_UTC = 2;

    auto [UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);

    double expected_UT1_TAI = -1;
    double expected_UTC_GPS = 17;
    double expected_UT1_GPS = 18;
    double expected_TT_UTC = 34.184;
    double expected_GPS_UTC = -17;

    _assert(fabs(expected_UT1_TAI - UT1_TAI)< 1e-10);
    _assert(fabs(expected_UTC_GPS - UTC_GPS)< 1e-10);
    _assert(fabs(expected_UT1_GPS - UT1_GPS)< 1e-10);
    _assert(fabs(expected_TT_UTC - TT_UTC)< 1e-10);
    _assert(fabs(expected_GPS_UTC - GPS_UTC)< 1e-10);
    return 0;
}

int i1_AzElPa_01(){

    Matrix& s = zeros(3);
    s(1) = 1; s(2) = 2; s(3) = 3;

    auto [Az, El, dAds, dEds] = AzElPa(s);

    double expected_Az = 0.463647609000806;
    double expected_El = 0.930274014115472;
    Matrix& expected_dAds = zeros(3);
    expected_dAds(1) = 0.4; expected_dAds(2) = -0.2; expected_dAds(3) = 0.0;
    Matrix& expected_dEds = zeros(3);
    expected_dEds(1) = -0.095831484749991; expected_dEds(2) = -0.191662969499982; expected_dEds(3) = 0.159719141249985;

    _assert(fabs(expected_Az - Az)< 1e-10);
    _assert(fabs(expected_El - El)< 1e-10);
    _assert(m_equals(expected_dAds, dAds, 1e-10));
    _assert(m_equals(expected_dEds, dEds, 1e-10));
    return 0;
}

int i1_IERS_01(){    
    double Mjd_UTC = AuxParam.Mjd_UTC;
    char interp = 'l';

    auto [x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC] = IERS(Mjd_UTC, interp);

    double expected_x_pole = -5.59378724204105e-07;
    double expected_y_pole = 2.33559834147171e-06;
    double expected_UT1_UTC = 0.325747632958789;
    double expected_LOD = 0.0027269897187419;
    double expected_dpsi = -1.16882953161739e-07;
    double expected_deps = -2.47835061986699e-08;
    double expected_dx_pole = -8.43027359620098e-10;
    double expected_dy_pole = -1.56811369104134e-09;
    double expected_TAI_UTC = 29;

    _assert(fabs(expected_x_pole - x_pole)< 1e-15);
    _assert(fabs(expected_y_pole - y_pole)< 1e-15);
    _assert(fabs(expected_UT1_UTC - UT1_UTC)< 1e-10);
    _assert(fabs(expected_LOD - LOD)< 1e-10);
    _assert(fabs(expected_dpsi - dpsi)< 1e-15);
    _assert(fabs(expected_deps - deps)< 1e-15);
    _assert(fabs(expected_dx_pole - dx_pole)< 1e-20);
    _assert(fabs(expected_dy_pole - dy_pole)< 1e-15);
    _assert(fabs(expected_TAI_UTC - TAI_UTC)< 1e-10);
    return 0;
}

int i1_Legendre_01(){
    int n = 2;
    int m = 3;
    double fi = 1.5;

    auto [pnm, dpnm] = Legendre(n, m, fi);

    Matrix& expected_pnm = zeros(n+1,m+1);
    expected_pnm(1,1) = 1; expected_pnm(1,2) = 0; expected_pnm(1,3) = 0; expected_pnm(1,4) = 0;
    expected_pnm(2,1) = 1.72771199709346; expected_pnm(2,2) = 0.122520427273707; expected_pnm(2,3) = 0; expected_pnm(2,4) = 0;
    expected_pnm(3,1) = 2.21928488408494; expected_pnm(3,2) = 0.273277720516261; expected_pnm(3,3) = 0.00968972350089721; expected_pnm(3,4) = 0;
    Matrix& expected_dpnm = zeros(n+1,m+1);
    expected_dpnm(1,1) = 0; expected_dpnm(1,2) = 0; expected_dpnm(1,3) = 0; expected_dpnm(1,4) = 0;
    expected_dpnm(2,1) = 0.122520427273707; expected_dpnm(2,2) = -1.72771199709346; expected_dpnm(2,3) = 0; expected_dpnm(2,4) = 0;
    expected_dpnm(3,1) = 0.473330896510772; expected_dpnm(3,2) = -3.83422445220383; expected_dpnm(3,3) = -0.273277720516261; expected_dpnm(3,4) = 0;

    _assert(m_equals(expected_pnm, pnm, 1e-10));
    _assert(m_equals(expected_dpnm, dpnm, 1e-10));
    return 0;
}

int i1_NutAngles_01(){
    double Mjd_TT = 1.5;

    auto [dpsi, deps] = NutAngles(Mjd_TT);

    double expected_dpsi = 2.72726579299296e-05;
    double expected_deps = 3.93565910019266e-05;

    _assert(fabs(expected_dpsi - dpsi)< 1e-10);
    _assert(fabs(expected_deps - deps)< 1e-10);
    return 0;
}

int i1_TimeUpdate_01(){
    Matrix& P = zeros(2,2);
    P(1,1) = 1; P(1,2) = 2;
    P(2,1) = 3; P(2,2) = 4;

    Matrix& Phi = zeros(2,2);
    Phi(1,1) = 5; Phi(1,2) = 6;
    Phi(2,1) = 7; Phi(2,2) = 8;

    double Qdt = 0.5;
    
    Matrix& R = TimeUpdate(P, Phi, Qdt);

    Matrix& expected = zeros(2,2);
    expected(1,1) = 319.5; expected(1,2) = 433.5;
    expected(2,1) = 431.5; expected(2,2) = 585.5;


    _assert(m_equals(expected, R, 1e-10));
    return 0;
}

int i2_AccelHarmonic_01(){
    Matrix& r = zeros(3);
    r(1)=7000e3; r(2)=0; r(3)=0;

    Matrix& E = eye(3);

    int n_max=2;
    int m_max=2;
    
    Matrix& R = AccelHarmonic(r,E,n_max,m_max);

    Matrix& expected = zeros(3);
    expected(1) = -8.14576607065686;
    expected(2) = -3.66267894892037e-05;
    expected(3) = -5.84508413583961e-09;

    _assert(m_equals(expected, R, 1e-10));
    return 0;
}

int i2_EqnEquinox_01(){
    double Mjd_TT=1.5;
    
    double R = EqnEquinox(Mjd_TT);

    double expected = 2.50186988674803e-05;

    _assert(fabs(R - expected) < 1e-10);
    return 0;
}

int i2_LTC_01(){
    double lon=1.5;
    double lat=1;
    
    Matrix& R = LTC(lon,lat);

    Matrix& expected = zeros(3,3);
    expected(1,1)= -0.997494986604054; expected(1,2)= 0.0707372016677029; expected(1,3)= 0;
    expected(2,1)= -0.0595233027498767; expected(2,2)= -0.839363088718653; expected(2,3)= 0.54030230586814; 
    expected(3,1)= 0.0382194731717195; expected(3,2)= 0.53894884135408; expected(3,3)= 0.841470984807897; 

    _assert(m_equals(R, expected, 1e-10));
    return 0;
}

int i2_NutMatrix_01(){
    double Mjd_TT=1.5;
    
    Matrix& R = NutMatrix(Mjd_TT);

    Matrix& expected = zeros(3,3);
    expected(1,1)= 0.999999999628101; expected(1,2)= -2.50186988643789e-05; expected(1,3)= -1.08564532657801e-05;
    expected(2,1)= 2.50182715720118e-05; expected(2,2)= 0.999999998912567; expected(2,3)= -3.93567267966133e-05; 
    expected(3,1)= 1.08574379080704e-05; expected(3,2)= 3.93564551722236e-05; expected(3,3)= 0.999999999166593; 

    _assert(m_equals(R, expected, 1e-10));
    return 0;
}

int i2_PoleMatrix_01(){
    double xp=2;
    double yp=3;
    
    Matrix& R = PoleMatrix(xp, yp);

    Matrix& expected = zeros(3,3);
    expected(1,1)= -0.416146836547142; expected(1,2)= 0.128320060202457; expected(1,3)= -0.900197629735517;
    expected(2,1)= 0; expected(2,2)= -0.989992496600445; expected(2,3)= -0.141120008059867; 
    expected(3,1)= -0.909297426825682; expected(3,2)= -0.058726644927621; expected(3,3)= 0.411982245665683; 

    _assert(m_equals(R, expected, 1e-10));
    return 0;
}

int i2_PrecMatrix_01(){
    double Mjd_1=2;
    double Mjd_2=3;
    
    Matrix& R = PrecMatrix(Mjd_1, Mjd_2);

    Matrix& expected = zeros(3,3);
    expected(1,1)= 0.999999999999778; expected(1,2)= -6.11707327974946e-07; expected(1,3)= -2.66201482252295e-07;
    expected(2,1)= 6.11707327974946e-07; expected(2,2)= 0.999999999999813; expected(2,3)= -8.14186990889656e-14; 
    expected(3,1)= 2.66201482252295e-07; expected(3,2)= -8.1418698322574e-14; expected(3,3)= 0.999999999999965; 

    _assert(m_equals(R, expected, 1e-10));
    return 0;
}

int i2_gmst_01(){
    double Mjd_UT1=0.5;
    
    double R = gmst(Mjd_UT1);

    double expected = 4.12340219740033;

    _assert(fabs(R - expected) < 1e-10);
    return 0;
}

int i2_JPL_Eph_DE430_01(){
    double Mjd_TDB = 60348;
    
    auto [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, 
        r_Neptune,r_Pluto,r_Moon,r_Sun] = JPL_Eph_DE430(Mjd_TDB);

    Matrix& expected_r_Mercury = zeros(3);
    expected_r_Mercury(1)=112623958311.779;
    expected_r_Mercury(2)=-150868623524.381;
    expected_r_Mercury(3)=-71779751443.9542;

    _assert(m_equals(r_Mercury, expected_r_Mercury,10e-4));

    Matrix& expected_r_Venus = zeros(3);
    expected_r_Venus(1)=67979665801.0159;
    expected_r_Venus(2)=-182040469650.972;
    expected_r_Venus(3)=-77754531631.0774;

    _assert(m_equals(r_Venus, expected_r_Venus,10e-4));

    Matrix& expected_r_Earth = zeros(3);
    expected_r_Earth(1)=-111448106096.032;
    expected_r_Earth(2)=89490608943.6678;
    expected_r_Earth(3)=38828100871.8833;

    _assert(m_equals(r_Earth, expected_r_Earth,10e-4));

    Matrix& expected_r_Mars = zeros(3);
    expected_r_Mars(1)=148466860011.68;
    expected_r_Mars(2)=-281603620116.961;
    expected_r_Mars(3)=-127931017643.326;

    _assert(m_equals(r_Mars, expected_r_Mars,10e-4));

    Matrix& expected_r_Jupiter = zeros(3);
    expected_r_Jupiter(1)=600849652579.016;
    expected_r_Jupiter(2)=431907465367.823;
    expected_r_Jupiter(3)=172747882009.303;

    _assert(m_equals(r_Jupiter, expected_r_Jupiter,10e-3));

    Matrix& expected_r_Saturn = zeros(3);
    expected_r_Saturn(1)=1466091955279.12;
    expected_r_Saturn(2)=-555189637093.051;
    expected_r_Saturn(3)=-289531886725.713;

    _assert(m_equals(r_Saturn, expected_r_Saturn,10e-4));

    Matrix& expected_r_Uranus = zeros(3);
    expected_r_Uranus(1)=1928309175333.92;
    expected_r_Uranus(2)=2027905259796.94;
    expected_r_Uranus(3)=862836636210.81;

    _assert(m_equals(r_Uranus, expected_r_Uranus,10e-3));

    Matrix& expected_r_Neptune = zeros(3);
    expected_r_Neptune(1)=4575619355154.4;
    expected_r_Neptune(2)=-280382588638.193;
    expected_r_Neptune(3)=-228103255917.356;

    _assert(m_equals(r_Neptune, expected_r_Neptune,10e-3));

    Matrix& expected_r_Pluto = zeros(3);
    expected_r_Pluto(1)=2700803532076.77;
    expected_r_Pluto(2)=-4144525986799.46;
    expected_r_Pluto(3)=-2084446068792.87;

    _assert(m_equals(r_Pluto, expected_r_Pluto,10e-2));

    Matrix& expected_r_Moon = zeros(3);
    expected_r_Moon(1)=130639413.73261;
    expected_r_Moon(2)=-298652884.800457;
    expected_r_Moon(3)=-164607636.963072;

    _assert(m_equals(r_Moon, expected_r_Moon,10e-5));

    Matrix& expected_r_Sun = zeros(3);
    expected_r_Sun(1)=110284675128.512;
    expected_r_Sun(2)=-89937859346.5398;
    expected_r_Sun(3)=-38988031316.0913;

    _assert(m_equals(r_Sun, expected_r_Sun,10e-4));

    return 0;
}

int i3_gast_01(){
    double Mjd_UT1=1.5;
    
    double R = gast(Mjd_UT1);

    double expected = 4.14063000738129;

    _assert(fabs(R - expected) < 1e-10);
    return 0;
}

int i3_G_AccelHarmonic_01(){
    Matrix& r = zeros(3);
    r(1) = 1000; r(2) = 2000; r(3) = 3000;

    Matrix& U = eye(3);

    int n_max = 4;
    int m_max = 4;
    
    Matrix& R = G_AccelHarmonic(r,U,n_max,m_max);

    Matrix& expected = zeros(3,3);
    expected(1,1) = 1207384439706; expected(1,2) = -1018762589128.5; expected(1,3) = -1723266165859.38;
    expected(2,1) = -1018762610136.75; expected(2,2) = 1190241885267.38; expected(2,3) = -1784910727629.88;
    expected(3,1) = -1723265773982.25; expected(3,2) = -1784910207689.25; expected(3,3) = -2397626334203.5;

    _assert(m_equals(R, expected, 1e1));
    return 0;
}

int i3_GHAMatrix_01(){
    double Mjd_UT1 = 1.5;
    
    Matrix& R = GHAMatrix(Mjd_UT1);

    Matrix& expected = zeros(3,3);
    expected(1,1) = -0.541112094250387; expected(1,2) = -0.840950475031651; expected(1,3) = 0;
    expected(2,1) = 0.840950475031651; expected(2,2) = -0.541112094250387; expected(2,3) = 0;
    expected(3,1) = 0; expected(3,2) = 0; expected(3,3) = 1;

    _assert(m_equals(R, expected, 1e-10));
    return 0;
}

int i3_MeasUpdate_01(){
    Matrix& x = zeros(2);
    x(1)=1000; x(2)=10;

    double z = 1010.0;
    double g = 1000.0;
    double s = 5.0;

    Matrix& G = zeros(2);
    G(1) = 1; G(2) = 0;

    Matrix& P = zeros(2,2);
    P(1,1) = 25; P(1,2) = 0;
    P(2,1) = 0; P(2,2) = 4;

    int n = 2;

    Matrix& K = zeros(1);
    
    tie(K, x, P) = MeasUpdate(x, z, g, s, G, P, n);

    Matrix& expected_K = zeros(2);
    expected_K(1) = 0.5; expected_K(2) = 0; 

    _assert(m_equals(K, expected_K, 1e-10));

    Matrix& expected_x = zeros(2);
    expected_x(1) = 1005; expected_x(2) = 10; 

    _assert(m_equals(x, expected_x, 1e-10));

    Matrix& expected_P = zeros(2,2);
    expected_P(1,1) = 12.5; expected_P(1,2) = 0; 
    expected_P(2,1) = 0; expected_P(2,2) = 4; 

    _assert(m_equals(P, expected_P, 1e-10));
    return 0;
}

int i3_Accel_01(){
    double x = 123456;
    Matrix& Y = zeros(6);
    Y(1) = 1000000;
    Y(2) = 2000000;
    Y(3) = 3000000;
    Y(4) = 4000000;
    Y(5) = 5000000;
    Y(6) = 6000000;
    
    Matrix& R = Accel(x,Y);

    Matrix& expected = zeros(6);
    expected(1) = 4000000;
    expected(2) = 5000000;
    expected(3) = 6000000;
    expected(4) = -5.70209501208335;
    expected(5) = -17.9604015906774;
    expected(6) = -32.0686615657293;

    _assert(m_equals(R, expected, 1e-8));
    return 0;
}

int i3_VarEqn_01(){
    double x = 0;
    Matrix& yPhi = zeros(42);
    yPhi(1) =       5542555.93722861;
    yPhi(2) =        3213514.8673492;
    yPhi(3) =       3990892.97587685;
    yPhi(4) =       5394.06842166351;
    yPhi(5) =      -2365.21337882342;
    yPhi(6) =      -7061.84554200295;
    yPhi(7) =                      1;
    yPhi(8) =                      0;
    yPhi(9) =                      0;
    yPhi(10) =                      0;
    yPhi(11) =                      0;
    yPhi(12) =                      0;
    yPhi(13) =                      0;
    yPhi(14) =                      1;
    yPhi(15) =                      0;
    yPhi(16) =                      0;
    yPhi(17) =                      0;
    yPhi(18) =                      0;
    yPhi(19) =                      0;
    yPhi(20) =                      0;
    yPhi(21) =                      1;
    yPhi(22) =                      0;
    yPhi(23) =                      0;
    yPhi(24) =                      0;
    yPhi(25) =                      0;
    yPhi(26) =                      0;
    yPhi(27) =                      0;
    yPhi(28) =                      1;
    yPhi(29) =                      0;
    yPhi(30) =                      0;
    yPhi(31) =                      0;
    yPhi(32) =                      0;
    yPhi(33) =                      0;
    yPhi(34) =                      0;
    yPhi(35) =                      1;
    yPhi(36) =                      0;
    yPhi(37) =                      0;
    yPhi(38) =                      0;
    yPhi(39) =                      0;
    yPhi(40) =                      0;
    yPhi(41) =                      0;
    yPhi(42) =                      1;
    
    AuxParam.Mjd_TT = 49746.1108586111;
    Matrix& R = VarEqn(x,yPhi);

    Matrix& expected = zeros(42);
    expected(1) =      5394.06842166351;
    expected(2) =      -2365.21337882342;
    expected(3) =      -7061.84554200295;
    expected(4) =       -5.1348367854085;
    expected(5) =      -2.97717622353621;
    expected(6) =      -3.70591776714204;
    expected(7) =                      0;
    expected(8) =                      0;
    expected(9) =                      0;
    expected(10) =   5.70032035795975e-07;
    expected(11) =   8.67651593239316e-07;
    expected(12) =   1.08169354007259e-06;
    expected(13) =                      0;
    expected(14) =                      0;
    expected(15) =                      0;
    expected(16) =   8.67651590574781e-07;
    expected(17) =  -4.23359106882515e-07;
    expected(18) =   6.27183702306411e-07;
    expected(19) =                      0;
    expected(20) =                      0;
    expected(21) =                      0;
    expected(22) =   1.08169353651988e-06;
    expected(23) =   6.27183704082768e-07;
    expected(24) =  -1.46672928913461e-07;
    expected(25) =                      1;
    expected(26) =                      0;
    expected(27) =                      0;
    expected(28) =                      0;
    expected(29) =                      0;
    expected(30) =                      0;
    expected(31) =                      0;
    expected(32) =                      1;
    expected(33) =                      0;
    expected(34) =                      0;
    expected(35) =                      0;
    expected(36) =                      0;
    expected(37) =                      0;
    expected(38) =                      0;
    expected(39) =                      1;
    expected(40) =                      0;
    expected(41) =                      0;
    expected(42) =                      0;

    _assert(m_equals(R, expected, 1e-8));
    return 0;
}

int i4_DEInteg_01(){
    double t = 0.0;
    double tout = -134.999991953373;
    double relerr = 1e-13;
    double abserr = 1e-6;
    int n_eqn = 6;
    Matrix& y = zeros(6);
    y(1) = 6221397.62857869;
    y(2) = 2867713.77965738;
    y(3) = 3006155.98509949;
    y(4) = 4645.04725161806;
    y(5) = -2752.21591588204;
    y(6) = -7507.99940987031;
    
    Matrix& R = DEInteg(Accel, t, tout, relerr, abserr, n_eqn, y);

    Matrix& expected = zeros(6);
    expected(1) = 5542555.89427451;
    expected(2) = 3213514.83814162;
    expected(3) = 3990892.92789074;
    expected(4) = 5394.06894044389;
    expected(5) = -2365.21290574021;
    expected(6) = -7061.8448137347;

    _assert(m_equals(R, expected, 1e-5));
    return 0;
}

int extra_Geodetic_01(){
    Matrix& r = zeros(3);
    r(1) = 4510732.0; r(2) = 4510732.0; r(3) = 4510732.0;
    auto [lon, lat, h] = Geodetic(r);
    _assert(fabs(lon - 0.785398163397448)<1e-10);
    _assert(fabs(lat - 0.61806355298986)<1e-10);
    _assert(fabs(h - 1441826.98849733)<1e-5);
    return 0;
}

int extra_angl_01(){
    Matrix& vec1 = zeros(3);
    vec1(1) = 1; vec1(2) = 2; vec1(3) = 3;
    Matrix& vec2 = zeros(3);
    vec2(1) = 4; vec2(2) = 5; vec2(3) = 6;
    double theta = angl(vec1,vec2);
    _assert(fabs(theta - 0.225726128552734)<1e-10);
    return 0;
}

int extra_elements_01(){
    Matrix& y = zeros(6);
    y(1) = 1; y(2) = 2; y(3) = 3;
    y(4) = 4; y(5) = 5; y(6) = 6;
    auto [p, a, e, i, Omega, omega, M] = elements(y);
    _assert(fabs(p - 1.35474011564823e-13)<1e-20);
    _assert(fabs(a - 1.87082869338765)<1e-10);
    _assert(fabs(e - 0.999999999999964)<1e-10);
    _assert(fabs(i - 1.99133066207886)<1e-10);
    _assert(fabs(Omega - 3.6052402625906)<1e-10);
    _assert(fabs(omega - 5.21086941752228)<1e-10);
    _assert(fabs(M - 3.14159030993265)<1e-10);
    return 0;
}

int extra_unit_01(){
    Matrix& vec = zeros(3);
    vec(1) = 3; vec(2) = 4; vec(3) = 0;
    Matrix& R = unit(vec);
    Matrix& expected = zeros(3);
    expected(1) = 0.6; expected(2) = 0.8; expected(3) = 0.0;
    _assert(m_equals(R,expected,1e-10));
    return 0;
}

int extra_hgibbs_01(){
    Matrix& r1 = zeros(3);
    r1(1) = 1000000; r1(2) = 2; r1(3) = 3;
    Matrix& r2 = zeros(3);
    r2(1) = 4000000; r2(2) = 5; r2(3) = 6;
    Matrix& r3 = zeros(3);
    r3(1) = 7000000; r3(2) = 8; r3(3) = 9;
    double Mjd1 = 10.5;
    double Mjd2 = 11.5;
    double Mjd3 = 12.5;
    auto [v2, theta,theta1,copa, error] = hgibbs (r1,r2,r3,Mjd1,Mjd2,Mjd3);
    Matrix& expected_v2 = zeros(3);
    expected_v2(1) = -2811318.55296047; expected_v2(2) = -5.67287456520078; expected_v2(3) = -8.53443057744108;
    _assert(m_equals(v2,expected_v2,1e-5));
    _assert(fabs(theta - 1.67709242756207e-06)<1e-15);
    _assert(fabs(theta1 - 2.39811494043933e-07)<1e-15);
    _assert(fabs(copa - 2.64697796016969e-22)<1e-20);
    _assert(error=="          ok");
    
    return 0;
}

int extra_gibbs_01(){
    Matrix& r1 = zeros(3);
    r1(1) = 7000000; r1(2) = 2; r1(3) = 3;
    Matrix& r2 = zeros(3);
    r2(1) = 7000000; r2(2) = 5; r2(3) = 6;
    Matrix& r3 = zeros(3);
    r3(1) = 7000000; r3(2) = 0; r3(3) = 0;
    auto [v2, theta,theta1,copa, error] = gibbs (r1,r2,r3);
    Matrix& expected_v2 = zeros(3);
    expected_v2(1) = 0.0; expected_v2(2) = 0.0167075643908466; expected_v2(3) = 0.0134726032730966;
    _assert(m_equals(v2,expected_v2,1e-10));
    _assert(fabs(theta - 6.06021267404052e-07)<1e-15);
    _assert(fabs(theta1 - 1.11579751739049e-06)<1e-15);
    _assert(fabs(copa - -5.48729485426625e-08)<1e15);
    _assert(error=="          ok");
    
    return 0;
}

int extra_anglesg_01(){
    double az1 = 1.0559084894933;
    double az2 = 1.36310214580757;
    double az3 = 1.97615602688759;
    double el1 = 0.282624656433946;
    double el2 = 0.453434794338875;
    double el3 = 0.586427138011591;
    double Mjd1 = 49746.1101504629;
    double Mjd2 = 49746.1112847221;
    double Mjd3 = 49746.1125347223;
    Matrix& Rs1 = zeros(3);
    Rs1(1) = -5512567.84003607; Rs1(2) = -2196994.44666933; Rs1(3) = 2330804.96614689;
    Matrix& Rs2 = Rs1*1;
    Matrix& Rs3 = Rs1*1;
    auto [r2, v2] = anglesg ( az1,az2,az3,el1,el2,el3,Mjd1,Mjd2,Mjd3,Rs1,Rs2,Rs3 );
    Matrix& expected_r2 = zeros(3);
    expected_r2(1) = 6221397.63831071; expected_r2(2) = 2867713.82724908; expected_r2(3) = 3006156.00237787;
    _assert(m_equals(r2,expected_r2,1e-7));
    Matrix& expected_v2 = zeros(3);
    expected_v2(1) = 4645.04722163732; expected_v2(2) = -2752.21634315407; expected_v2(3) = -7507.99964528836;
    _assert(m_equals(v2,expected_v2,1e-7));
    
    return 0;
}


int all_tests()
{
    // Matrix
    _verify(m_sum_01);
    _verify(m_sub_01);
	_verify(m_prod_01);
	_verify(m_asign_01);
    _verify(m_div_01);
    _verify(m_sum_02);
	_verify(m_sub_02);
	_verify(m_prod_02);
	_verify(m_div_02);
	_verify(m_transpose_01);
	_verify(m_extract_vector_01);
	_verify(m_union_vector_01);
    _verify(m_assign_row_01);
    _verify(m_assign_column_01);
    _verify(m_extract_row_01);
    _verify(m_extract_column_01);
    _verify(m_zeros_01);
	_verify(m_zeros_02);
	_verify(m_eye_01);
	_verify(m_norm_01);
    _verify(m_dot_01);
    _verify(m_cross_01);
    _verify(m_inv_01);

    // Iteration 1
    _verify(i1_R_x_01);
    _verify(i1_R_y_01);
    _verify(i1_R_z_01);
    _verify(i1_AccelPointMass_01);
    _verify(i1_Cheb3D_01);
    _verify(i1_EccAnom_01);
    _verify(i1_Frac_01);
    _verify(i1_MeanObliquity_01);
    _verify(i1_Mjday_01);
    _verify(i1_Mjday_TDB_01);
    _verify(i1_Position_01);
    _verify(i1_sign__01);
    _verify(i1_timediff_01);
    _verify(i1_AzElPa_01);
    _verify(i1_IERS_01);
    _verify(i1_Legendre_01);
    _verify(i1_NutAngles_01);
    _verify(i1_TimeUpdate_01);

    // Iteration 2
    _verify(i2_AccelHarmonic_01);
    _verify(i2_EqnEquinox_01);
    _verify(i2_LTC_01);
    _verify(i2_NutMatrix_01);
    _verify(i2_PoleMatrix_01);
    _verify(i2_PrecMatrix_01);
    _verify(i2_gmst_01);
    _verify(i2_JPL_Eph_DE430_01);

    // Iteration 3
    _verify(i3_gast_01);
    _verify(i3_G_AccelHarmonic_01);
    _verify(i3_GHAMatrix_01);
    _verify(i3_MeasUpdate_01);
    _verify(i3_Accel_01);
    _verify(i3_VarEqn_01);

    // Iteration 4
    _verify(i4_DEInteg_01);

    // Extra
    _verify(extra_Geodetic_01);
    _verify(extra_angl_01);
    _verify(extra_elements_01);
    _verify(extra_unit_01);
    _verify(extra_hgibbs_01);
    _verify(extra_gibbs_01);
    _verify(extra_anglesg_01);

    return 0;
}


int main()
{
    cout<<"Running tests..."<<endl;

    eop19620101(21413);
    GGM03S(181);
    DE430Coeff(2285, 1020);
    initializeAuxParam();

    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
