#include "..\include\matrix.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include "..\include\AccelPointMass.hpp"
#include "..\include\Cheb3D.hpp"
#include <cstdio>
#include <cmath>

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

int m_equals(Matrix A, Matrix B, double p) {
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
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = 2; C(1,2) =  2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = 8; C(2,2) = -3; C(2,3) = 1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = -2; C(3,3) = 0; C(3,4) = 7;
	
	Matrix R = A + B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_sub_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) = 0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = -2; C(1,2) = 2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = -6; C(2,2) = 1; C(2,3) = -1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = 4; C(3,3) = 0; C(3,4) = 3;
	
	Matrix R = A - B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_prod_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(c, f);
	B(1,1) = 2; B(1,2) = 7; B(1,3) = 0;
	B(2,1) = 0; B(2,2) = -2; B(2,3) = -3;
	B(3,1) = 0; B(3,2) = 1; B(3,3) = 0;
	B(4,1) = 0; B(4,2) = 0; B(4,3) = 2;
	
	Matrix C(f, f);
	C(1,1) = 0; C(1,2) = 4; C(1,3) = -6;
	C(2,1) = 2; C(2,2) = 9; C(2,3) = 3;
	C(3,1) = 0; C(3,2) = -2; C(3,3) = 7;
	
	Matrix R = A * B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_asign_01(){
	Matrix A(2, 2);
    A(1,1) = 1; A(1,2) = 2;
    A(2,1) = 3; A(2,2) = 4;

    Matrix B(2, 2);

    B = A;

	_assert(m_equals(A, B, 1e-10));

    A(1,1) = 99;

    _assert(B(1,1) == 1 && A(1,1) == 99);

    return 0;
}

int m_div_01(){
    Matrix A(2,2);
	A(1,1) = 1; A(1,2) = 2;
	A(2,1) = 3; A(2,2) = 4;

    Matrix B(2,2);
	B(1,1) = 5; B(1,2) = 6;
	B(2,1) = 7; B(2,2) = 8;

    Matrix R = A / B;

    Matrix C(2,2);
	C(1,1) = 3; C(1,2) = -2;
	C(2,1) = 2; C(2,2) = -1;

    _assert(m_equals(C, R, 1e-10));
    return 0;
}

int m_sum_02() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 3; B(1,2) = 5; B(1,3) = 11; B(1,4) = 3;
	B(2,1) = 4; B(2,2) = 2; B(2,3) = 3; B(2,4) = 3;
	B(3,1) = 3; B(3,2) = 4; B(3,3) = 3; B(3,4) = 8;
	
	Matrix R = A + 3;
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_sub_02() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = -2; B(1,2) = 0; B(1,3) = 6; B(1,4) = -2;
	B(2,1) = -1; B(2,2) = -3; B(2,3) = -2; B(2,4) = -2;
	B(3,1) = -2; B(3,2) = -1; B(3,3) = -2; B(3,4) = 3;
	
	Matrix R = A - 2;
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_prod_02() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 0; B(1,2) = 4; B(1,3) = 16; B(1,4) = 0;
	B(2,1) = 2; B(2,2) = -2; B(2,3) = 0; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = 2; B(3,3) = 0; B(3,4) = 10;
	
	Matrix R = A * 2;
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_div_02() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 0; B(1,2) = 4; B(1,3) = 16; B(1,4) = 0;
	B(2,1) = 2; B(2,2) = -2; B(2,3) = 0; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = 2; B(3,3) = 0; B(3,4) = 10;
	
	Matrix R = B / 2;
    
    _assert(m_equals(A, R, 1e-10));
    
    return 0;
}

int m_transpose_01() {
    int f = 3;
    int c = 4;
    
    Matrix A(f, c);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3; A(1,4) = 4;
    A(2,1) = 5; A(2,2) = 6; A(2,3) = 7; A(2,4) = 8;
    A(3,1) = 9; A(3,2) = 10; A(3,3) = 11; A(3,4) = 12;

    Matrix T(c, f);
    T(1,1) = 1; T(1,2) = 5; T(1,3) = 9;
    T(2,1) = 2; T(2,2) = 6; T(2,3) = 10;
    T(3,1) = 3; T(3,2) = 7; T(3,3) = 11;
    T(4,1) = 4; T(4,2) = 8; T(4,3) = 12;

    Matrix R = A.transpose();

    _assert(m_equals(T, R, 1e-10));

    return 0;
}

int m_extract_vector_01() {
    Matrix A(5);
    A(1) = 1; A(2) = 2; A(3) = 3; A(4) = 4; A(5) = 5;

    Matrix sub_vector = A.extract_vector(2,4);

    Matrix expected(3);
    expected(1) = 2; expected(2) = 3; expected(3) = 4;

    _assert(m_equals(sub_vector, expected, 1e-10));
    return 0;
}

int m_union_vector_01(){
	Matrix m1(3);
    Matrix m2(4);

    for (int i = 1; i <= 3; i++) {
        m1(i) = i;
    }
    for (int i = 1; i <= 4; i++) {
        m2(i) = i + 3;
    }

    Matrix result = m1.union_vector(m2);

	Matrix expected(7);
	for (int i = 1; i <= 7; i++) {
        expected(i) = i;
    }
	
	_assert(m_equals(result, expected, 1e-10));
    return 0;
}

int m_extract_row_01() {
    Matrix A(3, 3);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;
    A(3,1) = 7; A(3,2) = 8; A(3,3) = 9;

    Matrix row = A.extract_row(2);

    Matrix expected(3);
    expected(1) = 4; expected(2) = 5; expected(3) = 6;

    _assert(m_equals(row, expected, 1e-10));
    return 0;
}

int m_extract_column_01() {
    Matrix A(3, 3);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;
    A(3,1) = 7; A(3,2) = 8; A(3,3) = 9;

    Matrix col = A.extract_column(3);

    Matrix expected(3);
    expected(1) = 3; expected(2) = 6; expected(3) = 9;

    _assert(m_equals(col, expected, 1e-10));
    return 0;
}

int m_assign_row_01() {
    Matrix A(3, 3);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;
    A(3,1) = 7; A(3,2) = 8; A(3,3) = 9;

    Matrix row(3);
    row(1) = -1; row(2) = -2; row(3) = -3;

    A.assign_row(2, row);

    Matrix expected(3, 3);
    expected(1,1) = 1; expected(1,2) = 2;  expected(1,3) = 3;
    expected(2,1) = -1; expected(2,2) = -2; expected(2,3) = -3;
    expected(3,1) = 7; expected(3,2) = 8;  expected(3,3) = 9;

    _assert(m_equals(A, expected, 1e-10));
    return 0;
}

int m_assign_column_01() {
    Matrix A(3, 3);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;
    A(3,1) = 7; A(3,2) = 8; A(3,3) = 9;

    Matrix col(3);
    col(1) = -1; col(2) = -2; col(3) = -3;

    A.assign_column(2, col);

    Matrix expected(3, 3);
    expected(1,1) = 1; expected(1,2) = -1; expected(1,3) = 3;
    expected(2,1) = 4; expected(2,2) = -2; expected(2,3) = 6;
    expected(3,1) = 7; expected(3,2) = -3; expected(3,3) = 9;

    _assert(m_equals(A, expected, 1e-10));
    return 0;
}

int m_zeros_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 0; A(1,3) = 0; A(1,4) = 0;
	A(2,1) = 0; A(2,2) = 0; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 0; A(3,4) = 0;
	
	Matrix B = zeros(3, 4);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_zeros_02() {	
	Matrix A(4);
	A(1,1) = 0; A(2) = 0; A(3) = 0; A(4) = 0;
	
	Matrix B = zeros(4);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_eye_01() {		
	Matrix A(3,3);
	A(1,1) = 1; A(1,2) = 0; A(1,3) = 0;
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 1;
	
	Matrix B = eye(3);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_norm_01() {		
	Matrix A(3);
	A(1) = 12; A(2) = 3; A(3) = 4;
	
	double result = norm(A);

	double expected_norm = 13;
    
    _assert(fabs(result - expected_norm) < 1e-10);
    
    return 0;
}

int m_dot_01(){
    Matrix A(3);
	A(1) = 1; A(2) = 2; A(3) = 3;

    Matrix B(3);
	B(1) = 4; B(2) = 5; B(3) = 6;
	
	_assert(fabs(dot(A,B)-32.0)<1e-10);
    return 0;
}

int m_cross_01(){
    Matrix A(3);
	A(1) = 1; A(2) = 2; A(3) = 3;

    Matrix B(3);
	B(1) = 4; B(2) = 5; B(3) = 6;

    Matrix expected(3);
	expected(1) = -3; expected(2) = 6; expected(3) = -3;

    Matrix R = cross(A,B);
	
	_assert(m_equals(expected, R, 1e-10));
    return 0;
}

int m_inv_01(){
    Matrix A(3,3);
	A(1,1) = 2; A(1,2) = 1; A(1,3) = 1;
    A(2,1) = 1; A(2,2) = 1; A(2,3) = 1;
    A(3,1) = 1; A(3,2) = 1; A(3,3) = 2;

    Matrix R = inv(A);

    Matrix B(3,3);
	B(1,1) = 1; B(1,2) = -1; B(1,3) = 0;
    B(2,1) = -1; B(2,2) = 3; B(2,3) = -1;
    B(3,1) = 0; B(3,2) = -1; B(3,3) = 1;
    
    _assert(m_equals(B, R, 1e-10));
    return 0;
}


int i1_R_x_01(){
    Matrix R = R_x(3);
    
    Matrix A(3,3);
    A(1,1) = 1; A(1,2) = 0; A(1,3) = 0;
    A(2,1) = 0; A(2,2) = -0.989992496600445; A(2,3) = 0.141120008059867;
    A(3,1) = 0; A(3,2) = -0.141120008059867; A(3,3) = -0.989992496600445;

    _assert(m_equals(A, R, 1e-10));
    return 0;
}

int i1_R_y_01(){
    Matrix R = R_y(3);
    
    Matrix A(3,3);
    A(1,1) = -0.989992496600445; A(1,2) = 0; A(1,3) = -0.141120008059867;
    A(2,1) = 0; A(2,2) = 1; A(2,3) = 0;
    A(3,1) = 0.141120008059867; A(3,2) = 0; A(3,3) = -0.989992496600445;

    _assert(m_equals(A, R, 1e-10));
    return 0;
}

int i1_R_z_01(){
    Matrix R = R_z(3);
    
    Matrix A(3,3);
    A(1,1) = -0.989992496600445; A(1,2) = 0.141120008059867; A(1,3) = 0;
    A(2,1) = -0.141120008059867; A(2,2) = -0.989992496600445; A(2,3) = 0;
    A(3,1) = 0; A(3,2) = 0; A(3,3) = 1;

    _assert(m_equals(A, R, 1e-10));
    return 0;
}

int i1_AccelPointMass_01(){
    Matrix r(3);
    r(1)=1; r(2)=2; r(3)=3;

    Matrix s(3);
    s(1)=4; s(2)=5; s(3)=6;

    Matrix R = AccelPointMass(r,s,5);

    Matrix expected(3);
    expected(1)=0.0773165667868213; expected(2)=0.0699165293543773; expected(3)=0.0625164919219332;

    _assert(m_equals(expected, R, 1e-10));
    return 0;
}

int i1_Cheb3D_01(){
    int N = 4;
    double Ta = 0;
    double Tb = 1;
    double t = 0.5;

    Matrix Cx(N);
    Cx(1) = 1; Cx(2) = -0.5; Cx(3) = 0.3; Cx(4) = 0.1;
    Matrix Cy(N);
    Cy(1) = 0; Cy(2) = 0.7; Cy(3) = -0.2; Cy(4) = 0.05;
    Matrix Cz(N);
    Cz(1) = 0.5; Cz(2) = -0.1; Cz(3) = 0.2; Cz(4) = -0.05;

    Matrix R = Cheb3D(t, N, Ta, Tb, Cx, Cy, Cz);

    Matrix expected(3);
    expected(1) = 0.7; expected(2) = 0.2; expected(3) = 0.3;

    _assert(m_equals(expected, R, 1e-10));
    return 0;
}


int all_tests()
{
    //Matrix
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

    //Iteration 1
    _verify(i1_R_x_01);
    _verify(i1_R_y_01);
    _verify(i1_R_z_01);
    _verify(i1_AccelPointMass_01);
    _verify(i1_Cheb3D_01);

    return 0;
}


int main()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
