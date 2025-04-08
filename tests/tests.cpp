#include "..\include\matrix.h"
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


int all_tests()
{
    _verify(m_sum_01);
    _verify(m_sub_01);
	_verify(m_prod_01);
	_verify(m_asign_01);
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
