#ifndef _MATRIX_
#define _MATRIX_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

class Matrix
{
public:
	int n_row, n_column;
	double **data;

	// Parameterized constructor
	Matrix(const int v_size);
	Matrix(const int n_row, const int n_column);

	// Member operators
	double &operator()(const int n);
	double &operator()(const int row, const int column);
	Matrix &operator+(Matrix &m);
	Matrix &operator-(Matrix &m);
	Matrix &operator*(Matrix &m);
	Matrix &operator=(Matrix &m);
	Matrix& operator / (Matrix &m); //multiplicar por la inversa de la matriz

	Matrix& operator + (const double k); //suma escalar de cada elemento
	Matrix& operator - (const double k); //resta escalar de cada elemento
	Matrix& operator * (const double k); //producto escalar de cada elemento
	Matrix& operator / (const double k); //división escalar de cada elemento

	// Methods
	Matrix &transpose();
	Matrix& extract_vector(const int from, const int to); //devuelve el subvector desde el elemento 'from' hasta 'to'
	Matrix& union_vector(Matrix &v); //une el vector, lo concatena
	Matrix& extract_row(const int row); //devuelve la fila indicada (forma de vector)
	Matrix& extract_column(const int col); //devuelve la columna indicada (en forma de vector?)
	Matrix& assign_row(const int row, Matrix &v); //asigna el vector a la fila indicada (deben coincidir las longitudes)
	Matrix& assign_column(const int col, Matrix &v); //asigna el vector a la columna indicada (deben coincidir las longitudes, el vector debe estar en forma vector o vertical?)

	// Non-member operators
	friend ostream &operator<<(ostream &o, Matrix &m);
};

// Operator overloading
ostream &operator<<(ostream &o, Matrix &m);

// Methods
Matrix &zeros(const int n_row, const int n_column);
Matrix& zeros(const int v_size); //vector de ceros
Matrix& eye(const int n); //matriz identidad de tamaño n
double norm(Matrix &v); //norma euclidiana del vector
double dot(Matrix &v1, Matrix &v2); //producto escalar de vectores
Matrix& cross(Matrix &v1, Matrix &v2); //producto vectorial de vectores de longitud 3
Matrix &inv(Matrix &m); // devuelve la matriz inversa (si es cuadrada y con det!=0)

#endif