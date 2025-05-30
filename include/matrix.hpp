/**
 * @file Matrix.hpp
 * @brief Declaration of the Matrix class and related functions for basic matrix operations.
 *
 * This file contains the definition of the Matrix class, which supports
 * dynamic allocation of 2D matrices, arithmetic operators, and linear algebra methods.
 */

#ifndef _MATRIX_
#define _MATRIX_

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>

using namespace std;

/**
 * @class Matrix
 * @brief Class for basic matrix operations and storage.
 *
 * This class represents a matrix with dynamic memory allocation,
 * supports arithmetic operations with matrices and scalars,
 * and provides methods for common linear algebra operations.
 */
class Matrix
{
public:
	int n_row; /**< Number of rows in the matrix */
	int n_column; /**< Number of columns in the matrix */
	double **data; /**< 2D array to hold matrix elements */

	// Parameterized constructors
	/**
     * @brief Default constructor, creates an empty matrix.
     */
	Matrix();

	/**
     * @brief Constructs a vector with given size (1 row x v_size columns).
     * @param v_size Size of the vector.
     */
	Matrix(const int v_size);

	/**
     * @brief Constructs a matrix with specified number of rows and columns.
     * @param n_row Number of rows.
     * @param n_column Number of columns.
     */
	Matrix(const int n_row, const int n_column);

	// Member operators

	/**
     * @brief Access element by linear index (vector style, 1-based).
     * @param n Linear index (starting at 1).
     * @return Reference to the matrix element.
     */
	double& operator()(const int n);

	/**
     * @brief Access element by row and column indices (1-based).
     * @param row Row index (starting at 1).
     * @param column Column index (starting at 1).
     * @return Reference to the matrix element.
     */
	double& operator()(const int row, const int column);

	/**
     * @brief Matrix addition.
     * @param m Matrix to add.
     * @return Reference to the resulting matrix.
     */
	Matrix& operator+(Matrix& m);

	/**
     * @brief Matrix subtraction.
     * @param m Matrix to subtract.
     * @return Reference to the resulting matrix.
     */
	Matrix& operator-(Matrix& m);

	/**
     * @brief Matrix multiplication.
     * @param m Matrix to multiply by.
     * @return Reference to the resulting matrix.
     */
	Matrix& operator*(Matrix& m);

	/**
     * @brief Matrix assignment.
     * @param m Matrix to assign from.
     * @return Reference to this matrix.
     */
	Matrix& operator=(Matrix& m);

	/**
     * @brief Matrix division (multiplication by the inverse).
     * @param m Matrix to divide by.
     * @return Reference to the resulting matrix.
     */
	Matrix& operator / (Matrix& m); //multiplicar por la inversa de la matriz

	/**
     * @brief Adds a scalar to every element of the matrix.
     * @param k Scalar value to add.
     * @return Reference to the resulting matrix.
     */
	Matrix& operator + (const double k); //suma escalar de cada elemento

	/**
     * @brief Subtracts a scalar from every element of the matrix.
     * @param k Scalar value to subtract.
     * @return Reference to the resulting matrix.
     */
	Matrix& operator - (const double k); //resta escalar de cada elemento

	/**
     * @brief Multiplies every element of the matrix by a scalar.
     * @param k Scalar multiplier.
     * @return Reference to the resulting matrix.
     */
	Matrix& operator * (const double k); //producto escalar de cada elemento

	/**
     * @brief Divides every element of the matrix by a scalar.
     * @param k Scalar divisor.
     * @return Reference to the resulting matrix.
     */
	Matrix& operator / (const double k); //división escalar de cada elemento

	// Methods

	/**
     * @brief Transposes the matrix.
     * @return Reference to the transposed matrix.
     */
	Matrix& transpose();

	/**
     * @brief Extracts a subvector from the matrix elements.
     * @param from Starting index (1-based).
     * @param to Ending index (1-based).
     * @return Reference to the extracted vector.
     */
	Matrix& extract_vector(const int from, const int to); //devuelve el subvector desde el elemento 'from' hasta 'to'

	/**
     * @brief Concatenates another vector to this vector.
     * @param v Vector to concatenate.
     * @return Reference to the combined vector.
     */
	Matrix& union_vector(Matrix& v); //une el vector, lo concatena

	/**
     * @brief Extracts a specific row as a vector.
     * @param row Row index (1-based).
     * @return Reference to the extracted row vector.
     */
	Matrix& extract_row(const int row); //devuelve la fila indicada (forma de vector)

	/**
     * @brief Extracts a specific column as a vector.
     * @param col Column index (1-based).
     * @return Reference to the extracted column vector.
     */
	Matrix& extract_column(const int col); //devuelve la columna indicada (en forma de vector)

	/**
     * @brief Assigns a vector to a specific row.
     * @param row Row index (1-based).
     * @param v Vector to assign (must match row length).
     * @return Reference to the updated matrix.
     */
	Matrix& assign_row(const int row, Matrix& v); //asigna el vector a la fila indicada (deben coincidir las longitudes)

	/**
     * @brief Assigns a vector to a specific column.
     * @param col Column index (1-based).
     * @param v Vector to assign (must match column length).
     * @return Reference to the updated matrix.
     */
	Matrix& assign_column(const int col, Matrix& v); //asigna el vector a la columna indicada (deben coincidir las longitudes, el vector debe estar en forma vector o vertical?)

	// Non-member operators

	/**
     * @brief Outputs the matrix to an output stream.
     * @param o Output stream.
     * @param m Matrix to output.
     * @return Reference to the output stream.
     */
	friend ostream &operator<<(ostream &o, Matrix& m);
};

// Operator overloading

/**
 * @brief Outputs the matrix to an output stream.
 * @param o Output stream.
 * @param m Matrix to output.
 * @return Reference to the output stream.
 */
ostream &operator<<(ostream &o, Matrix& m);

// Methods

/**
 * @brief Creates a zero matrix with specified rows and columns.
 * @param n_row Number of rows.
 * @param n_column Number of columns.
 * @return Reference to the zero matrix.
 */
Matrix& zeros(const int n_row, const int n_column);

/**
 * @brief Creates a zero vector of specified size.
 * @param v_size Size of the vector.
 * @return Reference to the zero vector.
 */
Matrix& zeros(const int v_size); //vector de ceros

/**
 * @brief Creates an identity matrix of size n x n.
 * @param n Size of the identity matrix.
 * @return Reference to the identity matrix.
 */
Matrix& eye(const int n); //matriz identidad de tamaño n

/**
 * @brief Computes the Euclidean norm of a vector.
 * @param v Vector whose norm is computed.
 * @return Euclidean norm of the vector.
 */
double norm(Matrix& v); //norma euclidiana del vector

/**
 * @brief Computes the dot product of two vectors.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return Dot product result.
 */
double dot(Matrix& v1, Matrix& v2); //producto escalar de vectores

/**
 * @brief Computes the cross product of two 3D vectors.
 * @param v1 First vector (length 3).
 * @param v2 Second vector (length 3).
 * @return Reference to the cross product vector.
 */
Matrix& cross(Matrix& v1, Matrix& v2); //producto vectorial de vectores de longitud 3

/**
 * @brief Computes the inverse of a matrix (if square and invertible).
 * @param m Matrix to invert.
 * @return Reference to the inverse matrix.
 */
Matrix& inv(Matrix& m); // devuelve la matriz inversa (si es cuadrada y con det!=0)

#endif