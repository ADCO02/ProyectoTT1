#include "..\include\matrix.h"

// constructores

// vector
Matrix::Matrix(const int v_size)
{
	if (v_size <= 0)
	{
		cout << "Matrix create: error in v_size\n";
		exit(EXIT_FAILURE);
	}

	this->n_row = 1;
	this->n_column = v_size;
	this->data = (double **)malloc(n_row * sizeof(double *));

	if (this->data == NULL)
	{
		cout << "Matrix create: error in data\n";
		exit(EXIT_FAILURE);
	}

	this->data[0] = (double *)calloc(v_size, sizeof(double));
}

// matriz
Matrix::Matrix(const int n_row, const int n_column)
{
	if (n_row <= 0 || n_column <= 0)
	{
		cout << "Matrix create: error in n_row/n_column\n";
		exit(EXIT_FAILURE);
	}

	this->n_row = n_row;
	this->n_column = n_column;
	this->data = (double **)malloc(n_row * sizeof(double *));

	if (this->data == NULL)
	{
		cout << "Matrix create: error in data\n";
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < n_row; i++)
	{
		this->data[i] = (double *)malloc(n_column * sizeof(double));
	}
}

// operartor
double &Matrix::operator()(const int n)
{
	if (n <= 0 || n > this->n_row * this->n_column)
	{
		cout << "Matrix get: error in row/column\n";
		exit(EXIT_FAILURE);
	}

	return this->data[(n - 1) / this->n_column][(n - 1) % this->n_column];
}

double &Matrix::operator()(const int row, const int column)
{
	if (row <= 0 || row > this->n_row || column <= 0 || column > this->n_column)
	{
		cout << "Matrix get: error in row/column\n";
		exit(EXIT_FAILURE);
	}

	return this->data[row - 1][column - 1];
}

Matrix &Matrix::operator+(Matrix &m)
{
	if (this->n_row != m.n_row || this->n_column != m.n_column)
	{
		cout << "Matrix sum: error in n_row/n_column\n";
		exit(EXIT_FAILURE);
	}

	Matrix *m_aux = new Matrix(this->n_row, this->n_column);

	for (int i = 1; i <= this->n_row; i++)
	{
		for (int j = 1; j <= this->n_column; j++)
		{
			(*m_aux)(i, j) = (*this)(i, j) + m(i, j);
		}
	}

	return *m_aux;
}

Matrix &Matrix::operator-(Matrix &m)
{
	if (this->n_row != m.n_row || this->n_column != m.n_column)
	{
		cout << "Matrix sub: error in n_row/n_column\n";
		exit(EXIT_FAILURE);
	}

	Matrix *m_aux = new Matrix(this->n_row, this->n_column);

	for (int i = 1; i <= this->n_row; i++)
	{
		for (int j = 1; j <= this->n_column; j++)
		{
			(*m_aux)(i, j) = (*this)(i, j) - m(i, j);
		}
	}

	return *m_aux;
}

Matrix &Matrix::operator*(Matrix &m)
{
	if (this->n_column != m.n_row)
	{
		cout << "Matrix sub: error in n_row/n_column\n";
		exit(EXIT_FAILURE);
	}

	Matrix *m_aux = new Matrix(this->n_row, m.n_column);

	for (int i = 1; i <= this->n_row; i++)
	{
		for (int j = 1; j <= m.n_column; j++)
		{
			int suma = 0;
			for (int s = 1; s <= this->n_column; s++)
			{
				suma += (*this)(i, s) * m(s, j);
			}
			(*m_aux)(i, j) = suma;
		}
	}

	return *m_aux;
}

Matrix &Matrix::operator=(Matrix &m)
{
	if (this == &m)
		return *this;

	for (int i = 0; i < this->n_row; ++i)
	{
		free(this->data[i]);
	}
	free(this->data);

	this->n_row = m.n_row;
	this->n_column = m.n_column;
	this->data = (double **)malloc(n_row * sizeof(double *));

	if (this->data == NULL)
	{
		cout << "Matrix create: error in data\n";
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < n_row; i++)
	{
		this->data[i] = (double *)malloc(n_column * sizeof(double));
	}

	for (int i = 0; i < n_row; ++i)
	{
		for (int j = 0; j < n_column; ++j)
		{
			this->data[i][j] = m.data[i][j];
		}
	}

	return *this;
}

Matrix &Matrix::operator+(const double k)
{
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);

	for (int i = 1; i <= this->n_row; i++)
	{
		for (int j = 1; j <= this->n_column; j++)
		{
			(*m_aux)(i, j) = (*this)(i, j) + k;
		}
	}

	return *m_aux;
}

Matrix &Matrix::operator-(const double k)
{
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);

	for (int i = 1; i <= this->n_row; i++)
	{
		for (int j = 1; j <= this->n_column; j++)
		{
			(*m_aux)(i, j) = (*this)(i, j) - k;
		}
	}

	return *m_aux;
}

Matrix &Matrix::operator*(const double k)
{
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);

	for (int i = 1; i <= this->n_row; i++)
	{
		for (int j = 1; j <= this->n_column; j++)
		{
			(*m_aux)(i, j) = (*this)(i, j) * k;
		}
	}

	return *m_aux;
}

Matrix &Matrix::operator/(const double k)
{
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);

	for (int i = 1; i <= this->n_row; i++)
	{
		for (int j = 1; j <= this->n_column; j++)
		{
			(*m_aux)(i, j) = (*this)(i, j) / k;
		}
	}

	return *m_aux;
}

Matrix &Matrix::transpose()
{
	Matrix *m_aux = new Matrix(this->n_column, this->n_row);

	for (int i = 1; i <= this->n_row; i++)
	{
		for (int j = 1; j <= this->n_column; j++)
		{
			(*m_aux)(j, i) = (*this)(i, j);
		}
	}

	return *m_aux;
}


/*Matrix& Matrix::inv() {
    if (this->n_row != this->n_column) {
        cout << "Matrix inv: not a square matrix\n";
        exit(EXIT_FAILURE);
    }

    int n = this->n_row;
    Matrix* A = new Matrix(n, 2 * n);

    // Construir [A | I]
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            (*A)(i, j) = (*this)(i, j);
			(*A)(i, n + j) = (i == j) ? 1.0 : 0.0;
        }
    }

    // Gauss-Jordan
    for (int i = 1; i <= n; i++) {
		//elegimos el pivote de mayor valor
		int pivote = i;

		for(int j=i;j<=n;j++){
			if(abs((*A)(j,i))>abs((*A)(pivote,i))){
				pivote=j;
			}
		}
        double diag = (*A)(pivote, i);
        if (fabs(diag) < 1e-10) {
            cout << "Matrix inv: singular matrix (det = 0)\n";
            exit(EXIT_FAILURE);
        }

		//intercambiar la fila pivote si es necesario
		if(pivote!=i){
			Matrix& vaux = (*A).extract_row(pivote);
			(*A).assign_row(pivote,(*A).extract_row(i));
			(*A).assign_row(i,vaux);
			//liberar vaux??
		}

        // Normalizar fila i
		(*A).assign_row(i,(*A).extract_row(i)/diag);
        for (int j = i; j <= 2 * n; j++) {
            (*A)(i, j) /= diag;
        }

        // Eliminar otras filas
        for (int k = 1; k <= n; k++) {
            if (k == i) continue;
            double factor = (*A)(k, i);
            for (int j = 1; j <= 2 * n; j++) {
                (*A)(k, j) -= factor * (*A)(i, j);
            }
        }
    }

    // Extraer parte derecha (inversa)
    Matrix* inverse = new Matrix(n, n);
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            (*inverse)(i, j) = (*A)(i, n + j);
        }
    }

    delete A; // liberamos memoria auxiliar

    return *inverse;
}*/

Matrix& Matrix::extract_vector(const int from, const int to) {
	if (this->n_row != 1) {
        cout << "Matrix extract_vector: matrix must be a vector\n";
        exit(EXIT_FAILURE);
    }

    if (from <= 0 || to > this->n_row * this->n_column || from > to) {
        cout << "Matrix extract_vector: invalid from/to indices\n";
        exit(EXIT_FAILURE);
    }

    int sub_size = to - from + 1;
    Matrix* subvector = new Matrix(sub_size);

    int index = 0;
    for (int i = from; i <= to; ++i) {
        (*subvector)(++index) = (*this)(i);
    }

    return *subvector;
}

Matrix& Matrix::union_vector(Matrix &v) {
    if (this->n_row != 1) {
        cout << "Matrix union_vector: first matrix must be a vector\n";
        exit(EXIT_FAILURE);
    }
    if (v.n_row != 1) {
        cout << "Matrix union_vector: second matrix must be a vector\n";
        exit(EXIT_FAILURE);
    }

    int total_size = this->n_column + v.n_column;
    Matrix* result = new Matrix(total_size);

    int index = 0;
    for (int i = 1; i <= this->n_row * this->n_column; ++i) {
        (*result)(++index) = (*this)(i);
    }
    for (int i = 1; i <= v.n_row * v.n_column; ++i) {
        (*result)(++index) = v(i);
    }

    return *result;
}


Matrix& Matrix::extract_row(const int row) {
    if (row < 1 || row > this->n_row) {
        cout << "Matrix extract_row: row out of bounds\n";
        exit(EXIT_FAILURE);
    }

    Matrix* v = new Matrix(this->n_column);

    for (int j = 1; j <= this->n_column; ++j) {
        (*v)(j) = (*this)(row, j);
    }

    return *v;
}

Matrix& Matrix::extract_column(const int col) {
    if (col < 1 || col > this->n_column) {
        cout << "Matrix extract_column: column out of bounds\n";
        exit(EXIT_FAILURE);
    }

    Matrix* v = new Matrix(this->n_row);

    for (int i = 1; i <= this->n_row; ++i) {
        (*v)(i) = (*this)(i, col);
    }

    return *v;
}

Matrix& Matrix::assign_row(const int row, Matrix& v) {
    if (row < 1 || row > this->n_row) {
        cout << "Matrix assign_row: row out of bounds\n";
        exit(EXIT_FAILURE);
    }

    if (v.n_row != 1 || v.n_column != this->n_column) {
        cout << "Matrix assign_row: vector size mismatch\n";
        exit(EXIT_FAILURE);
    }

    for (int j = 1; j <= this->n_column; ++j) {
        (*this)(row, j) = v(j);
    }

    return *this;
}

Matrix& Matrix::assign_column(const int col, Matrix& v) {
    if (col < 1 || col > this->n_column) {
        cout << "Matrix assign_column: column out of bounds\n";
        exit(EXIT_FAILURE);
    }

    if (v.n_row != 1 || v.n_column != this->n_row) {
        cout << "Matrix assign_column: vector size mismatch\n";
        exit(EXIT_FAILURE);
    }

    for (int i = 1; i <= this->n_row; ++i) {
        (*this)(i, col) = v(i);
    }

    return *this;
}

ostream &operator<<(ostream &o, Matrix &m)
{
	for (int i = 1; i <= m.n_row; i++)
	{
		for (int j = 1; j <= m.n_column; j++)
			printf("%5.20lf ", m(i, j));
		o << "\n";
	}

	return o;
}

Matrix &zeros(const int n_row, const int n_column)
{
	Matrix *m_aux = new Matrix(n_row, n_column);

	for (int i = 1; i <= n_row; i++)
	{
		for (int j = 1; j <= n_column; j++)
		{
			(*m_aux)(i, j) = 0;
		}
	}

	return (*m_aux);
}

Matrix &zeros(const int v_size)
{
	Matrix *m_aux = new Matrix(v_size);

	for (int i = 1; i <= v_size; i++)
	{
		(*m_aux)(i) = 0;
	}

	return (*m_aux);
}

Matrix& eye(const int n)
{
	Matrix *m_aux = new Matrix(n, n);

	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= n; j++)
		{
			(*m_aux)(i, j) = (i==j) ? 1 : 0;
		}
	}

	return (*m_aux);
}

double norm(Matrix &v){
	if (v.n_row != 1) {
        cout << "Matrix norm: matrix must be a vector\n";
        exit(EXIT_FAILURE);
    }

	double result = 0;
	for (int i = 1; i <= v.n_row; i++)
	{
		result += v(i)*v(i);
	}

	return sqrt(result);
	
}