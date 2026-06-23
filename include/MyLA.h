//Load once
#pragma once

//Includes
#include <vector>
#include <stdexcept>
#include <cmath>

//-----------------------------------------------------------------
//-----------------------------------------------------------------
// Name Space
//-----------------------------------------------------------------
//-----------------------------------------------------------------


//Namespace creation
namespace myla {

    //Matrix Class
    class Matrix {
    private:
        //Members
        std::vector<double> data;
        size_t rows, cols;

    public:
        //Constructor, which will create a "flattened" matrix, of size m by n, and initialize all values to 0.
        Matrix(size_t m, size_t n)
        : data(m * n, 0.0), rows(m), cols(n) {

            //Make sure, matrix dimensions are valid
            if (m == 0 || n == 0) {
                throw std::invalid_argument("Tried to create a matrix with a dimension being 0");
            }
        }

        //Constructor, which allows us to manually initialize the values of data
        Matrix(size_t m, size_t n, const std::vector<double>& values)
        : data(values), rows(m), cols(n) {
            //Make sure, matrix dimensions are valid
            if (m == 0 || n == 0) {
                throw std::invalid_argument("Tried to create a matrix with a dimension being 0");
            }

            //Make sure the amount of data provided actual matches the dimensions of the matrix
            if (values.size() != (m * n)) {
                throw std::invalid_argument("Tried to initialize more/less data then dimensions allow");
            }
        }


        //Function to return the dimensions, not allowing direct access
        size_t m() const {
            return rows;
        }

        //Function to return the dimensions, not allowing direct access
        size_t n() const {
            return cols;
        }


        //-----------------------------------------------------------------
        //-----------------------------------------------------------------
        // Operators
        //-----------------------------------------------------------------
        //-----------------------------------------------------------------

        //Operator for (), allows us to directly say A(i,j)
        double& operator()(size_t i, size_t j) {
            //Checking if, index is "out of bounds" in the sense that it's not a valid representation of what we actually want
            if (i >= rows || j >= cols) {
                throw std::out_of_range("Tried to access index out of bounds");
            }

            //return proper entry
            return data[i * cols + j];
        }

        //Operator for (), allows us to directly say A(i,j), read only
        const double& operator()(size_t i, size_t j) const {
            if (i >= rows || j >= cols) {
                throw std::out_of_range("Tried to access index out of bounds");
            }

            //return proper entry
            return data[i * cols + j];
        }

        //Operator for ==, allows us to compare self
        bool operator==(const Matrix &B) const {
            //Check dims
            if (rows != B.m() || cols != B.n()) {
                return false;
            }

            //Loop through, compare each entry, first error, return false
            for (size_t i = 0; i < rows; i++) {
                for (size_t j = 0; j < cols; j++) {
                    if ((*this)(i,j) != B(i,j)) {
                        return false;
                    }
                }
            }

            return true;
        }



        //-----------------------------------------------------------------
        //-----------------------------------------------------------------
        // Self Operations
        //-----------------------------------------------------------------
        //-----------------------------------------------------------------

        //Return tr(a) the sum of all diagonal entries
        double trace() const {
            //Verify A is square
            if (rows != cols) {
                throw std::invalid_argument("Matrix is not square.");
            }

            //initializing trace var
            double tr = 0;

            //Summing the diagonal entries of A
            for (size_t i = 0; i < rows; i++) {
                tr += (*this)(i,i);
            }

            return tr;
        }

        //Check if self is square
        bool isSquare() const {
            //return boolean answer
            return rows == cols;
        }


        //Return the norm, ||v||
        double norm() const {
            //Check if n x 1
            if (cols != 1) {
                throw std::invalid_argument("Tried to take norm, matrix is not n x 1. Use normf() instead");
            }

            double sumsquares = 0;

            //Sum all the square components of the vector,
            for (size_t i = 0; i < rows; i++) {
                sumsquares += ((*this)(i,0) * (*this)(i,0));
            }

            return std::sqrt(sumsquares);
        }

        //Normalize a vector, A/||A||, returns a normalized copy of A
        Matrix normalized() const {
            //Check if n x 1
            if (cols != 1) {
                throw std::invalid_argument("Tried to normalize, matrix is not n x 1");
            }

            //Grab norm
            double nrm = norm();

            //vNM = vector normalized
            if (nrm == 0){
                throw std::invalid_argument("Tried to divide by zero | normalize()");
            }

            return scale(1/nrm);
        }

        //Normalize a vector, A/||A|| >> NOTE: THIS WILL MUTATE A! to return a copy, use normalized()
        void normalize() {
            *this = normalized();
        }



        //Return Frobenius norm ||A||_F
        double normf() const {

            double sumsquares = 0;

            //Sum the square of every entry in the matrix
            for (size_t i = 0; i < rows; i++) {
                for (size_t j = 0; j < cols; j++) {
                    sumsquares += ((*this)(i,j) * (*this)(i,j));
                }
            }

            return std::sqrt(sumsquares);
        }


        //Scale by r, and return scaled copy
        Matrix scale(const double r) const {

            //Initializing rA
            Matrix rA(rows, cols);

            //Scaling each entry of A by r, thus computing rA
            for (size_t i = 0; i < rows; i++) {
                for (size_t j = 0; j < cols; j++) {
                    rA(i,j) = (r * (*this)(i,j));
                }
            }

            return rA;
        }

        //Return identity matrix
        Matrix identity() const {

            //If not square, throw error
            if (!isSquare()) {
                throw std::invalid_argument("Tried to return identity of non-square matrix");
            }

            //Initializing identity matrix
            Matrix I(rows, cols);

            //Populate diagonal with 1's
            for (size_t i = 0; i < rows; i++) {
                I(i,i) = 1;
            }
            return I;

        }

        //Read the free function version for comments
        Matrix transpose() const {
            Matrix AT(cols, rows);
            for (size_t i = 0; i < rows; i++) {
                for (size_t j = 0; j < cols; j++) {
                    AT(j,i) = (*this)(i,j);
                }
            }

            return AT;
        }


        // Check if matrix is symmetric
         bool isSymmetric() const {
            //Check the obvious first
            if (!isSquare()) {
                return false;
            }

            //Return bool
            return (*this) == transpose();
         }

        //Set col, sets the indices of a col, to given values
        void setCol(const size_t j, const double value) {

            //looping down the jth column, assigning value
            for (size_t i = 0; i < rows; i++) {
                (*this)(i,j) = value;
            }

        }



        //Get Col, returns the entire column at position n
        Matrix getCol(const size_t n) const {
            //Check if index exists
            if (n >= cols) {
                throw std::invalid_argument("Tried to get column, column index out of bounds");
            }

            //Init col
            Matrix col = {rows, 1};

            //Loop through the col at pos. n, assign values to col
            for (size_t i = 0; i < rows; i++) {
                col(i,0) = (*this)(i,n);
            }

            //Return col
            return col;
        }


        //Get row, returns the entire row at position m
        Matrix getRow(const size_t m) const {
            //Check if index exists
            if (m >= rows) {
                throw std::invalid_argument("Tried to get row, row index out of bounds");
            }

            //Init row
            Matrix row = {1, cols};

            //Loop through the row at pos. n, assign values to row
            for (size_t j = 0; j < cols; j++) {
                row(0,j) = (*this)(m,j);
            }

            //Return row
            return row;
        }

        //-----------------------------------------------------------------
        //-----------------------------------------------------------------
        // Elementary Row Operations
        //-----------------------------------------------------------------
        //-----------------------------------------------------------------

        //interchange two rows
        void rowSwap(const size_t curLoc, const size_t swapLoc) {
            //Check if interchange position exists
            if (swapLoc >= rows || curLoc >= rows) {
                throw std::invalid_argument("Tried to interchange rows, indices is out of bounds");
            }

            //Return if row indices are equal
            if (curLoc == swapLoc) {
                return;
            }

            //Overwriting the current row, with the new row, and vice versa
            for (size_t j = 0; j < cols; j++) {
                double temp = (*this)(curLoc, j);
                (*this)(curLoc, j) = (*this)(swapLoc,j);
                (*this)(swapLoc,j) = temp;
            }

        }

        //Scale a row >> NOTE: THIS WILL MUTATE A! To scale a matrix/vector and return a copy, use scale(r)
        void rowScale(const size_t row, const double r) {
            //Check if row pos exists
            if (row >= rows) {
                throw std::invalid_argument("Tried to scale row, row index is out of bounds");
            }

            //Scaling row
            for (size_t j = 0; j < cols; j++) {
                (*this)(row,j) *= r;
            }

        }

        //Replace one row by the sum of itself and a scalar multiple of another row>> NOTE: THIS WILL MUTATE A! To add, rows and return a row, use myla::add(getRow(0), getRow(1));
        void rowAddScaled(const size_t firstRow, const double r, const size_t secondRow) {
            if (firstRow >= rows || secondRow >= rows) {
                throw std::invalid_argument("Tried to add rows, row index is out of bounds");
            }

            //Check if r is not 0
            if (r == 0) {
                return;
            }


            for (size_t j = 0; j < cols; j++) {
                (*this)(firstRow,j) += ( (r * (*this)(secondRow,j)));

            }
        }


    };


    //-----------------------------------------------------------------
    //-----------------------------------------------------------------
    // Primitive Operations
    //-----------------------------------------------------------------
    //-----------------------------------------------------------------

    //Print the matrix to the console
    void print(const Matrix &A);

    //Add two matrices A and B
    Matrix add(const Matrix &A, const Matrix &B);

    //Subtract two matrices A and B
    Matrix sub(const Matrix &A, const Matrix &B);

    //Multiply two matrices A and B
    Matrix multiply(const Matrix&A, const Matrix &B);

    //Transpose a matrix A
    Matrix transpose(const Matrix &A);

    //Compute dot product
    double dot(const Matrix &A, const Matrix &B);

    //Compute outer product >> NOTE: B must be 1 x n. I.e, if dealing with n x 1 matrices, simply pass yourSecondVec.transpose()
    Matrix outerProduct(const Matrix &A, const Matrix &B);

    //-----------------------------------------------------------------
    //-----------------------------------------------------------------
    // Basic Operations
    //-----------------------------------------------------------------
    //-----------------------------------------------------------------

    //Project A onto B
    Matrix proj(const Matrix &A, const Matrix &B);

    //Augment, Construct [ A | b ]
    Matrix augment(const Matrix &A, const Matrix&B);

    //Returns the determinant of a matrix -> Will allocate memory to compute, for inplace, use det_e
    double det(const Matrix &A);

    //-----------------------------------------------------------------
    //-----------------------------------------------------------------
    // Algorithms
    //-----------------------------------------------------------------
    //-----------------------------------------------------------------

    //Struct that stores LU, Permutation Vector (n x 1 matrix), and swapCount
    struct LUDecomp {
        Matrix LU;
        Matrix P;
        int swapCount = 0;
    };

    //Computes LU in place, returns packed LU factorization, and P >> NOTE: THIS WILL MUTATE A, USE LU_Safe IF NECESSARY
    LUDecomp LU(Matrix &A);

    //Copies A to B, then calls LU, keeping A intact
    LUDecomp LU_Safe(const Matrix &A);

    //Solving a Lower triangular system, via forward sub
    Matrix forSub(const Matrix &A, const Matrix &P, Matrix &B);

    //Solving an Upper triangular system, via back sub
    Matrix backSub(const Matrix &A, const Matrix &Y);

    //Solves Ax=b using an LU factorization with forSub/backSub >> NOTE: THIS WILL MUTATE A, use solve() to pass a copy
    Matrix solve_e(Matrix &A, Matrix &B);

    //Solves Ax=b using an LU factorization with forSub/backSub >> NOTE: This will allocate a copy of A.
    //To solve without allocating a copy of A, use solve_e()
    Matrix solve(const Matrix &A, Matrix &B);

    //Struct that stores the packed QR matrix, along with the beta values from the houseHolder computation;
    struct QRDecomp {
        Matrix QR;
        Matrix Betas;
    };


    //Factorizes QR in place, returns packed QR factorization >> NOTE: THIS WILL MUTATE A
    QRDecomp QR(Matrix &A);

    //Mutates a sub column of A into a householder vector, computes and sets diagonal R entry, and computes and returns beta >> NOTE: Not really intended for use outside QR
    double houseHolder(Matrix &A, size_t row, size_t col);

    //Compute QTB, i.e. Rx = QtB, compute the QtB component, using the householder vectors
    Matrix qTb(const Matrix &A, const Matrix &B);

    //Find a least squares solution to an overdetermined full rank design matrix X, with observation vector y
    Matrix leastSquares(const Matrix &X, const Matrix &y);


}





