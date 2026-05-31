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




    };


    //-----------------------------------------------------------------
    //-----------------------------------------------------------------
    // Primitive Operations
    //-----------------------------------------------------------------
    //-----------------------------------------------------------------

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

    //-----------------------------------------------------------------
    //-----------------------------------------------------------------
    // Basic Operations
    //-----------------------------------------------------------------
    //-----------------------------------------------------------------

    //Project A onto B
    Matrix proj(const Matrix &A, const Matrix &B);






}





