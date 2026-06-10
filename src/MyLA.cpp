//Includes
#include <MyLA/MyLA.h>
#include <iostream>


//-----------------------------------------------------------------
//-----------------------------------------------------------------
// Name Space
//-----------------------------------------------------------------
//-----------------------------------------------------------------

namespace myla {

    //-----------------------------------------------------------------
    //-----------------------------------------------------------------
    // Primitive Operations
    //-----------------------------------------------------------------
    //-----------------------------------------------------------------

    void print(const Matrix &A) {
        //Loop through matrix, printing each entry out
        for (size_t i = 0; i < A.m(); i++) {
            for (size_t j = 0; j < A.n(); j++) {
                std::cout << A(i, j) << "\t";
            }
            std::cout << "\n";
        }
    }



    Matrix add(const Matrix &A, const Matrix &B) {

        //Check dimension of each matrix
        if (A.n() != B.n() || A.m() != B.m()) {
            throw std::invalid_argument("Incompatible matrix dimensions.");
        }

        //Initializing C, to the dimension of A
        Matrix C(A.m(), A.n());

        //Looping through, populating C, with the addition of each entry
        for (size_t i = 0; i < A.m(); i++) {
            for (size_t j = 0; j < A.n(); j++) {
                C(i,j) = A(i,j) + B(i,j);
            }
        }


        return C;
    }

    Matrix sub(const Matrix &A, const Matrix &B) {

        //Check dimension of each matrix
        if (A.n() != B.n() || A.m() != B.m()) {
            throw std::invalid_argument("Incompatible matrix dimensions.");
        }

        //Initializing C, to the dimension of A
        Matrix C(A.m(), A.n());

        //Looping through, populating C, with the addition of each entry
        for (size_t i = 0; i < A.m(); i++) {
            for (size_t j = 0; j < A.n(); j++) {
                C(i,j) = A(i,j) - B(i,j);
            }
        }


        return C;
    }

    Matrix multiply(const Matrix&A, const Matrix&B) {

        //Checking dimensions are valid
        if (A.n() != B.m()) {
            throw std::invalid_argument("Incompatible matrix dimensions.");
        }

        // //Initializing matrix C, with dimensions R^(MxP)
        Matrix C(A.m(), B.n());

        //Multiplying via row-scaling (more explicit documentation is available @ https://github.com/Seth-Mazer/MiniCPP-Programs/blob/main/MatrixMultiplication/MatrixMultiplcation.cpp)
        //For each row in A
        for (size_t i = 0, rowA = A.m(); i < rowA; i++) {
            // Scale row q of B by A(i,q)
            for (size_t q = 0, rowB = B.m(); q < rowB; q++) {
                // Accumulate the scaled row into row i of C
                for (size_t j = 0; j < B.n(); j++) {
                    C(i,j) += A(i,q) * B(q,j);
                }
            }
        }

        return C;
    }

    Matrix transpose(const Matrix &A) {

        //Initializing AT, flipping from m x n, to n x m
        Matrix AT(A.n(), A.m());

        //A[i][j] = A^T[j][i] ==> we are simply assigning the "opposite" values
        for (size_t i = 0; i < A.m(); i++) {
            for (size_t j = 0; j <A.n(); j++) {
                AT(j,i) = A(i,j);
            }
        }

        return AT;
    }

    double dot(const Matrix &A, const Matrix &B) {

        //Check if n x 1 matrices were passed
        if (A.n() != 1 || B.n() != 1) {
            throw std::invalid_argument("Dot product not defined | tried to dot non-vectors");
        }

        if (A.m() != B.m()) {
            throw std::invalid_argument("Dot product not defined | vectors do not share the same dimension");
        }

        //Init dp at 0
        double dp = 0;

        //Multiplying the shared index of each vector together, and summing
        for (size_t i = 0; i < A.m(); i++) {
            dp += A(i,0) * B(i,0);
        }

        return dp;
    }

    Matrix outerProduct(const Matrix& A, const Matrix& B) {

        //Check if n x 1 matrices were passed
        if (A.n() != 1 || B.m() != 1) {
            throw std::invalid_argument("Outer product not defined | vectors are not in m x 1 and 1 x n form");
        }

        //Init OP
        Matrix oP(A.m(), B.n());

        for (size_t i = 0; i < A.m(); i++) {
            for (size_t j = 0; j < B.n(); j++) {
                oP(i,j) = A(i,0) *  B(0, j);
            }
        }

        return oP;

    }


    //-----------------------------------------------------------------
    //-----------------------------------------------------------------
    // Basic Operations
    //-----------------------------------------------------------------
    //-----------------------------------------------------------------

    Matrix proj(const Matrix &A, const Matrix &B) {
        // Calculating (uTv / vTv) * v
        double scalar = (dot(A,B) / dot(B,B));
        return B.scale(scalar);
    }

    Matrix augment(const Matrix &A, const Matrix &B) {
        //Check if b.m() = A.m()
        if (B.m() != A.m()) {
            throw std::invalid_argument("Tried to construct [A | b], b.m() != A.m()");
        }

        //Initializing AB (really A|B), with proper dims
        Matrix AB = {A.m(), (A.n() + B.n())};

        //Looping through AB
        for (size_t i = 0; i < AB.m(); i++) {
            for (size_t j = 0; j < AB.n(); j++) {
                //If we are "passed" A's cols, start entering from B
                if (j >= A.n()) {
                    //Offset j by the current A col index, to properly access Bs "starting point"
                    AB(i,j) = B(i,j - A.n());
                } else {
                    AB(i,j) = A(i,j);
                }
            }
        }

        return AB;
    }









}

