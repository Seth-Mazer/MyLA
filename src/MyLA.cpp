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


    double det(const Matrix &A) {
        //Call LU_Safe(), to return packed LU
        LUDecomp AD = LU_Safe(A);

        //Determine sign
        const double sign = (AD.swapCount % 2 == 0) ? 1.0 : -1.0;

        //Multiply down diagonal
        double det = 1.0;
        for (size_t j = 0; j < AD.LU.n(); j++) {
            det *= AD.LU(j,j);
        }

        return det * sign;
    }


    //-----------------------------------------------------------------
    //-----------------------------------------------------------------
    // Algorithms
    //-----------------------------------------------------------------
    //-----------------------------------------------------------------

    LUDecomp LU(Matrix &A) {

        //Init a matrix P, to keep track of permutations.
        Matrix P(A.m(), 1);

        //Init swapCount, to keep track of how many swaps occur, for det sign
        int swapCount = 0;

        //Init a constraint for how many diagonal entries are in A
        const size_t maxDiag = std::min(A.m(), A.n());

        //Epsilon

        //Main loop to iterate down pivot positions
        for (size_t j = 0; j < maxDiag; j++) {

            //Set up partial pivot, here, we are setting our current indices to be the largest pivot,
            //before iterating down the current column to potentially interchange
            size_t pivotRow = j;
            double largestPivot = std::fabs(A(j,j));

            //Search down the current column, comparing our current largest pivot to each entry beneath the col
            //If larger entry found, update largestPivot, and pivotRow.
            for (size_t q = j + 1; q < A.m(); q++) {
                if (std::fabs(A(q,j)) > largestPivot) {
                    largestPivot = std::fabs(A(q,j));
                    pivotRow = q;
                }
            }

            //Check if singular, if so throw error
            if (constexpr double epsilon = 1e-12; largestPivot < epsilon) {
                throw std::runtime_error("Matrix is singular, largestPivot < epsilon(1e-12)");
            }


            //Interchanging rows, incrementing swap count, and tracking row swaps with P
            if (j != pivotRow) {
                A.rowSwap(j,pivotRow);
                swapCount++;
            }
            P(j,0) = static_cast<double>(pivotRow);


            //Trailing Block Computation, here we essentially do R < R - mu, where m is our multipliers,
            //and u is the portion of our pivot row, to the immediate right of our pivot, essentially subtracting the outer matrix block
            //from the current sub matrix, then we do this process iteratively, hence the trailing
            for (size_t i = j + 1; i < A.m(); i++) {

                //Computing multipliers beneath each pivot
                //Here we are also transforming each entry under the pivot to the multiplier
                //Aka we are essentially creating L in place
                A(i,j) /= A(j,j);


                for (size_t col = j + 1; col < A.n(); col++) {
                    A(i,col) -= A(i,j) * A(j,col);
                }
            }

        }

        //Return the packed LU factorization, and P
        return {A, P, swapCount};
    }


    LUDecomp LU_Safe(const Matrix &A) {
        //Copying over A, to B, to run LU on B, keeping A intact
        Matrix B = A;
        return LU(B);
    }


    Matrix forSub(const Matrix &A, const Matrix &P, Matrix &B) {

        //Permutating B, i.e. b -> Pb
        for (size_t i = 0; i < B.m(); i++) {
            B.rowSwap(i, static_cast<size_t>(P(i,0)));
        }

        //init y
        Matrix Y(B.m(),1);

        //Solving yi + 1 = bi - L_ij(y_j) (omitted division, as its assumed the diagonal is 1)
        for (size_t i = 0; i < A.m(); i++) {
            double sum = 0.0;

            //Computing the sum portion
            for (size_t j = 0; j < i; j++) {
                sum += A(i,j) * Y(j,0);
            }

            //Computing each y_i
            Y(i,0) = B(i,0) - sum;

        }

        //Returning Y
        return Y;
    }


    Matrix backSub(const Matrix &A, const Matrix &Y) {

        //Init x
        Matrix X(Y.m(), 1);

        //Starting from the bottom, and looping up
        for (size_t i = A.n(); i-- > 0;) {
            double sum = 0.0;

            //Computing the sum portion
            for (size_t j = i + 1; j < A.n(); j++) {
                sum += A(i,j) * X(j,0);
            }

            //Computing each x_i
            X(i,0) = (Y(i,0) - sum)/ A(i,i);

        }

        //Returning X
        return X;
    }

    Matrix solve_e(Matrix &A, Matrix &B) {

        //Here I'm limiting solve() to square systems for now
        //I haven't yet implemented handling for underdetermined or overdetermined systems

        if (!A.isSquare())
            throw std::invalid_argument("Tried to solve non-square system | solve_e()");

        if (B.n() != 1)
            throw std::invalid_argument("Tried to solve Ax = b, where b is not a vector | solve_e()");

        if (B.m() != A.m())
            throw std::invalid_argument("Tried to solve Ax = b, dimension mismatch | solve_e()");

        //Call LU, and return X
        LUDecomp packed = LU(A);
        return backSub(packed.LU, forSub(packed.LU, packed.P, B));
    }


    Matrix solve(const Matrix& A, Matrix& B) {
        //Copy A to C.
        Matrix C = A;
        return solve_e(C, B);
    }





}

