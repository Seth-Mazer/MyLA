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
            throw std::invalid_argument("Incompatible matrix dimensions. | add()");
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
            throw std::invalid_argument("Incompatible matrix dimensions. | sub()");
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

    Matrix multiply(const Matrix &A, const Matrix &B) {

        //Checking dimensions are valid
        if (A.n() != B.m()) {
            throw std::invalid_argument("Incompatible matrix dimensions. | multiply()");
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
            throw std::invalid_argument("Dot product not defined | tried to dot non-vectors | dot()");
        }

        if (A.m() != B.m()) {
            throw std::invalid_argument("Dot product not defined | vectors do not share the same dimension | dot()");
        }

        //Init dp at 0
        double dp = 0;

        //Multiplying the shared index of each vector together, and summing
        for (size_t i = 0; i < A.m(); i++) {
            dp += A(i,0) * B(i,0);
        }

        return dp;
    }

    Matrix outerProduct(const Matrix &A, const Matrix &B) {

        //Check if n x 1 matrices were passed
        if (A.n() != 1 || B.m() != 1) {
            throw std::invalid_argument("Outer product not defined | vectors are not in m x 1 and 1 x n form | outerProduct()");
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
            throw std::invalid_argument("Tried to construct [A | b], b.m() != A.m() | augment()");
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
                throw std::runtime_error("Matrix is singular, largestPivot < epsilon(1e-12) | LU()");
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
        Matrix X(A.n(), 1);

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


    Matrix solve(const Matrix &A, Matrix &B) {
        //Copy A to C.
        Matrix C = A;
        return solve_e(C, B);
    }


    double houseHolder(Matrix &A, size_t row, size_t col) {

        //Calculating the sign of the first component of x, then calculating alpha, +-||x||
        const double sign = (A(row,col) >= 0) ? 1.0 : -1.0;

        //Calculating the norm of the active column we are in
        double normSq = 0.0;
        for (size_t i = row; i < A.m(); i++) {
            normSq += A(i,col)*A(i,col);
        }

        //Calculating ||x||
        const double norm = std::sqrt(normSq);

        //Updating the first coord basis position of A
        A(row,col) += sign*norm;

        //Calculating the diagonal entry of R, which is just -sign*norm
        //Then determining if we might result in a near-singular matrix, due to rDiag being near zero, if so, throw
        //to avoid division by near zero problems
        double rDiag = -sign*norm;
        if (constexpr double epsilon = 1e-12; std::fabs(rDiag) < epsilon) {
            throw std::runtime_error("Matrix is near-rank deficient, a diagonal value of R is < epsilon(1e-12) | houseHolder()");
        }


        //Grabbing the first entry of the active column
        double Azero = A(row,col);

        //Scaling the active column
        for (size_t i = row; i < A.m(); i++) {
            A(i,col) *= (1.0/Azero);
        }


        //Computing the dot here 'manually', because we cannot pass a column vector to dot() here
        //Note this dP != normSq
        double dot = 0.0;

        //Offsetting i, to account for the diagonal of A, being the diagonal entries of R
        //This does not matter as, since we normalized the first entry, we can treat it as one,
        //and just add 1 to the dot, as 1^2 is just 1, we can add 1 to the sum of squares
        for (size_t i = row + 1; i < A.m(); i++) {
            dot += A(i, col) * A(i, col);
        }

        //Setting the head of the householder vector to be the diagonal entry of R
        //This is so when we start at j + 1, in the QR function to avoid destruction, we just
        //set the diagonal here, and implicitly treat the head of the houseHolder vector as 1
        A(row, col) = rDiag;

        //Computing and returning beta, which is just 2/(1 + xTx), note we add one to account for the offset
        return (2.0 / (1 + dot));
    }


    QRDecomp QR(Matrix &A) {

        //Create a maximum index, prevent indexing out of bounds
        const size_t maxDiag = std::min(A.m(),A.n());

        //Initializing a vector to store beta values, restricted to the maxDiag
        Matrix Betas(maxDiag,1);

        //Looping down the columns of A
            for (size_t j = 0; j < maxDiag; j++) {
                //Calling houseHolder(), to compute and store the houseHolder vectors **in** A,
                //and return beta
                double beta = houseHolder(A,j,j);
                Betas(j,0) = beta;

                //Looping through the sub matrix, to compute vT*A where A is the sub matrix,
                //and vT is the stored houseHolder vectors already in A
                for (size_t k = j + 1; k < A.n(); k++) {

                    //Computing vTa and setting the first value of the dot product to A(j,k)
                    //because we cant dot the first index of v, as its the diagonal of R, but since its one,
                    //we can just set it as itself
                    double vTa = A(j,k);
                    for (size_t i = j + 1; i < A.m(); i++) {
                        vTa += A(i,j) * A(i,k);
                    }


                    //We do the same thing for the first index of A as well, for analogous reasons
                    A(j,k) -= beta * vTa;

                    //Here we are applying an implicitly constructed householder matrix, to each column of A,
                    //in turn creating R, while also storing the householder vectors in the lower triangular portion
                    for (size_t i = j + 1; i < A.m(); i++) {
                        A(i,k) -= beta * A(i,j) * vTa;
                    }

                }

            }

        //Returning the packed QR
        return {A, Betas};
    }


    Matrix qTb(const Matrix &A, const Matrix &betas, const Matrix &b) {

        //Init a new matrix, qtb, because we're gonna assume the caller doesn't want their original data overwritten
        Matrix qtb = b;

        //Loop to construct qtb
        for (size_t j = 0; j < A.n(); j++) {

            //Here we initialize vTb, to the first indice of b, we do this because, the head
            //of the householder vector is implicitly 1. However, recall that in the packed factorization,
            //the diagonals belong to R
            double vTb = qtb(j,0);

            //Computing the rest of vTb, below the implicit head of 1
            for (size_t i = j + 1; i < A.m(); i++) {
                vTb += A(i,j) * qtb(i,0);
            }

            //Update the first index of qtb for analgous reasons to vtb
            qtb(j,0) -= betas(j,0) * vTb;

            //Loop to construct QtB
            for (size_t i = j + 1; i < A.m(); i++) {
                //Applying the reflectors to qtb
                qtb(i,0) -= (betas(j,0) * vTb * A(i,j));
            }
        }

       //return qtb
        return qtb;
    }

//Rx = qTb, thus its really Ux = y, so backsub, where X is the design matrix.
    Matrix leastSquares(const Matrix &A, const Matrix &y) {

        //Check if matrix is square or tall, and not wide, to prevent rank deficient errors in the backsub portion
        //Near-rank deficiency is handled in the householder function, thus tall matrices are allowed
        if (A.m() < A.n()) {
            throw std::invalid_argument(
            "\nWide matrix Passed, Refuse to Compute QR backSub | leastSquares() \n"
            "Try to make the matrix a square, or tall");
        }

        //Allocate a copy, to not destroy original matrix, then factor
        Matrix X = A;
        QRDecomp Xqr = QR(X);

        //init and Compute qtb
        Matrix qtb = qTb(Xqr.QR, Xqr.Betas, y);

        //Compute Rx = QtB, via backsub
        Matrix xhat = backSub(Xqr.QR, qtb);

        //Return LS solution
        return xhat;

    }



}

