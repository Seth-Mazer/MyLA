# MyLA

A numerical linear algebra library written in  C++, with zero external dependencies.

Many libraries abstract away the math we use. 
MyLA was built to answer one question: What is actually under the abstraction?

---

#  Features

Here is what MyLA can do so far:

- Packed Householder QR Factorization of A
- Computation of Ax = B with LU Factorization via Outer Product and Partial Pivoting
- Matrix Addition / Subtraction
- Matrix Multiplication
- Matrix Augmentation
- Elementary Row Operations
- Transpose, Determinant, Trace, Norms, Normalization, Scaling
- Symmetry Checking
- Identity Matrix Generation
- Inner Products / Outer Products
- Vector Projections
- And more utility functions (e.g. getCol, isSquare, etc.)

---

# Design

MyLA uses a flattened row-major storage, chosen specifically for cache locality.
There are certain design choices made to maximize numerical stability, and minimize unnecessary memory allocations.
For example, in both LU and QR, MyLA overwrites A with packed factorizations, thus there is no additional L, U, Q, R or Householder vectors, stored in memory.
Everything needed out of, or from, the factorizations can be extracted, or implicitly calculated on the fly.

---
# Example #1: Solving Ax = B, via LU

```cpp
#include <MyLA/MyLA.h>
int main() {
    
    
    // Initializing a 3x3 Matrix, and populating its entries
    // We could also omit the data, however, then it would be initialized as zero matrix    
    myla::Matrix A(3, 3, {
        1, 2, 1,
        3, 2, 4,
        4, 1, 3,
    });


    // Let's now initialize our b vector
    myla::Matrix b(3,1, {
        5,
        13,
        11,
    });
    
    
    // Now, simply call solve() to solve the system!
    myla::Matrix x = solve(A,b);
    
    
    // Now we can Print x, our solution vector
    myla::print(x);
    
    // Prints
    // 1
    // 1
    // 2
    
    
    // Note, if we wanted to exercise, more "precise" control, we could instead do:
    // myla::LUDecomp LU_of_A = LU(A);
    // myla::Matrix Y = forSub(LU_of_A.LU, LU_of_A.P, b);
    // myla::Matrix x = backSub(LU_of_A.LU, Y);
    // myla::print(x);
    // This would yield the same results as solve(), without allocating a copy of A, as solve() passes a copy,
    // we could also do solve_e(), which is what solve() calls with the copy of A.
    
}
```
# Example #2: A Packed QR Factorization

```cpp
#include <MyLA/MyLA.h>
int main() {


    // Let's initialize another matrix A;
    myla::Matrix A(3, 3, {
        12, -51, 4,
        6, 167, -68,
        -4, 24, -41
    });


    // Let's now create a QR Decomp object, to store the results of the QR Factorization
    myla::QRDecomp QR_of_A = QR(A);

    // Now, we have a packed QR, along with a vector of Beta values, we will print both.
    myla::print(QR_of_A.QR);

    // Prints
    // -14	-21	14
    // 0.230769	-175	70
    // -0.153846	0.0555556	35

    myla::print(QR_of_A.Betas);

    // Prints
    //1.85714
    // 1.99385
    // 2

}
```



# Example #3: Matrix Multiplication

```cpp
#include <MyLA/MyLA.h>
int main() {
    
    
    // Lets initialize another matrix A;
    myla::Matrix A(3, 3, {
        1, 2, 3,
        4, 5, 6,
        7, 8, 9
    });
    
    
    // Lets now Initialize B, with the result of multiplying A by its transpose
    myla::Matrix B = myla::multiply(A, A.transpose());
    

    // Then print B to the console
    myla::print(B);
    
    // Prints
    // 14 32 50
    // 32 77 122
    // 50 122 194

}
```


# Build

---

Clone the repository:

```bash
git clone https://github.com/Seth-Mazer/MyLA.git
cd MyLA
```

Configure the build and compile:

```bash
cmake -S . -B build
cmake --build build --config Release
```

# Roadmap

---

## Version 0
- Primitive operations 

This version is intended as the absolute bare-bones release of MyLA, essentially, what operations do we need, in order to build future operations? For example, in order to build a function projecting a vector onto another vector, one would have to build the dot product before hand. Version 0 consists of most of the "elementary functions," if you will, of linear algebra, that must exist for the more complicated procedures.

---

## Version 1 (We are now here!)

- General Usability

To put it simply, version 1 aims to be able to do mostly "everything" you would need to do in practice, or, say, if someone did not want to do their linear algebra homework. This future release would include some of the more obvious things such as

- LU Factorization, and RREF, to solve equations in the form of Ax = b, and the like.
- QR Factorization
- Least Squares problems, and some Simple machine learning concepts
- Affine Geometry, such as, are these points within the convex hull S?

And more Quality of Life features like computing powers of matrices via exponentiation squaring.

For this version, one should think, "I could put together a project, and ship it. Will it be the most optimized fast solver in the world? No. But it will be fully functional, and it will work. That is the point."

---
## Version 2

- Numerical Accuracy
- Optimization

For this version, in the far future, this is where the accuracy of computational linear algebra will be evaluated, and improved. 
A large focus will be on numerical accuracy, such as stability, conditioning improvements,  rank determination, etc. 
For version 2, one should expect less new operations, but much, much more confidence and accuracy in computation. 

Along with more optimization, as I intend to overhaul the memory layout for improved performance

---


## Why did I build this?

Essentially, I never liked calling blackbox functions. I remember using pandas, and one of the first tutorials said to write blahblah.fit(). I was like, what does that do? As in, like why are we saying that? Whats the math? Etc.
Fast forward, I am now building my own fit.


## License

MIT License. See [LICENSE](LICENSE) for details.
