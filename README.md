# MyLA

A lightweight linear algebra library written in C++.

MyLA was built from scratch as an applied mathematics / numerical computation project. The library currently implements primitive matrix and vector operations using a custom Matrix class.
Many libraries abstract away the math we use.  MyLA was built for the purpose of answering one question, what is actually under the abstraction?

#  Features

Here is what MyLA can do so far:

- Matrix Addition / Subtraction
- Matrix Multiplication
- Matrix Augmentation
- Elementary Row Operations
- Transpose, Trace, Norms, Normalization, Scaling
- Symmetry Checking
- Identity Matrix Generation
- Inner Products / Outer Products
- Vector Projections
- And more utility functions (e.g. getCol, isSquare, etc.)


# Example

```cpp
#include <MyLA/MyLA.h>

int main() {
    // Initialzing a 3x3 Matrix, and populating its entries
    // We could also omit the data, then it would be initialized as zero matrix
    myla::Matrix A(3, 3, {
        1, 2, 3,
        4, 5, 6,
        7, 8, 9
    });
    
    // Lets now Initialize B, with the result of multiplying A by its transpose, then print B to the console
    myla::Matrix B = myla::multiply(A, A.transpose());
    myla::print(B);

}
```
# Example #2

```cpp
#include <MyLA/MyLA.h>
#include <iostream>

int main() {
    // Lets create another matrix A
    myla::Matrix A(3, 3, {
        5, 7, 1,
        1, 4, 6,
        8, 2, 9
    });
    
    // And an n x 1 (vector) matrix b
    myla::Matrix b(3, 1, {
    5,
    6,
    7
    });

    // Let's Augment A with b
    myla::Matrix Ab = myla::augment(A, b);
    // Now Ab is:
    // 5       7       1       5
    // 1       4       6       6
    // 8       2       9       7

    // Now, let's swap row 1 and 2, bringing 2 to the top.
    Ab.rowSwap(0,1);
    // Now Ab is:
    // 1       4       6       6
    // 5       7       1       5
    // 8       2       9       7

    // Let's check if Ab is square, if not, print Ab's dimensions
    if (Ab.isSquare()) {
        std::cout << "Ab is square";
    } else {
        std::cout << "Ab is not square, it's a " << Ab.m() << " by " << Ab.n() << " matrix";
    }
    
    // We know Ab.isSquare() will return false, as we started with a 3x3 matrix
    // and augmented it with a 3x1 matrix. Thus we will print "Ab is not square" along with its dimensions
    
}
```



# Build
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

## Version 0
- Primitive operations 

This version is intended as the absolute bare-bones release of MyLA, essentially, what operations do we need, in order to build future operations? For example, in order to build a function projecting a vector onto another vector, one would have to build the dot product before hand. Version 0 consists of most of the "elementary functions," if you will, of linear algebra, that must exist for the more complicated procedures.

## Version 1

- General Usability

To put it simply, version 1 aims to be able to do mostly "everything" you would need to do in practice, or, say, if someone did not want to do their linear algebra homework. This future release would include some of the more obvious things such as

- LU Factorization, and RREF, to solve equations in the form of Ax = b, and the like.
- QR Factorization
- Least Squares problems, and some Simple machine learning concepts
- Affine Geometry, such as, are these points within the convex hull S?

And more Quality of Life features like computing powers of matrices via exponentiation squaring.

For this version, one should think, "I could put together a project, and ship it. Will it be the most optimized fast solver in the world? No. But it will be fully functional, and it will work. That is the point."

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
