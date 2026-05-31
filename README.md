# MyLA

A lightweight linear algebra library written in C++.

MyLA was built from scratch as an applied mathematics / numerical computation project. The library currently implements primitive matrix and vector operations using a custom Matrix class.
Many libraries abstract away the math we use.  MyLA was built for the purpose of answering one question, what is actually under the abstraction?

#  Features

Here is what MyLA can do so far:

- Matrix Addition / Subtraction
- Matrix Multiplication
- Matrix Augmentation
- Transpose, Trace, Norms, Scaling
- Symmetry Checking
- Identity Matrix Generation
- Dot Products
- Vector Projections
- And many more utility functions


# Example

```cpp
#include <iostream>
#include <MyLA/MyLA.h>

int main() {
    myla::Matrix A(3, 3, {
        1, 2, 3,
        4, 5, 6,
        7, 8, 9
    });

    myla::Matrix B = myla::multiply(A, A.transpose());
    myla::print(B);

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
