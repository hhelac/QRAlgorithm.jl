# QR Algorithm for Matrix Eigenvalues (with Hessenberg Reduction)

A Julia implementation of the **QR Algorithm** used to compute the eigenvalues of square matrices. To optimize efficiency, the matrix is first reduced to its **Upper Hessenberg form** using orthogonal Householder transformations, making subsequent QR iterations significantly faster.

## Features
- **Hessenberg Reduction:** Uses Householder reflections to reduce computational complexity from $O(n^3)$ per iteration to $O(n^2)$.
- **Shift & Deflation Mechanisms:** Implemented to accelerate convergence and efficiently isolate computed eigenvalues.
- **Accuracy Verification:** Results are rigorously tested and validated against reference linear algebra frameworks.
- **Language:** Written completely in **Julia**, leveraging its high-performance numerical computing capabilities.

## About
This project was developed as a research seminar paper for the *Numerical Algorithms* course under the mentorship of Prof. Dr. Željko Jurić at the University of Sarajevo, Faculty of Electrical Engineering (February 2026).

*The full theoretical paper and mathematical analysis are available in the repository as a PDF.*
