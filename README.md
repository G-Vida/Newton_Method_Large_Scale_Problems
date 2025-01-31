# Newton Method with Finite differences (from scratch)
## Project for "Numerical methods for Large Scale Problems" (Course 2022/2023)


Collaboration with: https://github.com/GabPir

<br>

The main issue of optimization is the minimization of functions. There exist several techniques and many variants of them, like the popular Newton’s method in the original form and a modified form using the finite differences for the approximation of the gradient vector and the Hessian matrix.

<br>

### Pt1: Low-dimensional case
The goal is to minimize Rosenbrock's function, starting from two generic initial points:
<br>

![immagine](https://github.com/user-attachments/assets/2dac2735-cc07-4683-9894-cb0b02bf0567)

As a minimization method, the Newton method was applied, where the descent direction is iteratively computed by solving the linear system:

![immagine](https://github.com/user-attachments/assets/f4cb1f94-adcc-4daf-8890-df875e796b60)

The Newton method with exact resolution was compared with the one using approximations of the gradient and Hessian matrix, as well as a matrix-free implementation to avoid storing the Hessian. The results (in terms of gradient norm) were:

![immagine](https://github.com/user-attachments/assets/b41cc44a-f93b-44e0-a549-d9844dd2fb69)
![immagine](https://github.com/user-attachments/assets/df5f0012-b2ee-4431-b8fa-e0d9e3260747)

<br><br>


###Pt2: High-dimensional case
<br>
The following functions were considered to test the algorithms for non-trivial problems:

![immagine](https://github.com/user-attachments/assets/786d5ab8-fe86-41ce-932c-06f7d022893f)

The results obtained were:
<br>
![immagine](https://github.com/user-attachments/assets/8392bd44-46ca-484d-bc95-caf7c383922b)
<br>





### Additional aspects:
- In large-scale problems, we usually can’t directly solve the linear system and so we can use an iterative method like the Conjugate gradient method to find an approximated solution (pcg on Matlab). The iterative solver often requires many steps to achieve a good result and so to reduce the number of iterations it can be implemented the preconditioning of the system coefficient matrix, like the incomplete Cholesky factorization.
- Another big issue that comes into play when dealing with large-scale problems is the computation and storage of the gradient vector and the Hessian matrix. For this reason, it might be required to find approximations that call only for function evaluations at the given points. In the solver, all the finite differences formulas and the Matrix free implementation are implemented.
- In the case of nonlinear functions, we have frequently to deal with functions that aren’t convex in all the domains. To handle this, the Hessian matrix at a given point can be suitably modified to become a positive definitive one if it isn’t jet.
- Once the Newton descent direction is fixed, we would like to compute an optimal value of the step length 𝛼 that gives us a good decrease at each iteration. The Armijo condition will at each step guarantee a sufficient decrease along the descent direction.









