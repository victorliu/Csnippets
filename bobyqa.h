/*
This subroutine seeks the least value of a function of many variables,
by applying a trust region method that forms quadratic models by
interpolation. There is usually some freedom in the interpolation
conditions, which is taken up by minimizing the Frobenius norm of
the change to the second derivative of the model, beginning with the
zero matrix. The values of the variables are constrained by upper and
lower bounds. The arguments of the subroutine are as follows.

N must be set to the number of variables and must be at least two.
NPT is the number of interpolation conditions. Its value must be in
  the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not
  recommended.
Initial values of the variables must be set in X(1),X(2),...,X(N). They
  will be changed to the values that give the least calculated F.
For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper
  bounds, respectively, on X(I). The construction of quadratic models
  requires XL(I) to be strictly less than XU(I) for each I. Further,
  the contribution to a model from changes to the I-th variable is
  damaged severely by rounding errors if XU(I)-XL(I) is too small.
RHOBEG and RHOEND must be set to the initial and final values of a trust
  region radius, so both must be positive with RHOEND <= RHOBEG.
  Typically, RHOBEG should be about one tenth of the greatest
  expected change to a variable, while RHOEND should indicate the
  accuracy that is required in the final values of the variables. An
  error return occurs if any of the differences XU(I)-XL(I), I=1,...,N,
  is less than 2*RHOBEG.
The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
  amount of printing. Specifically, there is no output if IPRINT=0 and
  there is output only at the return if IPRINT=1. Otherwise, each new
  value of RHO is printed, with the best vector of variables so far and
  the corresponding value of the objective function. Further, each new
  value of F with its variables are output if IPRINT=3.
MAXFUN must be set to an upper bound on the number of calls of CALFUN.
The array W will be used for working space. Its length must be at least
  (NPT+5)*(NPT+N)+3*N*(N+5)/2.

SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set
F to the value of the objective function for the current values of the
variables X(1),X(2),...,X(N), which are generated automatically in a
way that satisfies the bounds given in XL and XU.
*/

void BOBYQA(
	int n, int npt,
	double *x, double *xl, double *xu,
	double (*func)(void *data, int n, const double *x), void *data,
	double rhobeg, double rhoend, int iprint, int maxfun,
	double *w
);
