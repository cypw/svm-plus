
#ifndef _SOLVER_NU_H_
#define _SOLVER_NU_H_

#include "common.h"
#include "Solver.h"
#include "Kernel.h"

//
// Solver for nu-svm classification and regression
//
// additional constraint: e^T \alpha = constant
//
class Solver_NU : public Solver
{
public:
  Solver_NU() {}
  void Solve(int l, const QMatrix& Q, const double *p, const schar *y,
	     double *alpha, double Cp, double Cn, double eps,
	     SolutionInfo* si, int shrinking)
  {
    this->si = si;
    Solver::Solve(l,Q,p,y,alpha,Cp,Cn,eps,si,shrinking);
  }
private:
  SolutionInfo *si;
  int select_working_set(int &i, int &j);
  double calculate_rho();
  bool be_shrunk(int i, double Gmax1, double Gmax2, double Gmax3, double Gmax4);
  void do_shrinking();
};

#endif
