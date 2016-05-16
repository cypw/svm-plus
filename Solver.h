
#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "common.h"
#include "Kernel.h"
#include "solve_linear_system.h"

// An SMO and conjugate SMO algorithms for SVM
// Solves:
//
//	min 0.5(\alpha^T Q \alpha) + p^T \alpha
//
//		y^T \alpha = \delta
//		y_i = +1 or -1
//		0 <= alpha_i <= Cp for y_i = 1
//		0 <= alpha_i <= Cn for y_i = -1
//
// Given:
//
//	Q, p, y, Cp, Cn, and an initial feasible point \alpha
//	l is the size of vectors and matrices
//	eps is the stopping tolerance
//
// solution will be put in \alpha, objective value will be put in obj
//
class Solver {
public:
  Solver(int optimizer_=-1) {
    if(optimizer_==-1)
      conjugate = false;
    else {
      conjugate = true;
      max_depth =  optimizer_; 
      A = new double*[max_depth+1];
      for(int i=0; i<=max_depth; i++)
	A[i] = new double[max_depth+1];
      b = new double[max_depth+1];
    }
    fprintf(stdout,"conjugate=%d\n",conjugate);
    fflush(stdout);

    sizeof_double = sizeof(double);
    sizeof_char = sizeof(char);
    sizeof_int = sizeof(int);
  };
  virtual ~Solver() {
    if(conjugate) {
      for(int i=0; i<=max_depth; i++)
	delete[] A[i];
      delete[] A;
      delete[] b;
    }
  };

  void Solve(int l, const QMatrix& Q, const double *p_, const schar *y_,
	     double *alpha_, double Cp, double Cn, double eps,
	     SolutionInfo* si, int shrinking);

  void Solve_cg(int l, const QMatrix& Q, const double *p_, const schar *y_,
		double *alpha_, double Cp, double Cn, double eps,
		SolutionInfo* si, int shrinking); 

protected:
  int active_size;
  schar *y;
  double *G;		// gradient of objective function
  double **G_cg;
  enum { LOWER_BOUND, UPPER_BOUND, FREE };
  char *alpha_status;	// LOWER_BOUND, UPPER_BOUND, FREE
  char **alpha_status_cg;
  double *alpha;
  double **alpha_cg;
  const QMatrix *Q;
  const Qfloat *QD;
  double eps;
  double Cp,Cn;
  double *p;
  int *active_set;           // maps permuted indices to the original ones
  double *G_bar;		// gradient, if we treat free variables as 0
  double **G_bar_cg;
  int l;
  bool unshrink;	// XXX
  bool conjugate;
  int *work;
  int max_depth;

  double get_C(int i)
  {
    return (y[i] > 0)? Cp : Cn;
  }
  void update_alpha_status(int i)
  {
    if(alpha[i] >= get_C(i)-1e-12)
      alpha_status[i] = UPPER_BOUND;
    else if(alpha[i] <= 1e-12)
      alpha_status[i] = LOWER_BOUND;
    else alpha_status[i] = FREE;
  }
  void update_alpha_status_cg(int i, int depth)
  {
    if(alpha_cg[depth][i] >= get_C(i)-1e-12)
      alpha_status_cg[depth][i] = UPPER_BOUND;
    else if(alpha_cg[depth][i] <= 1e-12)
      alpha_status_cg[depth][i] = LOWER_BOUND;
    else alpha_status_cg[depth][i] = FREE;
  }
  bool is_upper_bound(int i) { return alpha_status[i] == UPPER_BOUND; }
  bool is_upper_bound_cg(int i, int depth) { return alpha_status_cg[depth][i] == UPPER_BOUND; }
  bool is_lower_bound(int i) { return alpha_status[i] == LOWER_BOUND; }
  bool is_lower_bound_cg(int i, int depth) { return alpha_status_cg[depth][i] == LOWER_BOUND; }
  bool is_free(int i) { return alpha_status[i] == FREE; }
  bool is_free_cg(int i, int depth) { return alpha_status_cg[depth][i] == FREE; }
  void swap_index(int i, int j);
  void reconstruct_gradient();
  int generate_direction(double *u, int n, bool new_working_set);
  virtual int wss_first_order(double* u, double& lambda_star);
  virtual int select_working_set_hmg(double* u, double& lambda_star, double& gain);
  virtual int select_working_set_hmg2(double* u, double& lambda_star, double& gain);
  virtual int select_working_set(int &i, int &j);
  virtual double calculate_rho();
  virtual double calculate_rho_cg();
  virtual void do_shrinking();
  virtual bool do_shrinking_cg();

private:
  bool bumped();
  bool be_shrunk(int i, double Gmax1, double Gmax2);	
  bool be_shrunk_cg(int i, double Gmax1, double Gmax2);	
  void generate_direction3(double *u, bool new_working_set);
  void generate_direction4(double *u, bool new_working_set);
  int generate_direction_general(double *u, int n, bool new_working_set);
  int select_working_set_first_order();
  void compute_step_first_order(double* u, double& lambda_star);
  int select_working_set_incrementally(double* u, double& lambda_star, double& gain);
  bool feasible_direction(double *u);
  int iter;
  int curr_depth;
  double **A;
  double *b;
  int sizeof_double;
  int sizeof_int;
  int sizeof_char;
  bool *active;
};

#endif
