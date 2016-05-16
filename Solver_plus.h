
#ifndef _SOLVER_PLUS_H_
#define _SOLVER_PLUS_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include "common.h"
#include "Kernel.h"

// An SMO and conjugate SMO algorithm for SVM+
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
class Solver_plus{
 public:
  Solver_plus(int optimizer_) {
    if(optimizer_==-1)
      conjugate = false;
    else {
      conjugate = true;
      max_depth =  optimizer_; 

      A = new double*[max_depth+5];
      for(int i=0; i<=max_depth+4; i++)
	A[i] = new double[max_depth+5];
      b = new double[max_depth+5];
    }

    sizeof_double = sizeof(double);
    sizeof_char = sizeof(char);
    sizeof_int = sizeof(int);
  };
  virtual ~Solver_plus() {
    if(conjugate) {
      for(int i=0; i<=max_depth+4; i++)
	delete[] A[i];
      delete[] A;
      delete[] b;
    }
  };

  void Solve_plus(int l, const QMatrix& Q, const QMatrix& Q_star, const QMatrix& Q_star_beta, const schar *y_,
		  double *alpha_, double *beta_, double Cp, double Cn, double tau, double eps,
		  SolutionInfo* si, int shrinking);
  void Solve_plus_cg(int l, const QMatrix& Q, const QMatrix& Q_star, const QMatrix& Q_star_beta, const schar *y_,
		     double *alpha_, double *beta_, double Cp, double Cn, double tau, double eps,
		     SolutionInfo* si, int shrinking);
 protected:
  int active_size;
  int active_size_beta;
  schar *y;
  double *G;		// gradient of objective function
  double **G_cg;
  double *g;
  double **g_cg;
  double *g_beta;
  double **g_beta_cg;
  double *g_init;
  double *g_beta_init;
  int *work;
  enum { LOWER_BOUND, UPPER_BOUND, FREE };
  char *alpha_status;	// LOWER_BOUND, UPPER_BOUND, FREE
  char *beta_status;	// LOWER_BOUND, FREE
  char **alpha_status_cg;
  char **beta_status_cg;
  double *alpha;
  double *beta;
  double **alpha_cg;
  double **beta_cg;
  const QMatrix *Q;
  const QMatrix *Q_star;
  const QMatrix *Q_star_beta;
  const Qfloat *QD;
  const Qfloat *QD_star;
  const Qfloat *QD_star_beta;
  double eps;
  double Cp,Cn;
  double tau;
  double *p;
  int *active_set;
  int *active_set_beta;
  int *true_act_set;
  int *true_act_set_beta;
  double *G_bar;		// gradient, if we treat free variables as 0
  int l;
  bool unshrink;	// XXX

  double get_C(int i)
  {
    return (y[i] > 0)? Cp : Cn;
  }
  void update_alpha_status(int i)
  {
    if(alpha[i] <= 1e-8)
      alpha_status[i] = LOWER_BOUND;
    else alpha_status[i] = FREE;
  }
  void update_beta_status(int i)
  {
    if(beta[i] <= 1e-8)
      beta_status[i] = LOWER_BOUND;
    else beta_status[i] = FREE;
  }
  void update_alpha_status_cg(int i, int depth)
  {
    if(alpha_cg[depth][i] <= 1e-8)
      alpha_status_cg[depth][i] = LOWER_BOUND;
    else alpha_status_cg[depth][i] = FREE;
  }
  void update_beta_status_cg(int i, int depth)
  {
    if(beta_cg[depth][i] <= 1e-8)
      beta_status_cg[depth][i] = LOWER_BOUND;
    else beta_status_cg[depth][i] = FREE;
  }
  bool is_lower_bound(int i) { return alpha_status[i] == LOWER_BOUND; }
  bool is_lower_bound_beta(int i) { return beta_status[i] == LOWER_BOUND; }
  bool is_free(int i) { return alpha_status[i] == FREE; }
  bool is_free_beta(int i) { return beta_status[i] == FREE; }
  bool is_lower_bound_cg(int i, int depth) { return alpha_status_cg[depth][i] == LOWER_BOUND; }
  bool is_lower_bound_beta_cg(int i, int depth) { return beta_status_cg[depth][i] == LOWER_BOUND; }
  bool is_free_cg(int i, int depth) { return alpha_status_cg[depth][i] == FREE; }
  bool is_free_beta_cg(int i, int depth) { return beta_status_cg[depth][i] == FREE; }
  void swap_index_alpha(int i, int j);
  void swap_index_beta(int i, int j);
  void swap_index_alpha_cg(int i, int j);
  void swap_index_beta_cg(int i, int j);
  void reconstruct_gradient_plus();
  virtual int select_working_set_plus(int& set_type, int& i, int& j, int& k, int iter);
  virtual void calculate_rho_plus(double& rho, double& rho_star);
  virtual void calculate_rho_plus_cg(double& rho, double& rho_star);
  virtual void do_shrinking_plus();
  virtual bool do_shrinking_plus_cg();
  virtual int select_working_set(double* u, double& lambda_star);

 private:
 bool be_shrunk_alpha(int i, double max_B1, double max_A1, double max_A2,  double min_B1B2, double min_A1A3, double min_A2A4);
   bool be_shrunk_beta(int i, double max_B1, double max_A1, double max_A2,  double min_B1B2, double min_A1A3, double min_A2A4);
   bool be_shrunk_alpha_cg(int i, double max_B1, double max_A1, double max_A2,  double min_B1B2, double min_A1A3, double min_A2A4);
   bool be_shrunk_beta_cg(int i, double max_B1, double max_A1, double max_A2,  double min_B1B2, double min_A1A3, double min_A2A4);
   void generate_direction3(double *u, bool new_working_set);
   void generate_direction4(double *u, bool new_working_set);	
   int generate_direction_general(double *u, int n, bool new_working_set);
   void generate_direction4y(double *u, bool new_working_set);
   void generate_direction5y(double *u, bool new_working_set);	
   int generate_direction_y_general(double *u, int n, bool new_working_set); 
   int generate_direction(double *u, int n, bool new_working_set);
   int generate_direction_y(double *u, int n, bool new_working_set);
  void reconstruct_gradient_plus_cg();
  int wss_first_order(double *best_u, double& lambda_star);
  int select_working_set_plus_hmg(double *best_u, double& lambda_star, double& gain);
  int select_working_set_plus_hmg2(double *best_u, double& lambda_star, double& gain);
  int select_working_set_plus_incrementally(double *best_u, double& lambda_star, double& gain);
  bool compute_gain(int new_working_set_size, double &gain, double *curr_u, double &lambda, bool new_working_set);
  bool feasible_direction(double *u);
  double **A;
  double *b;
  int sizeof_double;
  int sizeof_int;
  int sizeof_char;
  bool conjugate;
  int max_depth;
  int curr_depth;
  int prev_depth;
  int working_set_size;
  bool *active;
  char working_set_type;
  enum {BETAS=0, ALPHAS=1, ALPHAS_BETAS=2, ALPHAS_DIFF_SIGN=3};
};


#endif

