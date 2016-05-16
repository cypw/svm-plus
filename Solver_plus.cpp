#include "common.h"
#include "Solver_plus.h"
#include "solve_linear_system.h"

static int* tmp_work;
static double* tmp_u;

int Solver_plus::generate_direction(double *u, int n, bool new_working_set)
{
  switch(n) {
  case 3:
    generate_direction3(u, new_working_set);
    return 0;
  case 4:
    generate_direction4(u, new_working_set);
    return 0;
  default:
    return generate_direction_general(u, n, new_working_set);
  }
}

int Solver_plus::generate_direction_y(double *u, int n, bool new_working_set)
{
  switch(n) {
  case 4:
    generate_direction4y(u, new_working_set);
    return 0;
  case 5:
    generate_direction5y(u, new_working_set);
    return 0;
  default:
    return generate_direction_y_general(u, n, new_working_set);
  }
}

void Solver_plus::generate_direction3(double *u, bool new_working_set)
{  
  int working2 = work[2];
  int working_new;

  double *G0 = G_cg[0];
  double *G1 = G_cg[1];
  double *g0 = g_cg[0];
  double *g1 = g_cg[1];
  double *g0beta = g_beta_cg[0];
  double *g1beta = g_beta_cg[1];
  int y_i;

  static double a0, a1;
  double a2;

  if(new_working_set) {
    int working0 = work[0];
    int working1 = work[1];

    if(working0<l) {
      y_i = y[working0];
      a0 = y_i*G0[working0]+g0[working0]/tau-y_i*G1[working0]-g1[working0]/tau;
    }
    else {
      working_new = working0 - l;
      a0 = g0beta[working_new]/tau-g1beta[working_new]/tau;
    }

    if(working1<l) {
      y_i = y[working1];
      a1 = y_i*G0[working1]+g0[working1]/tau-y_i*G1[working1]-g1[working1]/tau;
    }
    else {
      working_new = working1-l;
      a1 = g0beta[working_new]/tau-g1beta[working_new]/tau;
    }
  }

  if(working2<l) {
    y_i = y[working2];
    a2 = y_i*G0[working2]+g0[working2]/tau-y_i*G1[working2]-g1[working2]/tau;
  }
  else {
    working_new = working2-l;
    a2 = g0beta[working_new]/tau-g1beta[working_new]/tau;
  }

  u[0] = a2 - a1;
  u[1] = a0 - a2;
  u[2] = a1 - a0;
}

void Solver_plus::generate_direction4(double *u, bool new_working_set)
{
  int working3 = work[3];
  int working_new;
  int y_i;

  double *G0 = G_cg[0];
  double *G1 = G_cg[1];
  double *G2 = G_cg[2];

  double *g0 = g_cg[0];
  double *g1 = g_cg[1];
  double *g2 = g_cg[2];

  double *g0beta = g_beta_cg[0];
  double *g1beta = g_beta_cg[1];
  double *g2beta = g_beta_cg[2];
 
  double G0old;
  double a14, a24;
  static double a11, a12, a13, a21, a22, a23;
  static double d1, d2, d4;
  double d3, d5, d6;

  if(new_working_set) {
    int working0 = work[0];
    int working1 = work[1];
    int working2 = work[2];
    if(working0<l) {
      y_i = y[working0];
      G0old = y_i*G1[working0] + g1[working0]/tau;
      a11 = y_i*G0[working0] + g0[working0]/tau - G0old;
      a21 = G0old - y_i*G2[working0] - g2[working0]/tau;
    }
    else {
      working_new = working0 - l;
      G0old = g1beta[working_new]/tau;
      a11 = g0beta[working_new]/tau - G0old;
      a21 = G0old - g2beta[working_new]/tau;
    }

    if(working1<l) {
      y_i = y[working1];
      G0old = y_i*G1[working1] + g1[working1]/tau;
      a12 = y_i*G0[working1] + g0[working1]/tau - G0old;
      a22 = G0old - y_i*G2[working1] - g2[working1]/tau;
    }
    else {
      working_new = working1 - l;
      G0old = g1beta[working_new]/tau;
      a12 = g0beta[working_new]/tau - G0old;
      a22 = G0old - g2beta[working_new]/tau;
    }

    if(working2<l) {
      y_i = y[working2];
      G0old = y_i*G1[working2] + g1[working2]/tau;
      a13 = y_i*G0[working2] + g0[working2]/tau - G0old;
      a23 = G0old - y_i*G2[working2] - g2[working2]/tau;
    }
    else {
      working_new = working2 - l;
      G0old = g1beta[working_new]/tau;
      a13 = g0beta[working_new]/tau - G0old;
      a23 = G0old - g2beta[working_new]/tau;
    }

    d1 = a11*a22-a12*a21;
    d2 = a11*a23-a13*a21;
    d4 = a12*a23-a13*a22;
  }

  if(working3<l) {
    G0old = y[working3]*G1[working3] + g1[working3]/tau;
    a14 = y[working3]*G0[working3] + g0[working3]/tau - G0old;
    a24 = G0old - y[working3]*G2[working3] - g2[working3]/tau;
  }
  else {
    working_new = working3 - l;
    G0old = g1beta[working_new]/tau;
    a14 = g0beta[working_new]/tau - G0old;
    a24 = G0old - g2beta[working_new]/tau;
  }

  d3 = a11*a24-a21*a14;
  d5 = a12*a24-a14*a22;
  d6 = a13*a24-a14*a23;

  u[0] = d6-d5+d4;
  u[1] = -d6+d3-d2;
  u[2] = d5-d3+d1;
  u[3] = -d4+d2-d1;
}

int Solver_plus::generate_direction_general(double *u, int n, bool new_working_set)
{
  int i,j,result,worknew;
  int n_minus_one = n-1;
  int n_minus_two = n-2;
  int i_minus_one;
  double *A0 = A[0];

  // compute A matrix
  if(new_working_set) {
    for(j=0; j<n_minus_one; j++)
      A0[j] = 1;
  }

  double *Ai, *G, *G_old, *g, *g_old, *g_beta, *g_beta_old;
  int workj, start;

  if(new_working_set)
    start = 0;
  else
    start = n_minus_two;

  for(i=1; i<n_minus_one; i++) {
    i_minus_one = i - 1;
    Ai = A[i];
    G = G_cg[i_minus_one];
    G_old = G_cg[i];
    g = g_cg[i_minus_one];
    g_old = g_cg[i];
    g_beta = g_beta_cg[i_minus_one];
    g_beta_old = g_beta_cg[i];

    for(j=start; j<n_minus_one; j++) {
      workj = work[j+1];
      if(workj<l)
	Ai[j] = y[workj]*G[workj]+g[workj]/tau-y[workj]*G_old[workj]-g_old[workj]/tau;
      else {
	worknew = workj-l;
	Ai[j] = g_beta[worknew]/tau-g_beta_old[worknew]/tau;
      }
    }
  }

  // compute b vector
  if(new_working_set) {
    int work0 = work[0];
    b[0] = -1;
    if(work0<l)
      for(i=1; i<=n_minus_two; i++) {
	i_minus_one = i-1;
	b[i] = -y[work0]*G_cg[i_minus_one][work0]-g_cg[i_minus_one][work0]/tau+y[work0]*G_cg[i][work0]+g_cg[i][work0]/tau;
      }
    else {
      worknew = work0 - l; 
      for(i=1; i<=n_minus_two; i++)
	b[i] = -g_beta_cg[i-1][worknew]/tau+g_beta_cg[i][worknew]/tau;
    }
  }

  result = solve_linear_system(A,b,n_minus_one);
  if(result)
    return 1;
  u[0] = 1;
  memcpy(&u[1], b, sizeof_double*(n_minus_one));

  return 0;  
}

void Solver_plus::generate_direction4y(double *u, bool new_working_set)
{
  int working3 = work[3];
  int y_i;
  int working_new;
 
  static int a11, a12, a13;
  int a14 = y[working3];

  double *G0 = G_cg[0];
  double *G1 = G_cg[1];
  double *g0 = g_cg[0];
  double *g1 = g_cg[1];
  double *g0beta = g_beta_cg[0];
  double *g1beta = g_beta_cg[1];

  static double a21, a22, a23, d1, d2, d4;
  double a24, d3, d5, d6;

  if(new_working_set) {
    int working0 = work[0];
    int working1 = work[1];
    int working2 = work[2];
    a11 = y[working0];
    a12 = y[working1];
    a13 = y[working2];

    if(working0<l) {
      y_i = y[working0];
      a21 = y_i*G0[working0]+g0[working0]/tau-y_i*G1[working0]-g1[working0]/tau;
    }
    else {
      working_new = working0 - l;
      a21 = g0beta[working_new]/tau-g1beta[working_new]/tau;
    }

    if(working1<l) {
      y_i = y[working1];
      a22 = y_i*G0[working1]+g0[working1]/tau-y_i*G1[working1]-g1[working1]/tau;
    }
    else {
      working_new = working1 - l;
      a22 = g0beta[working_new]/tau-g1beta[working_new]/tau;
    }
  
    if(working2<l) {
      y_i = y[working2];
      a23 = y_i*G0[working2]+g0[working2]/tau-y_i*G1[working2]-g1[working2]/tau;
    }
    else {
      working_new = working2 - l;
      a23 = g0beta[working_new]/tau-g1beta[working_new]/tau;
    }

    d1 = a11*a22-a12*a21;
    d2 = a11*a23-a13*a21;
    d4 = a12*a23-a13*a22;
  }

  if(working3<l) {
    y_i = y[working3];
    a24 = y_i*G0[working3]+g0[working3]/tau-y_i*G1[working3]-g1[working3]/tau;
  }
  else {
    working_new = working3 - l;
    a24 = g0beta[working_new]/tau-g1beta[working_new]/tau;
  }

  
  d3 = a11*a24-a21*a14;
  d5 = a12*a24-a14*a22;
  d6 = a13*a24-a14*a23;

  u[0] = d6-d5+d4;
  u[1] = -d6+d3-d2;
  u[2] = d5-d3+d1;
  u[3] = -d4+d2-d1;
}

void Solver_plus::generate_direction5y(double *u, bool new_working_set)
{
  int working4 = work[4];
  int working_new;
  int y_i;

  int y5 = y[working4];
  static int y1, y2, y3, y4;

  double *G0=G_cg[0];
  double *G1=G_cg[1];
  double *G2=G_cg[2];

  double *g0=g_cg[0];
  double *g1=g_cg[1];
  double *g2=g_cg[2];

  double *g0beta=g_beta_cg[0];
  double *g1beta=g_beta_cg[1];
  double *g2beta=g_beta_cg[2];

  double G0old, a15, a25;  
  static double a11, a12, a13, a14, a21, a22, a23, a24; 
  static double d12, d13, d14, d23, d24, d34, d123, d124, d134, d234;

  if(new_working_set) {
    int working0 = work[0];
    int working1 = work[1];
    int working2 = work[2];
    int working3 = work[3];
    y1 = y[working0];
    y2 = y[working1];
    y3 = y[working2];
    y4 = y[working3];

    if(working0<l) {
      y_i = y[working0];
      G0old = y_i*G1[working0] + g1[working0]/tau;
      a11 = y_i*G0[working0] + g0[working0]/tau - G0old;
      a21 = G0old - y_i*G2[working0] - g2[working0]/tau;
    }
    else {
      working_new = working0 - l;
      G0old = g1beta[working_new]/tau;
      a11 = g0beta[working_new]/tau - G0old;
      a21 = G0old - g2beta[working_new]/tau;
    }

    if(working1<l) {
      y_i = y[working1];
      G0old = y_i*G1[working1] + g1[working1]/tau;
      a12 = y_i*G0[working1] + g0[working1]/tau - G0old;
      a22 = G0old - y_i*G2[working1] - g2[working1]/tau;
    }
    else {
      working_new = working1 - l;
      G0old = g1beta[working_new]/tau;
      a12 = g0beta[working_new]/tau - G0old;
      a22 = G0old - g2beta[working_new]/tau;
    }

    if(working2<l) {
      y_i = y[working2];
      G0old = y_i*G1[working2] + g1[working2]/tau;
      a13 = y_i*G0[working2] + g0[working2]/tau - G0old;
      a23 = G0old - y_i*G2[working2] - g2[working2]/tau;
    }
    else {
      working_new = working2 - l;
      G0old = g1beta[working_new]/tau;
      a13 = g0beta[working_new]/tau - G0old;
      a23 = G0old - g2beta[working_new]/tau;
    }

    if(working3<l) {
      y_i = y[working3];
      G0old = y_i*G1[working3] + g1[working3]/tau;
      a14 = y_i*G0[working3] + g0[working3]/tau - G0old;
      a24 = G0old - y_i*G2[working3] - g2[working3]/tau;
    }
    else {
      working_new = working3 - l;
      G0old = g1beta[working_new]/tau;
      a14 = g0beta[working_new]/tau - G0old;
      a24 = G0old - g2beta[working_new]/tau;
    }

    d12 = a11*a22-a21*a12;
    d13 = a11*a23-a13*a21;
    d14 = a11*a24-a14*a21;
    d23 = a12*a23-a13*a22;
    d24 = a12*a24-a14*a22;
    d34 = a13*a24-a14*a23;
    d123 = y1*d23-y2*d13+y3*d12;
    d124 = y1*d24-y2*d14+y4*d12;
    d134 = y1*d34-y3*d14+y4*d13;
    d234 = y2*d34-y3*d24+y4*d23;
  }

  if(working4<l) {
    y_i = y[working4];
    G0old = y_i*G1[working4] + g1[working4]/tau;
    a15 = y_i*G0[working4] + g0[working4]/tau - G0old;
    a25 = G0old - y_i*G2[working4] - g2[working4]/tau;
  }
  else {
    working_new = working4 - l;
    G0old = g1beta[working_new]/tau;
    a15 = g0beta[working_new]/tau - G0old;
    a25 = G0old - g2beta[working_new]/tau;
  }

  
  double d15 = a11*a25-a15*a21;
  double d25 = a12*a25-a15*a22;
  double d35 = a13*a25-a23*a15;
  double d45 = a14*a25-a15*a24;

  double d125 = y1*d25-y2*d15+y5*d12;
  double d135 = y1*d35-y3*d15+y5*d13;
  double d145 = y1*d45-y4*d15+y5*d14;
  double d235 = y2*d35-y3*d25+y5*d23;
  double d345 = y3*d45-y4*d35+y5*d34;
  double d245 = y2*d45-y4*d25+y5*d24;

  u[0] = d345 - d245 + d235 - d234;
  u[1] = -d345 + d145 - d135 + d134;
  u[2] = d245 - d145 + d125 - d124;
  u[3] = -d235 + d135 - d125 + d123;
  u[4] = d234 - d134 + d124 - d123;
}

int Solver_plus::generate_direction_y_general(double *u, int n, bool new_working_set)
{
  int i,j,result,worknew;
  double *A0 = A[0], *A1 = A[1];
  int n_minus_one=n-1, n_minus_two=n-2;
  int i_minus_one, i_minus_two, start;

  if(new_working_set)
    start = 0;
  else
    start = n_minus_two;

  // compute A matrix
  if(new_working_set)
    for(i=0; i<n_minus_one; i++)
      A0[i] = 1;
  for(i=start; i<n_minus_one; i++)
    A1[i] = y[work[i+1]];

  double *Ai, *G, *G_old, *g, *g_old, *g_beta, *g_beta_old;
  int workj, y_j, y_i;

  for(i=2; i<n_minus_one; i++) {
    i_minus_one = i - 1;
    i_minus_two = i - 2;
    Ai = A[i];
    G = G_cg[i_minus_two];
    G_old = G_cg[i_minus_one];
    g = g_cg[i_minus_two];
    g_old = g_cg[i_minus_one];
    g_beta = g_beta_cg[i_minus_two];
    g_beta_old = g_beta_cg[i_minus_one];

    for(j=start; j<n_minus_one; j++) {
      workj = work[j+1];
      if(workj<l) {
        y_j = y[workj];
	Ai[j] = y_j*G[workj]+g[workj]/tau-y_j*G_old[workj]-g_old[workj]/tau;
      }
      else {
	worknew = workj-l;
	Ai[j] = g_beta[worknew]/tau-g_beta_old[worknew]/tau;
      }
    }
  }

  // compute b vector
  if(new_working_set) {
    int work0 = work[0];
    b[0] = -1;
    b[1] = -y[work0];

    if(work0<l)
      for(i=2; i<=n_minus_two; i++) {
	i_minus_one = i - 1;
	i_minus_two = i - 2;
	y_i = y[work0];
	b[i] = -y_i*G_cg[i_minus_two][work0]-g_cg[i_minus_two][work0]/tau+y_i*G_cg[i_minus_one][work0]+g_cg[i_minus_one][work0]/tau;
      }
    else {
      worknew = work0 - l; 
      for(i=2; i<=n_minus_two; i++)
	b[i] = -g_beta_cg[i-2][worknew]/tau+g_beta_cg[i-1][worknew]/tau;
    }
  }

  result = solve_linear_system(A,b,n_minus_one);
  if(result)
    return 1;
  u[0] = 1;
  memcpy(&u[1], b, sizeof_double*(n_minus_one));
 
  return 0; 
}

void Solver_plus::swap_index_alpha(int i, int j)
{
  Q->swap_index(i,j);
  Q_star->swap_index(i,j);
  swap(y[i],y[j]);
  swap(alpha_status[i],alpha_status[j]);
  swap(alpha[i],alpha[j]);
  swap(true_act_set[active_set[i]],true_act_set[active_set[j]]);
  swap(active_set[i],active_set[j]);
  swap(G[i],G[j]);
  swap(g[i],g[j]);
  swap(g_init[i],g_init[j]);
}

void Solver_plus::swap_index_beta(int i, int j)
{
  Q_star_beta->swap_index(i,j);
  swap(beta_status[i],beta_status[j]);
  swap(beta[i],beta[j]);
  swap(true_act_set_beta[active_set_beta[i]],true_act_set_beta[active_set_beta[j]]);
  swap(active_set_beta[i],active_set_beta[j]);
  swap(g_beta[i],g_beta[j]);
  swap(g_beta_init[i],g_beta_init[j]);
}

void Solver_plus::swap_index_alpha_cg(int i, int j)
{
  int k;
  Q->swap_index(i,j);
  Q_star->swap_index(i,j);
  swap(y[i],y[j]);
  swap(true_act_set[active_set[i]],true_act_set[active_set[j]]);
  swap(active_set[i],active_set[j]);
  swap(g_init[i],g_init[j]);

  for(k=0; k < curr_depth+1; k++) {
    swap(alpha_status_cg[k][i],alpha_status_cg[k][j]);
    swap(alpha_cg[k][i],alpha_cg[k][j]);   
    swap(G_cg[k][i],G_cg[k][j]);
    swap(g_cg[k][i],g_cg[k][j]);  
  }

  for(k=0; k<working_set_size; k++) {
    if(i == work[k]) 
      work[k] = -1;
    else 
      if(j == work[k])
	work[k] = i;  
  }
  active[i] = active[j];
  active[j] = false;
}

void Solver_plus::swap_index_beta_cg(int i, int j)
{
  Q_star_beta->swap_index(i,j);
  swap(g_beta_init[i],g_beta_init[j]);
  swap(true_act_set_beta[active_set_beta[i]],true_act_set_beta[active_set_beta[j]]);
  swap(active_set_beta[i],active_set_beta[j]);

  int k;
  for(k=0; k < curr_depth+1; k++) {
    swap(beta_status_cg[k][i],beta_status_cg[k][j]);
    swap(beta_cg[k][i],beta_cg[k][j]);  
    swap(g_beta_cg[k][i],g_beta_cg[k][j]);
  }

  int il = i+l;
  int jl = j+l;

  for(k=0; k<working_set_size; k++) {
    if(il == work[k]) 
      work[k] = -1;
    else 
      if(jl == work[k])
	work[k] = il;  
  }
  active[il] = active[jl];
  active[jl] = false;
}

bool Solver_plus::do_shrinking_plus_cg()
{
  int i, j, y_i;
  double g_i, alpha_i, deriv_alpha_i;
  double max_B1=-1e20, min_B1B2=1e20, max_A2=-1e20, min_A2A4=1e20, max_A1=-1e20, min_A1A3=1e20;
  bool done_shrinking;
  double *alpha, *G, *g, *g_beta, *beta;

  // compute all maxima and minima related to alphas
  alpha = alpha_cg[0];
  beta = beta_cg[0];
  g = g_cg[0];
  G = G_cg[0];
  g_beta = g_beta_cg[0];
  for(i=0; i<active_size; i++) {
    alpha_i = alpha[i];
    g_i = g[i];
    y_i = y[i];
    deriv_alpha_i = y_i*G[i]+g_i/tau;

    // max A2
    if(alpha_i>1e-8 && y_i==-1 && deriv_alpha_i>max_A2) 
      max_A2 = deriv_alpha_i;

    // min A2A4
    if(y_i==-1 && deriv_alpha_i<min_A2A4) 
      min_A2A4 = deriv_alpha_i;
    
    // max A1
    if(alpha_i>1e-8 && y_i==1 && deriv_alpha_i>max_A1) 
      max_A1 = deriv_alpha_i;
    
    // min A1A3max_A2, min_A2A4, max_A1, min_A1A3
    if(y_i==1 && deriv_alpha_i<min_A1A3) 
      min_A1A3 = deriv_alpha_i;
  }

  // compute all maxima and minima related to betas
  for(i=0; i<active_size_beta; i++) {
    g_i = g_beta[i];

    // max B1
    if(beta[i]>1e-8 && g_i>max_B1) 
      max_B1 = g_i;

    // min B1B2
    if(g_i<min_B1B2) 
      min_B1B2 = g_i;
  }

  max_B1 /= tau;
  min_B1B2 /= tau;

  if(unshrink == false && max_B1-min_B1B2 < eps*10 &&
     max_A2-min_A2A4 < eps*10 && max_A1 - min_A1A3 < eps*10 &&
     2*max_B1+2-min_A1A3-min_A2A4 < eps*10 && max_A1+max_A2-2*min_B1B2-2 < eps*10) {
    unshrink = true;
    reconstruct_gradient_plus_cg();
    active_size = l;
     active_size_beta = l;
  }

  if(active_size_beta > 2) {
    for(i=0;i<active_size_beta;i++) {
      if(active_size_beta <= 2)
	break;
      if (be_shrunk_beta_cg(i, max_B1, max_A1, max_A2, min_B1B2, min_A1A3, min_A2A4)) {
	active_size_beta--;
	//        fprintf(stdout,"shrinking beta %d\n",i);
        //fflush(stdout);
        if(active_size_beta == i) {
          if(active[i+l]) {
            active[i+l] = false;
            for(j=0; j<working_set_size; j++)
	      if(work[j]==i+l) {
		work[j] = -1;
		break;
	      }
	  }
	}
        else {
	  while (active_size_beta > i) {
	    if (!be_shrunk_beta_cg(active_size_beta, max_B1, max_A1, max_A2, min_B1B2, min_A1A3, min_A2A4)) {
	      swap_index_beta_cg(i,active_size_beta);
	      done_shrinking = true;
	      break;
	    }
	    else {
	      if(active[active_size_beta+l]) {
		active[active_size_beta+l] = false;
		for(j=0; j<working_set_size; j++)
		  if(work[j]==active_size_beta+l) {
		    work[j] = -1;
		    break;
		  }
	      }
	      active_size_beta--;
	      if(active_size_beta <= 2)
		break;
	    }
	  }
	}
      }
    }
  }

  if(active_size>2) {
    for(i=0;i<active_size;i++) {
      if(active_size <= 2)
	break;
      if (be_shrunk_alpha_cg(i, max_B1, max_A1, max_A2, min_B1B2, min_A1A3, min_A2A4)) {
	active_size--;
	if(active_size == i) {
          if(active[i]) {
            active[i] = false;
            for(j=0; j<working_set_size; j++)
	      if(work[j]==i) {
		work[j] = -1;
		break;
	      }
	  }
	}
        else {
	  while (active_size > i) {
	    if (!be_shrunk_alpha_cg(active_size, max_B1, max_A1, max_A2, min_B1B2, min_A1A3, min_A2A4)) {
	      swap_index_alpha_cg(i,active_size);
	      done_shrinking = true;
	      break;
	    }
	    else {
	      if(active[active_size]) {
		active[active_size] = false;
		for(j=0; j<working_set_size; j++)
		  if(work[j]==active_size) {
		    work[j] = -1;
		    break;
		  }
	      }
	      active_size--;
	    }
	  }
	}
      }
    }
  }

  for(i=0; i<working_set_size; i++)
    if(work[i]==-1) {
      working_set_size--;
      while(working_set_size > i) {
	if(work[working_set_size] != -1) {
	  work[i] = work[working_set_size];
	  work[working_set_size] = -1;
	  break;
	}
	working_set_size--;
      }
    }

  switch(working_set_type) {
  case ALPHAS:
  case BETAS: 
    if(working_set_size <= 1) { 
      working_set_size = 0;
      curr_depth = 0;
    }
    break;    
  case ALPHAS_BETAS:
    if(working_set_size <= 2) { 
      working_set_size = 0;
      curr_depth = 0;
    }
    break;
  case ALPHAS_DIFF_SIGN:
    if(working_set_size <= 3) {
      working_set_size = 0;
      curr_depth = 0;
    }
  }
    
  return done_shrinking;
}

void Solver_plus::reconstruct_gradient_plus_cg()
{
  int i, j, k, true_i, act_set_i;
  double *alpha_k, *beta_k, *G_k, *g_k, *g_beta_k;

  //  fprintf(stdout,"reconstructing gradient\n");
  //fflush(stdout);

  if(active_size < l) {
    for(k=0; k<curr_depth+1; k++) {
      alpha_k = alpha_cg[k];
      beta_k = beta_cg[k];
      G_k = G_cg[k];
      g_k = g_cg[k];
      for(i=active_size;i<l;i++) {
	const Qfloat *Q_i = Q->get_Q(i,l);
	const Qfloat *Q_i_star = Q_star->get_Q(i,l);
       
	true_i = active_set[i];
	act_set_i = true_act_set_beta[true_i];

	const Qfloat *Q_i_star_beta = Q_star_beta->get_Q(act_set_i,l);
	G_k[i] = 0;
	g_k[i] = g_init[i];
	for(j=0;j<l;j++)
	  if(alpha_k[j]>1e-8) {
	    G_k[i] += alpha_k[j] * y[j] * Q_i[j];
	    g_k[i] += alpha_k[j] * Q_i_star[j];
	  }
	for(j=0;j<l;j++)
	  if(beta_k[j]>1e-8) 
	    g_k[i] += beta_k[j] * Q_i_star_beta[j];
      }
    }
  }

  if(active_size_beta < l) {
    for(k=0; k<curr_depth+1; k++) {
      alpha_k = alpha_cg[k];
      beta_k = beta_cg[k];
      g_beta_k = g_beta_cg[k];
      for(i=active_size_beta;i<l;i++) {
	const Qfloat *Q_i_star_beta = Q_star_beta->get_Q(i,l);
      
	true_i = active_set_beta[i];
	act_set_i = true_act_set[true_i];
	const Qfloat *Q_i_star = Q_star->get_Q(act_set_i,l);
       
	g_beta_k[i] = g_beta_init[i];
    
	for(j=0;j<l;j++)
	  if(beta_k[j]>1e-8) 
	    g_beta_k[i] += beta_k[j] * Q_i_star_beta[j];

	for(j=0;j<l;j++)
	  if(alpha_k[j]>1e-8) 
	    g_beta_k[i] += alpha_k[j] * Q_i_star[j];
      }
    }
  }
}

int Solver_plus::wss_first_order(double *best_u, double& lambda_star)
{
  int i, y_i;
  double g_i, alpha_i, deriv_alpha_i;
  double max_B1=-1e20, min_B1B2=1e20, max_A2=-1e20, min_A2A4=1e20, max_A1=-1e20, min_A1A3=1e20;
  double *alpha, *G, *g, *g_beta, *beta;
  int best_B1=-1, best_B1B2=-1, best_A2=-1, best_A2A4=-1, best_A1=-1, best_A1A3=-1;
  int type_selected[3], selected_indices[3][3], alpha_beta_type=-1;

  for(i=0; i<working_set_size; i++) {
    active[work[i]] = false;
    work[i] = -1;
  }
  working_set_size = 0;
  curr_depth = 0;

  // compute all maxima and minima related to alphas
  alpha = alpha_cg[0];
  beta = beta_cg[0];
  g = g_cg[0];
  G = G_cg[0];
  g_beta = g_beta_cg[0];
  for(i=0; i<active_size; i++) {
    alpha_i = alpha[i];
    g_i = g[i];
    y_i = y[i];
    deriv_alpha_i = y_i*G[i]+g_i/tau;

    // max A2
    if(alpha_i>1e-8 && y_i==-1 && deriv_alpha_i>max_A2) {
      max_A2 = deriv_alpha_i;
      best_A2 = i;
    }

    // min A2A4
    if(y_i==-1 && deriv_alpha_i<min_A2A4) {
      min_A2A4 = deriv_alpha_i;
      best_A2A4 = i;
    }

    // max A1
    if(alpha_i>1e-8 && y_i==1 && deriv_alpha_i>max_A1) {
      max_A1 = deriv_alpha_i;
      best_A1 = i;
    }
    
    // min A1A3max_A2, min_A2A4, max_A1, min_A1A3
    if(y_i==1 && deriv_alpha_i<min_A1A3) {
      min_A1A3 = deriv_alpha_i;
      best_A1A3 = i;
    }
  }

  // compute all maxima and minima related to betas
  for(i=0; i<active_size_beta; i++) {
    g_i = g_beta[i];

    // max B1
    if(beta[i]>1e-8 && g_i>max_B1) {
      max_B1 = g_i;
      best_B1 = i;
    }

    // min B1B2
    if(g_i<min_B1B2) {
      min_B1B2 = g_i;
      best_B1B2 = i;
    }
  }

  double gap[3];

  for(i=0; i <3; i++)
    gap[i] = -1;

  max_B1 /= tau;
  min_B1B2 /= tau;

  // select maximal violating pairs
  if(max_B1-min_B1B2 < eps)
    type_selected[0] = 0;
  else {
    type_selected[0] = 1;
    selected_indices[0][0] = best_B1;
    selected_indices[0][1] = best_B1B2;
    gap[0] = max_B1-min_B1B2;
  }

  if(max_A2-min_A2A4 < eps && max_A1 - min_A1A3 < eps)
    type_selected[1] = 0;
  else {
    type_selected[1] = 1;
    if(max_A2-min_A2A4 >  max_A1 - min_A1A3) {
      selected_indices[1][0] = best_A2;
      selected_indices[1][1] = best_A2A4;
      gap[1] = max_A2-min_A2A4;
    }
    else {
      selected_indices[1][0] = best_A1;
      selected_indices[1][1] = best_A1A3;
      gap[1] = max_A1-min_A1A3;
    }
  }

  if(2*max_B1+2-min_A1A3-min_A2A4 < eps && max_A1+max_A2-2*min_B1B2-2 < eps)
    type_selected[2] = 0;
  else {
    type_selected[2] = 1;
    if(2*max_B1+2-min_A1A3-min_A2A4 > max_A1+max_A2-2*min_B1B2-2) {
      selected_indices[2][0] = best_A1A3;
      selected_indices[2][1] = best_A2A4;
      selected_indices[2][2] = best_B1;
      gap[2] = 2*max_B1+2-min_A1A3-min_A2A4;
      alpha_beta_type = 0;
    }
    else {
      selected_indices[2][0] = best_A1;
      selected_indices[2][1] = best_A2;
      selected_indices[2][2] = best_B1B2;
      gap[2] = max_A1+max_A2-2*min_B1B2-2;
      alpha_beta_type = 1;
    }
  }

  if(type_selected[0]+type_selected[1]+type_selected[2] == 0) 
    return 1;

  if(gap[2]>=gap[1] && gap[2]>=gap[0]) {
    working_set_type = ALPHAS_BETAS;
    working_set_size = 3;
    work[0] = selected_indices[2][0];
    work[1] = selected_indices[2][1];
    work[2] = selected_indices[2][2]+l;
    active[work[0]] = true;
    active[work[1]] = true;
    active[work[2]] = true;
  }
  else
    if(gap[1]>=gap[0]) {
      working_set_type = ALPHAS;
      working_set_size = 2;
      work[0] = selected_indices[1][0];
      work[1] = selected_indices[1][1];
      active[work[0]] = true;
      active[work[1]] = true;
    }
    else {
      working_set_type = BETAS;
      working_set_size = 2;
      work[0] = selected_indices[0][0]+l;
      work[1] = selected_indices[0][1]+l;
      active[work[0]] = true;
      active[work[1]] = true;
    }

  Qfloat *Q_i, *Q_i_star, *Q_i_star_beta, *Q_j_star_beta, *Q_j, *Q_j_star, *Q_k_star, *Q_k_star_beta;
  int work0 = work[0], work1 = work[1], work2, true_k, act_set_k;

  switch(working_set_type) {
  case BETAS:
    best_u[0] = -1;
    best_u[1] = 1;
    work0 -= l;
    work1 -= l;
    Q_i_star_beta = Q_star_beta->get_Q(work0,active_size_beta);
    Q_j_star_beta = Q_star_beta->get_Q(work1,active_size_beta);     
    lambda_star = (g_beta[work0]-g_beta[work1])/(Q_i_star_beta[work0]+Q_j_star_beta[work1]-2*Q_i_star_beta[work1]);
    lambda_star = min(lambda_star, beta[work0]);
    break;

  case ALPHAS:
    best_u[0] = -1;
    best_u[1] = 1;
    Q_i = Q->get_Q(work0,active_size);
    Q_j = Q->get_Q(work1,active_size);
    Q_i_star = Q_star->get_Q(work0,active_size);
    Q_j_star = Q_star->get_Q(work1,active_size);
    lambda_star = y[work0]*G[work0]-y[work1]*G[work1]+(g[work0]-g[work1])/tau;
    lambda_star /= Q_i[work0]+Q_j[work1]-2*Q_i[work1]+(Q_i_star[work0]+Q_j_star[work1]-2*Q_i_star[work1])/tau;
    lambda_star = min(lambda_star, alpha[work[0]]);
    break;

  case ALPHAS_BETAS:
    work2 = work[2]-l;
    Q_i = Q->get_Q(work0,active_size);
    Q_j = Q->get_Q(work1,active_size);
    Q_i_star = Q_star->get_Q(work0,active_size);
    Q_j_star = Q_star->get_Q(work1,active_size);
    Q_k_star_beta = Q_star_beta->get_Q(work2,active_size_beta);
    true_k = active_set_beta[work2];
    act_set_k = true_act_set[true_k];
    Q_k_star = Q_star->get_Q(act_set_k, active_size);
    lambda_star = y[work0]*G[work0]+y[work1]*G[work1]-2+(g[work0]+g[work1]-2*g_beta[work2])/tau;
    lambda_star /= Q_i[work0]+Q_j[work1]-2*Q_i[work1]+(Q_i_star[work0]+Q_j_star[work1]+2*Q_i_star[work1]-4*Q_k_star[work0]-4*Q_k_star[work1]+4*Q_k_star_beta[work2])/tau;

    if(alpha_beta_type == 0) {
      best_u[0] = 1;
      best_u[1] = 1;
      best_u[2] = -2;
      lambda_star = min(-lambda_star, beta[work2]/2);
    }
    else {
      best_u[0] = -1;
      best_u[1] = -1;
      best_u[2] = 2;
      lambda_star = min(lambda_star, min(alpha[work0],alpha[work1]));
    }
  }

  return 0;
}

bool Solver_plus::feasible_direction(double *u)
{
  int i, work_i;
  int new_working_set_size = working_set_size+1;
  double *alpha = alpha_cg[0];
  double *beta = beta_cg[0]; 

  for(i=0; i<new_working_set_size; i++) {
    work_i = work[i];
    if(work_i<l) {
      if(alpha[work_i]<=1e-8 && u[i]<0)
        return false;
    }
    else {
      if(beta[work_i-l]<=1e-8 && u[i]<0)
        return false;
    }
  }

  return true;
}

bool Solver_plus::compute_gain(int new_working_set_size, double &gain, double *best_u, double &lambda, bool new_working_set)
{
  bool result, negligible, feasible_u, feasible_minus_u;
  double nominator, nominator1, nominator2, denominator, ui;
  int i, j, k, worki, workj, true_j, y_i;
  Qfloat *Q_i, *Q_i_star;
  double u[new_working_set_size], minus_u[new_working_set_size];
  double *G, *g, *g_beta, *alpha, *beta; 
  double tmp1, tmp2;
  double best_u_i, best_u_j;
  double u_i;

  G = G_cg[0];
  g = g_cg[0];
  g_beta = g_beta_cg[0];
  alpha = alpha_cg[0];
  beta = beta_cg[0];

  if(working_set_type == ALPHAS_BETAS || working_set_type == ALPHAS_DIFF_SIGN) 
    result = generate_direction_y(u,new_working_set_size, new_working_set);
  else
    result = generate_direction(u,new_working_set_size, new_working_set);
  if(result) 
    return false;

  for(i=0; i<new_working_set_size; i++) {
    u_i = u[i];
    if(u_i == 0) 
      return false;
    minus_u[i] = -u_i;
  }

  feasible_u = feasible_direction(u);
  feasible_minus_u = feasible_direction(minus_u);

  if(!feasible_u && !feasible_minus_u) 
    return false;  

  // choose which of {u,minus_u} is a descent direction
  nominator1 = 0;
  nominator2 = 0;
  nominator = 0;
  for(i = 0; i<new_working_set_size; i++) {
    worki = work[i];
    if(worki<l) {
      u_i = u[i];
      tmp1 = (y[worki]*G[worki]-1)*u_i;
      tmp2 = g[worki]*u_i;
      nominator += tmp1+tmp2/tau;
      nominator1 += tmp1;
      nominator2 += tmp2;
    }
    else {
      tmp1 = g_beta[worki-l]*u[i];
      nominator += tmp1/tau;
      nominator2 += tmp1;
    }
  }

  if(fabs(nominator1)<eps && fabs(nominator2)<eps)
    return false;

  if(fabs(nominator)<eps)
    return false;

  if(nominator<0) 
    if(feasible_u) 
      memcpy(best_u, u, sizeof_double*new_working_set_size);
    else
      return false;       
  else
    if(feasible_minus_u) {
      nominator *= -1;
      memcpy(best_u, minus_u, sizeof_double*new_working_set_size);    
    }
    else
      return false;

  // compute denominator
  denominator = 0;
  for(i = 0; i<new_working_set_size; i++) {
    worki = work[i];
    best_u_i = best_u[i];
    y_i = y[worki];
    if(worki<l) {
      Q_i = Q->get_Q(worki,active_size);
      Q_i_star = Q_star->get_Q(worki,l); //active_size);
      denominator += best_u_i*best_u_i*(QD[worki]+QD_star[worki]/tau);
      for(j = i+1; j<new_working_set_size; j++) {
	workj = work[j];
        best_u_j = best_u[j];
	if(workj<l)
	  denominator += 2*best_u_i*best_u_j*(y_i*y[workj]*Q_i[workj]+Q_i_star[workj]/tau);
	else {
	  k = active_set_beta[workj-l];
	  true_j = true_act_set[k];
	  denominator += 2*best_u_i*best_u_j*Q_i_star[true_j]/tau;
	}
      }
    }
    else {
      worki -= l;
      Q_i_star = Q_star_beta->get_Q(worki,l); //active_size_beta);
      denominator += best_u_i*best_u_i*QD_star_beta[worki]/tau;
      for(j = i+1; j<new_working_set_size; j++) {
	workj = work[j];
        best_u_j = best_u[j];
	if(workj<l) {
	  k = active_set[workj];
	  true_j = true_act_set_beta[k];
	  denominator += 2*best_u_i*best_u_j*Q_i_star[true_j]/tau;
	}
	else
	  denominator += 2*best_u_i*best_u_j*Q_i_star[workj-l]/tau;
      }
    }
  }

  if(denominator<=0)
    denominator=TAU;
  
  lambda = -nominator/denominator;
  for(i=0; i<new_working_set_size; i++) {
    ui = best_u[i];
    worki = work[i];
    if(ui>0)
      if(worki<l)  
	lambda = max(lambda, -alpha[worki]/ui);
      else
	lambda = max(lambda, -beta[worki-l]/ui);
    else 
      if(worki<l)  
	lambda = min(lambda, -alpha[worki]/ui);
      else
	lambda = min(lambda, -beta[worki-l]/ui);
  }

  if(fabs(lambda)<1e-3) 
    return false;
  negligible = false;
  for(i=0; i<new_working_set_size; i++)       
    if(fabs(lambda*best_u[i])<1e-3) {
      negligible = true;
      break;
    }
  if(negligible) 
    return false;

  gain = nominator*nominator/(2*denominator);
  return true;
}

int Solver_plus::select_working_set_plus_incrementally(double *best_u, double& lambda_star, double& best_gain)
{
  int k, best_k, new_working_set_size = working_set_size+1, n_betas=0, n_pos=0, n_neg=0;
  double gain, lambda;
  double *curr_u; 
  int work_k, work0, y_k, kl;
  int new_working_set_size_bytes = sizeof_double*new_working_set_size;
  bool chose_alpha = false, chose_beta = false;

  curr_u = new double[new_working_set_size];
  best_gain = 0;

  if(working_set_type == ALPHAS_BETAS) {
    for(k=0; k<working_set_size; k++) {
      work_k = work[k];
      if(work_k<l) {
	if(y[work_k]==-1)
	  n_neg++;
	else
	  n_pos++;
      }
      else
	n_betas++;
    }
  }   

  if(working_set_type == ALPHAS || (working_set_type == ALPHAS_BETAS && (n_betas>0 || n_pos+n_neg >= 3)) || working_set_type == ALPHAS_DIFF_SIGN) {  

    work0 = work[0];
    for(k = 0; k < active_size; k++) {

      y_k = y[k];
      if(working_set_type==ALPHAS && y_k != y[work0])
        continue;

      if(working_set_type==ALPHAS_BETAS) {
        if(n_pos == 0 && y_k == -1)
          continue;

        if(n_neg == 0 && y_k == 1)
          continue;

        if(n_betas == 0) {
          if(n_pos == 1 && y_k == -1)
            continue;
          if(n_neg == 1 && y_k == 1)
            continue;
	}
      }

      if(working_set_type==ALPHAS_DIFF_SIGN) {
        if(n_pos == 1 && y_k == -1)
          continue;
        if(n_pos == -1 && y_k == 1)
          continue;
      }

      // check that k is not in the working set
      if(active[k])
	continue;

      // generate direction
      work[working_set_size] = k;

      if(!compute_gain(new_working_set_size, gain, curr_u, lambda, k==0))
        continue;      
      
      if(gain > best_gain) {
        chose_alpha = true;
        best_gain = gain;
        best_k = k;
        memcpy(best_u, curr_u, new_working_set_size_bytes);
        lambda_star = lambda;
      }
    }
  }

  if(working_set_type == BETAS || (working_set_type == ALPHAS_BETAS && n_pos > 0 && n_neg > 0) || working_set_type == ALPHAS_DIFF_SIGN) 
    for(k = 0; k < active_size_beta; k++) {

      kl=k+l;
      // check that k is not in the working set
      if(active[kl])
	continue;

      // generate direction
      work[working_set_size] = kl;
     
      if(!compute_gain(new_working_set_size, gain, curr_u, lambda, best_gain==0 && k==0))
        continue;       

      if(gain > best_gain) {
        chose_beta = true;
        chose_alpha = false;
        best_gain = gain;
        best_k = kl;
        memcpy(best_u, curr_u, new_working_set_size_bytes);
        lambda_star = lambda;
      }
    }

  if(best_gain == 0) {
    delete[] curr_u;
    return 1;
  }

  work[working_set_size] = best_k;
  active[best_k] = true;
  working_set_size++;

  if(working_set_type == ALPHAS_DIFF_SIGN && chose_beta)
    working_set_type = ALPHAS_BETAS;
  else
    if(working_set_type == ALPHAS_BETAS && chose_alpha && (n_betas == 0))
      working_set_type = ALPHAS_DIFF_SIGN;

  delete[] curr_u;
  return 0;
}

int Solver_plus::select_working_set_plus_hmg(double *best_u, double& lambda_star, double& incr_gain)
{
  int i;

  if(prev_depth <  max_depth) 
    return select_working_set_plus_incrementally(best_u, lambda_star, incr_gain);

  // throw away the oldest active example and select the working set incrementally
  int working_set_size_minus_one = working_set_size-1;
  active[work[0]] = false;
  for(i=0; i<working_set_size_minus_one; i++)
    work[i] = work[i+1];
  work[working_set_size-1]=-1;
  working_set_size--;

  return select_working_set_plus_incrementally(best_u, lambda_star, incr_gain);
}

// return 1 if already optimal, return 0 otherwise
int Solver_plus::select_working_set_plus_hmg2(double *best_u, double& lambda_star, double& absolute_best_z)
{
  int i, j, best_B1=-1, best_B1B2=-1, best_A2=-1, best_A2A4=-1, best_A1=-1, best_A1A3=-1, i_ind, j_ind, k_ind;
  int type_selected[3], selected_indices[3][3];
  double max_B1=-1e20, min_B1B2=1e20, max_A2=-1e20, min_A2A4=1e20, max_A1=-1e20, min_A1A3=1e20;
  double alpha_i, G_i, g_i, g_j, y_i, y_j, deriv_alpha_i, first_order_criterion;
  double max_z[3], z, nominator, nominator_base, denominator, denominator_base, j_deriv, Q_star_ii;
  int best_z_index[3], true_k, act_set_k;
  Qfloat *Q_i, *Q_i_star, *Q_k_star, *Q_k_star_beta;
  double *alpha, *G, *g, *g_beta, *beta, lambda[3];
  double nominator_base1, nominator_base2, j_deriv1, j_deriv2, nominator1, nominator2;
  double max_A2_1, max_A2_2, min_A2A4_1, min_A2A4_2, max_A1_1, max_A1_2, min_A1A3_1, min_A1A3_2;
  double deriv_alpha_i1, deriv_alpha_i2;


  absolute_best_z = 0;
  // first-order working set selection 

  for(i=0; i<working_set_size; i++) {
    active[work[i]] = false;
    work[i] = -1;
  }
  working_set_size = 0;
  curr_depth = 0;

  // compute all maxima and minima related to alphas
  alpha = alpha_cg[0];
  beta = beta_cg[0];
  g = g_cg[0];
  G = G_cg[0];
  g_beta = g_beta_cg[0];
  for(i=0; i<active_size; i++) {
    alpha_i = alpha[i];
    G_i = G[i];
    g_i = g[i];
    y_i = y[i];
    deriv_alpha_i1 = y_i*G[i];
    deriv_alpha_i2 = g_i;
    deriv_alpha_i = deriv_alpha_i1+deriv_alpha_i2/tau;

    // max A2
    if(alpha_i>1e-8 && y_i==-1 && deriv_alpha_i>max_A2) {
      max_A2 = deriv_alpha_i;
      best_A2 = i;
      max_A2_1 = deriv_alpha_i1;
      max_A2_2 = deriv_alpha_i2;
    }

    // min A2A4
    if(y_i==-1 && deriv_alpha_i<min_A2A4) {
      min_A2A4 = deriv_alpha_i;
      best_A2A4 = i;
      min_A2A4_1 = deriv_alpha_i1;
      min_A2A4_2 = deriv_alpha_i2;
    }

    // max A1
    if(alpha_i>1e-8 && y_i==1 && deriv_alpha_i>max_A1) {
      max_A1 = deriv_alpha_i;
      best_A1 = i;
      max_A1_1 = deriv_alpha_i1;
      max_A1_2 = deriv_alpha_i2;
    }
    
    // min A1A3max_A2, min_A2A4, max_A1, min_A1A3
    if(y_i==1 && deriv_alpha_i<min_A1A3) {
      min_A1A3 = deriv_alpha_i;
      min_A1A3_1 = deriv_alpha_i1;
      min_A1A3_2 = deriv_alpha_i2;      
      best_A1A3 = i;
    }
  }

  // compute all maxima and minima related to betas
  for(i=0; i<active_size_beta; i++) {
    g_i = g_beta[i];

    // max B1
    if(beta[i]>1e-8 && g_i>max_B1) {
      max_B1 = g_i;
      best_B1 = i;
    }

    // min B1B2
    if(g_i<min_B1B2) {
      min_B1B2 = g_i;
      best_B1B2 = i;
    }
  }

  max_B1 /= tau;
  min_B1B2 /= tau;

  // select maximal violating pairs
  if(max_B1-min_B1B2 < eps)
    type_selected[0] = 0;
  else {
    type_selected[0] = 1;
    selected_indices[0][0] = best_B1;
    selected_indices[0][1] = best_B1B2;
  }

  if(((max_A2 - min_A2A4 < eps) || ((max_A2_1 - min_A2A4_1 < eps) && (max_A2_2 - min_A2A4_2 < eps))) && 
     ((max_A1 - min_A1A3 < eps) || ((max_A1_1 - min_A1A3_1 < eps) && (max_A1_2 - min_A1A3_2 < eps))))
    type_selected[1] = 0;
  else {
    if((max_A2 - min_A2A4 > max_A1 - min_A1A3) && ((max_A2_1 - min_A2A4_1 >= eps) || (max_A2_2 - min_A2A4_2 >= eps))) {
      type_selected[1] = 1;
      selected_indices[1][0] = best_A2;
      selected_indices[1][1] = best_A2A4;
    }
    else {
      if((max_A2 - min_A2A4 <= max_A1 - min_A1A3) && ((max_A1_1 - min_A1A3_1 >= eps) || (max_A1_2 - min_A1A3_2 >= eps))) {
        type_selected[1] = 1;
        selected_indices[1][0] = best_A1;
        selected_indices[1][1] = best_A1A3;
      }
      else
        type_selected[1] = 0;
    }
  }

  if(((2*max_B1+2-min_A1A3-min_A2A4 < eps) || ((2-min_A1A3_1-min_A2A4_1<eps) && (2*max_B1*tau-min_A1A3_2-min_A2A4_2<eps))) && 
     ((max_A1+max_A2-2*min_B1B2-2 < eps) || ((max_A1_1+max_A2_1-2<eps) && (max_A1_2+max_A2_2-2*min_B1B2*tau<eps))))
    type_selected[2] = 0;
  else {
    if((2*max_B1+2-min_A1A3-min_A2A4 > max_A1+max_A2-2*min_B1B2-2) && ((2-min_A1A3_1-min_A2A4_1>=eps) || (2*max_B1*tau-min_A1A3_2-min_A2A4_2>=eps))) {
      type_selected[2] = 1;
      selected_indices[2][0] = best_A1A3;
      selected_indices[2][1] = best_A2A4;
      selected_indices[2][2] = best_B1;
    }
    else {
      if((2*max_B1+2-min_A1A3-min_A2A4 <= max_A1+max_A2-2*min_B1B2-2) && ((max_A1_1+max_A2_1-2>=eps) || (max_A1_2+max_A2_2-2*min_B1B2*tau>=eps))) {
        type_selected[2] = 1;
        selected_indices[2][0] = best_A1;
        selected_indices[2][1] = best_A2;
        selected_indices[2][2] = best_B1B2;
      }
      else
        type_selected[2] = 0;
    }
  }

  if(type_selected[0]+type_selected[1]+type_selected[2] == 0) 
    return 1;

  for(i=0; i<3; i++)
    max_z[i] = -1e20;    

  // second-order working set selection
  if(type_selected[0] == 1) {
     i_ind = selected_indices[0][0];
     g_i = g_beta[i_ind];
     Q_i_star = Q_star_beta->get_Q(i_ind,active_size_beta);
     Q_star_ii = Q_i_star[i_ind];
     for(j = 0; j < active_size_beta; j++) {
        g_j = g_beta[j];
	if(eps+g_j/tau < g_i/tau) {
	  nominator = g_i-g_j;
          denominator = Q_star_ii+QD_star_beta[j]-2*Q_i_star[j]; 
	  z = nominator*nominator/(2*tau*denominator);
	  if(z > max_z[0]) {
	    max_z[0] = z;
	    best_z_index[0] = j;
            lambda[0] = nominator/denominator;
	  }
	}
      }
  }

  if(type_selected[1] == 1) {
     i_ind = selected_indices[1][0];
     y_i = y[i_ind];
     Q_i = Q->get_Q(i_ind,active_size);
     Q_i_star = Q_star->get_Q(i_ind,active_size);
     nominator_base = y_i*G[i_ind] + g[i_ind]/tau;
     nominator_base1 = y_i*G[i_ind];
     nominator_base2 = g[i_ind];
     denominator_base = 2*(Q_i[i_ind]+Q_i_star[i_ind]/tau);
   
     for(j = 0; j < active_size; j++) { 
        y_j = y[j];
        j_deriv = y_j*G[j]+g[j]/tau;
        j_deriv1 = y_j*G[j];
        j_deriv2 = g[j];
        
	if(y_j == y_i && j_deriv+eps < nominator_base && ((j_deriv1+eps < nominator_base1) || (j_deriv2+eps < nominator_base2))) {        
          j_deriv = j_deriv1+j_deriv2/tau;
	  nominator = nominator_base-j_deriv;
          denominator = denominator_base + 2*(QD[j]-2*Q_i[j]+(QD_star[j]-2*Q_i_star[j])/tau);
	  z = nominator*nominator/denominator;
	  if(z > max_z[1]) {
	    max_z[1] = z;
	    best_z_index[1] = j;
            lambda[1] = nominator/(denominator/2);
	  }
	}
    }
  }

  if(type_selected[2] == 1) {
    i_ind = selected_indices[2][0];
    j_ind = selected_indices[2][1];
    k_ind = selected_indices[2][2];
    Q_i = Q->get_Q(i_ind,active_size);
    Q_i_star = Q_star->get_Q(i_ind,active_size);
    Q_k_star_beta = Q_star_beta->get_Q(k_ind,active_size_beta);

    true_k = active_set_beta[k_ind];
    act_set_k = true_act_set[true_k];
    Q_k_star = Q_star->get_Q(act_set_k, active_size);

    nominator_base = y[i_ind]*G[i_ind]-2+(g[i_ind]-2*g_beta[k_ind])/tau;
    nominator_base1 = y[i_ind]*G[i_ind]-2;
    nominator_base2 = g[i_ind]-2*g_beta[k_ind];
    denominator_base = 2*(Q_i[i_ind]+(Q_i_star[i_ind]-4*Q_k_star[i_ind]+4*Q_k_star_beta[k_ind])/tau);
    first_order_criterion = nominator_base+y[j_ind]*G[j_ind]+g[j_ind]/tau;
    for(j = 0; j < active_size; j++) {
      if(y[j] == -1) {        
        nominator1 = nominator_base1+y[j]*G[j];
        nominator2 = nominator_base2+g[j];
        nominator = nominator_base+y[j]*G[j]+g[j]/tau;
        if((first_order_criterion < 0 && nominator<-eps && ((nominator1 < -eps) || (nominator2 < -eps))) ||
           (first_order_criterion > 0 && alpha[j] > 1e-8 && nominator > eps && ((nominator1 > eps) || (nominator2 > eps)))) {
          denominator = denominator_base + 2*(QD[j]-2*Q_i[j]+(QD_star[j]+2*Q_i_star[j]-4*Q_k_star[j])/tau);
	  z = nominator*nominator/denominator;
	  if(z > max_z[2]) {
	    max_z[2] = z;
	    best_z_index[2] = j;
            lambda[2] = nominator/(denominator/2);
	  }
        }
      }
    }
  }

  // choose the best type
  absolute_best_z = -1;
  for(i=0; i<3; i++) {
    if((type_selected[i] == 1) && (max_z[i] > absolute_best_z)) {  
      absolute_best_z = max_z[i];
      working_set_type = (char)i;
      work[0] = selected_indices[i][0];
      work[1] = best_z_index[i];
      if(i==0) {
        work[0] += l;
        work[1] += l;
      }
      if(i == 2)
	work[2] = selected_indices[i][2]+l;
      lambda_star = lambda[i];
    }
  }

  active[work[0]] = true;
  active[work[1]] = true;
  if(working_set_type==2) {
    active[work[2]] = true;
    working_set_size = 3;
  }
  else
    working_set_size = 2;

  if(absolute_best_z == -1) {
    working_set_size = -1;
    return 1;
  }

  switch(working_set_type) {
  case BETAS:
    best_u[0] = -1;
    best_u[1] = 1;
    lambda_star = min(lambda_star, beta[work[0]-l]);
    break;
  case ALPHAS:
    best_u[0] = -1;
    best_u[1] = 1;
    lambda_star = min(lambda_star, alpha[work[0]]);
    break;
  case ALPHAS_BETAS:
    if(lambda_star>0) {
      best_u[0] = -1;
      best_u[1] = -1;
      best_u[2] = 2;
      lambda_star = min(lambda_star, min(alpha[work[0]],alpha[work[1]]));
    }
    else {
      lambda_star = -lambda_star;
      best_u[0] = 1;
      best_u[1] = 1;
      best_u[2] = -2;
      lambda_star = min(lambda_star, beta[work[2]-l]/2);
    }
  }
  
  return 0;
}

void Solver_plus::reconstruct_gradient_plus()
{
  int i, j, true_i, act_set_i;

  if(active_size < l) {
    for(i=active_size;i<l;i++) {
      const Qfloat *Q_i = Q->get_Q(i,l);
      const Qfloat *Q_i_star = Q_star->get_Q(i,l);
       
      true_i = active_set[i];
      act_set_i = true_act_set_beta[true_i];

      const Qfloat *Q_i_star_beta = Q_star_beta->get_Q(act_set_i,l);
      G[i] = 0;
      g[i] = g_init[i];
      for(j=0;j<l;j++)
	if(alpha[j]>1e-8) {
	  G[i] += alpha[j] * y[j] * Q_i[j];
	  g[i] += alpha[j] * Q_i_star[j];
	}
      for(j=0;j<l;j++)
	if(beta[j]>1e-8) 
	  g[i] += beta[j] * Q_i_star_beta[j];
    }
  }

  if(active_size_beta < l) {
    for(i=active_size_beta;i<l;i++) {
      const Qfloat *Q_i_star_beta = Q_star_beta->get_Q(i,l);
      
      true_i = active_set_beta[i];
      act_set_i = true_act_set[true_i];
      const Qfloat *Q_i_star = Q_star->get_Q(act_set_i,l);
       
      g_beta[i] = g_beta_init[i];
    
      for(j=0;j<l;j++)
	if(beta[j]>1e-8) 
	  g_beta[i] += beta[j] * Q_i_star_beta[j];

      for(j=0;j<l;j++)
	if(alpha[j]>1e-8) 
	  g_beta[i] += alpha[j] * Q_i_star[j];
    }
  }
}

int Solver_plus::select_working_set(double* u, double& lambda_star)
{
  static double gain_hmg, gain2, tmp_lambda;
  static int tmp_depth, tmp_working_set_size, i, result;
  static char tmp_working_set_type;

  if(select_working_set_plus_hmg(u, lambda_star, gain_hmg)!=0) { 
    // cannot find conjugate direction, try to go in the  direction of gradient 
    memcpy(tmp_work, work, sizeof_int*working_set_size);
    tmp_depth = curr_depth;
    tmp_working_set_size = working_set_size;
    tmp_working_set_type = working_set_type;
    for(i=0; i<working_set_size; i++) {
      active[work[i]] = false;
      work[i] = -1;
    }
    curr_depth = 0;
    working_set_size = 0;
    result = select_working_set_plus_hmg2(u, lambda_star, gain2);
    if(result == 1) {
      curr_depth = tmp_depth;
      working_set_size = tmp_working_set_size;
      working_set_type = tmp_working_set_type;
      memcpy(work, tmp_work, sizeof_int*tmp_working_set_size);
      for(i=0; i<working_set_size; i++)
	active[work[i]] = true;
    }
    return result;
  }
  else {
    // store current direction, step size and the working set
    memcpy(tmp_work, work, sizeof_int*working_set_size);
    memcpy(tmp_u, u, sizeof_double*working_set_size);
    tmp_depth = curr_depth;
    tmp_working_set_size = working_set_size;
    tmp_working_set_type = working_set_type;
    tmp_lambda = lambda_star;
    for(i=0; i<working_set_size; i++) {
      active[work[i]] = false;
      work[i] = -1;
    }
    curr_depth = 0;
    working_set_size = 0;
    if(select_working_set_plus_hmg2(u, lambda_star, gain2) == 0) {
      if(gain_hmg>gain2) {
	for(i=0; i<working_set_size; i++)
	  active[work[i]] = false;
	curr_depth = tmp_depth;
	working_set_size = tmp_working_set_size;
	working_set_type = tmp_working_set_type;
	lambda_star = tmp_lambda;
	memcpy(work, tmp_work, sizeof_int*tmp_working_set_size);
	memcpy(u, tmp_u, sizeof_double*working_set_size);
	for(i=0; i<working_set_size; i++)
	  active[work[i]] = true;
      }
    }
    else {
      curr_depth = tmp_depth;
      working_set_size = tmp_working_set_size;
      working_set_type = tmp_working_set_type;
      lambda_star = tmp_lambda;
      memcpy(u, tmp_u, sizeof_double*working_set_size);
      memcpy(work, tmp_work, sizeof_int*tmp_working_set_size);
      for(i=0; i<working_set_size; i++)
	active[work[i]] = true;
    }
    return 0;
  }
}

void Solver_plus::Solve_plus_cg(int l, const QMatrix& Q, const QMatrix& Q_star, const QMatrix& Q_star_beta, const schar *y_,
				double *alpha_, double *beta_, double Cp, double Cn, double tau_, double eps,
				SolutionInfo* si, int shrinking)
{
  int i,j;

  this->l = l;
  this->Q = &Q;
  this->Q_star = &Q_star;
  this->Q_star_beta = &Q_star_beta;
  QD = Q.get_QD();
  QD_star = Q_star.get_QD();
  QD_star_beta = Q_star_beta.get_QD();
  clone(alpha,alpha_,l);
  clone(beta,beta_,l);
  this->Cp = Cp;
  this->Cn = Cn;
  this->eps = eps;
  tau = tau_;
  unshrink = false;
  prev_depth = -1;
  curr_depth = 0;
  working_set_size = 0;

  int l2 = 2*l;
  double *G,*g,*g_beta;
  bool done_shrinking = false;

  y = new schar[l2];
  memcpy(y,y_, sizeof_char*l);
  for(i = l; i<l2; i++)
    y[i] = 0;

  work = new int[max_depth+3];
  for(i=0; i<max_depth+3; i++)
    work[i] = -1;

  active = new bool[l2];
  for(i=0; i<l2; i++)
    active[i] = false;

  // initialize alpha's 
  alpha_cg = new double*[max_depth+1];
  for(i = 0; i < max_depth + 1; i++) 
    clone(alpha_cg[i],alpha_,l);

  // initialize beta's 
  beta_cg = new double*[max_depth+1];
  for(i = 0; i < max_depth + 1; i++) 
    clone(beta_cg[i],beta_,l);

  // initialize alpha_status
  alpha_status_cg = new char*[max_depth+1];
  for(i=0; i < max_depth + 1; i++) {
    alpha_status_cg[i] = new char[l];
    for(j=0;j<l;j++) {
      update_alpha_status_cg(j,i);
    }
  }

  // initialize beta_status
  beta_status_cg = new char*[max_depth+1];
  for(i=0; i < max_depth + 1; i++) {
    beta_status_cg[i] = new char[l];
    for(j=0;j<l;j++) {
      update_beta_status_cg(j,i);
    }
  }

  // initialize gradient
  G_cg =  new double*[max_depth+1];
  g_cg =  new double*[max_depth+1];
  g_beta_cg =  new double*[max_depth+1];

  G_cg[0] = new double[l];
  g_cg[0] = new double[l];
  g_beta_cg[0] = new double[l];
  g_init = new double[l];
  g_beta_init = new double[l];

  G = G_cg[0];
  g = g_cg[0];
  g_beta = g_beta_cg[0];

  for(i=0; i<l; i++) {
    G[i] = 0;
    g[i] = 0;
    g_init[i] =0;
  }

  for(i=0;i<l;i++) {
    const Qfloat *Q_i_star = Q_star.get_Q(i,l);
	  
    for(j=0; j<l; j++) {
      g[j] -= Cp*Q_i_star[j];
      g_init[j] -= Cp*Q_i_star[j];
    }

    if(!is_lower_bound_cg(i,0)) {  
      const Qfloat *Q_i = Q.get_Q(i,l);
      double alpha_i = alpha[i];
      double y_i = y[i];
      for(j=0;j<l;j++) {
	G[j] += alpha_i*y_i*Q_i[j];
	g[j] += alpha_i*Q_i_star[j];
      }
    }
	    
    if(!is_lower_bound_beta_cg(i,0)) {
      double beta_i = beta[i];
      for(j=0;j<l;j++) 
	g[j] += beta_i*Q_i_star[j];    
    }
  }

  clone(g_beta, g, l);
  clone(g_beta_init, g_init, l);
  for(i=1; i<max_depth+1; i++) {
    clone(G_cg[i], G, l);
    clone(g_cg[i], g, l);
    clone(g_beta_cg[i], g_beta, l);
  }

  active_set = new int[l];
  active_set_beta = new int[l];
  true_act_set = new int[l];
  true_act_set_beta = new int[l];
  for(i=0; i<l; i++) {
    active_set[i] = i;
    active_set_beta[i] = i;
    true_act_set[i] = i;
    true_act_set_beta[i] = i;
  }
  active_size = l;
  active_size_beta = l;

  int iter=0, counter = min(l,1000)+1, next_depth, n_conjugate=0, worki, work0;
  bool corner = false;
  double *tmp_alpha, *tmp_beta, *tmp_G, *tmp_g, *tmp_g_beta;
  char *tmp_status, *tmp_status_beta;
  double *u, lambda_star, *alpha, *beta;
  Qfloat *Q_i_star_beta, *Q_i, *Q_i_star, *Q_k_star, *Q_k_star_beta;
  double *alpha_old, *beta_old, diff_i, diff_i_y, diff_k;
  int r, true_i, act_set_i,  true_k, act_set_k;
  double gain2;

  tmp_work = new int[max_depth+3];
  u =  new double[max_depth+3];
  tmp_u =  new double[max_depth+3];

  while(1) { 

    if(--counter == 0) {
      counter = min(l,1000);
      if(shrinking) 
        done_shrinking = do_shrinking_plus_cg();
    }

    lambda_star = 0;

    if(corner) {
      if(wss_first_order(u, lambda_star)!=0) {
 	reconstruct_gradient_plus_cg();
 	active_size = l;
	active_size_beta = l;
	if(wss_first_order(u, lambda_star)!=0) 
	  break;
	else
	  counter = 1;
      }
    }
    else 
      if(curr_depth > 0) {
       	if(select_working_set(u, lambda_star)!=0) { 
	  // reconstruct the whole gradient
	  reconstruct_gradient_plus_cg();
	  // reset active set size 
	  active_size = l;
	  active_size_beta = l;
	  if(select_working_set(u, lambda_star)==0)
	    counter = 1;
	  else 
	    break;
	}
      }
      else {
	if(select_working_set_plus_hmg2(u, lambda_star, gain2)!=0) {
	  reconstruct_gradient_plus_cg();
 	  active_size = l;
	  active_size_beta = l;
          if(select_working_set_plus_hmg2(u, lambda_star, gain2)==0)
	    counter = 1;
          else
	    break;  
	}
      }

    iter++;
    fprintf(stdout,"iter=%d\n",iter);
    fflush(stdout);
    prev_depth = curr_depth;

    // shift old alpha's, G's and G_bar's
    n_conjugate++;
    next_depth = min(curr_depth+1,max_depth);
    tmp_alpha = alpha_cg[next_depth];
    tmp_beta = beta_cg[next_depth];
    tmp_G = G_cg[next_depth];
    tmp_g = g_cg[next_depth];
    tmp_g_beta = g_beta_cg[next_depth];
    tmp_status = alpha_status_cg[next_depth];
    tmp_status_beta = beta_status_cg[next_depth];

    for(i=next_depth; i>0; i--) {
      alpha_cg[i] = alpha_cg[i-1];
      beta_cg[i] = beta_cg[i-1];
      G_cg[i] = G_cg[i-1];
      g_cg[i] = g_cg[i-1];
      g_beta_cg[i] = g_beta_cg[i-1];
      alpha_status_cg[i] = alpha_status_cg[i-1];
      beta_status_cg[i] = beta_status_cg[i-1];
    }

    if(!done_shrinking || (curr_depth == max_depth && curr_depth == prev_depth)) {
      memcpy(tmp_alpha, alpha_cg[0], active_size*sizeof_double);
      memcpy(tmp_beta, beta_cg[0], active_size_beta*sizeof_double);
      memcpy(tmp_G, G_cg[0], active_size*sizeof_double);
      memcpy(tmp_g, g_cg[0], active_size*sizeof_double);
      memcpy(tmp_g_beta, g_beta_cg[0], active_size_beta*sizeof_double);
      memcpy(tmp_status, alpha_status_cg[0], active_size*sizeof_char);
      memcpy(tmp_status_beta, beta_status_cg[0], active_size_beta*sizeof_char);
      done_shrinking = false;
    }
    else {
      memcpy(tmp_alpha, alpha_cg[0], l*sizeof_double);
      memcpy(tmp_beta, beta_cg[0], l*sizeof_double);
      memcpy(tmp_G, G_cg[0], l*sizeof_double);
      memcpy(tmp_g, g_cg[0], l*sizeof_double);
      memcpy(tmp_g_beta, g_beta_cg[0], l*sizeof_double);
      memcpy(tmp_status, alpha_status_cg[0], l*sizeof_char);
      memcpy(tmp_status_beta, beta_status_cg[0], l*sizeof_char);
    }

    alpha_cg[0] = tmp_alpha;
    beta_cg[0] = tmp_beta;
    G_cg[0] = tmp_G;
    g_cg[0] = tmp_g;
    g_beta_cg[0] = tmp_g_beta;
    alpha_status_cg[0] = tmp_status;
    beta_status_cg[0] = tmp_status_beta;
    curr_depth = next_depth;  

    alpha = alpha_cg[0];
    beta = beta_cg[0];
    g_beta = g_beta_cg[0];
    g = g_cg[0];
    G = G_cg[0];
    alpha_old = alpha_cg[1];
    beta_old = beta_cg[1];

    for(i=0; i<working_set_size; i++) {
      worki = work[i];
      if(worki<l) {
	alpha[worki] += lambda_star*u[i];
	update_alpha_status_cg(worki,0);
      }
      else  {
	beta[worki-l] += lambda_star*u[i];
	update_beta_status_cg(worki-l,0);
      }
    }
    fflush(stdout);

    // check if we are at the corner
    corner = true;
    for(i=0; i<working_set_size; i++) {
      worki = work[i];
      if(worki<l) {
	if(alpha[worki] >= 1e-8) {
	  corner = false;
	  break;
	}
      }
      else
	if(beta[worki-l] >= 1e-8) {
	  corner = false;
	  break;
	}
    }

    // update gradients
    switch(working_set_type) {
    case BETAS:
      for(i=0; i<working_set_size; i++) {
	work0 = work[i]-l;
	Q_i_star_beta = Q_star_beta.get_Q(work0,active_size_beta);
	diff_i = beta[work0]-beta_old[work0];
	for(r=0; r<active_size_beta; r++)   
	  g_beta[r] += diff_i*Q_i_star_beta[r]; 

	true_i = active_set_beta[work0];
	act_set_i = true_act_set[true_i];
	Q_i_star = Q_star.get_Q(act_set_i,active_size);
    
	for(r=0; r<active_size; r++) 
	  g[r] += diff_i*Q_i_star[r];
      }
      break;

    case ALPHAS:
    case ALPHAS_DIFF_SIGN:
      for(i=0; i<working_set_size; i++) {
	work0 = work[i];
	Q_i = Q.get_Q(work0,active_size);
	Q_i_star = Q_star.get_Q(work0,active_size); 
	diff_i = alpha[work0]-alpha_old[work0];
	diff_i_y = diff_i * y[work0];
     
	for (r=0; r<active_size; r++) {
	  G[r] += diff_i_y*Q_i[r];
	  g[r] += diff_i*Q_i_star[r]; 
	}

	true_i = active_set[work0];
	act_set_i = true_act_set_beta[true_i];
	Q_i_star_beta = Q_star_beta.get_Q(act_set_i,active_size_beta);

	for (r=0; r<active_size_beta; r++) 
	  g_beta[r] += diff_i*Q_i_star_beta[r]; 
      }
      break;

    case ALPHAS_BETAS:   
      for(i=0; i<working_set_size; i++) 
	if(work[i]<l) {
	  work0 = work[i];
	  Q_i = Q.get_Q(work0,active_size);
	  Q_i_star = Q_star.get_Q(work0,active_size);
	  diff_i = alpha[work0]-alpha_old[work0];
	  diff_i_y = diff_i * y[work0];

	  for (r=0; r<active_size; r++) {
	    G[r] += diff_i_y*Q_i[r];
	    g[r] += diff_i*Q_i_star[r]; 
	  }

	  true_i = active_set[work0];
	  act_set_i = true_act_set_beta[true_i];
	  Q_i_star_beta = Q_star_beta.get_Q(act_set_i,active_size_beta);

	  for(r=0; r<active_size_beta; r++)
	    g_beta[r] += diff_i*Q_i_star_beta[r];

	  update_alpha_status_cg(work0,0);
	}
	else {
	  work0 = work[i] - l;
	  Q_k_star_beta = Q_star_beta.get_Q(work0,active_size_beta);
	  true_k = active_set_beta[work0];
	  act_set_k = true_act_set[true_k];
	  Q_k_star = Q_star.get_Q(act_set_k, active_size);
	  diff_k = beta[work0]-beta_old[work0];
	  for (r=0; r<active_size; r++) 
	    g[r] += diff_k*Q_k_star[r]; 
	  for(r=0; r<active_size_beta; r++)
	    g_beta[r] += diff_k*Q_k_star_beta[r];
	} 
    }
  }
  calculate_rho_plus_cg(si->rho,si->rho_star);

  // put back the solution
  alpha = alpha_cg[0];
  beta = beta_cg[0];
  for(i=0;i<l;i++) {
    alpha_[active_set[i]] = alpha[i];
    beta_[active_set_beta[i]] = beta[i];
  }

  info("\noptimization finished, #iter = %d\n",iter);

  si->rho *= -1;

/*
  for(i=0; i<max_depth+1; i++) {
    delete[] alpha_cg[i];
    delete[] beta_cg[i];
    delete[] alpha_status_cg[i];
    delete[] beta_status_cg[i];
    delete[] G_cg[i];
    delete[] g_cg[i];
    delete[] g_beta_cg[i];
  }
  
  delete[] alpha_cg;
  delete[] beta_cg;
  delete[] alpha_status_cg;
  delete[] beta_status_cg;
  delete[] G_cg;
  delete[] g_cg;
  delete[] g_beta_cg;
  delete[] g_init;
  delete[] g_beta_init;
  delete[] active;
  delete[] y;
  delete[] work;
  delete[] u;
  delete[] tmp_work;
  delete[] tmp_u;
*/
}

void Solver_plus::Solve_plus(int l, const QMatrix& Q, const QMatrix& Q_star, const QMatrix& Q_star_beta, const schar *y_,
			     double *alpha_, double *beta_, double Cp, double Cn, double tau_, double eps,
			     SolutionInfo* si, int shrinking)
{
  int i,j;
  struct tms init_time, fin_time;
  long t_ini, t_fin;
  double  net_time;

  this->l = l;
  this->Q = &Q;
  this->Q_star = &Q_star;
  this->Q_star_beta = &Q_star_beta;
  QD = Q.get_QD();
  QD_star = Q_star.get_QD();
  QD_star_beta = Q_star_beta.get_QD();
  clone(y, y_,l);
  clone(alpha,alpha_,l);
  clone(beta,beta_,l);
  this->Cp = Cp;
  this->Cn = Cn;
  this->eps = eps;
  tau = tau_;
  unshrink = false;

  alpha_status = new char[l];
  for(i=0;i<l;i++)
    update_alpha_status(i);
  beta_status = new char[l];
  for(i=0;i<l;i++)
    update_beta_status(i);

  t_ini = times( &init_time);

  // initialize gradient
  G = new double[l];
  g = new double[l];
  g_beta = new double[l];
  g_init = new double[l];
  g_beta_init = new double[l];

  for(i=0; i<l; i++) {
    G[i] = 0;
    g[i] = 0;
    g_init[i] =0;
  }

  for(i=0;i<l;i++) {
    const Qfloat *Q_i_star = Q_star.get_Q(i,l);
	  
    for(j=0; j<l; j++) {
      g[j] -= Cp*Q_i_star[j];
      g_init[j] -= Cp*Q_i_star[j];
    }

    if(!is_lower_bound(i)) {  
      const Qfloat *Q_i = Q.get_Q(i,l);
      double alpha_i = alpha[i];
      double y_i = y[i];
      for(j=0;j<l;j++) {
	G[j] += alpha_i*y_i*Q_i[j];
	g[j] += alpha_i*Q_i_star[j];
      }
    }
	    
    if(!is_lower_bound_beta(i)) {
      double beta_i = beta[i];
      for(j=0;j<l;j++) 
	g[j] += beta_i*Q_i_star[j];    
    }
  }
  for(i=0; i<l; i++) {
    g_beta[i] = g[i];
    g_beta_init[i] = g_init[i];
  }

  active_set = new int[l];
  active_set_beta = new int[l];
  true_act_set = new int[l];
  true_act_set_beta = new int[l];
  for(int i=0; i<l; i++) {
    active_set[i] = i;
    active_set_beta[i] = i;
    true_act_set[i] = i;
    true_act_set_beta[i] = i;
  }
  active_size = l;
  active_size_beta = l;

  int counter = min(l,1000)+1;

  // optimization step
  int iter = 0, y_i, y_j;
  Qfloat *Q_i, *Q_j, *Q_i_star, *Q_j_star, *Q_k_star, *Q_i_star_beta, *Q_j_star_beta, *Q_k_star_beta;
  double Delta, beta_i_old, beta_j_old, alpha_i_old, alpha_j_old, beta_k_old, nominator, denominator, min_alpha, alpha_change;
  double diff_i, diff_j, diff_k, beta_i, beta_k, alpha_i, diff_i_y, diff_j_y;
  int true_i, true_j, true_k, act_set_i, act_set_j, act_set_k;
  
  while(iter<1e7) {

    int i,j,k,set_type,r;

    if(--counter == 0) {
      counter = min(l,1000);
      if(shrinking) 
        do_shrinking_plus();
    }

    if(select_working_set_plus(set_type, i,j,k, iter) != 0) {
      
      // reconstruct the whole gradient
      reconstruct_gradient_plus();

      // reset active set size and check
      active_size = l;
      active_size_beta = l;
	  
      if(select_working_set_plus(set_type, i,j,k, iter) != 0)
	break;
      else
	counter = 1;	// do shrinking next iteration
    }
			
    ++iter;
	 
    switch(set_type) {

    case BETA_I_BETA_J: 
      Q_i_star_beta = Q_star_beta.get_Q(i,active_size_beta);
      Q_j_star_beta = Q_star_beta.get_Q(j,active_size_beta);
      beta_i_old = beta[i];
      beta_j_old = beta[j];
      Delta = beta_i_old + beta_j_old;
      beta[i] += (g_beta[j]-g_beta[i])/(Q_i_star_beta[i]+Q_j_star_beta[j]-2*Q_i_star_beta[j]);
      beta_i = beta[i];
      if (beta_i < 0)
	beta[i] = 0;
      if (beta_i > Delta)
	beta[i] = Delta;
      beta[j] = Delta - beta[i];
     
      diff_i = beta[i]-beta_i_old;
      diff_j = beta[j]-beta_j_old;
      for(r=0; r<active_size_beta; r++) 
	g_beta[r] += diff_i*Q_i_star_beta[r]+diff_j*Q_j_star_beta[r]; 

      true_i = active_set_beta[i];
      act_set_i = true_act_set[true_i];
      true_j = active_set_beta[j];
      act_set_j = true_act_set[true_j];
      Q_i_star = Q_star.get_Q(act_set_i,active_size);
      Q_j_star = Q_star.get_Q(act_set_j,active_size);

      for(r=0; r<active_size; r++) 
        g[r] += diff_i*Q_i_star[r]+diff_j*Q_j_star[r];

      update_beta_status(i);
      update_beta_status(j);
      //      fprintf(stdout,"beta_i_old=%f beta_i=%f beta_j_old=%f beta_j=%f\n",beta_i_old,beta[i],beta_j_old,beta[j]);
      // fflush(stdout);
      break;

    case ALPHA_I_ALPHA_J:
      Q_i = Q.get_Q(i,active_size);
      Q_j = Q.get_Q(j,active_size);
      Q_i_star = Q_star.get_Q(i,active_size);
      Q_j_star = Q_star.get_Q(j,active_size);
      alpha_i_old = alpha[i];
      alpha_j_old = alpha[j];
      y_i = y[i];
      y_j = y[j];
      Delta = alpha_i_old + alpha_j_old;
      nominator = y_j*G[j]-y_i*G[i]+(g[j]-g[i])/tau;
      denominator = Q_i[i]+Q_j[j]-2*Q_i[j]+(Q_i_star[i]+Q_j_star[j]-2*Q_i_star[j])/tau;
      alpha[i] += nominator/denominator;
      alpha_i = alpha[i];
      if (alpha_i < 0)
	alpha[i] = 0;
      if (alpha_i > Delta)
	alpha[i] = Delta;
      alpha[j] = Delta - alpha[i];

      diff_i = alpha[i]-alpha_i_old;
      diff_j = alpha[j]-alpha_j_old;
      diff_i_y = diff_i * y_i;
      diff_j_y = diff_j * y_j;      
      for (r=0; r<active_size; r++) {
	G[r] += diff_i_y*Q_i[r]+diff_j_y*Q_j[r];
	g[r] += diff_i*Q_i_star[r]+diff_j*Q_j_star[r]; 
      }

      true_i = active_set[i];
      act_set_i = true_act_set_beta[true_i];
      true_j = active_set[j];
      act_set_j = true_act_set_beta[true_j];
      Q_i_star_beta = Q_star_beta.get_Q(act_set_i,active_size_beta);
      Q_j_star_beta = Q_star_beta.get_Q(act_set_j,active_size_beta);

      for (r=0; r<active_size_beta; r++) 
	g_beta[r] += diff_i*Q_i_star_beta[r]+diff_j*Q_j_star_beta[r]; 
      
      update_alpha_status(i);
      update_alpha_status(j);
      break;

    case ALPHA_I_ALPHA_J_BETA_K:
      Q_i = Q.get_Q(i,active_size);
      Q_j = Q.get_Q(j,active_size);
      Q_i_star = Q_star.get_Q(i,active_size);
      Q_j_star = Q_star.get_Q(j,active_size);
      Q_k_star_beta = Q_star_beta.get_Q(k,active_size_beta);
 
      true_k = active_set_beta[k];
      act_set_k = true_act_set[true_k];
      Q_k_star = Q_star.get_Q(act_set_k, active_size);

      alpha_i_old = alpha[i];
      alpha_j_old = alpha[j];
      beta_k_old = beta[k];
      y_i = y[i];
      y_j = y[j];
      if(alpha_i_old < alpha_j_old)
	min_alpha = alpha_i_old;
      else
	min_alpha = alpha_j_old;
      Delta = beta_k_old + 2*min_alpha;
      nominator = y[i]*G[i]+y[j]*G[j]-2+(g[i]+g[j]-2*g_beta[k])/tau;
      denominator = Q_i[i]+Q_j[j]-2*Q_i[j]+(Q_i_star[i]+Q_j_star[j]+2*Q_i_star[j]-4*Q_k_star[i]-4*Q_k_star[j]+4*Q_k_star_beta[k])/tau;
      beta[k] += 2*nominator/denominator;
      beta_k = beta[k];
      if (beta_k < 0)
	beta[k] = 0;
      if (beta_k > Delta)
	beta[k] = Delta;
      alpha_change = (beta_k_old-beta[k])/2;
      alpha[i] += alpha_change;
      alpha[j] += alpha_change;

      diff_i = alpha[i]-alpha_i_old;
      diff_j = alpha[j]-alpha_j_old;
      diff_k = beta[k]-beta_k_old;
      diff_i_y = diff_i * y_i;
      diff_j_y = diff_j * y_j;
    
      for (r=0; r<active_size; r++) {
	G[r] += diff_i_y*Q_i[r]+diff_j_y*Q_j[r];
	g[r] += diff_i*Q_i_star[r]+diff_j*Q_j_star[r]+diff_k*Q_k_star[r]; 
      }

      true_i = active_set[i];
      act_set_i = true_act_set_beta[true_i];
      true_j = active_set[j];
      act_set_j = true_act_set_beta[true_j];
      Q_i_star_beta = Q_star_beta.get_Q(act_set_i,active_size_beta);
      Q_j_star_beta = Q_star_beta.get_Q(act_set_j,active_size_beta);

      for(r=0; r<active_size_beta; r++)
        g_beta[r] += diff_i*Q_i_star_beta[r]+diff_j*Q_j_star_beta[r]+diff_k*Q_k_star_beta[r];

      update_alpha_status(i);
      update_alpha_status(j);
      update_beta_status(k);
	      
      break;
    }
  }

  calculate_rho_plus(si->rho,si->rho_star);

  t_fin = times( &fin_time);
  net_time = (double) (fin_time.tms_utime - init_time.tms_utime)/HZ;
  fprintf(stdout,"Solver time = %3.3e\n", net_time);
  fflush(stdout);

  // put back the solution
  for(i=0;i<l;i++) {
    alpha_[active_set[i]] = alpha[i];
    beta_[active_set_beta[i]] = beta[i];
  }
		
  si->upper_bound_p = Cp;
  si->upper_bound_n = Cn;

  info("\noptimization finished, #iter = %d\n",iter);

  si->rho *= -1;

  delete[] G;
  delete[] g;
  delete[] g_init;
  delete[] g_beta;
  delete[] g_beta_init;
  delete[] alpha_status;
  delete[] beta_status;
  delete[] alpha;
  delete[] beta;
}

void Solver_plus::calculate_rho_plus(double& rho, double& rho_star)
{
  int i, pos_size=0, neg_size=0;
  double pos_sum = 0, neg_sum = 0;

  for(i=0; i < active_size; i++)
    if(alpha[i]>1e-8) {
      if(y[i]==1) {
        pos_size++;
        pos_sum += 1-G[i]-g[i]/tau; 
      }
      else {
        neg_size++;
        neg_sum += -1-G[i]+g[i]/tau;
      }
    }

  if(pos_size != 0)
    pos_sum /= pos_size;

  if(neg_size != 0)
    neg_sum /= neg_size;

  rho = (pos_sum+neg_sum)/2;
  rho_star = pos_sum - rho; 
   
}

void Solver_plus::calculate_rho_plus_cg(double& rho, double& rho_star)
{
  int i, pos_size=0, neg_size=0;
  double pos_sum = 0, neg_sum = 0;
  double *alpha, *G, *g;

  alpha = alpha_cg[0];
  G = G_cg[0];
  g = g_cg[0];
  
  for(i=0; i < active_size; i++)
    if(alpha[i]>1e-8) {
      if(y[i]==1) {
        pos_size++;
        pos_sum += 1-G[i]-g[i]/tau; 
      }
      else {
        neg_size++;
        neg_sum += -1-G[i]+g[i]/tau;
      }
    }

  if(pos_size != 0)
    pos_sum /= pos_size;

  if(neg_size != 0)
    neg_sum /= neg_size;

  rho = (pos_sum+neg_sum)/2;
  rho_star = pos_sum - rho; 
   
}


// return 1 if already optimal, return 0 otherwise
int Solver_plus::select_working_set_plus(int& set_type, int& i_out, int& j_out, int& k_out, int iter)
{
  double gap[3];

  for(int i=0; i <3; i++)
    gap[i] = -1;

  int i, j, best_B1=-1, best_B1B2=-1, best_A2=-1, best_A2A4=-1, best_A1=-1, best_A1A3=-1, i_ind, j_ind, k_ind;
  int type_selected[3], selected_indices[3][3];
  double max_B1=-1e20, min_B1B2=1e20, max_A2=-1e20, min_A2A4=1e20, max_A1=-1e20, min_A1A3=1e20;
  double alpha_i, g_i, g_j, y_i, y_j, deriv_alpha_i, first_order_criterion;
  double max_z[3], z, absolute_best_z, nominator, nominator_base, denominator_base, j_deriv, tau2, Q_star_ii;
  int best_z_index[3], true_k, act_set_k;
  Qfloat *Q_i, *Q_i_star, *Q_k_star, *Q_k_star_beta;
  double deriv_alpha_i1, deriv_alpha_i2, nominator1, nominator2, nominator_base1, nominator_base2, j_deriv1, j_deriv2;  
  double max_A2_1, max_A2_2, min_A2A4_1, min_A2A4_2, max_A1_1, max_A1_2, min_A1A3_1, min_A1A3_2;

  // first-order working set selection 

  // compute all maxima and minima related to alphas
  for(i=0; i<active_size; i++) {
    alpha_i = alpha[i];
    g_i = g[i];
    y_i = y[i];
    deriv_alpha_i1 = y_i*G[i];
    deriv_alpha_i2 = g_i;
    deriv_alpha_i = deriv_alpha_i1+deriv_alpha_i2/tau;

    // max A2
    if(alpha_i>1e-8 && y_i==-1 && deriv_alpha_i>max_A2) {
      max_A2 = deriv_alpha_i;
      best_A2 = i;
      max_A2_1 = deriv_alpha_i1;
      max_A2_2 = deriv_alpha_i2;      
    }

    // min A2A4
    if(y_i==-1 && deriv_alpha_i<min_A2A4) {
      min_A2A4 = deriv_alpha_i;
      best_A2A4 = i;
      min_A2A4_1 = deriv_alpha_i1;
      min_A2A4_2 = deriv_alpha_i2;      
    }
    
    // max A1
    if(alpha_i>1e-8 && y_i==1 && deriv_alpha_i>max_A1) {
      max_A1 = deriv_alpha_i;
      best_A1 = i;
      max_A1_1 = deriv_alpha_i1;
      max_A1_2 = deriv_alpha_i2;      
    }
    
    // min A1A3
    if(y_i==1 && deriv_alpha_i<min_A1A3) {
      min_A1A3 = deriv_alpha_i;
      best_A1A3 = i;
      min_A1A3_1 = deriv_alpha_i1;
      min_A1A3_2 = deriv_alpha_i2;      
    } 
  }

  // compute all maxima and minima related to betas
  for(i=0; i<active_size_beta; i++) {
    g_i = g_beta[i];

    // max B1
    if(beta[i]>1e-8 && g_i>max_B1) {
      max_B1 = g_i;
      best_B1 = i;
    }

    // min B1B2
    if(g_i<min_B1B2) {
      min_B1B2 = g_i;
      best_B1B2 = i;
    }
  }
 
  max_B1 /= tau;
  min_B1B2 /= tau;

  // select maximal violating pairs
  if(max_B1-min_B1B2 < eps)
    type_selected[0] = 0;
  else {
    type_selected[0] = 1;
    selected_indices[0][0] = best_B1;
    selected_indices[0][1] = best_B1B2;
    gap[0] = max_B1-min_B1B2;
  }

 if(((max_A2 - min_A2A4 < eps) || ((max_A2_1 - min_A2A4_1 < eps) && (max_A2_2 - min_A2A4_2 < eps))) &&
     ((max_A1 - min_A1A3 < eps) || ((max_A1_1 - min_A1A3_1 < eps) && (max_A1_2 - min_A1A3_2 < eps))))
    type_selected[1] = 0;
  else {
    if((max_A2 - min_A2A4 > max_A1 - min_A1A3) && ((max_A2_1 - min_A2A4_1 >= eps) || (max_A2_2 - min_A2A4_2 >= eps))) {
      type_selected[1] = 1;
      selected_indices[1][0] = best_A2;
      selected_indices[1][1] = best_A2A4;
    }
    else {
      if((max_A2 - min_A2A4 <= max_A1 - min_A1A3) && ((max_A1_1 - min_A1A3_1 >= eps) || (max_A1_2 - min_A1A3_2 >= eps))) {
        type_selected[1] = 1;
        selected_indices[1][0] = best_A1;
        selected_indices[1][1] = best_A1A3;
      }
      else
        type_selected[1] = 0;
    }
  }

   if(((2*max_B1+2-min_A1A3-min_A2A4 < eps) || ((2-min_A1A3_1-min_A2A4_1<eps) && (2*max_B1*tau-min_A1A3_2-min_A2A4_2<eps))) &&
     ((max_A1+max_A2-2*min_B1B2-2 < eps) || ((max_A1_1+max_A2_1-2<eps) && (max_A1_2+max_A2_2-2*min_B1B2*tau<eps))))
    type_selected[2] = 0;
  else {
    if((2*max_B1+2-min_A1A3-min_A2A4 > max_A1+max_A2-2*min_B1B2-2) && ((2-min_A1A3_1-min_A2A4_1>=eps) || (2*max_B1*tau-min_A1A3_2-min_A2A4_2>=eps))) {
      type_selected[2] = 1;
      selected_indices[2][0] = best_A1A3;
      selected_indices[2][1] = best_A2A4;
      selected_indices[2][2] = best_B1;
    }
    else {
      if((2*max_B1+2-min_A1A3-min_A2A4 <= max_A1+max_A2-2*min_B1B2-2) && ((max_A1_1+max_A2_1-2>=eps) || (max_A1_2+max_A2_2-2*min_B1B2*tau>=eps))) {
        type_selected[2] = 1;
        selected_indices[2][0] = best_A1;
        selected_indices[2][1] = best_A2;
        selected_indices[2][2] = best_B1B2;
      }
      else
        type_selected[2] = 0;
    }
  }

  if(type_selected[0]+type_selected[1]+type_selected[2] == 0) 
    return 1;

  for(i=0; i<3; i++)
    max_z[i] = -1e20;    

  // second-order working set selection
  if(type_selected[0] == 1) {
     i_ind = selected_indices[0][0];
     g_i = g_beta[i_ind];
     Q_i_star = Q_star_beta->get_Q(i_ind,active_size_beta);
     Q_star_ii = Q_i_star[i_ind];
     tau2 = 2*tau;
     for(j = 0; j < active_size_beta; j++) {
        g_j = g_beta[j];
	if(eps+g_j/tau < g_i/tau) {
	  nominator = g_i-g_j;
	  z = nominator*nominator/(tau2*(Q_star_ii+QD_star_beta[j]-2*Q_i_star[j]));
	  if(z > max_z[0]) {
	    max_z[0] = z;
	    best_z_index[0] = j;
	  }
	}
      }
  }

  if(type_selected[1] == 1) {
     i_ind = selected_indices[1][0];
     y_i = y[i_ind];
     Q_i = Q->get_Q(i_ind,active_size);
     Q_i_star = Q_star->get_Q(i_ind,active_size);
     nominator_base = y_i*G[i_ind] + g[i_ind]/tau;
      nominator_base1 = y_i*G[i_ind];
     nominator_base2 = g[i_ind];
     denominator_base = 2*(Q_i[i_ind]+Q_i_star[i_ind]/tau);
   
     for(j = 0; j < active_size; j++) { 
        y_j = y[j];
        j_deriv = y_j*G[j]+g[j]/tau;
         j_deriv1 = y_j*G[j];
        j_deriv2 = g[j];

        if(y_j == y_i && j_deriv+eps < nominator_base && ((j_deriv1+eps < nominator_base1) || (j_deriv2+eps < nominator_base2))) {
	  nominator = nominator_base-j_deriv;
	  z = nominator*nominator/(denominator_base + 2*(QD[j]-2*Q_i[j]+(QD_star[j]-2*Q_i_star[j])/tau));
	  if(z > max_z[1]) {
	    max_z[1] = z;
	    best_z_index[1] = j;
	  }
	}
    }
  }

  if(type_selected[2] == 1) {
    i_ind = selected_indices[2][0];
    j_ind = selected_indices[2][1];
    k_ind = selected_indices[2][2];
    Q_i = Q->get_Q(i_ind,active_size);
    Q_i_star = Q_star->get_Q(i_ind,active_size);
    Q_k_star_beta = Q_star_beta->get_Q(k_ind,active_size_beta);

    true_k = active_set_beta[k_ind];
    act_set_k = true_act_set[true_k];
    Q_k_star = Q_star->get_Q(act_set_k, active_size);

    nominator_base = y[i_ind]*G[i_ind]-2+(g[i_ind]-2*g_beta[k_ind])/tau;
      nominator_base1 = y[i_ind]*G[i_ind]-2;
    nominator_base2 = g[i_ind]-2*g_beta[k_ind];
    denominator_base = 2*(Q_i[i_ind]+(Q_i_star[i_ind]-4*Q_k_star[i_ind]+4*Q_k_star_beta[k_ind])/tau);
    first_order_criterion = nominator_base+y[j_ind]*G[j_ind]+g[j_ind]/tau;
    for(j = 0; j < active_size; j++) 
      if(y[j] == -1) {
         nominator1 = nominator_base1+y[j]*G[j];
        nominator2 = nominator_base2+g[j];
        nominator = nominator_base+y[j]*G[j]+g[j]/tau;
        if((first_order_criterion < 0 && nominator<-eps && ((nominator1 < -eps) || (nominator2 < -eps))) ||
           (first_order_criterion > 0 && alpha[j] > 1e-8 && nominator > eps && ((nominator1 > eps) || (nominator2 > eps)))) {
	  z = nominator*nominator/(denominator_base + 2*(QD[j]-2*Q_i[j]+(QD_star[j]+2*Q_i_star[j]-4*Q_k_star[j])/tau));
	  if(z > max_z[2]) {
	    max_z[2] = z;
	    best_z_index[2] = j;
	  }
        }
      }
  }

  // choose the best type
  absolute_best_z = -1;
  for(i=0; i<3; i++)
    if((type_selected[i] == 1) && (max_z[i] > absolute_best_z)) {  
      absolute_best_z = max_z[i];
      set_type = i;
      i_out = selected_indices[i][0];
      j_out = best_z_index[i];
      if(i == 2)
	k_out = selected_indices[i][2];
    }

  return 0;
}

bool Solver_plus::be_shrunk_alpha(int i, double max_B1, double max_A1, double max_A2, double min_B1B2, double min_A1A3, double min_A2A4)
{
  int y_i = y[i];
  double deriv_i = y_i*G[i]+g[i]/tau;

  if(alpha[i]<=1e-8) {
    if(y_i == 1 && deriv_i <= max_A1+eps)
      return false;
    if(y_i == -1 && deriv_i <= max_A2+eps)
      return false;
    return deriv_i+min_A1A3+eps>2*max_B1+2;
  }
  else {
    if(y_i == 1) 
      return max_A1-deriv_i<eps && deriv_i-min_A1A3<eps && 2*max_B1+2-deriv_i-min_A2A4<eps && deriv_i+max_A2-2*min_B1B2-2<eps; 
    else 
      return max_A2-deriv_i<eps && deriv_i-min_A2A4<eps && 2*max_B1+2-min_A1A3-deriv_i<eps && max_A1+deriv_i-2*min_B1B2-2<eps;
  }

}

bool Solver_plus::be_shrunk_beta(int i, double max_B1, double max_A1, double max_A2, double min_B1B2, double min_A1A3, double min_A2A4)
{
  double g_beta_i = g_beta[i]/tau;

  if(beta[i]<=1e-8) 
    return (g_beta_i +eps> max_B1 && 2*g_beta_i+2 + eps > max_A1 + max_A2);
  else 
    return (g_beta_i-min_B1B2<eps && max_B1-g_beta_i<eps && 
            2*g_beta_i+2-min_A1A3-min_A2A4<eps && max_A1+max_A2-2*g_beta_i-2<eps);

}

bool Solver_plus::be_shrunk_alpha_cg(int i, double max_B1, double max_A1, double max_A2, double min_B1B2, double min_A1A3, double min_A2A4)
{
  if(active[i])
    return false;

  int y_i = y[i];
  double deriv_i = y_i*G_cg[0][i]+g_cg[0][i]/tau;

  if(alpha_cg[0][i]<=1e-8) {
    if(y_i == 1 && deriv_i <= max_A1+eps)
      return false;
    if(y_i == -1 && deriv_i <= max_A2+eps)
      return false;
    return deriv_i+min_A1A3+eps>2*max_B1+2;
  }
  else {
    if(y_i == 1) 
      return max_A1-deriv_i<eps && deriv_i-min_A1A3<eps && 2*max_B1+2-deriv_i-min_A2A4<eps && deriv_i+max_A2-2*min_B1B2-2<eps; 
    else 
      return max_A2-deriv_i<eps && deriv_i-min_A2A4<eps && 2*max_B1+2-min_A1A3-deriv_i<eps && max_A1+deriv_i-2*min_B1B2-2<eps;
  }
}

bool Solver_plus::be_shrunk_beta_cg(int i, double max_B1, double max_A1, double max_A2, double min_B1B2, double min_A1A3, double min_A2A4)
{
  if(active[i+l])
    return false;

  double g_beta_i = g_beta_cg[0][i]/tau;

  if(beta_cg[0][i]<=1e-8) 
    return (g_beta_i +eps> max_B1 && 2*g_beta_i+2 + eps > max_A1 + max_A2);
  else 
    return (g_beta_i-min_B1B2<eps && max_B1-g_beta_i<eps && 
            2*g_beta_i+2-min_A1A3-min_A2A4<eps && max_A1+max_A2-2*g_beta_i-2<eps);
}

void Solver_plus::do_shrinking_plus() 
{
  int i, y_i;
  double g_i, alpha_i, deriv_alpha_i;
  double max_B1=-1e20, min_B1B2=1e20, max_A2=-1e20, min_A2A4=1e20, max_A1=-1e20, min_A1A3=1e20;

  // compute all maxima and minima related to alphas
  for(i=0; i<active_size; i++) {
    alpha_i = alpha[i];
    g_i = g[i];
    y_i = y[i];
    deriv_alpha_i = y_i*G[i]+g_i/tau;

    // max A2
    if(alpha_i>1e-8 && y_i==-1 && deriv_alpha_i>max_A2) 
      max_A2 = deriv_alpha_i;

    // min A2A4
    if(y_i==-1 && deriv_alpha_i<min_A2A4) 
      min_A2A4 = deriv_alpha_i;
    
    // max A1
    if(alpha_i>1e-8 && y_i==1 && deriv_alpha_i>max_A1) 
      max_A1 = deriv_alpha_i;
    
    // min A1A3max_A2, min_A2A4, max_A1, min_A1A3
    if(y_i==1 && deriv_alpha_i<min_A1A3) 
      min_A1A3 = deriv_alpha_i;
  }

  // compute all maxima and minima related to betas
  for(i=0; i<active_size_beta; i++) {
    g_i = g_beta[i];

    // max B1
    if(beta[i]>1e-8 && g_i>max_B1) 
      max_B1 = g_i;

    // min B1B2
    if(g_i<min_B1B2) 
      min_B1B2 = g_i;
  }

  if(unshrink == false && max_B1-min_B1B2 < eps*10 &&
     max_A2-min_A2A4 < eps*10 && max_A1 - min_A1A3 < eps*10 &&
     2*max_B1+2-min_A1A3-min_A2A4 < eps*10 && max_A1+max_A2-2*min_B1B2-2 < eps*10) {
    unshrink = true;
    reconstruct_gradient_plus();
    active_size = l;
    active_size_beta = l;
  }

  if(active_size_beta > 2) {
    for(i=0;i<active_size_beta;i++) {
      if(active_size_beta <= 2)
	break;
      if (be_shrunk_beta(i, max_B1, max_A1, max_A2, min_B1B2, min_A1A3, min_A2A4)) {
	active_size_beta--;
	while (active_size_beta > i) {
	  if (!be_shrunk_beta(active_size_beta, max_B1, max_A1, max_A2, min_B1B2, min_A1A3, min_A2A4)) {
	    swap_index_beta(i,active_size_beta);
	    break;
	  }
	  active_size_beta--;
	  if(active_size_beta <= 2)
	    break;
	}
      }
    }
  }

  for(i=0;i<active_size;i++) {
    if (be_shrunk_alpha(i, max_B1, max_A1, max_A2, min_B1B2, min_A1A3, min_A2A4)) {
      active_size--;
      while (active_size > i) {
	if (!be_shrunk_alpha(active_size, max_B1, max_A1, max_A2, min_B1B2, min_A1A3, min_A2A4)) {
	  swap_index_alpha(i,active_size);
	  break;
	}
	active_size--;
      }
    }
  }
}

