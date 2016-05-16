#include "common.h"
#include "solve_linear_system.h"
#include "Solver.h"
#include "Kernel.h"

bool Solver::bumped()
{
  for(int i=0; i<curr_depth+2; i++)
    if(is_upper_bound_cg(work[i],0) || is_lower_bound_cg(work[i],0))
      return true;
  return false;
}

void Solver::swap_index(int i, int j)
{
  Q->swap_index(i,j);
  swap(y[i],y[j]);
  swap(p[i],p[j]);
  swap(active_set[i],active_set[j]);

  if(!conjugate) {   
    swap(G[i],G[j]);
    swap(alpha_status[i],alpha_status[j]);
    swap(alpha[i],alpha[j]);
    swap(G_bar[i],G_bar[j]);
  }
  else {
    int k;
    for(k=0; k < curr_depth+1; k++) {
      swap(G_cg[k][i], G_cg[k][j]);
      swap(alpha_status_cg[k][i],alpha_status_cg[k][j]);
      swap(alpha_cg[k][i],alpha_cg[k][j]);
      swap(G_bar_cg[k][i],G_bar_cg[k][j]);
    }
    for(k=0; k<curr_depth+2; k++) {
      if(i == work[k]) 
        work[k] = -1;
      else 
        if(j == work[k])
	  work[k] = i;  
    }
    active[i] = active[j];
    active[j] = false;
  }
}

void Solver::reconstruct_gradient()
{
  // reconstruct inactive elements of G and G_old from G_bar, G_bar_old and free variables
  int i,j;
  int nr_free = 0;

  if(!conjugate) {
     
    if(active_size == l) 
      return;

    for(j=active_size;j<l;j++) 
      G[j] = G_bar[j] + p[j];
    for(j=0;j<active_size;j++) {
      if(is_free(j))
	nr_free++;
    }

    if(2*nr_free < active_size)
      info("\nWarning: using -h 0 may be faster\n");

    if (nr_free*l > 2*active_size*(l-active_size)) 
      for(i=active_size;i<l;i++) {
	const Qfloat *Q_i = Q->get_Q(i,active_size);
	for(j=0;j<active_size;j++) 
	  if(is_free(j))
	    G[i] += alpha[j] * Q_i[j];
      }
    else 
      for(i=0;i<active_size;i++)
	if(is_free(i)) {
	  const Qfloat *Q_i = Q->get_Q(i,l);
	  double alpha_i = alpha[i];
	
	  if(is_free(i))
	    for(j=active_size;j<l;j++) 
	      G[j] += alpha_i * Q_i[j];
	}
  }
  else {

    if(active_size == l) 
      return;

    //fprintf(stdout,"reconstructing gradient\n");
    //        fflush(stdout);

    for(int k=0; k<curr_depth+1; k++) {
      for(j=active_size;j<l;j++) 
	G_cg[k][j] = G_bar_cg[k][j] + p[j]; 

      nr_free = 0;
      for(j=0;j<active_size;j++) 
	if(is_free_cg(j,k))
	  nr_free++;

      if (nr_free*l > 2*active_size*(l-active_size)) 
	for(i=active_size;i<l;i++) {
	  const Qfloat *Q_i = Q->get_Q(i,active_size);
	  for(j=0;j<active_size;j++) 
	    if(is_free_cg(j,k))
	      G_cg[k][i] += alpha_cg[k][j] * Q_i[j];
	}
      
      else 
	for(i=0;i<active_size;i++)
	  if(is_free_cg(i,k)) {
	    const Qfloat *Q_i = Q->get_Q(i,l);
	    double alpha_i = alpha_cg[k][i];
	
	    if(is_free_cg(i,k))
	      for(j=active_size;j<l;j++) 
		G_cg[k][j] += alpha_i * Q_i[j];
	  }
    }
  }
}

void Solver::Solve_cg(int l, const QMatrix& Q, const double *p_, const schar *y_,
		      double *alpha_, double Cp, double Cn, double eps,
		      SolutionInfo* si, int shrinking)
{
  int i, j;
  double delta_alpha_i, C_i;

  this->l = l;
  this->Q = &Q;
  QD=Q.get_QD();
  clone(p, p_,l);
  clone(y, y_,l);
  this->Cp = Cp;
  this->Cn = Cn;
  this->eps = eps;
  unshrink = false;
  curr_depth = -1;
  active_size = l;

  work = new int[max_depth+2];
  for(i=0; i<max_depth+2; i++)
    work[i] = -1;

  active = new bool[l];
  for(i=0; i<l; i++)
    active[i] = false;

  // initialize alpha's 
  alpha_cg = new double*[max_depth+1];
  for(i = 0; i < max_depth + 1; i++) 
    clone(alpha_cg[i],alpha_,l);
  
  // initialize alpha_status
  alpha_status_cg = new char*[max_depth+1];
  for(i=0; i < max_depth + 1; i++) {
    alpha_status_cg[i] = new char[l];
    for(j=0;j<l;j++) {
      update_alpha_status_cg(j,i);
    }
  }

  // initialize active set (for shrinking)
  active_set = new int[l];
  for(i=0;i<l;i++)
    active_set[i] = i;
 
  // initialize gradient
  G_cg = new double*[max_depth+1];
  G_bar_cg = new double*[max_depth+1];
  G_cg[0] = new double[l];
  G_bar_cg[0] = new double[l]; 

  for(i=0; i<l; i++) {
    G_cg[0][i] = p[i];
    G_bar_cg[0][i] = 0;
  }

  for(i=0;i<l;i++)
    if(!is_lower_bound_cg(i,0)) {
      const Qfloat *Q_i = Q.get_Q(i,l);
      double alpha_i = alpha_cg[0][i];
      int j;
      for(j=0;j<l;j++)
	G_cg[0][j] += alpha_i*Q_i[j];
      if(is_upper_bound_cg(i,0)) {
        C_i = get_C(i);
	for(j=0;j<l;j++) 
	  G_bar_cg[0][j] += C_i * Q_i[j];
      }
    }

  for(i=1; i<max_depth+1; i++) {
    clone(G_cg[i], G_cg[0], l);
    clone(G_bar_cg[i], G_bar_cg[0], l);
  }
   
  // optimization step

  iter = 0;
  int counter = min(l,1000)+1;
  double *u, *tmp_u;
  double lambda_star, tmp_lambda, gain_hmg, gain2;
  bool corner=false, upper_bound_i, bumped=false;
  double alpha_i;
  Qfloat *Q_i;
  double *tmp_alpha, *tmp_G, *tmp_G_bar;
  char *tmp_status;
  int work_i, prev_depth;
  int *tmp_work, tmp_depth, tmp_working_set_size, working_set_size;
  double old_alpha[2];
  bool old_upper_bound[2], done_shrinking = false;
  int last_updated_G_bar=-1;
  int n_conjugate = 0, next_depth;

  u = new double[max_depth+2];
  tmp_u = new double[max_depth+2];
  tmp_work = new int[max_depth+2];

  while(1) {

    // do shrinking
    ///////////////////////////////////////
    if(--counter == 0) {
      counter = min(l,1000);
      if(shrinking) 
	done_shrinking = do_shrinking_cg();
    }
    ///////////////////////////////////////

    // do working set selection
    ////////////////////////////////////////////////////////////     
   
    prev_depth = curr_depth;

    if(corner || iter == 0) {
      if(wss_first_order(u, lambda_star)!=0) {
	reconstruct_gradient();
	active_size = l;
	if(wss_first_order(u, lambda_star)!=0) 
	  break;
	else
	  counter = 1;
      }
    }
    else {
      if(curr_depth >= 0 && max_depth > 0) {
	if(select_working_set_hmg(u, lambda_star, gain_hmg)!=0) {
	  memcpy(tmp_work, work, sizeof_int*(curr_depth+2));
	  tmp_depth = curr_depth;
          for(i=0; i<curr_depth+2; i++) {
            active[work[i]] = false;
            work[i] = -1;
	  }
          curr_depth = -1;
	  if(select_working_set_hmg2(u, lambda_star, gain2)!=0) {
	    memcpy(work, tmp_work, sizeof_int*(tmp_depth+2));
	    curr_depth = tmp_depth;
            for(i=0; i<curr_depth+2; i++) 
              active[work[i]] = true;
	    // reconstruct the whole gradient
	    reconstruct_gradient();
	    // reset active set size 
            active_size = l;
	    if(select_working_set_hmg(u, lambda_star, gain_hmg)==0)
	      counter = 1;
	    else {
	      if(select_working_set_hmg2(u, lambda_star, gain2)==0)
		counter = 1;	
	      else 
		break;
	    }
	  }
	}
        else {
          // try to improve the working set
         
	  // store current direction, step size and the working set
          working_set_size = curr_depth + 2;
	  memcpy(tmp_work, work, sizeof_int*working_set_size);
	  memcpy(tmp_u, u, sizeof_double*working_set_size);
	  tmp_depth = curr_depth;
	  tmp_working_set_size = working_set_size;
	  tmp_lambda = lambda_star;
	  for(i=0; i<working_set_size; i++) {
	    active[work[i]] = false;
	    work[i] = -1;
	  }
	  curr_depth = 0;
	  working_set_size = 0;
	  if(select_working_set_hmg2(u, lambda_star, gain2) == 0) {
	    if(gain_hmg>gain2) {
	      for(i=0; i<working_set_size; i++)
		active[work[i]] = false;
	      curr_depth = tmp_depth;
	      working_set_size = tmp_working_set_size;
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
	    lambda_star = tmp_lambda;
	    memcpy(u, tmp_u, sizeof_double*working_set_size);
	    memcpy(work, tmp_work, sizeof_int*tmp_working_set_size);
	    for(i=0; i<working_set_size; i++)
	      active[work[i]] = true;
	  }
          
	}
      }
      else 
	if(select_working_set_hmg2(u, lambda_star, gain2)!=0) {
	  reconstruct_gradient();
	  active_size = l;
          if(select_working_set_hmg2(u, lambda_star, gain2)==0)
	    counter = 1;
          else
	    break;  
	}
    }
    iter++;

    /////////////////////////////////////////////////////////////

    // work array - the chosen working set     
    // working_set_size - size of the working set

    // shift old alpha's, G's and G_bar's
    if(max_depth > 0) {
      n_conjugate++;
      next_depth = min(curr_depth+1,max_depth);
      tmp_alpha = alpha_cg[next_depth];
      tmp_G = G_cg[next_depth];
      tmp_G_bar = G_bar_cg[next_depth];
      tmp_status = alpha_status_cg[next_depth];

      for(i=next_depth; i>0; i--) {
	alpha_cg[i] = alpha_cg[i-1];
	G_cg[i] = G_cg[i-1];
	G_bar_cg[i] = G_bar_cg[i-1];
	alpha_status_cg[i] = alpha_status_cg[i-1];
      }
		
      if(!done_shrinking || (curr_depth == max_depth && curr_depth == prev_depth)) {
	memcpy(tmp_alpha, alpha_cg[0], active_size*sizeof_double);
	memcpy(tmp_G, G_cg[0], active_size*sizeof_double);
	if((iter-last_updated_G_bar <= max_depth) || (curr_depth > prev_depth))
	  memcpy(tmp_G_bar, G_bar_cg[0], l*sizeof_double);
	memcpy(tmp_status, alpha_status_cg[0], active_size*sizeof_char);
	done_shrinking = false;
      }
      else {
	memcpy(tmp_alpha, alpha_cg[0], l*sizeof_double);
	memcpy(tmp_G, G_cg[0], l*sizeof_double);
	memcpy(tmp_G_bar, G_bar_cg[0], l*sizeof_double);
	memcpy(tmp_status, alpha_status_cg[0], l*sizeof_char);
      }

      alpha_cg[0] = tmp_alpha;
      G_cg[0] = tmp_G;
      G_bar_cg[0] = tmp_G_bar;
      alpha_status_cg[0] = tmp_status;
    }
    else {
      old_alpha[0] = alpha_cg[0][work[0]];
      old_alpha[1] = alpha_cg[0][work[1]];
      old_upper_bound[0] = is_upper_bound_cg(work[0],0);
      old_upper_bound[1] = is_upper_bound_cg(work[1],0);
    }

    // update alpha's and alpha_status
    for(i=0; i<curr_depth+2; i++) {
      work_i = work[i];
      alpha_cg[0][work_i] += lambda_star*u[i];
      update_alpha_status_cg(work_i,0);
    }

    // check if we are at the corner
    corner = true;
    bumped = false;
    for(i=0; i<curr_depth+2; i++) {
      work_i = work[i];
      alpha_i = alpha_cg[0][work_i];
      C_i = get_C(work_i);
      if(alpha_i >= 1e-8*C_i && alpha_i <= C_i-1e-8) {
	corner = false;
	break;
      }
      else 
        bumped = true;
    }
     
    // update G and G_bar
    if(max_depth > 0)
      for(i=0; i<curr_depth+2; i++) {
	work_i = work[i];
	delta_alpha_i = alpha_cg[0][work_i] - alpha_cg[1][work_i];
	Q_i = Q.get_Q(work_i,active_size);
	for(j=0; j<active_size; j++)
	  G_cg[0][j] += Q_i[j]*delta_alpha_i;
	upper_bound_i = is_upper_bound_cg(work_i,0);
	if(upper_bound_i !=  is_upper_bound_cg(work_i,1)) {
          last_updated_G_bar = iter;
	  C_i = get_C(work_i);
          Q_i = Q.get_Q(work_i,l);
	  if(upper_bound_i)
	    for(j=0; j<l; j++)
	      G_bar_cg[0][j] += Q_i[j]*C_i; 
	  else
	    for(j=0; j<l; j++)
	      G_bar_cg[0][j] -= Q_i[j]*C_i; 
	}
      }
    else
      for(i=0; i<curr_depth+2; i++) {
	work_i = work[i];
	delta_alpha_i = alpha_cg[0][work_i] - old_alpha[i];
	Q_i = Q.get_Q(work_i,active_size);
	for(j=0; j<active_size; j++)
	  G_cg[0][j] += Q_i[j]*delta_alpha_i;
	upper_bound_i = is_upper_bound_cg(work_i,0);
	if(upper_bound_i !=  old_upper_bound[i]) {
          last_updated_G_bar = iter;
	  C_i = get_C(work_i);
          Q_i = Q.get_Q(work_i,l);
	  if(upper_bound_i)
	    for(j=0; j<l; j++)
	      G_bar_cg[0][j] += Q_i[j]*C_i; 
	  else
	    for(j=0; j<l; j++)
	      G_bar_cg[0][j] -= Q_i[j]*C_i; 
	}
      }
  }

  // calculate rho
  si->rho = calculate_rho_cg();

  // calculate objective value
  double v = 0;
  for(i=0;i<l;i++)
    v += alpha_cg[0][i] * (G_cg[0][i] + p[i]);

  si->obj = v/2;

  // put back the solution
  for(i=0;i<l;i++)
    alpha_[active_set[i]] = alpha_cg[0][i];

  si->upper_bound_p = Cp;
  si->upper_bound_n = Cn;

  info("\noptimization finished, #iter = %d\n",iter);

  delete[] p;
  delete[] y;
  delete[] u;
  delete[] tmp_u;
  delete[] work;
  delete[] tmp_work;
  delete[] active_set;
  delete[] active;
  /*
  for(i=0; i<max_depth+1; i++) {
    delete[] alpha_cg[i];
    delete[] alpha_status_cg[i]; 
    delete[] G_cg[i];
    delete[] G_bar_cg[i];
  }
  delete[] alpha_cg;
  delete[] alpha_status_cg;
  delete[] G_cg;
  delete[] G_bar_cg;
  */
}


void Solver::Solve(int l, const QMatrix& Q, const double *p_, const schar *y_,
		   double *alpha_, double Cp, double Cn, double eps,
		   SolutionInfo* si, int shrinking)
{
  this->l = l;
  this->Q = &Q;
  QD=Q.get_QD();
  clone(p, p_,l);
  clone(y, y_,l);
  clone(alpha,alpha_,l);
  this->Cp = Cp;
  this->Cn = Cn;
  this->eps = eps;
  unshrink = false;

  // initialize alpha_status
  {
    alpha_status = new char[l];
    for(int i=0;i<l;i++)
      update_alpha_status(i);
  }

  // initialize active set (for shrinking)
  {
    active_set = new int[l];
    for(int i=0;i<l;i++)
      active_set[i] = i;
    active_size = l;
  }

  // initialize gradient
  {
    G = new double[l];
    G_bar = new double[l];
    int i;
    for(i=0;i<l;i++)
      {
	G[i] = p[i];
	G_bar[i] = 0;
      }
    for(i=0;i<l;i++)
      if(!is_lower_bound(i))
	{
	  const Qfloat *Q_i = Q.get_Q(i,l);
	  double alpha_i = alpha[i];
	  int j;
	  for(j=0;j<l;j++)
	    G[j] += alpha_i*Q_i[j];
	  if(is_upper_bound(i))
	    for(j=0;j<l;j++)
	      G_bar[j] += get_C(i) * Q_i[j];
	}
  }

  // optimization step

  iter = 0;
  int counter = min(l,1000)+1;
  int nbumped = 0;

  while(1)
    {
      // show progress and do shrinking

      if(--counter == 0)
	{
	  counter = min(l,1000);
	  if(shrinking) do_shrinking();
	  info(".");
	}

      int i,j;
      if(select_working_set(i,j)!=0)
	{
	  // reconstruct the whole gradient
	  reconstruct_gradient();
	  // reset active set size and check
	  active_size = l;
	  info("*");
	  if(select_working_set(i,j)!=0)
	    break;
	  else
	    counter = 1;	// do shrinking next iteration
	}
		
      ++iter;

      // update alpha[i] and alpha[j], handle bounds carefully
		
      const Qfloat *Q_i = Q.get_Q(i,active_size);
      const Qfloat *Q_j = Q.get_Q(j,active_size);

      double C_i = get_C(i);
      double C_j = get_C(j);

      double old_alpha_i = alpha[i];
      double old_alpha_j = alpha[j];

      if(y[i]!=y[j])
	{
	  double quad_coef = Q_i[i]+Q_j[j]+2*Q_i[j];
	  if (quad_coef <= 0)
	    quad_coef = TAU;
	  double delta = (-G[i]-G[j])/quad_coef;
	  double diff = alpha[i] - alpha[j];
	  alpha[i] += delta;
	  alpha[j] += delta;
			
	  if(diff > 0)
	    {
	      if(alpha[j] < 0)
		{
		  alpha[j] = 0;
		  alpha[i] = diff;
		}
	    }
	  else
	    {
	      if(alpha[i] < 0)
		{
		  alpha[i] = 0;
		  alpha[j] = -diff;
		}
	    }
	  if(diff > C_i - C_j)
	    {
	      if(alpha[i] > C_i)
		{
		  alpha[i] = C_i;
		  alpha[j] = C_i - diff;
		}
	    }
	  else
	    {
	      if(alpha[j] > C_j)
		{
		  alpha[j] = C_j;
		  alpha[i] = C_j + diff;
		}
	    }
	}
      else
	{
	  double quad_coef = Q_i[i]+Q_j[j]-2*Q_i[j];
	  if (quad_coef <= 0)
	    quad_coef = TAU;
	  double delta = (G[i]-G[j])/quad_coef;
	  double sum = alpha[i] + alpha[j];
	  alpha[i] -= delta;
	  alpha[j] += delta;

	  if(sum > C_i)
	    {
	      if(alpha[i] > C_i)
		{
		  alpha[i] = C_i;
		  alpha[j] = sum - C_i;
		}
	    }
	  else
	    {
	      if(alpha[j] < 0)
		{
		  alpha[j] = 0;
		  alpha[i] = sum;
		}
	    }
	  if(sum > C_j)
	    {
	      if(alpha[j] > C_j)
		{
		  alpha[j] = C_j;
		  alpha[i] = sum - C_j;
		}
	    }
	  else
	    {
	      if(alpha[i] < 0)
		{
		  alpha[i] = 0;
		  alpha[j] = sum;
		}
	    }
	}

      // update G

      double delta_alpha_i = alpha[i] - old_alpha_i;
      double delta_alpha_j = alpha[j] - old_alpha_j;
		
      for(int k=0;k<active_size;k++)
	{
	  G[k] += Q_i[k]*delta_alpha_i + Q_j[k]*delta_alpha_j;
	}

      // update alpha_status and G_bar

      {
	bool ui = is_upper_bound(i);
	bool uj = is_upper_bound(j);
	update_alpha_status(i);
	update_alpha_status(j);
	int k;
	if(ui != is_upper_bound(i))
	  {
	    Q_i = Q.get_Q(i,l);
	    if(ui)
	      for(k=0;k<l;k++)
		G_bar[k] -= C_i * Q_i[k];
	    else
	      for(k=0;k<l;k++)
		G_bar[k] += C_i * Q_i[k];
	  }

	if(uj != is_upper_bound(j))
	  {
	    Q_j = Q.get_Q(j,l);
	    if(uj)
	      for(k=0;k<l;k++)
		G_bar[k] -= C_j * Q_j[k];
	    else
	      for(k=0;k<l;k++)
		G_bar[k] += C_j * Q_j[k];
	  }
      }
     
      double bumped_i = is_lower_bound(i) || is_upper_bound(i);
      double bumped_j = is_lower_bound(j) || is_upper_bound(j);

      if((bumped_i && !bumped_j) || (!bumped_i && bumped_j))
	nbumped++;
    }

  // calculate rho

  si->rho = calculate_rho();

  // calculate objective value
  {
    double v = 0;
    int i;
    for(i=0;i<l;i++)
      v += alpha[i] * (G[i] + p[i]);

    si->obj = v/2;
  }

  // put back the solution
  {
    for(int i=0;i<l;i++)
      alpha_[active_set[i]] = alpha[i];
  }

  // juggle everything back
  /*{
    for(int i=0;i<l;i++)
    while(active_set[i] != i)
    swap_index(i,active_set[i]);
    // or Q.swap_index(i,active_set[i]);
    }*/

  si->upper_bound_p = Cp;
  si->upper_bound_n = Cn;

  info("\noptimization finished, #iter = %d\n",iter);

  delete[] p;
  delete[] y;
  delete[] alpha;
  delete[] alpha_status;
  delete[] active_set;
  delete[] G;
  delete[] G_bar;
}

int Solver::generate_direction(double *u, int n, bool new_working_set)
{
  //   return generate_direction_general(u, n);

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

int Solver::generate_direction_general(double *u, int n, bool new_working_set)
{
  int i,j,result;
  int n_minus_one = n-1;
  int n_minus_two = n-2;
  int start, workj, work0;
  double *A0, *Ai, *G, *G_old;  

  if(new_working_set)
    start = 0;
  else
    start = n_minus_two;

  // compute A matrix
  A0 = A[0];
  for(i=start; i<n_minus_one; i++)
    A0[i] = y[work[i+1]];
  for(i=1; i<=n_minus_two; i++) {
    Ai = A[i];
    G_old = G_cg[i-1];
    G = G_cg[i]; 
    for(j=start; j<n_minus_one; j++) {
      workj = work[j+1];
      Ai[j] = G_old[workj]-G[workj];
    }
  }

  // compute b vector
  if(new_working_set) {
    work0 = work[0];
    b[0] = -y[work0];
    for(i=1; i<=n_minus_two; i++)
      b[i] = -G_cg[i-1][work0]+G_cg[i][work0];
  }

  result = solve_linear_system(A,b,n_minus_one);
  if(result)
    return 1;
  u[0] = 1;
  memcpy(&u[1], b, sizeof(double)*(n_minus_one));
 
  return 0;   
}

void Solver::generate_direction3(double *u, bool new_working_set)
{
  static double a_00, a_01, a_10, a_11;
  double a_02, a_12;
  double *G0, *G1;

  G0 = G_cg[0];
  G1 = G_cg[1];
  if(new_working_set) {
    int working0 = work[0];
    int working1 = work[1];
    a_00 = y[working0];
    a_01 = y[working1];
    a_10 = G0[working0]-G1[working0];
    a_11 = G0[working1]-G1[working1];
  }

  int working2 = work[2]; 
  a_02 = y[working2];
  a_12 = G0[working2]-G1[working2];

  u[0] = a_01*a_12-a_02*a_11;
  u[1] = a_02*a_10-a_00*a_12;
  u[2] = a_00*a_11-a_01*a_10;
}

void Solver::generate_direction4(double *u, bool new_working_set)
{
  static double a_00, a_01, a_02, a_10, a_11, a_12, a_20, a_21, a_22, m_01, m_02, m_12;
  double *G0, *G1, *G2;

  G0 = G_cg[0];
  G1 = G_cg[1];
  G2 = G_cg[2];

  int working0 = work[0];
  int working1 = work[1];
  int working2 = work[2];
  int working3 = work[3];

  double G_10 = G1[working0]; 
  double G_11 = G1[working1];
  double G_12 = G1[working2];
  double G_13 = G1[working3];

  if(new_working_set) {
    a_00 = y[working0];
    a_01 = y[working1];
    a_02 = y[working2];
    a_10 = G0[working0]-G_10;
    a_11 = G0[working1]-G_11;
    a_12 = G0[working2]-G_12;
    a_20 = G_10-G2[working0];
    a_21 = G_11-G2[working1];
    a_22 = G_12-G2[working2];

    m_01 = a_10*a_21-a_11*a_20;
    m_02 = a_10*a_22-a_12*a_20;
    m_12 = a_11*a_22-a_12*a_21;
  }
  
  double a_03= y[working3];
  double a_13= G0[working3]-G_13;
  double a_23= G_13-G2[working3];
  
  double m_03 = a_10*a_23-a_13*a_20;
  double m_13 = a_11*a_23-a_13*a_21;
  double m_23 = a_12*a_23-a_13*a_22;

  u[0] = a_01*m_23-a_02*m_13+a_03*m_12;
  u[1] = -a_00*m_23+a_02*m_03-a_03*m_02;
  u[2] = a_00*m_13-a_01*m_03+a_03*m_01;
  u[3] = -a_00*m_12+a_01*m_02-a_02*m_01;
}

// return 1 if already optimal, return 0 otherwise
int Solver::select_working_set(int &out_i, int &out_j)
{
  // return i,j such that
  // i: maximizes -y_i * grad(f)_i, i in I_up(\alpha)
  // j: minimizes the decrease of obj value
  //    (if quadratic coefficeint <= 0, replace it with tau)
  //    -y_j*grad(f)_j < -y_i*grad(f)_i, j in I_low(\alpha)
	
  double Gmax = -INF;
  double Gmax2 = -INF;
  int Gmax_idx = -1;
  int Gmin_idx = -1;
  double obj_diff_min = INF;

  for(int t=0;t<active_size;t++)
    if(y[t]==+1)	
      {
	if(!is_upper_bound(t))
	  if(-G[t] >= Gmax)
	    {
	      Gmax = -G[t];
	      Gmax_idx = t;
	    }
      }
    else
      {
	if(!is_lower_bound(t))
	  if(G[t] >= Gmax)
	    {
	      Gmax = G[t];
	      Gmax_idx = t;
	    }
      }

  int i = Gmax_idx;
  const Qfloat *Q_i = NULL;
  if(i != -1) // NULL Q_i not accessed: Gmax=-INF if i=-1
    Q_i = Q->get_Q(i,active_size);

  for(int j=0;j<active_size;j++)
    {
      if(y[j]==+1)
	{
	  if (!is_lower_bound(j))
	    {
	      double grad_diff=Gmax+G[j];
	      if (G[j] >= Gmax2)
		Gmax2 = G[j];
	      if (grad_diff > 0)
		{
		  double obj_diff; 
		  double quad_coef=Q_i[i]+QD[j]-2.0*y[i]*Q_i[j];
		  if (quad_coef > 0)
		    obj_diff = -(grad_diff*grad_diff)/quad_coef;
		  else
		    obj_diff = -(grad_diff*grad_diff)/TAU;

		  if (obj_diff <= obj_diff_min)
		    {
		      Gmin_idx=j;
		      obj_diff_min = obj_diff;
		    }
		}
	    }
	}
      else
	{
	  if (!is_upper_bound(j))
	    {
	      double grad_diff= Gmax-G[j];
	      if (-G[j] >= Gmax2)
		Gmax2 = -G[j];
	      if (grad_diff > 0)
		{
		  double obj_diff; 
		  double quad_coef=Q_i[i]+QD[j]+2.0*y[i]*Q_i[j];
		  if (quad_coef > 0)
		    obj_diff = -(grad_diff*grad_diff)/quad_coef;
		  else
		    obj_diff = -(grad_diff*grad_diff)/TAU;

		  if (obj_diff <= obj_diff_min)
		    {
		      Gmin_idx=j;
		      obj_diff_min = obj_diff;
		    }
		}
	    }
	}
    }

  if(Gmax+Gmax2 < eps)
    return 1;

  out_i = Gmax_idx;
  out_j = Gmin_idx;

  return 0;
}

int Solver::select_working_set_hmg2(double* u, double& lambda_star, double& best_gain)
{
  for(int i=0; i<curr_depth+2; i++) {
    active[work[i]] = false;
    work[i] = -1;
  }
  curr_depth = -1;

  if(select_working_set_first_order())
    return 1;

  // generic situation: use the MG selection
  int a, b, bb, out_i=-1, out_j=-1;
  double da, db;					// diagonal entries of Q
  double ga, gb;					// gradient in coordinates a and b
  double gt, gs;					// gradient in coordinate a+b or a-b
  double alpha_a, alpha_b;
  double progress;
  double nu;
  double Ca, Cb;
  double lambda;
  Qfloat* q;

  double best = 0.0;
  double g_best = 0.0;

  // try combinations with the last working set
  for(bb=0; bb<2; bb++) {
    b = work[bb];

    q = Q->get_Q(b, active_size);

    db = QD[b];
    Cb = get_C(b);
    gb = G_cg[0][b];
    alpha_b = alpha_cg[0][b];
 
    for(a=0; a<active_size; a++) {
      if(a == b) 
	continue;
      da = QD[a];
      Ca = get_C(a);
      ga = G_cg[0][a];
      alpha_a = alpha_cg[0][a];
      nu = q[a];

      if(y[a]*y[b] > 0.0) {
	lambda = da + db - 2.0*nu;
	gs = gt = (ga - gb) / lambda;
	if(gs > 0.0) {
	  if(alpha_a <= 1e-12 || alpha_b >= Cb-1e-12) 
	    continue;
	  if(gs < -alpha_a) 
	    gs = -alpha_a;
	  if(gs < alpha_b - Cb) 
	    gs = alpha_b - Cb;
	}
	else {
	  if(alpha_a >= Ca-1e-12 || alpha_b <= 1e-12) 
	    continue;
	  if(gs > Ca - alpha_a) 
	    gs = Ca - alpha_a;
	  if(gs > alpha_b) 
	    gs = alpha_b;
	}
	progress = gs * (2.0 * gt - gs) * lambda;
      }
      else {
	lambda = da + db + 2.0*nu;
	gs = gt = (ga + gb) / lambda;
	if(gs > 0.0) {
	  if(alpha_a <= 1e-12 || alpha_b <= 1e-12) 
	    continue;
	  if(gs < -alpha_a) 
	    gs = -alpha_a;
	  if(gs < -alpha_b) 
	    gs = -alpha_b;
	}
	else {
	  if(alpha_a >= Ca-1e-12 || alpha_b >= Cb-1e-12) 
	    continue;
	  if(gs > Ca - alpha_a) 
	    gs = Ca - alpha_a;
	  if(gs > Cb - alpha_b) 
	    gs = Cb - alpha_b;
	}
	progress = gs * (2.0 * gt - gs) * lambda;
      }

      // select the largest progress
      if(progress > best) {
	best = progress;
	g_best = gs;
	out_i = a;
	out_j = b;
      }
    }
  }

  // stopping condition
  if(fabs(g_best) < eps) 
    return 1;			// optimal

  best_gain = best;
  work[0] = out_i;
  work[1] = out_j;
  active[out_i] = true;
  active[out_j] = true;

  compute_step_first_order(u, lambda_star);
  curr_depth = 0;
  return 0;
}

bool Solver::feasible_direction(double *u)
{
  int i;
  double alpha_i;
  double *alpha0 = alpha_cg[0];

  for(i=0; i<curr_depth+3; i++) {
    alpha_i = alpha0[work[i]];
    if((u[i]>0 && alpha_i>=get_C(work[i])) || (u[i]<0 && alpha_i<=0))
      return false;
  }

  return true;
}

// select k given i and j
int Solver::select_working_set_incrementally(double* best_u, double& out_lambda_star, double& final_gain)
{
  int i,j,k,best_indx=-1;
  double gain, best_gain = 0;
  double nominator, denominator;
  double A[curr_depth+3], B[curr_depth+3], lambda_init, lambda_star;
  bool feasible_u, feasible_minus_u, negligible;
  double u[curr_depth+3], minus_u[curr_depth+3], *curr_u;
  Qfloat* Q_row[curr_depth+3];
  double *G0 = G_cg[0];
  Qfloat *Q_row_i;
  double *alpha0;

  for(i=0; i<curr_depth+2; i++)
    Q_row[i] = Q->get_Q(work[i],active_size);
  
  for(k = 0; k < active_size; k++) {

    // check that k is not in the working set
    if(active[k])
      continue;

    // generate direction
    work[curr_depth+2] = k;
    if(generate_direction(u,curr_depth+3, k==0))
      continue;

    // check if u or -u is feasible
    for(i=0; i<curr_depth+3; i++)
      minus_u[i] = -u[i];

    feasible_u = feasible_direction(u);
    feasible_minus_u = feasible_direction(minus_u);

    if(!feasible_u && !feasible_minus_u)
      continue;

    // choose which of {u,minus_u} is a descent direction
    nominator = 0;
    for(i = 0; i<curr_depth+3; i++)
      nominator += G0[work[i]]*u[i];
    
    if(nominator<0) 
      if(feasible_u)
	curr_u = u;
      else
        continue;       
    else
      if(feasible_minus_u) {
	nominator *= -1;
        curr_u = minus_u;     
      }
      else
	continue;

    // compute denominator
    denominator = 0;
    for(i = 0; i<curr_depth+3; i++) {
      denominator += curr_u[i]*curr_u[i]*QD[work[i]];
      Q_row_i = Q_row[i];
      for(j = i+1; j<curr_depth+3; j++)
        denominator += 2*curr_u[i]*curr_u[j]*Q_row_i[work[j]];
    }

    if(denominator<=0)
      denominator=TAU;
   
    lambda_init = -nominator/denominator;
    lambda_star = lambda_init;

    for(i=0; i<curr_depth+3; i++)
      if(curr_u[i]>0) {
	A[i] = 0;
	B[i] = get_C(work[i]);
      }
      else {
	A[i] = get_C(work[i]);
	B[i] = 0;
      }

    alpha0 = alpha_cg[0];
    for(i=0; i<curr_depth+3; i++)
      lambda_star = max(lambda_star,(A[i]-alpha0[work[i]])/curr_u[i]);
    for(i=0; i<curr_depth+3; i++)
      lambda_star = min(lambda_star,(B[i]-alpha0[work[i]])/curr_u[i]);

    if(fabs(lambda_star)<1e-3)
      continue;
    negligible = false;
    for(i=0; i<curr_depth+3; i++)       
      if(fabs(lambda_star*curr_u[i])<1e-3) {
        negligible = true;
        break;
      }
    if(negligible)
      continue;

    gain = denominator*lambda_star*(2*lambda_init-lambda_star); 

    if(gain > best_gain) {
      best_gain = gain;
      best_indx = k;
      memcpy(best_u, curr_u, sizeof_double*(curr_depth+3));
      out_lambda_star = lambda_star;
    }
  }

  final_gain = best_gain;

  if(best_indx != -1) {
    work[curr_depth+2] = best_indx;
    active[best_indx] = true;
    curr_depth++;
    return 0;
  }
  else {
    work[curr_depth+2] = -1;
    return 1; 
  }
}

// select working set incrementally
int Solver::select_working_set_hmg(double* u, double& lambda_star, double& best_gain)
{
  int i;

  if(curr_depth != max_depth) 
    return select_working_set_incrementally(u, lambda_star, best_gain);

  // throw away the oldest active example and select the working set incrementally
  active[work[0]] = false;
  for(i=0; i<=curr_depth; i++)
    work[i] = work[i+1];
  work[curr_depth+1]=-1;
  curr_depth--;

  return select_working_set_incrementally(u, lambda_star, best_gain);
}

int Solver::wss_first_order(double* u, double& lambda_star)
{
   for(int i=0; i<curr_depth+2; i++) {
		active[work[i]] = false;
		work[i] = -1;
	      }
	      curr_depth = -1;
  if(select_working_set_first_order())
    return 1;
  active[work[0]] = true;
  active[work[1]] = true;
  compute_step_first_order(u, lambda_star);
  curr_depth = 0;
  return 0;
}

// return 1 if already optimal, return 0 otherwise
int Solver::select_working_set_first_order()
{
  
  // return i,j which maximize -grad(f)^T d , under constraint
  // if alpha_i == C, d != +1
  // if alpha_i == 0, d != -1

  double Gmax1 = -INF;		// max { -grad(f)_i * d | y_i*d = +1 }
  int Gmax1_idx = -1;

  double Gmax2 = -INF;		// max { -grad(f)_i * d | y_i*d = -1 }// 	
  int Gmax2_idx = -1;
  double *G0 = G_cg[0];
  
  for(int i=0;i<active_size;i++) {
    if(y[i]==+1) {	// y = +1
      if(!is_upper_bound_cg(i,0)) {	// d = +1
	if(-G0[i] >= Gmax1) {
	  Gmax1 = -G0[i];
	  Gmax1_idx = i;
	}
      }
      if(!is_lower_bound_cg(i,0)) {	// d = -1
	if(G0[i] >= Gmax2) {
	  Gmax2 = G0[i];
	  Gmax2_idx = i;
	}
      }
    }
    else {		// y = -1
      if(!is_upper_bound_cg(i,0)) {	// d = +1
	if(-G0[i] >= Gmax2) {
	  Gmax2 = -G0[i];
	  Gmax2_idx = i;
	}
      }
      if(!is_lower_bound_cg(i,0)) {	// d = -1
	if(G0[i] >= Gmax1) {
	  Gmax1 = G0[i];
	  Gmax1_idx = i;
	}
      }
    }
  }

  if(Gmax1+Gmax2 < eps) 
    return 1;
	
  work[0] = Gmax1_idx;
  work[1] = Gmax2_idx;
 
  return 0;
}

void Solver::compute_step_first_order(double* u, double& lambda_star)
{
  int out_i = work[0];
  int out_j = work[1];
  
  Qfloat *Q_i = Q->get_Q(out_i, active_size);
  Qfloat *Q_j = Q->get_Q(out_j, active_size);
 
  double C_i = get_C(out_i);
  double C_j = get_C(out_j);
  double old_alpha_i = alpha_cg[0][out_i];
  double old_alpha_j = alpha_cg[0][out_j];
  double new_alpha_i, new_alpha_j;
  if(y[out_i]!=y[out_j]) {
    double quad_coef = Q_i[out_i]+Q_j[out_j]+2*Q_i[out_j];
    if (quad_coef <= 0)
      quad_coef = TAU;
    lambda_star = (-G_cg[0][out_i]-G_cg[0][out_j])/quad_coef;
    u[0] = 1;
    u[1] = 1;
    double diff = old_alpha_i - old_alpha_j;
    new_alpha_i = old_alpha_i + lambda_star;
    new_alpha_j = old_alpha_j + lambda_star;			
    if(diff > 0) {
      if(new_alpha_j < 0)
	lambda_star = -old_alpha_j;
    }
    else {
      if(new_alpha_i < 0)
	lambda_star = -old_alpha_i;
    }
    if(diff > C_i - C_j) {
      if(new_alpha_i > C_i)
	lambda_star = C_i-old_alpha_i;
    }
    else {
      if(new_alpha_j > C_j)
	lambda_star = C_j-old_alpha_j;
    }
  }
  else {
    double quad_coef = Q_i[out_i]+Q_j[out_j]-2*Q_i[out_j];
    if (quad_coef <= 0)
      quad_coef = TAU;
    lambda_star = (G_cg[0][out_i]-G_cg[0][out_j])/quad_coef;
    u[0] = -1;
    u[1] = 1;
    double sum = old_alpha_i + old_alpha_j;
    new_alpha_i = old_alpha_i - lambda_star;
    new_alpha_j = old_alpha_j + lambda_star;

    if(sum > C_i) {
      if(new_alpha_i > C_i)
	lambda_star = -C_i+old_alpha_i;
    }
    else {
      if(new_alpha_j < 0)
	lambda_star = -old_alpha_j;
    }
    if(sum > C_j) {
      if(new_alpha_j > C_j)
	lambda_star = C_j-old_alpha_j;  
    }
    else {
      if(new_alpha_i < 0)
	lambda_star = old_alpha_i;
    }
  }
}

bool Solver::be_shrunk(int i, double Gmax1, double Gmax2)
{
  if(is_upper_bound(i))
    {
      if(y[i]==+1)
	return(-G[i] > Gmax1);
      else
	return(-G[i] > Gmax2);
    }
  else if(is_lower_bound(i))
    {
      if(y[i]==+1)
	return(G[i] > Gmax2);
      else	
	return(G[i] > Gmax1);
    }
  else
    return(false);
}

bool Solver::be_shrunk_cg(int i, double Gmax1, double Gmax2)
{
  if(is_upper_bound_cg(i,0))
    {
      if(y[i]==+1)
	return(-G_cg[0][i] > Gmax1);
      else
	return(-G_cg[0][i] > Gmax2);
    }
  else if(is_lower_bound_cg(i,0))
    {
      if(y[i]==+1)
	return(G_cg[0][i] > Gmax2);
      else	
	return(G_cg[0][i] > Gmax1);
    }
  else
    return(false);
}


void Solver::do_shrinking()
{
  int i;
  double Gmax1 = -INF;		// max { -y_i * grad(f)_i | i in I_up(\alpha) }
  double Gmax2 = -INF;		// max { y_i * grad(f)_i | i in I_low(\alpha) }

  // find maximal violating pair first
  for(i=0;i<active_size;i++)
    {
      if(y[i]==+1)	
	{
	  if(!is_upper_bound(i))	
	    {
	      if(-G[i] >= Gmax1)
		Gmax1 = -G[i];
	    }
	  if(!is_lower_bound(i))	
	    {
	      if(G[i] >= Gmax2)
		Gmax2 = G[i];
	    }
	}
      else	
	{
	  if(!is_upper_bound(i))	
	    {
	      if(-G[i] >= Gmax2)
		Gmax2 = -G[i];
	    }
	  if(!is_lower_bound(i))	
	    {
	      if(G[i] >= Gmax1)
		Gmax1 = G[i];
	    }
	}
    }

  if(unshrink == false && Gmax1 + Gmax2 <= eps*10) 
    {
      unshrink = true;
      reconstruct_gradient();
      active_size = l;
      info("*");
    }

  for(i=0;i<active_size;i++) {
    if (be_shrunk(i, Gmax1, Gmax2))
      {
	active_size--;
	while (active_size > i)
	  {
	    if (!be_shrunk(active_size, Gmax1, Gmax2))
	      {
		swap_index(i,active_size);
		break;
	      }
	    active_size--;
	  }
      }
  }
}

bool Solver::do_shrinking_cg()
{
  int i,j;
  double Gmax1 = -INF;		// max { -y_i * grad(f)_i | i in I_up(\alpha) }
  double Gmax2 = -INF;		// max { y_i * grad(f)_i | i in I_low(\alpha) }
  bool done_shrinking = false;

  // find maximal violating pair first
  for(i=0;i<active_size;i++)
    if(y[i]==1) {
      if(!is_upper_bound_cg(i,0) && -G_cg[0][i] >= Gmax1) 
	  Gmax1 = -G_cg[0][i];
      if(!is_lower_bound_cg(i,0) && G_cg[0][i] >= Gmax2) 
	  Gmax2 = G_cg[0][i];
    }
    else {
      if(!is_upper_bound_cg(i,0) && -G_cg[0][i] >= Gmax2) 
	  Gmax2 = -G_cg[0][i];
      if(!is_lower_bound_cg(i,0) && G_cg[0][i] >= Gmax1) 
	  Gmax1 = G_cg[0][i];
    }

  if(unshrink == false && Gmax1 + Gmax2 <= eps*10) {
    unshrink = true;
    reconstruct_gradient();
    active_size = l;
    info("*");
  }

  for(i=0;i<active_size;i++) {
    if(be_shrunk_cg(i, Gmax1, Gmax2)) {
      active_size--;
      while(active_size > i) {
	if(!be_shrunk_cg(active_size, Gmax1, Gmax2)) {
	  swap_index(i,active_size);
          done_shrinking = true;
	  break;
	}
        else {
          if(active[active_size]) {
            active[active_size] = false;
            for(j=0; j<curr_depth+2; j++)
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
  // shrink working set
  if(done_shrinking) {
    for(i=0; i<curr_depth+2; i++) {
      active[work[i]] = false;
      work[i] = -1;
    }
  curr_depth = -1;
  }

  return done_shrinking;
}

double Solver::calculate_rho()
{
  double r;
  int nr_free = 0;
  double ub = INF, lb = -INF, sum_free = 0;
  for(int i=0;i<active_size;i++)
    {
      double yG = y[i]*G[i];

      if(is_upper_bound(i))
	{
	  if(y[i]==-1)
	    ub = min(ub,yG);
	  else
	    lb = max(lb,yG);
	}
      else if(is_lower_bound(i))
	{
	  if(y[i]==+1)
	    ub = min(ub,yG);
	  else
	    lb = max(lb,yG);
	}
      else
	{
	  ++nr_free;
	  sum_free += yG;
	}
    }

  if(nr_free>0)
    r = sum_free/nr_free;
  else
    r = (ub+lb)/2;

  return r;
}

double Solver::calculate_rho_cg()
{
  double r;
  int nr_free = 0;
  double ub = INF, lb = -INF, sum_free = 0;
  for(int i=0;i<active_size;i++)
    {
      double yG = y[i]*G_cg[0][i];

      if(is_upper_bound_cg(i,0))
	{
	  if(y[i]==-1)
	    ub = min(ub,yG);
	  else
	    lb = max(lb,yG);
	}
      else if(is_lower_bound_cg(i,0))
	{
	  if(y[i]==+1)
	    ub = min(ub,yG);
	  else
	    lb = max(lb,yG);
	}
      else
	{
	  ++nr_free;
	  sum_free += yG;
	}
    }

  if(nr_free>0)
    r = sum_free/nr_free;
  else
    r = (ub+lb)/2;

  return r;
}

