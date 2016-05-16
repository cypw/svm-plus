#include "common.h"
#include "Kernel.h"
#include "Solver.h"
#include "Solver_plus.h"
#include "Solver_NU.h"

//
// construct and solve various formulations
//
static void solve_c_svc(
			const svm_problem *prob, const svm_parameter* param,
			double *alpha, SolutionInfo* si, double Cp, double Cn)
{
  int l = prob->l;
  double *minus_ones = new double[l];
  schar *y = new schar[l];

  int i;

  for(i=0;i<l;i++)
    {
      alpha[i] = 0;
      minus_ones[i] = -1;
      if(prob->y[i] > 0) y[i] = +1; else y[i]=-1;
    }

  Solver s(param->optimizer);
  if(param->optimizer != -1)
    s.Solve_cg(l, SVC_Q(*prob,*param,y), minus_ones, y,
	       alpha, Cp, Cn, param->eps, si, param->shrinking);
  else
    s.Solve(l, SVC_Q(*prob,*param,y), minus_ones, y,
	    alpha, Cp, Cn, param->eps, si, param->shrinking);

  double sum_alpha=0;
  for(i=0;i<l;i++)
    sum_alpha += alpha[i];

  if (Cp==Cn)
    info("nu = %f\n", sum_alpha/(Cp*prob->l));

  for(i=0;i<l;i++)
    alpha[i] *= y[i];

  delete[] minus_ones;
  delete[] y;
}

static void solve_svm_plus(const svm_problem *prob, const svm_parameter* param,
			   double *alpha, double *beta, SolutionInfo* si, double Cp, double Cn)
{
  int l = prob->l;
  schar *y = new schar[l];
  schar *y_true = new schar[l];
  svm_parameter cheat_param, cheat_param2; 
  svm_problem cheat_problem, cheat_problem2;
  int i;
  int l_pos=0, l_neg=0;

  // initialize alphas and betas
  for(i=0;i<l;i++) {
      alpha[i] = 0;
      beta[i] = Cp; 
      y[i] = 1;
      if(prob->y[i] > 0) { 
        y_true[i] = +1; 
        l_pos++;
      }
      else {
        y_true[i]=-1;
        l_neg++;
      }
  }

  cheat_param = *param;
  cheat_param.kernel_type = 2;
  cheat_param.gamma = param->gamma_star;
  cheat_problem = *prob;
  cheat_problem.x = Malloc(struct svm_node*,prob->l);
  memcpy(cheat_problem.x, prob->x_star, l*sizeof(struct svm_node*));
  cheat_param2 = *param;
  cheat_param2.kernel_type = 2;
  cheat_param2.gamma = param->gamma_star;
  cheat_problem2 = *prob;
  cheat_problem2.x = Malloc(struct svm_node*,prob->l);
  memcpy(cheat_problem2.x, prob->x_star, l*sizeof(struct svm_node*));
  
  SVC_Q kernel1 = SVC_Q(*prob,*param,y);
  SVC_Q kernel2 = SVC_Q(cheat_problem, cheat_param, y);
  SVC_Q kernel3 = SVC_Q(cheat_problem2, cheat_param2, y);

  Solver_plus s(param->optimizer);
  
  if(param->optimizer != -1) 
    s.Solve_plus_cg(l, kernel1, kernel2, 
		    kernel3, y_true,
		    alpha, beta, Cp, Cn, param->tau, param->eps, si, param->shrinking);
  else
    s.Solve_plus(l, kernel1, kernel2, 
		 kernel3, y_true,
		 alpha, beta, Cp, Cn, param->tau, param->eps, si, param->shrinking);

  // produce the same output as SVM
  for(i=0;i<l;i++) { 
    if(alpha[i]<0)
      alpha[i] = 0;
    alpha[i] *= y_true[i];
  }

  delete[] y;
  delete[] y_true;
}

static void solve_nu_svc(
			 const svm_problem *prob, const svm_parameter *param,
			 double *alpha, SolutionInfo* si)
{
  int i;
  int l = prob->l;
  double nu = param->nu;

  schar *y = new schar[l];

  for(i=0;i<l;i++)
    if(prob->y[i]>0)
      y[i] = +1;
    else
      y[i] = -1;

  double sum_pos = nu*l/2;
  double sum_neg = nu*l/2;

  for(i=0;i<l;i++)
    if(y[i] == +1)
      {
	alpha[i] = min(1.0,sum_pos);
	sum_pos -= alpha[i];
      }
    else
      {
	alpha[i] = min(1.0,sum_neg);
	sum_neg -= alpha[i];
      }

  double *zeros = new double[l];

  for(i=0;i<l;i++)
    zeros[i] = 0;

  Solver_NU s;
  s.Solve(l, SVC_Q(*prob,*param,y), zeros, y,
	  alpha, 1.0, 1.0, param->eps, si,  param->shrinking);
  double r = si->r;

  info("C = %f\n",1/r);

  for(i=0;i<l;i++)
    alpha[i] *= y[i]/r;

  si->rho /= r;
  si->obj /= (r*r);
  si->upper_bound_p = 1/r;
  si->upper_bound_n = 1/r;

  delete[] y;
  delete[] zeros;
}

static void solve_one_class(
			    const svm_problem *prob, const svm_parameter *param,
			    double *alpha, SolutionInfo* si)
{
  int l = prob->l;
  double *zeros = new double[l];
  schar *ones = new schar[l];
  int i;

  int n = (int)(param->nu*prob->l);	// # of alpha's at upper bound

  for(i=0;i<n;i++)
    alpha[i] = 1;
  if(n<prob->l)
    alpha[n] = param->nu * prob->l - n;
  for(i=n+1;i<l;i++)
    alpha[i] = 0;

  for(i=0;i<l;i++)
    {
      zeros[i] = 0;
      ones[i] = 1;
    }

  Solver s;
  s.Solve(l, ONE_CLASS_Q(*prob,*param), zeros, ones,
	  alpha, 1.0, 1.0, param->eps, si, param->shrinking);

  delete[] zeros;
  delete[] ones;
}

static void solve_epsilon_svr(
			      const svm_problem *prob, const svm_parameter *param,
			      double *alpha, SolutionInfo* si)
{
  int l = prob->l;
  double *alpha2 = new double[2*l];
  double *linear_term = new double[2*l];
  schar *y = new schar[2*l];
  int i;

  for(i=0;i<l;i++)
    {
      alpha2[i] = 0;
      linear_term[i] = param->p - prob->y[i];
      y[i] = 1;

      alpha2[i+l] = 0;
      linear_term[i+l] = param->p + prob->y[i];
      y[i+l] = -1;
    }

  Solver s;
  s.Solve(2*l, SVR_Q(*prob,*param), linear_term, y,
	  alpha2, param->C, param->C, param->eps, si, param->shrinking);

  double sum_alpha = 0;
  for(i=0;i<l;i++)
    {
      alpha[i] = alpha2[i] - alpha2[i+l];
      sum_alpha += fabs(alpha[i]);
    }
  info("nu = %f\n",sum_alpha/(param->C*l));

  delete[] alpha2;
  delete[] linear_term;
  delete[] y;
}

static void solve_nu_svr(
			 const svm_problem *prob, const svm_parameter *param,
			 double *alpha, SolutionInfo* si)
{
  int l = prob->l;
  double C = param->C;
  double *alpha2 = new double[2*l];
  double *linear_term = new double[2*l];
  schar *y = new schar[2*l];
  int i;

  double sum = C * param->nu * l / 2;
  for(i=0;i<l;i++)
    {
      alpha2[i] = alpha2[i+l] = min(sum,C);
      sum -= alpha2[i];

      linear_term[i] = - prob->y[i];
      y[i] = 1;

      linear_term[i+l] = prob->y[i];
      y[i+l] = -1;
    }

  Solver_NU s;
  s.Solve(2*l, SVR_Q(*prob,*param), linear_term, y,
	  alpha2, C, C, param->eps, si, param->shrinking);

  info("epsilon = %f\n",-si->r);

  for(i=0;i<l;i++)
    alpha[i] = alpha2[i] - alpha2[i+l];

  delete[] alpha2;
  delete[] linear_term;
  delete[] y;
}

//
// decision_function
//
struct decision_function
{
  double *alpha;
  double *beta;
  double rho;	
  double rho_star;
};

decision_function svm_train_one(const svm_problem *prob, const svm_parameter *param,
				double Cp, double Cn)
{
  double *alpha = Malloc(double,prob->l);
  double *beta = NULL;

  if(param->svm_type == SVM_PLUS)
    beta = Malloc(double,prob->l);

  SolutionInfo si;
  switch(param->svm_type) {
  case C_SVC:
    solve_c_svc(prob,param,alpha,&si,Cp,Cn);
    break;
  case SVM_PLUS:
    solve_svm_plus(prob,param,alpha,beta,&si,Cp,Cn); 
    break;
  case NU_SVC:
    solve_nu_svc(prob,param,alpha,&si);
    break;
  case ONE_CLASS:
    solve_one_class(prob,param,alpha,&si);
    break;
  case EPSILON_SVR:
    solve_epsilon_svr(prob,param,alpha,&si);
    break;
  case NU_SVR:
    solve_nu_svr(prob,param,alpha,&si);
    break;
  }

  info("obj = %f, rho = %f\n",si.obj,si.rho);

  // output SVs

  int nSV = 0;
  int nSV_star = 0;
  int nBSV = 0;
  int nBSV_star = 0;
  for(int i=0;i<prob->l;i++) {
    if(fabs(alpha[i]) > 0) {
      ++nSV;
      if (param->svm_type != SVM_PLUS) {
	if(prob->y[i] > 0) {
	  if(fabs(alpha[i]) >= si.upper_bound_p)
	    ++nBSV;
	}
	else {
	  if(fabs(alpha[i]) >= si.upper_bound_n)
	    ++nBSV;
	}
      }
    }
    if (param->svm_type == SVM_PLUS) {
      if(fabs(beta[i]) > 0)
	++nSV_star;
    }
  }

  info("nSV = %d, nBSV = %d\n",nSV,nBSV);
  if (param->svm_type == SVM_PLUS)
    info("nSV_star = %d nBSV_star = %d \n",nSV_star,nBSV_star);

  decision_function f;
  f.alpha = alpha;
  f.rho = si.rho;
  if (param->svm_type == SVM_PLUS) {
    f.beta = beta;
    f.rho_star = si.rho_star;
  }

  return f;
}

//
// svm_model
//
struct svm_model
{
  svm_parameter param;	// parameter
  int nr_class;		// number of classes, = 2 in regression/one class svm
  int l;			// total #SV
  int l_star;		// total #SV_stars
  svm_node **SV;		// SVs (SV[l])
  svm_node **SV_star;     // SV_stars (SV_star[l_star])
  double **sv_coef;	// coefficients for SVs in decision functions (sv_coef[k-1][l])
  double **sv_coef_star;	// coefficients for SVs in correcting functions (sv_coef[k-1][l]) 
  double *rho;		// constants in decision functions (rho[k*(k-1)/2])
  double *rho_star;       // constants in correcting functions (rho[k*(k-1)/2])
  double *probA;		// pariwise probability information
  double *probB;

  // for classification only

  int *label;		// label of each class (label[k])
  int *nSV;		// number of SVs for each class (nSV[k])
  // nSV[0] + nSV[1] + ... + nSV[k-1] = l
  int *nSV_star;		// number of SV_stars for each class (nSV_star[k])
				// nSV_star[0] + nSV_star[1] + ... + nSV_star[k-1] = l_star
  // XXX
  int free_sv;		// 1 if svm_model is created by svm_load_model
  // 0 if svm_model is created by svm_train
};

// Platt's binary SVM Probablistic Output: an improvement from Lin et al.
void sigmoid_train(
		   int l, const double *dec_values, const double *labels, 
		   double& A, double& B)
{
  double prior1=0, prior0 = 0;
  int i;

  for (i=0;i<l;i++)
    if (labels[i] > 0) prior1+=1;
    else prior0+=1;
	
  int max_iter=100;	// Maximal number of iterations
  double min_step=1e-10;	// Minimal step taken in line search
  double sigma=1e-12;	// For numerically strict PD of Hessian
  double eps=1e-5;
  double hiTarget=(prior1+1.0)/(prior1+2.0);
  double loTarget=1/(prior0+2.0);
  double *t=Malloc(double,l);
  double fApB,p,q,h11,h22,h21,g1,g2,det,dA,dB,gd,stepsize;
  double newA,newB,newf,d1,d2;
  int iter; 
	
  // Initial Point and Initial Fun Value
  A=0.0; B=log((prior0+1.0)/(prior1+1.0));
  double fval = 0.0;

  for (i=0;i<l;i++)
    {
      if (labels[i]>0) t[i]=hiTarget;
      else t[i]=loTarget;
      fApB = dec_values[i]*A+B;
      if (fApB>=0)
	fval += t[i]*fApB + log(1+exp(-fApB));
      else
	fval += (t[i] - 1)*fApB +log(1+exp(fApB));
    }
  for (iter=0;iter<max_iter;iter++)
    {
      // Update Gradient and Hessian (use H' = H + sigma I)
      h11=sigma; // numerically ensures strict PD
      h22=sigma;
      h21=0.0;g1=0.0;g2=0.0;
      for (i=0;i<l;i++)
	{
	  fApB = dec_values[i]*A+B;
	  if (fApB >= 0)
	    {
	      p=exp(-fApB)/(1.0+exp(-fApB));
	      q=1.0/(1.0+exp(-fApB));
	    }
	  else
	    {
	      p=1.0/(1.0+exp(fApB));
	      q=exp(fApB)/(1.0+exp(fApB));
	    }
	  d2=p*q;
	  h11+=dec_values[i]*dec_values[i]*d2;
	  h22+=d2;
	  h21+=dec_values[i]*d2;
	  d1=t[i]-p;
	  g1+=dec_values[i]*d1;
	  g2+=d1;
	}

      // Stopping Criteria
      if (fabs(g1)<eps && fabs(g2)<eps)
	break;

      // Finding Newton direction: -inv(H') * g
      det=h11*h22-h21*h21;
      dA=-(h22*g1 - h21 * g2) / det;
      dB=-(-h21*g1+ h11 * g2) / det;
      gd=g1*dA+g2*dB;


      stepsize = 1;		// Line Search
      while (stepsize >= min_step)
	{
	  newA = A + stepsize * dA;
	  newB = B + stepsize * dB;

	  // New function value
	  newf = 0.0;
	  for (i=0;i<l;i++)
	    {
	      fApB = dec_values[i]*newA+newB;
	      if (fApB >= 0)
		newf += t[i]*fApB + log(1+exp(-fApB));
	      else
		newf += (t[i] - 1)*fApB +log(1+exp(fApB));
	    }
	  // Check sufficient decrease
	  if (newf<fval+0.0001*stepsize*gd)
	    {
	      A=newA;B=newB;fval=newf;
	      break;
	    }
	  else
	    stepsize = stepsize / 2.0;
	}

      if (stepsize < min_step)
	{
	  info("Line search fails in two-class probability estimates\n");
	  break;
	}
    }

  if (iter>=max_iter)
    info("Reaching maximal iterations in two-class probability estimates\n");
  free(t);
}

double sigmoid_predict(double decision_value, double A, double B)
{
  double fApB = decision_value*A+B;
  if (fApB >= 0)
    return exp(-fApB)/(1.0+exp(-fApB));
  else
    return 1.0/(1+exp(fApB)) ;
}

// Method 2 from the multiclass_prob paper by Wu, Lin, and Weng
void multiclass_probability(int k, double **r, double *p)
{
  int t,j;
  int iter = 0, max_iter=max(100,k);
  double **Q=Malloc(double *,k);
  double *Qp=Malloc(double,k);
  double pQp, eps=0.005/k;
	
  for (t=0;t<k;t++)
    {
      p[t]=1.0/k;  // Valid if k = 1
      Q[t]=Malloc(double,k);
      Q[t][t]=0;
      for (j=0;j<t;j++)
	{
	  Q[t][t]+=r[j][t]*r[j][t];
	  Q[t][j]=Q[j][t];
	}
      for (j=t+1;j<k;j++)
	{
	  Q[t][t]+=r[j][t]*r[j][t];
	  Q[t][j]=-r[j][t]*r[t][j];
	}
    }
  for (iter=0;iter<max_iter;iter++)
    {
      // stopping condition, recalculate QP,pQP for numerical accuracy
      pQp=0;
      for (t=0;t<k;t++)
	{
	  Qp[t]=0;
	  for (j=0;j<k;j++)
	    Qp[t]+=Q[t][j]*p[j];
	  pQp+=p[t]*Qp[t];
	}
      double max_error=0;
      for (t=0;t<k;t++)
	{
	  double error=fabs(Qp[t]-pQp);
	  if (error>max_error)
	    max_error=error;
	}
      if (max_error<eps) break;
		
      for (t=0;t<k;t++)
	{
	  double diff=(-Qp[t]+pQp)/Q[t][t];
	  p[t]+=diff;
	  pQp=(pQp+diff*(diff*Q[t][t]+2*Qp[t]))/(1+diff)/(1+diff);
	  for (j=0;j<k;j++)
	    {
	      Qp[j]=(Qp[j]+diff*Q[t][j])/(1+diff);
	      p[j]/=(1+diff);
	    }
	}
    }
  if (iter>=max_iter)
    info("Exceeds max_iter in multiclass_prob\n");
  for(t=0;t<k;t++) free(Q[t]);
  free(Q);
  free(Qp);
}

// Cross-validation decision values for probability estimates
void svm_binary_svc_probability(
				const svm_problem *prob, const svm_parameter *param,
				double Cp, double Cn, double& probA, double& probB)
{
  int i;
  int nr_fold = 5;
  int *perm = Malloc(int,prob->l);
  double *dec_values = Malloc(double,prob->l);

  // random shuffle
  for(i=0;i<prob->l;i++) perm[i]=i;
  for(i=0;i<prob->l;i++)
    {
      int j = i+rand()%(prob->l-i);
      swap(perm[i],perm[j]);
    }
  for(i=0;i<nr_fold;i++)
    {
      int begin = i*prob->l/nr_fold;
      int end = (i+1)*prob->l/nr_fold;
      int j,k;
      struct svm_problem subprob;

      subprob.l = prob->l-(end-begin);
      subprob.x = Malloc(struct svm_node*,subprob.l);
      subprob.y = Malloc(double,subprob.l);
			
      k=0;
      for(j=0;j<begin;j++)
	{
	  subprob.x[k] = prob->x[perm[j]];
	  subprob.y[k] = prob->y[perm[j]];
	  ++k;
	}
      for(j=end;j<prob->l;j++)
	{
	  subprob.x[k] = prob->x[perm[j]];
	  subprob.y[k] = prob->y[perm[j]];
	  ++k;
	}
      int p_count=0,n_count=0;
      for(j=0;j<k;j++)
	if(subprob.y[j]>0)
	  p_count++;
	else
	  n_count++;

      if(p_count==0 && n_count==0)
	for(j=begin;j<end;j++)
	  dec_values[perm[j]] = 0;
      else if(p_count > 0 && n_count == 0)
	for(j=begin;j<end;j++)
	  dec_values[perm[j]] = 1;
      else if(p_count == 0 && n_count > 0)
	for(j=begin;j<end;j++)
	  dec_values[perm[j]] = -1;
      else
	{
	  svm_parameter subparam = *param;
	  subparam.probability=0;
	  subparam.C=1.0;
	  subparam.nr_weight=2;
	  subparam.weight_label = Malloc(int,2);
	  subparam.weight = Malloc(double,2);
	  subparam.weight_label[0]=+1;
	  subparam.weight_label[1]=-1;
	  subparam.weight[0]=Cp;
	  subparam.weight[1]=Cn;
	  struct svm_model *submodel = svm_train(&subprob,&subparam);
	  for(j=begin;j<end;j++)
	    {
	      svm_predict_values(submodel,prob->x[perm[j]],&(dec_values[perm[j]])); 
	      // ensure +1 -1 order; reason not using CV subroutine
	      dec_values[perm[j]] *= submodel->label[0];
	    }		
	  svm_destroy_model(submodel);
	  svm_destroy_param(&subparam);
	}
      free(subprob.x);
      free(subprob.y);
    }		
  sigmoid_train(prob->l,dec_values,prob->y,probA,probB);
  free(dec_values);
  free(perm);
}

// Return parameter of a Laplace distribution 
double svm_svr_probability(
			   const svm_problem *prob, const svm_parameter *param)
{
  int i;
  int nr_fold = 5;
  double *ymv = Malloc(double,prob->l);
  double mae = 0;

  svm_parameter newparam = *param;
  newparam.probability = 0;
  svm_cross_validation(prob,&newparam,nr_fold,ymv);
  for(i=0;i<prob->l;i++)
    {
      ymv[i]=prob->y[i]-ymv[i];
      mae += fabs(ymv[i]);
    }		
  mae /= prob->l;
  double std=sqrt(2*mae*mae);
  int count=0;
  mae=0;
  for(i=0;i<prob->l;i++)
    if (fabs(ymv[i]) > 5*std) 
      count=count+1;
    else 
      mae+=fabs(ymv[i]);
  mae /= (prob->l-count);
  info("Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma= %g\n",mae);
  free(ymv);
  return mae;
}


// label: label name, start: begin of each class, count: #data of classes, perm: indices to the original data
// perm, length l, must be allocated before calling this subroutine
void svm_group_classes(const svm_problem *prob, int *nr_class_ret, int **label_ret, int **start_ret, int **count_ret, int *perm)
{
  int l = prob->l;
  int max_nr_class = 16;
  int nr_class = 0;
  int *label = Malloc(int,max_nr_class);
  int *count = Malloc(int,max_nr_class);
  int *data_label = Malloc(int,l);	
  int i;

  for(i=0;i<l;i++)
    {
      int this_label = (int)prob->y[i];
      int j;
      for(j=0;j<nr_class;j++)
	{
	  if(this_label == label[j])
	    {
	      ++count[j];
	      break;
	    }
	}
      data_label[i] = j;
      if(j == nr_class)
	{
	  if(nr_class == max_nr_class)
	    {
	      max_nr_class *= 2;
	      label = (int *)realloc(label,max_nr_class*sizeof(int));
	      count = (int *)realloc(count,max_nr_class*sizeof(int));
	    }
	  label[nr_class] = this_label;
	  count[nr_class] = 1;
	  ++nr_class;
	}
    }

  int *start = Malloc(int,nr_class);
  start[0] = 0;
  for(i=1;i<nr_class;i++)
    start[i] = start[i-1]+count[i-1];
  for(i=0;i<l;i++)
    {
      perm[start[data_label[i]]] = i;
      ++start[data_label[i]];
    }
  start[0] = 0;
  for(i=1;i<nr_class;i++)
    start[i] = start[i-1]+count[i-1];

  *nr_class_ret = nr_class;
  *label_ret = label;
  *start_ret = start;
  *count_ret = count;
  free(data_label);
}

//
// Interface functions
//
svm_model *svm_train(const svm_problem *prob, const svm_parameter *param)
{
  svm_model *model = Malloc(svm_model,1);
  model->param = *param;
  model->free_sv = 0;	// XXX

  if(param->svm_type == ONE_CLASS ||
     param->svm_type == EPSILON_SVR ||
     param->svm_type == NU_SVR)
    {
      // regression or one-class-svm
      model->nr_class = 2;
      model->label = NULL;
      model->nSV = NULL;
      model->probA = NULL; model->probB = NULL;
      model->sv_coef = Malloc(double *,1);

      if(param->probability && 
	 (param->svm_type == EPSILON_SVR ||
	  param->svm_type == NU_SVR))
	{
	  model->probA = Malloc(double,1);
	  model->probA[0] = svm_svr_probability(prob,param);
	}

      decision_function f = svm_train_one(prob,param,0,0);
      model->rho = Malloc(double,1);
      model->rho[0] = f.rho;

      int nSV = 0;
      int i;
      for(i=0;i<prob->l;i++)
	if(fabs(f.alpha[i]) > 0) ++nSV;
      model->l = nSV;
      model->SV = Malloc(svm_node *,nSV);
      model->sv_coef[0] = Malloc(double,nSV);
      int j = 0;
      for(i=0;i<prob->l;i++)
	if(fabs(f.alpha[i]) > 0)
	  {
	    model->SV[j] = prob->x[i];
	    model->sv_coef[0][j] = f.alpha[i];
	    ++j;
	  }		

      free(f.alpha);
    }
  else
    {
      // classification
      int l = prob->l;
      int nr_class;
      int *label = NULL;
      int *start = NULL;
      int *count = NULL;
      int *perm = Malloc(int,l);

      /*XXX*/

      // group training data of the same class
      svm_group_classes(prob,&nr_class,&label,&start,&count,perm);		
      svm_node **x = Malloc(svm_node *,l);
      svm_node **x_star = NULL;
      int i;
      for(i=0;i<l;i++)
	x[i] = prob->x[perm[i]];
      if(param->svm_type == SVM_PLUS) {
	x_star = Malloc(svm_node *,l);
	for(i=0;i<l;i++)
	  x_star[i] = prob->x_star[perm[i]];
      }

      // calculate weighted C

      double *weighted_C = Malloc(double, nr_class);
      for(i=0;i<nr_class;i++)
	weighted_C[i] = param->C;
     
      for(i=0;i<param->nr_weight;i++)
	{	
	  int j;
	  for(j=0;j<nr_class;j++)
	    if(param->weight_label[i] == label[j])
	      break;
	  if(j == nr_class)
	    fprintf(stderr,"warning: class label %d specified in weight is not found\n", param->weight_label[i]);
	  else
	    weighted_C[j] *= param->weight[i];
	}

      // train k*(k-1)/2 models
		
      bool *nonzero = Malloc(bool,l);
      bool *nonzero_star = Malloc(bool,l);
      for(i=0;i<l;i++) {
	nonzero[i] = false;
	nonzero_star[i] = false;
      }
      decision_function *f = Malloc(decision_function,nr_class*(nr_class-1)/2);

      double *probA=NULL,*probB=NULL;
      if (param->probability)
	{
	  probA=Malloc(double,nr_class*(nr_class-1)/2);
	  probB=Malloc(double,nr_class*(nr_class-1)/2);
	}

      int p = 0, p_star = 0;
      for(i=0;i<nr_class;i++)
	for(int j=i+1;j<nr_class;j++)
	  {
	    svm_problem sub_prob;
	    int si = start[i], sj = start[j];
	    int ci = count[i], cj = count[j];
	    sub_prob.l = ci+cj;
	    sub_prob.x = Malloc(svm_node *,sub_prob.l);
	    if(param->svm_type == SVM_PLUS)
	      sub_prob.x_star = Malloc(svm_node *,sub_prob.l);    
	    sub_prob.y = Malloc(double,sub_prob.l);
	    int k;
	    for(k=0;k<ci;k++)
	      {
		sub_prob.x[k] = x[si+k];
		sub_prob.y[k] = +1;
		if(param->svm_type == SVM_PLUS)
		  sub_prob.x_star[k] = x_star[si+k];
	      }
	    for(k=0;k<cj;k++)
	      {
		sub_prob.x[ci+k] = x[sj+k];
		sub_prob.y[ci+k] = -1;
		if(param->svm_type == SVM_PLUS)
		  sub_prob.x_star[ci+k] = x_star[sj+k];
	      }

	    if(param->probability)
	      svm_binary_svc_probability(&sub_prob,param,weighted_C[i],weighted_C[j],probA[p],probB[p]);

	    f[p] = svm_train_one(&sub_prob,param,weighted_C[i],weighted_C[j]);

	    //if(param->svm_type == SVM_PLUS)  
	    //  for(int i=0; i<prob->l; i++)
	    //	fprintf(stdout,"%f %f\n", f[p].alpha[i], f[p].beta[i]);
	    //  fflush(stdout);

	    for(k=0;k<ci;k++) {
	      if(!nonzero[si+k] && fabs(f[p].alpha[k]) > 0)
		nonzero[si+k] = true;
	      if(param->svm_type == SVM_PLUS)
		if(!nonzero_star[si+k] && f[p].beta[k] > 0) 
		  nonzero_star[si+k] = true;
	    }
	    for(k=0;k<cj;k++) {
	      if(!nonzero[sj+k] && fabs(f[p].alpha[ci+k]) > 0)
		nonzero[sj+k] = true;
	      if(param->svm_type == SVM_PLUS)
		if(!nonzero_star[sj+k] && f[p].beta[ci+k] > 0) 
		  nonzero_star[sj+k] = true;
	    }
	    free(sub_prob.x);
	    if(param->svm_type == SVM_PLUS)
	      free(sub_prob.x_star);
	    free(sub_prob.y);
	    ++p;
	  }

      // build output

      model->nr_class = nr_class;
		
      model->label = Malloc(int,nr_class);
      for(i=0;i<nr_class;i++)
	model->label[i] = label[i];
		
      model->rho = Malloc(double,nr_class*(nr_class-1)/2);
      for(i=0;i<nr_class*(nr_class-1)/2;i++)
	model->rho[i] = f[i].rho;

      if(param->svm_type == SVM_PLUS) {
	model->rho_star = Malloc(double,nr_class*(nr_class-1)/2);
	for(i=0;i<nr_class*(nr_class-1)/2;i++)
	  model->rho_star[i] = f[i].rho_star;
      }

      if(param->probability)
	{
	  model->probA = Malloc(double,nr_class*(nr_class-1)/2);
	  model->probB = Malloc(double,nr_class*(nr_class-1)/2);
	  for(i=0;i<nr_class*(nr_class-1)/2;i++)
	    {
	      model->probA[i] = probA[i];
	      model->probB[i] = probB[i];
	    }
	}
      else
	{
	  model->probA=NULL;
	  model->probB=NULL;
	}

      int total_sv = 0;
      int total_sv_star = 0;
      int *nz_count = Malloc(int,nr_class);
      int *nz_count_star = Malloc(int,nr_class);
      model->nSV = Malloc(int,nr_class);
      model->nSV_star = Malloc(int,nr_class);
      for(i=0;i<nr_class;i++)
	{
	  int nSV = 0;
	  int nSV_star = 0;
	  for(int j=0;j<count[i];j++) {
	    if(nonzero[start[i]+j])
	      {	
		++nSV;
		++total_sv;
	      }
	    if(nonzero_star[start[i]+j])
	      {	
		++nSV_star;
		++total_sv_star;
	      } 
	  }
	  model->nSV[i] = nSV;
	  nz_count[i] = nSV;
	  model->nSV_star[i] = nSV_star;
	  nz_count_star[i] = nSV_star;
	}
		
      info("Total nSV = %d\n",total_sv);
      info("Total nSV_star = %d\n",total_sv_star);

      model->l = total_sv;
      model->SV = Malloc(svm_node *,total_sv);
      p = 0;
      for(i=0;i<l;i++)
	if(nonzero[i]) model->SV[p++] = x[i];

      if(model->param.svm_type == SVM_PLUS) {
        model->SV_star = Malloc(svm_node *,total_sv_star);
        p_star = 0;
        for(i=0;i<l;i++)
	  if(nonzero_star[i]) model->SV_star[p_star++] = x_star[i]; 
      }

      int *nz_start = Malloc(int,nr_class);
      nz_start[0] = 0;
      for(i=1;i<nr_class;i++)
	nz_start[i] = nz_start[i-1]+nz_count[i-1];

      model->sv_coef = Malloc(double *,nr_class-1);
      for(i=0;i<nr_class-1;i++)
	model->sv_coef[i] = Malloc(double,total_sv);

      int *nz_start_star=NULL;
      if(model->param.svm_type == SVM_PLUS) {
        nz_start_star = Malloc(int,nr_class);
        nz_start_star[0] = 0;
        for(i=1;i<nr_class;i++)
	  nz_start_star[i] = nz_start_star[i-1]+nz_count_star[i-1];

        model->sv_coef_star = Malloc(double *,nr_class-1);
        for(i=0;i<nr_class-1;i++)
	  model->sv_coef_star[i] = Malloc(double,total_sv_star);
      }

      p = 0;
      for(i=0;i<nr_class;i++)
	for(int j=i+1;j<nr_class;j++)
	  {
	    // classifier (i,j): coefficients with
	    // i are in sv_coef[j-1][nz_start[i]...],
	    // j are in sv_coef[i][nz_start[j]...]

	    int si = start[i];
	    int sj = start[j];
	    int ci = count[i];
	    int cj = count[j];
				
	    int q = nz_start[i];
            int q_star = 0;
            if(model->param.svm_type == SVM_PLUS) 
              q_star = nz_start_star[i];
	    int k;
	    for(k=0;k<ci;k++) {
	      if(nonzero[si+k])
		model->sv_coef[j-1][q++] = f[p].alpha[k];
              if(model->param.svm_type == SVM_PLUS) 
                if(nonzero_star[si+k])
		  model->sv_coef_star[j-1][q_star++] = f[p].beta[k];
            }
	    q = nz_start[j];
            if(model->param.svm_type == SVM_PLUS) 
              q_star = nz_start_star[j];
	    for(k=0;k<cj;k++) {
	      if(nonzero[sj+k])
		model->sv_coef[i][q++] = f[p].alpha[ci+k];
              if(model->param.svm_type == SVM_PLUS) 
                if(nonzero_star[sj+k])
		  model->sv_coef_star[i][q_star++] = f[p].beta[ci+k];
            }
	    ++p;
	  }

      //  if(model->param.svm_type == SVM_PLUS)  
      //        for(int i=0; i<model->l; i++)
      //          fprintf(stdout,"%f %f\n", model->sv_coef[0][i], model->sv_coef_star[0][i]);
      //      fflush(stdout);	
	
      free(label);
      free(probA);
      free(probB);
      free(count);
      free(perm);
      free(start);
      free(x);
      free(weighted_C);
      free(nonzero);
      for(i=0;i<nr_class*(nr_class-1)/2;i++)
	free(f[i].alpha);
      free(f);
      free(nz_count);
      free(nz_start);
      if(model->param.svm_type == SVM_PLUS) {
        free(x_star);
        for(i=0;i<nr_class*(nr_class-1)/2;i++)
	  free(f[i].beta);
        free(nz_count_star);
        free(nz_start_star);
      }

    }
  return model;
}

// Stratified cross validation
void svm_cross_validation(const svm_problem *prob, const svm_parameter *param, int nr_fold, double *target)
{
  int i;
  int *fold_start = Malloc(int,nr_fold+1);
  int l = prob->l;
  int *perm = Malloc(int,l);
  int nr_class;

  // stratified cv may not give leave-one-out rate
  // Each class to l folds -> some folds may have zero elements
  if((param->svm_type == C_SVC ||
      param->svm_type == NU_SVC) && nr_fold < l)
    {
      int *start = NULL;
      int *label = NULL;
      int *count = NULL;
      svm_group_classes(prob,&nr_class,&label,&start,&count,perm);

      // random shuffle and then data grouped by fold using the array perm
      int *fold_count = Malloc(int,nr_fold);
      int c;
      int *index = Malloc(int,l);
      for(i=0;i<l;i++)
	index[i]=perm[i];
      for (c=0; c<nr_class; c++) 
	for(i=0;i<count[c];i++)
	  {
	    int j = i+rand()%(count[c]-i);
	    swap(index[start[c]+j],index[start[c]+i]);
	  }
      for(i=0;i<nr_fold;i++)
	{
	  fold_count[i] = 0;
	  for (c=0; c<nr_class;c++)
	    fold_count[i]+=(i+1)*count[c]/nr_fold-i*count[c]/nr_fold;
	}
      fold_start[0]=0;
      for (i=1;i<=nr_fold;i++)
	fold_start[i] = fold_start[i-1]+fold_count[i-1];
      for (c=0; c<nr_class;c++)
	for(i=0;i<nr_fold;i++)
	  {
	    int begin = start[c]+i*count[c]/nr_fold;
	    int end = start[c]+(i+1)*count[c]/nr_fold;
	    for(int j=begin;j<end;j++)
	      {
		perm[fold_start[i]] = index[j];
		fold_start[i]++;
	      }
	  }
      fold_start[0]=0;
      for (i=1;i<=nr_fold;i++)
	fold_start[i] = fold_start[i-1]+fold_count[i-1];
      free(start);	
      free(label);
      free(count);	
      free(index);
      free(fold_count);
    }
  else
    {
      for(i=0;i<l;i++) perm[i]=i;
      for(i=0;i<l;i++)
	{
	  int j = i+rand()%(l-i);
	  swap(perm[i],perm[j]);
	}
      for(i=0;i<=nr_fold;i++)
	fold_start[i]=i*l/nr_fold;
    }

  for(i=0;i<nr_fold;i++)
    {
      int begin = fold_start[i];
      int end = fold_start[i+1];
      int j,k;
      struct svm_problem subprob;

      subprob.l = l-(end-begin);
      subprob.x = Malloc(struct svm_node*,subprob.l);
      subprob.y = Malloc(double,subprob.l);
			
      k=0;
      for(j=0;j<begin;j++)
	{
	  subprob.x[k] = prob->x[perm[j]];
	  subprob.y[k] = prob->y[perm[j]];
	  ++k;
	}
      for(j=end;j<l;j++)
	{
	  subprob.x[k] = prob->x[perm[j]];
	  subprob.y[k] = prob->y[perm[j]];
	  ++k;
	}
      struct svm_model *submodel = svm_train(&subprob,param);
      if(param->probability && 
	 (param->svm_type == C_SVC || param->svm_type == NU_SVC))
	{
	  double *prob_estimates=Malloc(double,svm_get_nr_class(submodel));
	  for(j=begin;j<end;j++)
	    target[perm[j]] = svm_predict_probability(submodel,prob->x[perm[j]],prob_estimates);
	  free(prob_estimates);			
	}
      else
	for(j=begin;j<end;j++)
	  target[perm[j]] = svm_predict(submodel,prob->x[perm[j]]);
      svm_destroy_model(submodel);
      free(subprob.x);
      free(subprob.y);
    }		
  free(fold_start);
  free(perm);	
}


int svm_get_svm_type(const svm_model *model)
{
  return model->param.svm_type;
}

int svm_get_nr_class(const svm_model *model)
{
  return model->nr_class;
}

void svm_get_labels(const svm_model *model, int* label)
{
  if (model->label != NULL)
    for(int i=0;i<model->nr_class;i++)
      label[i] = model->label[i];
}

double svm_get_svr_probability(const svm_model *model)
{
  if ((model->param.svm_type == EPSILON_SVR || model->param.svm_type == NU_SVR) &&
      model->probA!=NULL)
    return model->probA[0];
  else
    {
      fprintf(stderr,"Model doesn't contain information for SVR probability inference\n");
      return 0;
    }
}

void svm_predict_values(const svm_model *model, const svm_node *x, double* dec_values)
{
  if(model->param.svm_type == ONE_CLASS ||
     model->param.svm_type == EPSILON_SVR ||
     model->param.svm_type == NU_SVR)
    {
      double *sv_coef = model->sv_coef[0];
      double sum = 0;
      for(int i=0;i<model->l;i++)
	sum += sv_coef[i] * Kernel::k_function(x,model->SV[i],model->param);
      sum -= model->rho[0];
      *dec_values = sum;
    }
  else
    {
      int i;
      int nr_class = model->nr_class;
      int l = model->l;
		
      double *kvalue = Malloc(double,l);
      for(i=0;i<l;i++)
	kvalue[i] = Kernel::k_function(x,model->SV[i],model->param);

      int *start = Malloc(int,nr_class);
      start[0] = 0;
      for(i=1;i<nr_class;i++)
	start[i] = start[i-1]+model->nSV[i-1];

      int p=0;
      for(i=0;i<nr_class;i++)
	for(int j=i+1;j<nr_class;j++)
	  {
	    double sum = 0;
	    int si = start[i];
	    int sj = start[j];
	    int ci = model->nSV[i];
	    int cj = model->nSV[j];
				
	    int k;
	    double *coef1 = model->sv_coef[j-1];
	    double *coef2 = model->sv_coef[i];
	    for(k=0;k<ci;k++)
	      sum += coef1[si+k] * kvalue[si+k];
	    for(k=0;k<cj;k++)
	      sum += coef2[sj+k] * kvalue[sj+k];
	    sum -= model->rho[p];
	    dec_values[p] = sum;
	    p++;
	  }

      free(kvalue);
      free(start);
    }
}

double svm_predict(const svm_model *model, const svm_node *x)
{
  if(model->param.svm_type == ONE_CLASS ||
     model->param.svm_type == EPSILON_SVR ||
     model->param.svm_type == NU_SVR)
    {
      double res;
      svm_predict_values(model, x, &res);
		
      if(model->param.svm_type == ONE_CLASS)
	return (res>0)?1:-1;
      else
	return res;
    }
  else
    {
      int i;
      int nr_class = model->nr_class;
      double *dec_values = Malloc(double, nr_class*(nr_class-1)/2);
      svm_predict_values(model, x, dec_values);
      double soft_classification = dec_values[0]*model->label[0];

      int *vote = Malloc(int,nr_class);
      for(i=0;i<nr_class;i++)
	vote[i] = 0;
      int pos=0;
      for(i=0;i<nr_class;i++)
	for(int j=i+1;j<nr_class;j++)
	  {
	    if(dec_values[pos++] > 0)
	      ++vote[i];
	    else
	      ++vote[j];
	  }

      int vote_max_idx = 0;
      for(i=1;i<nr_class;i++)
	if(vote[i] > vote[vote_max_idx])
	  vote_max_idx = i;
      free(vote);
      free(dec_values);
      return soft_classification; // model->label[vote_max_idx];
    }
}

double svm_predict_probability(
			       const svm_model *model, const svm_node *x, double *prob_estimates)
{
  if ((model->param.svm_type == C_SVC || model->param.svm_type == NU_SVC) &&
      model->probA!=NULL && model->probB!=NULL)
    {
      int i;
      int nr_class = model->nr_class;
      double *dec_values = Malloc(double, nr_class*(nr_class-1)/2);
      svm_predict_values(model, x, dec_values);

      double min_prob=1e-7;
      double **pairwise_prob=Malloc(double *,nr_class);
      for(i=0;i<nr_class;i++)
	pairwise_prob[i]=Malloc(double,nr_class);
      int k=0;
      for(i=0;i<nr_class;i++)
	for(int j=i+1;j<nr_class;j++)
	  {
	    pairwise_prob[i][j]=min(max(sigmoid_predict(dec_values[k],model->probA[k],model->probB[k]),min_prob),1-min_prob);
	    pairwise_prob[j][i]=1-pairwise_prob[i][j];
	    k++;
	  }
      multiclass_probability(nr_class,pairwise_prob,prob_estimates);

      int prob_max_idx = 0;
      for(i=1;i<nr_class;i++)
	if(prob_estimates[i] > prob_estimates[prob_max_idx])
	  prob_max_idx = i;
      for(i=0;i<nr_class;i++)
	free(pairwise_prob[i]);
      free(dec_values);
      free(pairwise_prob);	     
      return model->label[prob_max_idx];
    }
  else 
    return svm_predict(model, x);
}

const char *svm_type_table[] =
  {
    "c_svc","nu_svc","one_class","epsilon_svr","nu_svr","svmp_plus",NULL
  };

const char *kernel_type_table[]=
  {
    "linear","polynomial","rbf","sigmoid","precomputed",NULL
  };

int svm_save_model(const char *model_file_name, const svm_model *model)
{
  FILE *fp = fopen(model_file_name,"w");
  fprintf(stdout,"model_file_name=%s", model_file_name);
  fflush(stdout);
  if(fp==NULL) return -1;

  const svm_parameter& param = model->param;

  fprintf(fp,"svm_type %s\n", svm_type_table[param.svm_type]);
  fprintf(fp,"kernel_type %s\n", kernel_type_table[param.kernel_type]);

  if(param.kernel_type == POLY)
    fprintf(fp,"degree %d\n", param.degree);

  if(param.kernel_type == POLY || param.kernel_type == RBF || param.kernel_type == SIGMOID)
    fprintf(fp,"gamma %g\n", param.gamma);

  if(param.kernel_type == POLY || param.kernel_type == SIGMOID)
    fprintf(fp,"coef0 %g\n", param.coef0);

  int nr_class = model->nr_class;
  int l = model->l;
  fprintf(fp, "nr_class %d\n", nr_class);
  fprintf(fp, "total_sv %d\n",l);
  fprintf(stdout, "total_sv=%d",l);
  fflush(stdout);
	
  {
    fprintf(fp, "rho");
    for(int i=0;i<nr_class*(nr_class-1)/2;i++)
      fprintf(fp," %g",model->rho[i]);
    fprintf(fp, "\n");
  }
	
  if(model->label)
    {
      fprintf(fp, "label");
      for(int i=0;i<nr_class;i++)
	fprintf(fp," %d",model->label[i]);
      fprintf(fp, "\n");
    }

  if(model->probA) // regression has probA only
    {
      fprintf(fp, "probA");
      for(int i=0;i<nr_class*(nr_class-1)/2;i++)
	fprintf(fp," %g",model->probA[i]);
      fprintf(fp, "\n");
    }
  if(model->probB)
    {
      fprintf(fp, "probB");
      for(int i=0;i<nr_class*(nr_class-1)/2;i++)
	fprintf(fp," %g",model->probB[i]);
      fprintf(fp, "\n");
    }

  if(model->nSV)
    {
      fprintf(fp, "nr_sv");
      for(int i=0;i<nr_class;i++)
	fprintf(fp," %d",model->nSV[i]);
      fprintf(fp, "\n");
    }

  fprintf(fp, "SV\n");
  const double * const *sv_coef = model->sv_coef;
  const svm_node * const *SV = model->SV;

  for(int i=0;i<l;i++)
    {
      for(int j=0;j<nr_class-1;j++) {
	fprintf(fp, "%.16g ",sv_coef[j][i]);
        fflush(stdout);
      }
      const svm_node *p = SV[i];

      if(param.kernel_type == PRECOMPUTED)
	fprintf(fp,"0:%d ",(int)(p->value));
      else
	while(p->index != -1)
	  {
	    fprintf(fp,"%d:%.8g ",p->index,p->value);
	    p++;
	  }
      fprintf(fp, "\n");
    }
  exit(1);

  if (ferror(fp) != 0 || fclose(fp) != 0) return -1;
  else return 0;
}

static char *line = NULL;
static int max_line_len;

static char* readline(FILE *input)
{
  int len;

  if(fgets(line,max_line_len,input) == NULL)
    return NULL;

  while(strrchr(line,'\n') == NULL)
    {
      max_line_len *= 2;
      line = (char *) realloc(line,max_line_len);
      len = (int) strlen(line);
      if(fgets(line+len,max_line_len-len,input) == NULL)
	break;
    }
  return line;
}

svm_model *svm_load_model(const char *model_file_name)
{
  FILE *fp = fopen(model_file_name,"rb");
  if(fp==NULL) return NULL;
	
  // read parameters

  svm_model *model = Malloc(svm_model,1);
  svm_parameter& param = model->param;
  model->rho = NULL;
  model->probA = NULL;
  model->probB = NULL;
  model->label = NULL;
  model->nSV = NULL;

  char cmd[81];
  while(1)
    {
      fscanf(fp,"%80s",cmd);

      if(strcmp(cmd,"svm_type")==0)
	{
	  fscanf(fp,"%80s",cmd);
	  int i;
	  for(i=0;svm_type_table[i];i++)
	    {
	      if(strcmp(svm_type_table[i],cmd)==0)
		{
		  param.svm_type=i;
		  break;
		}
	    }
	  if(svm_type_table[i] == NULL)
	    {
	      fprintf(stderr,"unknown svm type.\n");
	      free(model->rho);
	      free(model->label);
	      free(model->nSV);
	      free(model);
	      return NULL;
	    }
	}
      else if(strcmp(cmd,"kernel_type")==0)
	{		
	  fscanf(fp,"%80s",cmd);
	  int i;
	  for(i=0;kernel_type_table[i];i++)
	    {
	      if(strcmp(kernel_type_table[i],cmd)==0)
		{
		  param.kernel_type=i;
		  break;
		}
	    }
	  if(kernel_type_table[i] == NULL)
	    {
	      fprintf(stderr,"unknown kernel function.\n");
	      free(model->rho);
	      free(model->label);
	      free(model->nSV);
	      free(model);
	      return NULL;
	    }
	}
      else if(strcmp(cmd,"degree")==0)
	fscanf(fp,"%d",&param.degree);
      else if(strcmp(cmd,"gamma")==0)
	fscanf(fp,"%lf",&param.gamma);
      else if(strcmp(cmd,"coef0")==0)
	fscanf(fp,"%lf",&param.coef0);
      else if(strcmp(cmd,"nr_class")==0)
	fscanf(fp,"%d",&model->nr_class);
      else if(strcmp(cmd,"total_sv")==0)
	fscanf(fp,"%d",&model->l);
      else if(strcmp(cmd,"rho")==0)
	{
	  int n = model->nr_class * (model->nr_class-1)/2;
	  model->rho = Malloc(double,n);
	  for(int i=0;i<n;i++)
	    fscanf(fp,"%lf",&model->rho[i]);
	}
      else if(strcmp(cmd,"label")==0)
	{
	  int n = model->nr_class;
	  model->label = Malloc(int,n);
	  for(int i=0;i<n;i++)
	    fscanf(fp,"%d",&model->label[i]);
	}
      else if(strcmp(cmd,"probA")==0)
	{
	  int n = model->nr_class * (model->nr_class-1)/2;
	  model->probA = Malloc(double,n);
	  for(int i=0;i<n;i++)
	    fscanf(fp,"%lf",&model->probA[i]);
	}
      else if(strcmp(cmd,"probB")==0)
	{
	  int n = model->nr_class * (model->nr_class-1)/2;
	  model->probB = Malloc(double,n);
	  for(int i=0;i<n;i++)
	    fscanf(fp,"%lf",&model->probB[i]);
	}
      else if(strcmp(cmd,"nr_sv")==0)
	{
	  int n = model->nr_class;
	  model->nSV = Malloc(int,n);
	  for(int i=0;i<n;i++)
	    fscanf(fp,"%d",&model->nSV[i]);
	}
      else if(strcmp(cmd,"SV")==0)
	{
	  while(1)
	    {
	      int c = getc(fp);
	      if(c==EOF || c=='\n') break;	
	    }
	  break;
	}
      else
	{
	  fprintf(stderr,"unknown text in model file: [%s]\n",cmd);
	  free(model->rho);
	  free(model->label);
	  free(model->nSV);
	  free(model);
	  return NULL;
	}
    }

  // read sv_coef and SV

  int elements = 0;
  long pos = ftell(fp);

  max_line_len = 1024;
  line = Malloc(char,max_line_len);
  char *p,*endptr,*idx,*val;

  while(readline(fp)!=NULL)
    {
      p = strtok(line,":");
      while(1)
	{
	  p = strtok(NULL,":");
	  if(p == NULL)
	    break;
	  ++elements;
	}
    }
  elements += model->l;

  fseek(fp,pos,SEEK_SET);

  int m = model->nr_class - 1;
  int l = model->l;
  model->sv_coef = Malloc(double *,m);
  int i;
  for(i=0;i<m;i++)
    model->sv_coef[i] = Malloc(double,l);
  model->SV = Malloc(svm_node*,l);
  svm_node *x_space = NULL;
  if(l>0) x_space = Malloc(svm_node,elements);

  int j=0;
  for(i=0;i<l;i++)
    {
      readline(fp);
      model->SV[i] = &x_space[j];

      p = strtok(line, " \t");
      model->sv_coef[0][i] = strtod(p,&endptr);
      for(int k=1;k<m;k++)
	{
	  p = strtok(NULL, " \t");
	  model->sv_coef[k][i] = strtod(p,&endptr);
	}

      while(1)
	{
	  idx = strtok(NULL, ":");
	  val = strtok(NULL, " \t");

	  if(val == NULL)
	    break;
	  x_space[j].index = (int) strtol(idx,&endptr,10);
	  x_space[j].value = strtod(val,&endptr);

	  ++j;
	}
      x_space[j++].index = -1;
    }
  free(line);

  if (ferror(fp) != 0 || fclose(fp) != 0)
    return NULL;

  model->free_sv = 1;	// XXX
  return model;
}

void svm_destroy_model(svm_model* model)
{
  
  int i;
  if(model->free_sv && model->l > 0) {
    free((void *)(model->SV[0]));
    /*
      if(model->param.svm_type == SVM_PLUS)
      free((void *)(model->SV_star[0]));
    */
  }
  for(i=0;i<model->nr_class-1;i++) 
    free(model->sv_coef[i]);
  /*
    if(model->param.svm_type == SVM_PLUS) {
    for(i=0;i<model->nr_class-1;i++) 
    free(model->sv_coef[i]);
    free(model->rho_star);
    free(model->SV_star);
    free(model->sv_coef_star);
    free(model->nSV_star);
    }
  */
  free(model->SV);
  free(model->sv_coef);
  free(model->rho);
  free(model->label);
  free(model->probA);
  free(model->probB);
  free(model->nSV);
  free(model);
}

void svm_destroy_param(svm_parameter* param)
{
  free(param->weight_label);
  free(param->weight);
}

const char *svm_check_parameter(const svm_problem *prob, const svm_parameter *param)
{
  // svm_type

  int svm_type = param->svm_type;
  if(svm_type != C_SVC &&
     svm_type != NU_SVC &&
     svm_type != ONE_CLASS &&
     svm_type != EPSILON_SVR &&
     svm_type != NU_SVR &&
     svm_type != SVM_PLUS)
    return "unknown svm type";
	
  // kernel_type, degree
	
  int kernel_type = param->kernel_type;
  if(kernel_type != LINEAR &&
     kernel_type != POLY &&
     kernel_type != RBF &&
     kernel_type != SIGMOID &&
     kernel_type != PRECOMPUTED)
    return "unknown kernel type";

  int kernel_type_star = param->kernel_type_star;
  if(kernel_type_star != LINEAR &&
     kernel_type_star != POLY &&
     kernel_type_star != RBF &&
     kernel_type_star != SIGMOID &&
     kernel_type_star != PRECOMPUTED)
    return "unknown kernel type for the correcting space";

  if(param->degree < 0)
    return "degree of polynomial kernel < 0";

  if(param->degree_star < 0)
    return "degree of polynomial kernel for the correcting space < 0";

  // cache_size,eps,C,tau,nu,p,shrinking

  if(param->cache_size <= 0)
    return "cache_size <= 0";

  if(param->eps <= 0)
    return "eps <= 0";

  if(svm_type == C_SVC ||
     svm_type == SVM_PLUS ||
     svm_type == EPSILON_SVR ||
     svm_type == NU_SVR)
    if(param->C <= 0)
      return "C <= 0";

  if(svm_type == SVM_PLUS)
    if(param->tau <= 0)
      return "tau <= 0";

  if(svm_type == NU_SVC ||
     svm_type == ONE_CLASS ||
     svm_type == NU_SVR)
    if(param->nu <= 0 || param->nu > 1)
      return "nu <= 0 or nu > 1";

  if(svm_type == EPSILON_SVR)
    if(param->p < 0)
      return "p < 0";

  if(param->shrinking != 0 &&
     param->shrinking != 1)
    return "shrinking != 0 and shrinking != 1";

  if(param->probability != 0 &&
     param->probability != 1)
    return "probability != 0 and probability != 1";

  if(param->probability == 1 &&
     svm_type == ONE_CLASS)
    return "one-class SVM probability output not supported yet";


  // check whether nu-svc is feasible
	
  if(svm_type == NU_SVC)
    {
      int l = prob->l;
      int max_nr_class = 16;
      int nr_class = 0;
      int *label = Malloc(int,max_nr_class);
      int *count = Malloc(int,max_nr_class);

      int i;
      for(i=0;i<l;i++)
	{
	  int this_label = (int)prob->y[i];
	  int j;
	  for(j=0;j<nr_class;j++)
	    if(this_label == label[j])
	      {
		++count[j];
		break;
	      }
	  if(j == nr_class)
	    {
	      if(nr_class == max_nr_class)
		{
		  max_nr_class *= 2;
		  label = (int *)realloc(label,max_nr_class*sizeof(int));
		  count = (int *)realloc(count,max_nr_class*sizeof(int));
		}
	      label[nr_class] = this_label;
	      count[nr_class] = 1;
	      ++nr_class;
	    }
	}
	
      for(i=0;i<nr_class;i++)
	{
	  int n1 = count[i];
	  for(int j=i+1;j<nr_class;j++)
	    {
	      int n2 = count[j];
	      if(param->nu*(n1+n2)/2 > min(n1,n2))
		{
		  free(label);
		  free(count);
		  return "specified nu is infeasible";
		}
	    }
	}
      free(label);
      free(count);
    }

  return NULL;
}

int svm_check_probability_model(const svm_model *model)
{
  return ((model->param.svm_type == C_SVC || model->param.svm_type == NU_SVC) &&
	  model->probA!=NULL && model->probB!=NULL) ||
    ((model->param.svm_type == EPSILON_SVR || model->param.svm_type == NU_SVR) &&
     model->probA!=NULL);
}
