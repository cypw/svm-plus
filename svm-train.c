#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include "common.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

void print_null(const char *s) {}

void exit_with_help()
{
  printf(
	 "Usage: svm-train [options] training_set_file [model_file]\n"
	 "options:\n"
	 "-s svm_type : set type of SVM (default 0)\n"
	 "	0 -- C-SVC\n"
	 "	1 -- nu-SVC\n"
	 "	2 -- one-class SVM\n"
	 "	3 -- epsilon-SVR\n"
	 "	4 -- nu-SVR\n"
	 "	5 -- SVM+\n"
         "-a n : optimization method \n"
         "    -1  -- Max Unconstrained Gain SMO (default)\n"
         "     0  -- Max Constrained Gain SMO (Glassmachers&Igel, JMLR2006)\n"
         "    k>0 -- Conjugate SMO of order k\n"
	 "-t kernel_type : set type of kernel function (default 2)\n"
	 "	0 -- linear: u'*v\n"
	 "	1 -- polynomial: (gamma*u'*v + coef0)^degree\n"
	 "	2 -- radial basis function: exp(-gamma*|u-v|^2)\n"
	 "	3 -- sigmoid: tanh(gamma*u'*v + coef0)\n"
	 "	4 -- precomputed kernel (kernel values in training_set_file)\n"
	 "-T kernel_type_star : set type of kernel function for the correcting space (default 2), for SVM+\n"
	 "	0 -- linear: u'*v\n"
	 "	1 -- polynomial: (gamma*u'*v + coef0)^degree\n"
	 "	2 -- radial basis function: exp(-gamma*|u-v|^2)\n"
	 "	3 -- sigmoid: tanh(gamma*u'*v + coef0)\n"
	 "	4 -- precomputed kernel (kernel values in training_set_file)\n"
	 "-f star_file : name of the file containing star examples. Necessary parameter for SVM+ \n"
	 "-d degree : set degree in kernel function (default 3)\n"
	 "-D degree_star : set degree_star in kernel function in the correcting space (default 3)\n" 
	 "-g gamma : set gamma in kernel function (default 1/number of features)\n"
	 "-G gamma_star : set gamma_star in kernel function in the correcting space (default 1/number of features in the  correcting space)\n"
	 "-r coef0 : set coef0 in kernel function (default 0)\n"
	 "-R coef0_star : set coef0_star in kernel function (default 0)\n"
	 "-c cost : set the parameter C of C-SVC, epsilon-SVR, nu-SVR and SVM+ (default 1)\n"
	 "-C tau : set the parameter tau in SVM+ (default 1)\n"
	 "-n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)\n"
	 "-p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)\n"
	 "-m cachesize : set cache memory size in MB (default 100)\n"
	 "-e epsilon : set tolerance of termination criterion (default 0.001)\n"
	 "-h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)\n"
	 "-b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)\n"
	 "-wi weight : set the parameter C of class i to weight*C, for C-SVC and SVM+ (default 1)\n"
	 "-v n: n-fold cross validation mode\n"
	 "-q : quiet mode (no outputs)\n"
	 );
  exit(1);
}

void exit_input_error(int line_num)
{
  fprintf(stderr,"Wrong input format at line %d\n", line_num);
  exit(1);
}

void parse_command_line(int argc, char **argv, char *input_file_name, char *input_file_name_star, char *model_file_name);
void read_problem(const char *filename, struct svm_problem *prob, struct svm_node **x_space, int *max_index);
void do_cross_validation();
void check_kernel_input(struct svm_problem prob, int max_index);
void check_compatibility(struct svm_problem prob, struct svm_problem prob_star);

struct svm_parameter param;		// set by parse_command_line
struct svm_problem prob;		// set by read_problem
struct svm_problem prob_star;		// set by read_problem
struct svm_model *model;
struct svm_node *x_space;
struct svm_node *x_space_star;
int cross_validation;
int nr_fold;

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

int main(int argc, char **argv)
{
  char input_file_name[1024], input_file_name_star[1024];
  char model_file_name[1024];
  const char *error_msg;
  int max_index;

  parse_command_line(argc, argv, input_file_name, input_file_name_star, model_file_name);
  read_problem(input_file_name,&prob, &x_space, &max_index);
  if(param.gamma == 0 && max_index > 0)
    param.gamma = 1.0/max_index;        
  if(param.kernel_type == PRECOMPUTED)
    check_kernel_input(prob, max_index);

  if(param.svm_type == SVM_PLUS) {          
    read_problem(input_file_name_star,&prob_star, &x_space_star, &max_index);
    if(param.gamma_star == 0 && max_index > 0)
      param.gamma_star = 1.0/max_index; 
    if(param.kernel_type_star == PRECOMPUTED)
      check_kernel_input(prob_star, max_index);
    check_compatibility(prob,prob_star);
    prob.x_star = prob_star.x;   /* merge two problems into a single one */
  }

  error_msg = svm_check_parameter(&prob,&param);

  if(error_msg) {
      fprintf(stderr,"Error: %s\n",error_msg);
      exit(1);
  }

  if(cross_validation) {
      do_cross_validation();
  }
  else {
      model = svm_train(&prob,&param);
      svm_save_model(model_file_name,model);
      svm_destroy_model(model);
  }
  svm_destroy_param(&param);
  free(prob.y);
  free(prob.x);
  free(prob.x_star);
  free(x_space);
  free(x_space_star);
  free(line);

  return 0;
}

void do_cross_validation()
{
  int i;
  int total_correct = 0;
  double total_error = 0;
  double sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;
  double *target = Malloc(double,prob.l);

  svm_cross_validation(&prob,&param,nr_fold,target);
  if(param.svm_type == EPSILON_SVR ||
     param.svm_type == NU_SVR)
    {
      for(i=0;i<prob.l;i++)
	{
	  double y = prob.y[i];
	  double v = target[i];
	  total_error += (v-y)*(v-y);
	  sumv += v;
	  sumy += y;
	  sumvv += v*v;
	  sumyy += y*y;
	  sumvy += v*y;
	}
      printf("Cross Validation Mean squared error = %g\n",total_error/prob.l);
      printf("Cross Validation Squared correlation coefficient = %g\n",
	     ((prob.l*sumvy-sumv*sumy)*(prob.l*sumvy-sumv*sumy))/
	     ((prob.l*sumvv-sumv*sumv)*(prob.l*sumyy-sumy*sumy))
	     );
    }
  else
    {
      for(i=0;i<prob.l;i++)
	if(target[i] == prob.y[i])
	  ++total_correct;
      printf("Cross Validation Accuracy = %g%%\n",100.0*total_correct/prob.l);
    }
  free(target);
}

void parse_command_line(int argc, char **argv, char *input_file_name, char *input_file_name_star, char *model_file_name)
{
  int i;

  // default values
  param.svm_type = C_SVC;
  param.kernel_type = RBF;
  param.kernel_type_star = RBF;
  param.degree = 3;
  param.degree_star = 3;
  param.gamma = 0;	// 1/k
  param.gamma_star = 0; // 1/k
  param.coef0 = 0;
  param.coef0_star = 0;
  param.nu = 0.5;
  param.cache_size = 100;
  param.C = 1;
  param.tau = 0;
  param.eps = 1e-3;
  param.p = 0.1;
  param.shrinking = 1;
  param.probability = 0;
  param.nr_weight = 0;
  param.weight_label = NULL;
  param.weight = NULL;
  param.optimizer = -1;
  cross_validation = 0;
  input_file_name_star[0] = '\0';

  // parse options
  for(i=1;i<argc;i++)
    {
      if(argv[i][0] != '-') break;
      if(++i>=argc)
	exit_with_help();
      switch(argv[i-1][1])
	{
	case 'a':
	  param.optimizer = atoi(argv[i]);
	  break;
	case 's':
	  param.svm_type = atoi(argv[i]);
	  break;
	case 't':
	  param.kernel_type = atoi(argv[i]);
	  break;
	case 'T':
	  param.kernel_type_star = atoi(argv[i]);
	  break;
	case 'd':
	  param.degree = atoi(argv[i]);
	  break;
	case 'D':
	  param.degree_star = atoi(argv[i]);
	  break;
	case 'g':
	  param.gamma = atof(argv[i]);
	  break;
	case 'G':
	  param.gamma_star = atof(argv[i]);
	  break;
	case 'r':
	  param.coef0 = atof(argv[i]);
	  break;
	case 'R':
	  param.coef0_star = atof(argv[i]);
	  break;
	case 'n':
	  param.nu = atof(argv[i]);
	  break;
	case 'm':
	  param.cache_size = atof(argv[i]);
	  break;
	case 'c':
	  param.C = atof(argv[i]);
	  break;
	case 'C':
	  param.tau = atof(argv[i]);
	  break;
	case 'e':
	  param.eps = atof(argv[i]);
	  break;
	case 'p':
	  param.p = atof(argv[i]);
	  break;
	case 'h':
	  param.shrinking = atoi(argv[i]);
	  break;
	case 'b':
	  param.probability = atoi(argv[i]);
	  break;
	case 'q':
	  svm_print_string = &print_null;
	  i--;
	  break;
	case 'v':
	  cross_validation = 1;
	  nr_fold = atoi(argv[i]);
	  if(nr_fold < 2)
	    {
	      fprintf(stderr,"n-fold cross validation: n must >= 2\n");
	      exit_with_help();
	    }
	  break;
	case 'w':
	  ++param.nr_weight;
	  param.weight_label = (int *)realloc(param.weight_label,sizeof(int)*param.nr_weight);
	  param.weight = (double *)realloc(param.weight,sizeof(double)*param.nr_weight);
	  param.weight_label[param.nr_weight-1] = atoi(&argv[i-1][2]);
	  param.weight[param.nr_weight-1] = atof(argv[i]);
	  break;
	case 'f': 
	  strcpy(input_file_name_star, argv[i]);
	  break;
	default:
	  fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
	  exit_with_help();
	}
    }

  // determine filenames
  if((param.svm_type == SVM_PLUS || param.svm_type == SVM_PLUS) && input_file_name_star[0] == '\0')
    exit_with_help();

  if(i>=argc)
    exit_with_help();

  strcpy(input_file_name, argv[i]);

  if(i<argc-1)
    strcpy(model_file_name,argv[i+1]);
  else
    {
      char *p = strrchr(argv[i],'/');
      if(p==NULL)
	p = argv[i];
      else
	++p;
      sprintf(model_file_name,"%s.model",p);
    }
}

// read in a problem (in svmlight format)

void read_problem(const char *filename, struct svm_problem *prob, struct svm_node **x_space, int *prob_max_index)
{
  int elements, max_index, inst_max_index, i, j;
  FILE *fp = fopen(filename,"r");
  char *endptr;
  char *idx, *val, *label;

  if(fp == NULL)
    {
      fprintf(stderr,"can't open input file %s\n",filename);
      exit(1);
    }

  prob->l = 0;
  elements = 0;

  max_line_len = 1024;
  line = Malloc(char,max_line_len);
  while(readline(fp)!=NULL)
    {
      char *p = strtok(line," \t"); // label

      // features
      while(1)
	{
	  p = strtok(NULL," \t");
	  if(p == NULL || *p == '\n') // check '\n' as ' ' may be after the last feature
	    break;
	  ++elements;
	}
      ++elements;
      ++(prob->l);
    }
  rewind(fp);

  prob->y = Malloc(double,prob->l);
  prob->x = Malloc(struct svm_node *,prob->l);
  *x_space = Malloc(struct svm_node,elements);

  max_index = 0;
  j=0;
  for(i=0;i<prob->l;i++)
    {
      inst_max_index = -1; // strtol gives 0 if wrong format, and precomputed kernel has <index> start from 0
      readline(fp);
      prob->x[i] = &((*x_space)[j]);
      label = strtok(line," \t");
      prob->y[i] = strtod(label,&endptr);
      if(endptr == label)
	exit_input_error(i+1);

      while(1)
	{
	  idx = strtok(NULL,":");
	  val = strtok(NULL," \t");

	  if(val == NULL)
	    break;

	  errno = 0;
	  (*x_space)[j].index = (int) strtol(idx,&endptr,10);
	  if(endptr == idx || errno != 0 || *endptr != '\0' || (*x_space)[j].index <= inst_max_index)
	    exit_input_error(i+1);
	  else
	    inst_max_index = (*x_space)[j].index;

	  errno = 0;
	  (*x_space)[j].value = strtod(val,&endptr);
	  if(endptr == val || errno != 0 || (*endptr != '\0' && !isspace(*endptr)))
	    exit_input_error(i+1);

	  ++j;
	}

      if(inst_max_index > max_index)
	max_index = inst_max_index;
      (*x_space)[j++].index = -1;
    }
  *prob_max_index = max_index;
        
  fclose(fp);     
}

void check_kernel_input(struct svm_problem prob, int max_index)
{
  int i;
  for(i=0;i<prob.l;i++)
    {
      if (prob.x[i][0].index != 0)
	{
	  fprintf(stderr,"Wrong input format: first column must be 0:sample_serial_number\n");
	  exit(1);
	}
      if ((int)prob.x[i][0].value <= 0 || (int)prob.x[i][0].value > max_index)
	{
	  fprintf(stderr,"Wrong input format: sample_serial_number out of range\n");
	  exit(1);
	}
    }
}

void check_compatibility(struct svm_problem prob, struct svm_problem prob_star)
{
  int i;

  if(prob.l != prob_star.l) {
    fprintf(stderr,"Different number of examples in X and X* space\n");
    exit(1);
  }

  for(i = 0; i < prob.l; i++)
    if(prob.y[i] != prob_star.y[i]) {
      fprintf(stderr,"Different labels in example %d\n",i);
      exit(1);
    }
}

