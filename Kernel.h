
#ifndef _KERNEL_H_
#define _KERNEL_H_

#include "common.h"
#include "Cache.h"

//
// Kernel evaluation
//
// the static method k_function is for doing single kernel evaluation
// the constructor of Kernel prepares to calculate the l*l kernel matrix
// the member function get_Q is for getting one column from the Q Matrix
//
class QMatrix {
public:
  virtual Qfloat *get_Q(int column, int len) const = 0;
  virtual Qfloat *get_QD() const = 0;
  virtual void swap_index(int i, int j) const = 0;
  virtual ~QMatrix() {}
};

class Kernel: public QMatrix {
public:
  Kernel(int l, svm_node * const * x, const svm_parameter& param);
  virtual ~Kernel();

  static double k_function(const svm_node *x, const svm_node *y,
			   const svm_parameter& param);
  virtual Qfloat *get_Q(int column, int len) const = 0;
  virtual Qfloat *get_QD() const = 0;
  virtual void swap_index(int i, int j) const	// no so const...
  {
    swap(x[i],x[j]);
    if(x_square) swap(x_square[i],x_square[j]);
  }
protected:

  double (Kernel::*kernel_function)(int i, int j) const;

private:
  const svm_node **x;
  double *x_square;

  // svm_parameter
  const int kernel_type;
  const int degree;
  const double gamma;
  const double coef0;

  static double dot(const svm_node *px, const svm_node *py);
  double kernel_linear(int i, int j) const
  {
    return dot(x[i],x[j]);
  }
  double kernel_poly(int i, int j) const
  {
    return powi(gamma*dot(x[i],x[j])+coef0,degree);
  }
  double kernel_rbf(int i, int j) const
  {
    return exp(-gamma*(x_square[i]+x_square[j]-2*dot(x[i],x[j])));
  }
  double kernel_sigmoid(int i, int j) const
  {
    return tanh(gamma*dot(x[i],x[j])+coef0);
  }
  double kernel_precomputed(int i, int j) const
  {
    return x[i][(int)(x[j][0].value)].value;
  }
};

//
// Q matrices for various formulations
//
class SVC_Q: public Kernel
{ 
public:
  SVC_Q(const svm_problem& prob, const svm_parameter& param, const schar *y_)
    :Kernel(prob.l, prob.x, param)
  {
    clone(y,y_,prob.l);
    cache = new Cache(prob.l,(long int)(param.cache_size*(1<<20)));
    QD = new Qfloat[prob.l];
    for(int i=0;i<prob.l;i++)
      QD[i]= (Qfloat)(this->*kernel_function)(i,i);
  }
	
  Qfloat *get_Q(int i, int len) const
  {
    Qfloat *data;
    int start, j;
    if((start = cache->get_data(i,&data,len)) < len)
      {
	for(j=start;j<len;j++)
	  data[j] = (Qfloat)(y[i]*y[j]*(this->*kernel_function)(i,j));
      }
    return data;
  }

  Qfloat *get_QD() const
  {
    return QD;
  }

  void swap_index(int i, int j) const
  {
    cache->swap_index(i,j);
    Kernel::swap_index(i,j);
    swap(y[i],y[j]);
    swap(QD[i],QD[j]);
  }

  ~SVC_Q()
  {
    delete[] y;
    delete cache;
    delete[] QD;
  }
private:
  schar *y;
  Cache *cache;
  Qfloat *QD;
};

class ONE_CLASS_Q: public Kernel
{
public:
  ONE_CLASS_Q(const svm_problem& prob, const svm_parameter& param)
    :Kernel(prob.l, prob.x, param)
  {
    cache = new Cache(prob.l,(long int)(param.cache_size*(1<<20)));
    QD = new Qfloat[prob.l];
    for(int i=0;i<prob.l;i++)
      QD[i]= (Qfloat)(this->*kernel_function)(i,i);
  }
	
  Qfloat *get_Q(int i, int len) const
  {
    Qfloat *data;
    int start, j;
    if((start = cache->get_data(i,&data,len)) < len)
      {
	for(j=start;j<len;j++)
	  data[j] = (Qfloat)(this->*kernel_function)(i,j);
      }
    return data;
  }

  Qfloat *get_QD() const
  {
    return QD;
  }

  void swap_index(int i, int j) const
  {
    cache->swap_index(i,j);
    Kernel::swap_index(i,j);
    swap(QD[i],QD[j]);
  }

  ~ONE_CLASS_Q()
  {
    delete cache;
    delete[] QD;
  }
private:
  Cache *cache;
  Qfloat *QD;
};

class SVR_Q: public Kernel
{ 
public:
  SVR_Q(const svm_problem& prob, const svm_parameter& param)
    :Kernel(prob.l, prob.x, param)
  {
    l = prob.l;
    cache = new Cache(l,(long int)(param.cache_size*(1<<20)));
    QD = new Qfloat[2*l];
    sign = new schar[2*l];
    index = new int[2*l];
    for(int k=0;k<l;k++)
      {
	sign[k] = 1;
	sign[k+l] = -1;
	index[k] = k;
	index[k+l] = k;
	QD[k]= (Qfloat)(this->*kernel_function)(k,k);
	QD[k+l]=QD[k];
      }
    buffer[0] = new Qfloat[2*l];
    buffer[1] = new Qfloat[2*l];
    next_buffer = 0;
  }

  void swap_index(int i, int j) const
  {
    swap(sign[i],sign[j]);
    swap(index[i],index[j]);
    swap(QD[i],QD[j]);
  }
	
  Qfloat *get_Q(int i, int len) const
  {
    Qfloat *data;
    int j, real_i = index[i];
    if(cache->get_data(real_i,&data,l) < l)
      {
	for(j=0;j<l;j++)
	  data[j] = (Qfloat)(this->*kernel_function)(real_i,j);
      }

    // reorder and copy
    Qfloat *buf = buffer[next_buffer];
    next_buffer = 1 - next_buffer;
    schar si = sign[i];
    for(j=0;j<len;j++)
      buf[j] = (Qfloat) si * (Qfloat) sign[j] * data[index[j]];
    return buf;
  }

  Qfloat *get_QD() const
  {
    return QD;
  }

  ~SVR_Q()
  {
    delete cache;
    delete[] sign;
    delete[] index;
    delete[] buffer[0];
    delete[] buffer[1];
    delete[] QD;
  }
private:
  int l;
  Cache *cache;
  schar *sign;
  int *index;
  mutable int next_buffer;
  Qfloat *buffer[2];
  Qfloat *QD;
};

#endif
