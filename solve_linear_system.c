
#include <math.h>
#include <stdio.h>
#include <memory.h>

#define TINY 1.0e-20

void lubksb(double **a, int n, int *indx, double b[]);
int ludcmp(double **a, int n, int *indx,  int *d);

int solve_linear_system(double **a, double *b, int n)
{
  int *indx,result;
  int d;

  indx = new int[n];
  result = ludcmp(a,n,indx,&d);

  if(result){
    delete[] indx;
    return 1;
  }
  lubksb(a,n,indx,b);

  delete[] indx;
  return 0;
}

void lubksb(double **a, int n, int *indx, double b[])
{
  int i,ii=-1,ip,j;
  double sum;

  for(i=0;i<n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    b[i] = sum;
    if(ii != -1)
      for(j=ii;j<=i-1;j++) 
	sum -= a[i][j]*b[j];
    else 
      if(sum) 
	ii=i;
    b[i]=sum;
  }
 
  for(i=n-1;i>=0;i--) {
    sum=b[i];
    for(j=i+1;j<n;j++) 
      sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

int ludcmp(double **a, int n, int *indx, int *d)
{
  int i,imax,j,k;
  double big,dum,sum,temp;
  double *vv;

  vv = new double[n];
  *d=1;
  for(i=0;i<n;i++) {
    big=0.0;
    for(j=0;j<n;j++)
      if((temp=fabs(a[i][j])) > big) 
        big=temp;
    if(big == 0.0) 
      return 1;
    vv[i]=1.0/big;
  }
  for(j=0;j<n;j++) {
    for(i=0;i<j;i++) {
      sum=a[i][j];
      for(k=0;k<i;k++) 
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for(i=j;i<n;i++) {
      sum=a[i][j];
      for(k=0;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if((dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if(j != imax) {
      for(k=0;k<n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;

    if(a[j][j] == 0.0) 
      a[j][j]=TINY;
    if(j != n-1) {
      dum=1.0/(a[j][j]);
      for(i=j+1;i<n;i++) 
        a[i][j] *= dum;
    }
  }
  
  delete[] vv;
  return 0;
}
#undef TINY

