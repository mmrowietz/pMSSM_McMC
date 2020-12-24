#ifndef __NUM_OUT_omg_p21a33
#define __NUM_OUT_omg_p21a33

#include<stdlib.h>
#include<string.h> 
#include<math.h>

#include"nType.h"

#define maxNp 20

extern  int    FError;

extern const int nin_omg_p21a33;
extern const int nout_omg_p21a33;
extern const int nprc_omg_p21a33;
extern const int nvar_omg_p21a33;
extern const int nfunc_omg_p21a33;

extern char * pinf_omg_p21a33(int nsub,int nprtcl,REAL* pmass,int*num);
extern int   pinfAux_omg_p21a33(int nsub, int nprtcl,int *spin2, int* color,int*neutral,int*ndf);
extern char * varName_omg_p21a33[];

extern double sqme_omg_p21a33(int nsub,double GG, REAL * momenta, REAL*cb_coeff,int * err);
extern int calcFunc_omg_p21a33(void);
extern double BWrange_omg_p21a33;
extern int twidth_omg_p21a33, gtwidth_omg_p21a33, gswidth_omg_p21a33;
extern double (*aWidth_omg_p21a33)(char *);
extern REAL va_omg_p21a33[];

extern  char * den_info_omg_p21a33(int nsub, int n, int * mass, int * width,int *pnum);


typedef  struct  { int pow; int nC; int * chains;} colorBasis;

extern colorBasis cb_omg_p21a33[];  

extern double (*aWidth_omg_p21a33)(char *);

#ifndef  __CALCHEP_INTERFACE__
#define  __CALCHEP_INTERFACE__
typedef struct CalcHEP_interface
{

  int forceUG;
  char * CALCHEP;

  int nvar;
  int nfunc;
  char ** varName;
  REAL * va;
  
  int nin;
  int nout;
  int nprc;
  char* (*pinf)(int, int , REAL*,int *);
  int  (*pinfAux)(int, int,int *,int*,int*,int*);
  char** polarized;
  int (*calcFunc)(void);
  double * BWrange;
  int    * twidth;    
  int *   gtwidth;
  int *   gswidth;
  double (**aWidth)(char *);

  double (*sqme)(int,double,REAL*,REAL*,int*);

  char * (*den_info)(int, int, int *, int*,int*);
  colorBasis *cb;  
} CalcHEP_interface;

extern int    jobInit(CalcHEP_interface * interface);
extern int    jobInState(int nProc, double P1, double P2, char* strf1, char*strf2);
extern int    jobCut2(char * par,double min,  double max);
extern int    jobCutMin(char * par,double min);
extern int    jobCutMax(char * par,double max);
extern void   jobCutDel(void);
extern int    jobHist(double min, char * par, double max);
extern void   jobHistDel(void);
extern double jobVegas(int nSess,int nCalls,int clear,int*err_,double*dI,double*chi2);
#endif

extern CalcHEP_interface interface_omg_p21a33;
extern CalcHEP_interface * PtrInterface_omg_p21a33;

extern void link_process(CalcHEP_interface * interface);

extern  int    OnlyTEQ0;

#define  DENOMINATOR_ERROR   2
#endif
