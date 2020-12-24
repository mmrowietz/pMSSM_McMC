/*******************************
*    CalcHEP  3.8.7*
*******************************/
#include<math.h>
#include<complex.h>
#include"num_out.h"
#include"num_in.h"
double BWrange_a32a32_a10a10p1=2.7;
int twidth_a32a32_a10a10p1=0;
int gtwidth_a32a32_a10a10p1=0;
int gswidth_a32a32_a10a10p1=0;
 REAL va_a32a32_a10a10p1[35]={0};
const int nin_a32a32_a10a10p1 = 2;

const int nout_a32a32_a10a10p1 = 3;

const int nprc_a32a32_a10a10p1 = 1;

char * pinf_a32a32_a10a10p1(int nsub,int nprtcl,REAL* pmass,int * num)
{
int n;
 static char *names[1][5] ={
{"~L1","~L1","A","L","L"}};
int const nvalue[1][5]={
{8,8,0,3,3}};
int const pcode[1][5]={
{-1000015,-1000015,22,-15,-15}};
if  (nsub<0 ||nsub>1||nprtcl<0||nprtcl>5) return NULL;
if(pmass)
{
  n=nvalue[nsub-1][nprtcl-1];
  if (n==0) *pmass=0; else *pmass=va_a32a32_a10a10p1[n];
  if (*pmass<0) (*pmass)=-(*pmass);
}
if(num)*num=pcode[nsub-1][nprtcl-1];
return names[nsub-1][nprtcl-1];
}
char * polarized_a32a32_a10a10p1[3]={"",",",","};
int pinfAux_a32a32_a10a10p1(int nsub,int nprtcl,int*spin2,int*color,int*neutral,int*ndf)
{
int const pcode[1][5]={
{-1000015,-1000015,22,-15,-15}};
int const Spin2[1][5]={
{0,0,2,1,1}};
int const Color[1][5]={
{1,1,1,1,1}};
int const Neutral[1][5]={
{0,0,1,0,0}};
int const NDF[1][5]={
{1,1,2,2,2}};
if(nsub<0 ||nsub>1||nprtcl<0||nprtcl>5) return 0;
if(spin2) *spin2=Spin2[nsub-1][nprtcl-1];
if(color) *color=Color[nsub-1][nprtcl-1];
if(neutral) *neutral=Neutral[nsub-1][nprtcl-1];
if(ndf) *ndf=NDF[nsub-1][nprtcl-1];
return pcode[nsub-1][nprtcl-1];
}
const int nvar_a32a32_a10a10p1 = 26;

const int nfunc_a32a32_a10a10p1 = 8;

char * varName_a32a32_a10a10p1[35]={"P(cms)"
,"MW"
,"EE"
,"Ml"
,"CW"
,"SW"
,"S2W"
,"cb"
,"MSl1"
,"Zl11"
,"Zl12"
,"MNE1"
,"MNE2"
,"MNE3"
,"MNE4"
,"Zn11"
,"Zn12"
,"Zn13"
,"Zn21"
,"Zn22"
,"Zn23"
,"Zn31"
,"Zn32"
,"Zn33"
,"Zn41"
,"Zn42"
,"Zn43"
,"B00230"
,"B00231"
,"B00234"
,"B00235"
,"B00238"
,"B00239"
,"B00242"
,"B00243"};

 char * den_info_a32a32_a10a10p1(int nsub,int n, int * mass, int * width, int*pnum)
{
 switch(nsub){
 case 1: switch(n){
    case 1: *mass=8; *width=0;  if(pnum) *pnum=31; return "\1\3";
    case 2: *mass=14; *width=0;  if(pnum) *pnum=25; return "\2\4";
    case 3: *mass=14; *width=0;  if(pnum) *pnum=25; return "\2\5";
    case 4: *mass=14; *width=0;  if(pnum) *pnum=25; return "\1\5";
    case 5: *mass=14; *width=0;  if(pnum) *pnum=25; return "\1\4";
    case 6: *mass=13; *width=0;  if(pnum) *pnum=24; return "\2\4";
    case 7: *mass=13; *width=0;  if(pnum) *pnum=24; return "\2\5";
    case 8: *mass=13; *width=0;  if(pnum) *pnum=24; return "\1\5";
    case 9: *mass=13; *width=0;  if(pnum) *pnum=24; return "\1\4";
    case 10: *mass=12; *width=0;  if(pnum) *pnum=23; return "\2\4";
    case 11: *mass=12; *width=0;  if(pnum) *pnum=23; return "\2\5";
    case 12: *mass=12; *width=0;  if(pnum) *pnum=23; return "\1\5";
    case 13: *mass=12; *width=0;  if(pnum) *pnum=23; return "\1\4";
    case 14: *mass=11; *width=0;  if(pnum) *pnum=22; return "\2\4";
    case 15: *mass=11; *width=0;  if(pnum) *pnum=22; return "\2\5";
    case 16: *mass=3; *width=0;  if(pnum) *pnum=9; return "\3\4";
    case 17: *mass=3; *width=0;  if(pnum) *pnum=9; return "\3\5";
    case 18: *mass=11; *width=0;  if(pnum) *pnum=22; return "\1\5";
    case 19: *mass=8; *width=0;  if(pnum) *pnum=31; return "\2\3";
    case 20: *mass=11; *width=0;  if(pnum) *pnum=22; return "\1\4";
    default:*mass=0; *width=0; if(pnum) *pnum=0; return NULL;
                  }
   default: *mass=0; *width=0; return NULL;
            }
}

CalcHEP_interface interface_a32a32_a10a10p1={ 0,
"/nfs/dust/cms/user/mrowietm/python_scan/pMSSM_McMC/packages/micromegas_5.2.4/CalcHEP_src/"
,26, 8, varName_a32a32_a10a10p1,va_a32a32_a10a10p1,2, 3, 1, &pinf_a32a32_a10a10p1, &pinfAux_a32a32_a10a10p1, polarized_a32a32_a10a10p1, &calcFunc_a32a32_a10a10p1, &BWrange_a32a32_a10a10p1,&twidth_a32a32_a10a10p1,&gtwidth_a32a32_a10a10p1,&gswidth_a32a32_a10a10p1, &aWidth_a32a32_a10a10p1, &sqme_a32a32_a10a10p1,&den_info_a32a32_a10a10p1,cb_a32a32_a10a10p1};

CalcHEP_interface * PtrInterface_a32a32_a10a10p1=&interface_a32a32_a10a10p1;
