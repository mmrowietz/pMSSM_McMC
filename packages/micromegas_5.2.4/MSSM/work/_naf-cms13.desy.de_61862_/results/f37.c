/*******************************
*    CalcHEP  3.8.7*
*******************************/
#include"num_out.h"
#include"num_in.h"
extern FNN F37_a32a32_a10a10p1;
static void C37(REAL*V,REAL * C)
{
REAL S[10];                                                                 
     
S[0]=V[8]*V[8];
S[1]=V[30]*V[30];
S[2]=V[29]*V[29];
S[3]=V[28]*V[28];
C[14]=+S[0]*(V[3]*(V[3]*(V[27]*(V[29]*(V[3]*(V[28]*(4*V[30]*V[3]+2*V[29]*
 V[12])+2*V[30]*V[27]*V[11])+4*V[30]*V[28]*S[0]+V[29]*V[27]*V[12]*V[11])+
 S[1]*(V[12]*(2*V[28]*V[3]+V[27]*V[11])))+S[3]*(V[11]*(V[12]*(S[1]+S[2])+2*
 V[30]*V[29]*V[3])))+S[0]*(V[27]*(V[12]*(V[28]*(S[1]+S[2]))+V[30]*V[29]*
 V[27]*V[11])+V[30]*V[29]*S[3]*V[11])));
S[4]=V[27]*V[27];
C[13]=+V[3]*(V[3]*(V[27]*(V[29]*(V[3]*(V[28]*(4*V[30]*V[3]+2*V[29]*V[12])+2*
 V[30]*V[27]*V[11])+4*V[30]*V[28]*S[0]+V[29]*V[27]*V[12]*V[11])+S[1]*(V[12]*
 (2*V[28]*V[3]+V[27]*V[11])))+S[3]*(V[11]*(V[12]*(S[1]+S[2])+2*V[30]*V[29]*
 V[3])))+S[0]*(V[11]*(V[29]*(V[30]*(2*(S[3]+S[4]))))));
C[12]=+V[3]*(V[29]*(V[30]*(V[11]*(2*(S[3]+S[4]))+4*V[28]*V[27]*V[3])));
C[11]=+V[3]*(V[3]*(V[27]*(V[29]*(V[3]*(V[28]*(-4*V[30]*V[3]-2*V[29]*V[12])-
 2*V[30]*V[27]*V[11])+8*V[30]*V[28]*S[0]-V[29]*V[27]*V[12]*V[11])+S[1]*(
 V[12]*(-2*V[28]*V[3]-V[27]*V[11])))+S[3]*(V[11]*(V[12]*(-S[1]-S[2])-2*
 V[30]*V[29]*V[3])))+S[0]*(3*(V[27]*(V[12]*(V[28]*(S[1]+S[2]))+V[30]*V[29]*
 V[27]*V[11])+V[30]*V[29]*S[3]*V[11])))+S[0]*(V[11]*(V[12]*(2*(S[1]*S[3]+
 S[2]*S[4]))));
C[10]=+V[3]*(V[27]*(V[28]*(V[12]*(S[1]+S[2])+4*V[30]*V[29]*V[3])+3*V[30]*
 V[29]*V[27]*V[11])+3*V[30]*V[29]*S[3]*V[11])+V[11]*(V[12]*(2*(S[1]*S[3]+
 S[2]*S[4])));
S[5]=V[3]*V[3];
C[9]=+V[27]*(V[29]*(V[28]*(V[30]*(2*(-S[0]-S[5])))+V[29]*V[27]*V[12]*
 V[11]))+S[1]*S[3]*V[12]*V[11];
C[8]=+V[3]*(V[27]*(V[28]*(V[12]*(4*(S[1]+S[2]))+12*V[30]*V[29]*V[3])+4*
 V[30]*V[29]*V[27]*V[11])+4*V[30]*V[29]*S[3]*V[11])+V[11]*(V[12]*(2*(S[1]*
 S[3]+S[2]*S[4])));
C[7]=+V[3]*(V[27]*(V[28]*(V[12]*(2*(S[1]+S[2]))+8*V[30]*V[29]*V[3])+2*V[30]*
 V[29]*V[27]*V[11])+2*V[30]*V[29]*S[3]*V[11]);
C[6]=+V[3]*(V[3]*(V[27]*(V[29]*(V[3]*(V[28]*(4*V[30]*V[3]+2*V[29]*V[12])+2*
 V[30]*V[27]*V[11])+4*V[30]*V[28]*S[0]+V[29]*V[27]*V[12]*V[11])+S[1]*(V[12]*
 (2*V[28]*V[3]+V[27]*V[11])))+S[3]*(V[11]*(V[12]*(S[1]+S[2])+2*V[30]*V[29]*
 V[3])))+S[0]*(V[27]*(V[12]*(V[28]*(S[1]+S[2]))+V[30]*V[29]*V[27]*V[11])+
 V[30]*V[29]*S[3]*V[11]));
C[5]=+V[3]*(V[27]*(V[12]*(V[28]*(S[1]+S[2]))-V[30]*V[29]*V[27]*V[11])-V[30]*
 V[29]*S[3]*V[11]);
C[4]=+V[27]*(V[29]*(V[3]*(V[28]*(2*V[30]*V[3]+V[29]*V[12])+V[30]*V[27]*
 V[11])+V[29]*V[27]*V[12]*V[11]-2*V[30]*V[28]*S[0])+S[1]*V[28]*V[12]*V[3])+
 S[3]*(V[11]*(V[30]*(V[30]*V[12]+V[29]*V[3])));
C[3]=+V[3]*(V[27]*(V[28]*(V[12]*(2*(S[1]+S[2]))+8*V[30]*V[29]*V[3])+2*V[30]*
 V[29]*V[27]*V[11])+2*V[30]*V[29]*S[3]*V[11])+4*V[30]*V[29]*V[28]*V[27]*
 S[0];
C[2]=+4*V[30]*V[29]*V[28]*V[27];
S[6]=V[7]*V[7]*V[7]*V[7];
S[7]=V[6]*V[6]*V[6]*V[6];
S[8]=V[1]*V[1]*V[1]*V[1];
C[1]=+S[6]*S[7]*S[8];
S[9]=V[2]*V[2]*V[2]*V[2]*V[2]*V[2];
C[0]=+32*S[9];
}
REAL F37_a32a32_a10a10p1(double GG, REAL*V, REAL*DP,REAL*Q0,COMPLEX*Q1,REAL*
 Q2,REAL*cb_coeff,int gsw,int gtw)
{
REAL N,D,R; COMPLEX Prop;
static REAL C[15];                                                          
     
if(!DP){C37(V,C); return 0;} 
  REAL N_p1p2_=1/DP[0];
N=-C[0];
D=+C[1];
R=+DP[3]*(DP[3]*(C[2]*(DP[0]-DP[1]+DP[2]-DP[3]+DP[4])+C[7])+DP[0]*(C[10]-
 C[2]*DP[2])+C[11]-C[9]*DP[1]-C[8]*DP[2]-C[3]*DP[4])+DP[0]*(C[5]*(DP[0]-
 DP[1]-DP[4])+C[12]*DP[2]-C[13])+DP[4]*(C[6]+C[5]*DP[1]-C[4]*DP[2])+C[6]*
 DP[2]-C[14];
R*=(N/D);
Prop=1*(gtw ? creal(Q1[1]):Q1[1])*(gsw ? creal(Q1[4]):Q1[4])*(gtw ? creal(
 Q1[8]):conj(Q1[8]))*(gtw ? creal(Q1[2]):conj(Q1[2]));
R*=creal(Prop);
 return R;
}
