/*******************************
*    CalcHEP  3.8.7*
*******************************/
#include"num_out.h"
#include"num_in.h"
extern FNN F176_a32a32_a10a10p1;
static void C176(REAL*V,REAL * C)
{
REAL S[10];                                                                 
     
S[0]=V[8]*V[8];
S[1]=V[32]*V[32];
S[2]=V[31]*V[31];
S[3]=V[28]*V[28];
S[4]=V[27]*V[27];
S[5]=V[3]*V[3];
C[17]=+S[5]*(V[27]*(V[31]*(V[3]*(V[28]*(-8*V[32]*V[3]-4*V[31]*V[13])-4*
 V[32]*V[27]*V[11])-4*V[32]*V[28]*S[0]-2*V[31]*V[27]*V[13]*V[11])+S[1]*(
 V[13]*(-4*V[28]*V[3]-2*V[27]*V[11])))+S[3]*(V[11]*(V[13]*(2*(-S[1]-S[2]))-
 4*V[32]*V[31]*V[3])))+S[0]*(V[11]*(V[13]*(2*(S[1]*S[3]+S[2]*S[4]))));
C[16]=+V[27]*(V[31]*(V[28]*(V[32]*(4*(-S[0]-S[5])))+2*V[31]*V[27]*V[13]*
 V[11]))+2*S[1]*S[3]*V[13]*V[11];
C[15]=+V[27]*(V[28]*(2*(V[3]*(V[13]*(S[1]+S[2])-V[32]*V[31]*V[3])-V[32]*
 V[31]*S[0]))+3*S[2]*V[27]*V[13]*V[11])+3*S[1]*S[3]*V[13]*V[11];
C[14]=+V[13]*(V[27]*(S[2]*(V[28]*V[3]+V[27]*V[11])+S[1]*V[28]*V[3])+S[1]*
 S[3]*V[11]);
C[13]=+S[5]*(V[27]*(V[31]*(V[3]*(V[28]*(-4*V[32]*V[3]-2*V[31]*V[13])-2*
 V[32]*V[27]*V[11])-2*V[32]*V[28]*S[0]-V[31]*V[27]*V[13]*V[11])+S[1]*(V[13]*
 (-2*V[28]*V[3]-V[27]*V[11])))+S[3]*(V[11]*(V[13]*(-S[1]-S[2])-2*V[32]*
 V[31]*V[3])))+S[0]*(V[11]*(V[13]*(S[1]*S[3]+S[2]*S[4])));
C[12]=+V[31]*(V[27]*(V[32]*(2*(V[3]*(V[27]*V[11]-V[28]*V[3])-V[28]*S[0]))+3*
 V[31]*V[27]*V[13]*V[11])+2*V[32]*S[3]*V[11]*V[3])+3*S[1]*S[3]*V[13]*V[11];
C[11]=+V[3]*(V[27]*(V[13]*(V[28]*(S[1]+S[2]))+V[32]*V[31]*V[27]*V[11])+
 V[32]*V[31]*S[3]*V[11])+V[11]*(V[13]*(2*(S[1]*S[3]+S[2]*S[4])));
C[10]=+V[11]*(V[31]*(S[4]*(V[32]*V[3]+V[31]*V[13])+V[32]*S[3]*V[3])+S[1]*
 S[3]*V[13]);
C[9]=+V[31]*(V[32]*(V[3]*(V[11]*(4*(S[3]+S[4]))+8*V[28]*V[27]*V[3])+4*V[28]*
 V[27]*S[0]));
C[8]=+V[31]*(V[32]*(V[3]*(V[11]*(2*(S[3]+S[4]))+4*V[28]*V[27]*V[3])+2*V[28]*
 V[27]*S[0]));
C[7]=+V[27]*(V[28]*(V[3]*(V[13]*(4*(S[1]+S[2]))+8*V[32]*V[31]*V[3])+4*V[32]*
 V[31]*S[0]));
C[6]=+6*V[32]*V[31]*V[28]*V[27];
C[5]=+V[27]*(V[28]*(V[3]*(V[13]*(2*(S[1]+S[2]))+4*V[32]*V[31]*V[3])+2*V[32]*
 V[31]*S[0]));
C[4]=+2*V[32]*V[31]*V[28]*V[27];
C[3]=+8*V[32]*V[31]*V[28]*V[27];
C[2]=+4*V[32]*V[31]*V[28]*V[27];
S[6]=V[7]*V[7]*V[7]*V[7];
S[7]=V[6]*V[6]*V[6]*V[6];
S[8]=V[1]*V[1]*V[1]*V[1];
C[1]=+S[6]*S[7]*S[8];
S[9]=V[2]*V[2]*V[2]*V[2]*V[2]*V[2];
C[0]=+32*S[9];
}
REAL F176_a32a32_a10a10p1(double GG, REAL*V, REAL*DP,REAL*Q0,COMPLEX*Q1,REAL*
 Q2,REAL*cb_coeff,int gsw,int gtw)
{
REAL N,D,R; COMPLEX Prop;
static REAL C[18];                                                          
     
if(!DP){C176(V,C); return 0;} 
  REAL N_p1p2_=1/DP[0];
N=-C[0];
D=+C[1];
R=+DP[0]*(DP[0]*(C[2]*(DP[3]-DP[0]+DP[4])+C[6]*(DP[1]+DP[2])+C[16])+DP[1]*(
 C[4]*(-DP[1]-DP[3])-C[15]-C[2]*DP[2]-C[6]*DP[4])+DP[2]*(C[4]*(-DP[2]-
 DP[4])-C[12]-C[6]*DP[3])+DP[3]*(C[9]-C[3]*DP[4])+C[17]+C[7]*DP[4])+DP[1]*(
 DP[4]*(C[4]*(DP[1]+DP[2])+C[2]*DP[3]-C[5])+DP[2]*(C[11]+C[4]*DP[3])+C[14]*
 DP[1]-C[13]-C[8]*DP[3])+DP[2]*(DP[3]*(C[4]*DP[2]-C[8]+C[2]*DP[4])+C[10]*
 DP[2]-C[13]-C[5]*DP[4]);
R*=(N/D);
Prop=1*(gtw ? creal(Q1[12]):Q1[12])*(gtw ? creal(Q1[2]):Q1[2])*(gtw ? creal(
 Q1[20]):conj(Q1[20]))*(gtw ? creal(Q1[7]):conj(Q1[7]));
R*=creal(Prop);
 return R;
}
