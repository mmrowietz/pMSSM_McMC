/*******************************
*    CalcHEP  3.8.7*
*******************************/
#include"num_out.h"
#include"num_in.h"
extern FNN F203_a32a32_a10a10p1;
static void C203(REAL*V,REAL * C)
{
REAL S[11];                                                                 
     
S[0]=V[8]*V[8];
S[1]=V[34]*V[34];
S[2]=V[33]*V[33];
S[3]=V[32]*V[32];
S[4]=V[8]*V[8]*V[8]*V[8];
C[20]=+S[0]*(V[3]*(V[3]*(V[31]*(V[33]*(V[3]*(V[32]*(4*V[34]*V[3]+2*V[33]*
 V[14])+2*V[34]*V[31]*V[13])+V[33]*V[31]*V[14]*V[13]-4*V[34]*V[32]*S[0])+
 S[1]*(V[14]*(2*V[32]*V[3]+V[31]*V[13])))+S[3]*(V[13]*(V[14]*(S[1]+S[2])+2*
 V[34]*V[33]*V[3])))+S[0]*(V[31]*(V[14]*(V[32]*(-S[1]-S[2]))-V[34]*V[33]*
 V[31]*V[13])-V[34]*V[33]*S[3]*V[13]))+2*V[34]*V[33]*V[32]*V[31]*S[4]);
S[5]=V[31]*V[31];
C[19]=+V[3]*(V[3]*(V[31]*(V[33]*(V[3]*(V[32]*(4*V[34]*V[3]+2*V[33]*V[14])+2*
 V[34]*V[31]*V[13])+V[33]*V[31]*V[14]*V[13]-4*V[34]*V[32]*S[0])+S[1]*(V[14]*
 (2*V[32]*V[3]+V[31]*V[13])))+S[3]*(V[13]*(V[14]*(S[1]+S[2])+2*V[34]*V[33]*
 V[3])))+S[0]*(V[13]*(V[33]*(V[34]*(2*(-S[3]-S[5]))))))+2*V[34]*V[33]*V[32]*
 V[31]*S[4];
C[18]=+V[3]*(V[31]*(V[14]*(V[32]*(S[1]+S[2]))-V[34]*V[33]*V[31]*V[13])-
 V[34]*V[33]*S[3]*V[13])-2*V[34]*V[33]*V[32]*V[31]*S[0];
C[17]=+S[0]*(V[31]*(V[32]*(V[3]*(V[14]*(S[1]+S[2])+2*V[34]*V[33]*V[3])-2*
 V[34]*V[33]*S[0])));
C[16]=+V[3]*(V[3]*(V[31]*(V[33]*(V[3]*(V[32]*(4*V[34]*V[3]+2*V[33]*V[14])+2*
 V[34]*V[31]*V[13])+V[33]*V[31]*V[14]*V[13]-2*V[34]*V[32]*S[0])+S[1]*(V[14]*
 (2*V[32]*V[3]+V[31]*V[13])))+S[3]*(V[13]*(V[14]*(S[1]+S[2])+2*V[34]*V[33]*
 V[3])))+S[0]*(V[13]*(V[33]*(V[34]*(-S[3]-S[5])))))+2*V[34]*V[33]*V[32]*
 V[31]*S[4];
C[15]=+V[3]*(V[31]*(V[32]*(V[14]*(S[1]+S[2])-2*V[34]*V[33]*V[3])-2*V[34]*
 V[33]*V[31]*V[13])-2*V[34]*V[33]*S[3]*V[13]);
C[14]=+V[3]*(V[3]*(V[31]*(V[33]*(V[3]*(V[32]*(-4*V[34]*V[3]-2*V[33]*V[14])-
 2*V[34]*V[31]*V[13])-4*V[34]*V[32]*S[0]-V[33]*V[31]*V[14]*V[13])+S[1]*(
 V[14]*(-2*V[32]*V[3]-V[31]*V[13])))+S[3]*(V[13]*(V[14]*(-S[1]-S[2])-2*
 V[34]*V[33]*V[3])))+S[0]*(V[31]*(V[14]*(V[32]*(-S[1]-S[2]))+V[34]*V[33]*
 V[31]*V[13])+V[34]*V[33]*S[3]*V[13]))+S[0]*(V[13]*(V[14]*(2*(S[1]*S[3]+
 S[2]*S[5])))+4*V[34]*V[33]*V[32]*V[31]*S[0]);
C[13]=+V[31]*(V[33]*(V[3]*(V[34]*V[31]*V[13]+V[33]*V[32]*V[14])+4*V[34]*
 V[32]*S[0]+2*V[33]*V[31]*V[14]*V[13])+S[1]*V[32]*V[14]*V[3])+S[3]*(V[13]*(
 V[34]*(2*V[34]*V[14]+V[33]*V[3])));
C[12]=+V[31]*(V[33]*(V[3]*(V[32]*(4*V[34]*V[3]+V[33]*V[14])+2*V[34]*V[31]*
 V[13])+4*V[34]*V[32]*S[0]+V[33]*V[31]*V[14]*V[13])+S[1]*V[32]*V[14]*V[3])+
 S[3]*(V[13]*(V[34]*(V[34]*V[14]+2*V[33]*V[3])));
S[6]=V[3]*V[3];
C[11]=+V[14]*(V[31]*(S[2]*(V[32]*V[3]+2*V[31]*V[13])+S[1]*V[32]*V[3])+2*
 S[1]*S[3]*V[13])-2*V[34]*V[33]*V[32]*V[31]*S[6];
C[10]=+V[33]*(V[34]*(V[3]*(V[13]*(2*(S[3]+S[5]))+4*V[32]*V[31]*V[3])+6*
 V[32]*V[31]*S[0]));
C[9]=+V[3]*(V[3]*(V[13]*(V[14]*(S[5]*(S[1]+S[2])+S[3]*(S[1]+S[2]))+V[3]*(
 V[33]*(V[34]*(2*(S[3]+S[5])))))+V[3]*(V[31]*(V[32]*(V[14]*(2*(S[1]+S[2]))+
 4*V[34]*V[33]*V[3]))))+S[0]*(V[31]*(V[14]*(V[32]*(S[1]+S[2]))-V[34]*V[33]*
 V[31]*V[13])-V[34]*V[33]*S[3]*V[13]));
C[8]=+V[3]*(V[31]*(V[32]*(V[14]*(3*(S[1]+S[2]))+4*V[34]*V[33]*V[3])-V[34]*
 V[33]*V[31]*V[13])-V[34]*V[33]*S[3]*V[13])-4*V[34]*V[33]*V[32]*V[31]*S[0];
C[7]=+V[3]*(V[33]*(V[34]*(V[13]*(S[3]+S[5])+2*V[32]*V[31]*V[3])));
C[6]=+V[3]*(V[31]*(V[32]*(V[14]*(2*(S[1]+S[2]))+4*V[34]*V[33]*V[3])+V[34]*
 V[33]*V[31]*V[13])+V[34]*V[33]*S[3]*V[13])+V[13]*(V[14]*(S[1]*S[3]+S[2]*
 S[5]));
C[5]=+V[3]*(V[31]*(V[32]*(V[14]*(2*(S[1]+S[2]))+8*V[34]*V[33]*V[3])+2*V[34]*
 V[33]*V[31]*V[13])+2*V[34]*V[33]*S[3]*V[13]);
C[4]=+V[31]*(V[32]*(V[3]*(V[14]*(2*(S[1]+S[2]))+4*V[34]*V[33]*V[3])-2*V[34]*
 V[33]*S[0]));
C[3]=+2*V[34]*V[33]*V[32]*V[31];
C[2]=+4*V[34]*V[33]*V[32]*V[31];
S[7]=V[7]*V[7]*V[7]*V[7];
S[8]=V[6]*V[6]*V[6]*V[6];
S[9]=V[1]*V[1]*V[1]*V[1];
C[1]=+S[7]*S[8]*S[9];
S[10]=V[2]*V[2]*V[2]*V[2]*V[2]*V[2];
C[0]=+32*S[10];
}
REAL F203_a32a32_a10a10p1(double GG, REAL*V, REAL*DP,REAL*Q0,COMPLEX*Q1,REAL*
 Q2,REAL*cb_coeff,int gsw,int gtw)
{
REAL N,D,R; COMPLEX Prop;
static REAL C[21];                                                          
     
if(!DP){C203(V,C); return 0;} 
  REAL N_p1p2_=1/DP[0];
N=+C[0];
D=+C[1];
R=+DP[0]*(C[3]*(DP[0]*(DP[1]-DP[0]+DP[2])+DP[3]*(DP[1]+DP[3])+DP[4]*(-DP[2]-
 DP[4]))+DP[4]*(C[2]*(DP[0]-DP[1]-DP[3])-C[8])+C[19]+C[18]*DP[0]+C[7]*DP[1]-
 C[15]*DP[2]-C[13]*DP[3])+DP[3]*(DP[4]*(C[3]*(DP[2]-DP[1])+C[2]*(DP[4]-
 DP[3])-C[5])+DP[2]*(C[11]-C[3]*DP[3])+C[12]*DP[1]-C[14]+C[10]*DP[3])+DP[4]*
 (DP[1]*(C[3]*DP[4]-C[7])+C[6]*DP[2]-C[9]+C[4]*DP[4])+C[20]+C[17]*DP[1]-
 C[16]*DP[2];
R*=(N/D);
Prop=1*(gtw ? creal(Q1[15]):Q1[15])*(gsw ? creal(Q1[4]):Q1[4])*(gtw ? creal(
 Q1[17]):conj(Q1[17]))*(gtw ? creal(Q1[2]):conj(Q1[2]));
R*=creal(Prop);
 return R;
}
