/*******************************
*    CalcHEP  3.8.7*
*******************************/
#include"num_out.h"
#include"num_in.h"
extern FNN F238_a32a32_a10a10p1;
static void C238(REAL*V,REAL * C)
{
REAL S[10];                                                                 
     
S[0]=V[34]*V[34];
S[1]=V[33]*V[33];
S[2]=V[30]*V[30];
S[3]=V[29]*V[29];
S[4]=V[8]*V[8];
S[5]=V[8]*V[8]*V[8]*V[8];
C[20]=+S[4]*(V[3]*(V[3]*(V[12]*(V[14]*(S[3]*(S[0]+S[1])+S[2]*(S[0]+S[1]))+
 V[3]*(V[33]*(V[34]*(2*(S[2]+S[3])))))+V[3]*(V[29]*(V[30]*(V[14]*(2*(S[0]+
 S[1]))+4*V[34]*V[33]*V[3]))))+S[4]*(V[29]*(V[14]*(V[30]*(S[0]+S[1]))-V[34]*
 V[33]*V[29]*V[12])-V[34]*V[33]*S[2]*V[12]))-2*V[34]*V[33]*V[30]*V[29]*
 S[5]);
C[19]=+V[3]*(V[3]*(V[29]*(V[33]*(V[3]*(V[30]*(4*V[34]*V[3]+2*V[33]*V[14])+2*
 V[34]*V[29]*V[12])+4*V[34]*V[30]*S[4]+V[33]*V[29]*V[14]*V[12])+S[0]*(V[14]*
 (2*V[30]*V[3]+V[29]*V[12])))+S[2]*(V[12]*(V[14]*(S[0]+S[1])+2*V[34]*V[33]*
 V[3])))+S[4]*(V[14]*(V[29]*(V[30]*(2*(S[0]+S[1]))))))-2*V[34]*V[33]*V[30]*
 V[29]*S[5];
C[18]=+V[3]*(V[29]*(V[30]*(V[14]*(S[0]+S[1])+4*V[34]*V[33]*V[3])+V[34]*
 V[33]*V[29]*V[12])+V[34]*V[33]*S[2]*V[12])+2*V[34]*V[33]*V[30]*V[29]*S[4];
C[17]=+V[3]*(V[3]*(V[29]*(V[33]*(V[3]*(V[30]*(4*V[34]*V[3]+2*V[33]*V[14])+2*
 V[34]*V[29]*V[12])+2*V[34]*V[30]*S[4]+V[33]*V[29]*V[14]*V[12])+S[0]*(V[14]*
 (2*V[30]*V[3]+V[29]*V[12])))+S[2]*(V[12]*(V[14]*(S[0]+S[1])+2*V[34]*V[33]*
 V[3])))+S[4]*(V[14]*(V[29]*(V[30]*(S[0]+S[1])))))-2*V[34]*V[33]*V[30]*
 V[29]*S[5];
C[16]=+V[3]*(V[29]*(V[30]*(V[14]*(2*(S[0]+S[1]))+6*V[34]*V[33]*V[3])+V[34]*
 V[33]*V[29]*V[12])+V[34]*V[33]*S[2]*V[12]);
C[15]=+S[4]*(V[33]*(V[34]*(V[3]*(V[12]*(S[2]+S[3])+2*V[30]*V[29]*V[3])+2*
 V[30]*V[29]*S[4])));
C[14]=+V[3]*(V[3]*(V[29]*(V[33]*(V[3]*(V[30]*(4*V[34]*V[3]+2*V[33]*V[14])+2*
 V[34]*V[29]*V[12])+4*V[34]*V[30]*S[4]+V[33]*V[29]*V[14]*V[12])+S[0]*(V[14]*
 (2*V[30]*V[3]+V[29]*V[12])))+S[2]*(V[12]*(V[14]*(S[0]+S[1])+2*V[34]*V[33]*
 V[3])))+S[4]*(V[29]*(V[14]*(V[30]*(S[0]+S[1]))+V[34]*V[33]*V[29]*V[12])+
 V[34]*V[33]*S[2]*V[12]));
C[13]=+V[3]*(V[29]*(V[30]*(V[14]*(S[0]+S[1])+8*V[34]*V[33]*V[3])+3*V[34]*
 V[33]*V[29]*V[12])+3*V[34]*V[33]*S[2]*V[12])+4*V[34]*V[33]*V[30]*V[29]*
 S[4];
C[12]=+V[3]*(V[29]*(V[30]*(V[14]*(S[0]+S[1])+4*V[34]*V[33]*V[3])+2*V[34]*
 V[33]*V[29]*V[12])+2*V[34]*V[33]*S[2]*V[12])+V[12]*(V[14]*(S[0]*S[2]+S[1]*
 S[3]));
C[11]=+V[3]*(V[29]*(V[30]*(V[14]*(S[0]+S[1])+2*V[34]*V[33]*V[3])));
C[10]=+V[33]*(V[34]*(V[3]*(V[12]*(2*(S[2]+S[3]))+4*V[30]*V[29]*V[3])+2*
 V[30]*V[29]*S[4]));
C[9]=+V[3]*(V[3]*(V[12]*(V[14]*(S[3]*(-S[0]-S[1])+S[2]*(-S[0]-S[1]))+V[3]*(
 V[33]*(V[34]*(2*(-S[2]-S[3])))))+V[3]*(V[29]*(V[30]*(V[14]*(2*(-S[0]-
 S[1]))-4*V[34]*V[33]*V[3]))))+S[4]*(V[29]*(V[14]*(V[30]*(3*(S[0]+S[1])))-
 V[34]*V[33]*V[29]*V[12])-V[34]*V[33]*S[2]*V[12]))+S[4]*(V[12]*(V[14]*(2*(
 S[0]*S[2]+S[1]*S[3])))-4*V[34]*V[33]*V[30]*V[29]*S[4]);
C[8]=+V[29]*(V[33]*(V[3]*(V[30]*(4*V[34]*V[3]+3*V[33]*V[14])+V[34]*V[29]*
 V[12])+2*V[33]*V[29]*V[14]*V[12]-4*V[34]*V[30]*S[4])+3*S[0]*V[30]*V[14]*
 V[3])+S[2]*(V[12]*(V[34]*(2*V[34]*V[14]+V[33]*V[3])));
C[7]=+V[3]*(V[29]*(V[30]*(V[14]*(4*(S[0]+S[1]))+6*V[34]*V[33]*V[3])+V[34]*
 V[33]*V[29]*V[12])+V[34]*V[33]*S[2]*V[12])+V[12]*(V[14]*(2*(S[0]*S[2]+S[1]*
 S[3])));
C[6]=+V[12]*(V[33]*(S[3]*(V[34]*V[3]+V[33]*V[14])+V[34]*S[2]*V[3])+S[0]*
 S[2]*V[14])-4*V[34]*V[33]*V[30]*V[29]*S[4];
C[5]=+V[3]*(2*(V[29]*(V[14]*(V[30]*(S[0]+S[1]))-V[34]*V[33]*V[29]*V[12])-
 V[34]*V[33]*S[2]*V[12]));
C[4]=+V[29]*(V[30]*(V[3]*(V[14]*(2*(S[0]+S[1]))+4*V[34]*V[33]*V[3])+6*V[34]*
 V[33]*S[4]));
C[3]=+2*V[34]*V[33]*V[30]*V[29];
C[2]=+4*V[34]*V[33]*V[30]*V[29];
S[6]=V[7]*V[7]*V[7]*V[7];
S[7]=V[6]*V[6]*V[6]*V[6];
S[8]=V[1]*V[1]*V[1]*V[1];
C[1]=+S[6]*S[7]*S[8];
S[9]=V[2]*V[2]*V[2]*V[2]*V[2]*V[2];
C[0]=+32*S[9];
}
REAL F238_a32a32_a10a10p1(double GG, REAL*V, REAL*DP,REAL*Q0,COMPLEX*Q1,REAL*
 Q2,REAL*cb_coeff,int gsw,int gtw)
{
REAL N,D,R; COMPLEX Prop;
static REAL C[21];                                                          
     
if(!DP){C238(V,C); return 0;} 
  REAL N_p1p2_=1/DP[0];
N=+C[0];
D=+C[1];
R=+DP[0]*(C[3]*(DP[0]*(DP[0]-DP[1]-DP[2])+DP[3]*(DP[1]+DP[3])+DP[4]*(-DP[2]-
 DP[4]))+DP[3]*(C[2]*(DP[2]-DP[0]+DP[4])-C[13])+C[19]+C[18]*DP[0]-C[16]*
 DP[1]-C[11]*DP[2]-C[8]*DP[4])+DP[3]*(DP[4]*(C[3]*(DP[2]-DP[1])+C[2]*(DP[4]-
 DP[3])+C[5])+DP[2]*(C[11]-C[3]*DP[3])+C[12]*DP[1]-C[14]+C[10]*DP[3])+DP[4]*
 (DP[1]*(C[7]+C[3]*DP[4])+C[6]*DP[2]-C[9]-C[4]*DP[4])+C[20]-C[17]*DP[1]+
 C[15]*DP[2];
R*=(N/D);
Prop=1*(gtw ? creal(Q1[16]):Q1[16])*(gsw ? creal(Q1[4]):Q1[4])*(gtw ? creal(
 Q1[20]):conj(Q1[20]))*(gtw ? creal(Q1[11]):conj(Q1[11]));
R*=creal(Prop);
 return R;
}
