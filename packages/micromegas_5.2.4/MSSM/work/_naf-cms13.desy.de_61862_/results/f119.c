/*******************************
*    CalcHEP  3.8.7*
*******************************/
#include"num_out.h"
#include"num_in.h"
extern FNN F119_a32a32_a10a10p1;
static void C119(REAL*V,REAL * C)
{
REAL S[14];                                                                 
     
S[0]=V[8]*V[8];
S[1]=V[12]*V[12];
S[2]=V[3]*V[3];
S[3]=V[29]*V[29];
S[4]=V[30]*V[30]*V[30];
S[5]=V[30]*V[30]*V[30]*V[30];
S[6]=V[30]*V[30];
S[7]=V[8]*V[8]*V[8]*V[8];
S[8]=V[8]*V[8]*V[8]*V[8]*V[8]*V[8];
C[21]=+V[3]*(V[3]*(V[29]*(V[29]*(V[30]*(V[3]*(V[3]*(V[30]*(12*S[0]-4*S[1]-8*
 S[2])-8*V[29]*V[12]*V[3])+12*V[29]*V[12]*S[0])+S[0]*(V[30]*(4*(S[1]+
 S[0]))))+S[1]*(S[3]*(4*S[0]-2*S[2])))+S[4]*(V[3]*(V[12]*(12*S[0]-8*
 S[2]))))+S[5]*(S[1]*(4*S[0]-2*S[2])))+S[7]*(V[12]*(V[29]*(V[30]*(2*(S[6]+
 S[3]))))))+4*S[6]*S[3]*S[8];
S[9]=V[29]*V[29]*V[29];
C[20]=+V[3]*(V[29]*(V[30]*(V[3]*(V[29]*(V[30]*(12*(S[0]+S[2])+4*S[1])+12*
 V[29]*V[12]*V[3])+12*S[6]*V[12]*V[3])+S[0]*(V[12]*(4*(S[6]+S[3]))))+4*S[9]*
 S[1]*V[3])+4*S[5]*S[1]*V[3])+8*S[6]*S[3]*S[7];
C[19]=+V[29]*(V[30]*(V[3]*(V[12]*(2*(S[6]+S[3]))+8*V[30]*V[29]*V[3])+4*
 V[30]*V[29]*S[0]));
C[18]=+S[2]*(V[12]*(V[29]*(V[29]*(V[12]*(2*S[6]+S[3])+4*V[30]*V[29]*V[3])+4*
 S[4]*V[3])+S[5]*V[12])+4*S[6]*S[3]*S[2])+4*S[6]*S[3]*S[7];
C[17]=+4*S[6]*S[3]*S[0];
C[16]=+V[3]*(V[29]*(V[30]*(V[3]*(V[29]*(V[30]*(2*S[1]+8*S[0]+12*S[2])+8*
 V[29]*V[12]*V[3])+8*S[6]*V[12]*V[3])+S[0]*(V[12]*(2*(S[6]+S[3]))))+S[9]*
 S[1]*V[3])+S[5]*S[1]*V[3])+8*S[6]*S[3]*S[7];
C[15]=+V[29]*(V[30]*(V[3]*(V[12]*(2*(S[6]+S[3]))+8*V[30]*V[29]*V[3])+12*
 V[30]*V[29]*S[0]));
C[14]=+V[29]*(V[30]*(V[3]*(V[12]*(2*(S[6]+S[3]))+4*V[30]*V[29]*V[3])+4*
 V[30]*V[29]*S[0]));
C[13]=+V[29]*(V[30]*(V[3]*(V[3]*(V[3]*(V[12]*(4*(-S[6]-S[3]))-8*V[30]*V[29]*
 V[3])+28*V[30]*V[29]*S[0])+S[0]*(V[12]*(16*(S[6]+S[3]))))+4*V[30]*V[29]*
 S[7])+4*S[9]*S[1]*S[0])+4*S[5]*S[1]*S[0];
C[12]=+V[29]*(V[29]*(V[30]*(V[3]*(28*V[30]*V[3]+16*V[29]*V[12])+8*V[30]*
 S[0])+4*S[3]*S[1])+16*S[4]*V[12]*V[3])+4*S[5]*S[1];
C[11]=+V[29]*(V[29]*(V[30]*(4*(V[3]*(V[30]*V[3]+V[29]*V[12])+V[30]*S[0]))+2*
 S[3]*S[1])+4*S[4]*V[12]*V[3])+2*S[5]*S[1];
C[10]=+V[29]*(V[29]*(V[30]*(V[3]*(24*V[30]*V[3]+10*V[29]*V[12])+8*V[30]*
 S[0])+2*S[3]*S[1])+10*S[4]*V[12]*V[3])+2*S[5]*S[1];
C[9]=+V[29]*(V[30]*(V[3]*(V[3]*(V[3]*(V[12]*(4*(S[6]+S[3]))+8*V[30]*V[29]*
 V[3])+12*V[30]*V[29]*S[0])+S[0]*(V[12]*(8*(S[6]+S[3]))))+4*V[30]*V[29]*
 S[7])+4*S[9]*S[1]*S[0])+4*S[5]*S[1]*S[0];
C[8]=+V[29]*(V[29]*(V[30]*(V[3]*(12*V[30]*V[3]+8*V[29]*V[12])+8*V[30]*S[0])+
 4*S[3]*S[1])+8*S[4]*V[12]*V[3])+4*S[5]*S[1];
C[7]=+V[12]*(2*(V[29]*(S[3]*(V[30]*V[3]+V[29]*V[12])+S[4]*V[3])+S[5]*
 V[12]))+4*S[6]*S[3]*S[0];
C[6]=+V[29]*(V[29]*(V[30]*(V[3]*(12*V[30]*V[3]+8*V[29]*V[12])+8*V[30]*S[0])+
 2*S[3]*S[1])+8*S[4]*V[12]*V[3])+2*S[5]*S[1];
C[5]=+8*S[6]*S[3];
C[4]=+V[12]*(V[29]*(S[3]*(12*V[30]*V[3]+4*V[29]*V[12])+12*S[4]*V[3])+4*S[5]*
 V[12])+24*S[6]*S[3]*S[2];
C[3]=+V[29]*(V[29]*(V[30]*(V[3]*(8*V[30]*V[3]+6*V[29]*V[12])+4*V[30]*S[0])+
 2*S[3]*S[1])+6*S[4]*V[12]*V[3])+2*S[5]*S[1];
C[2]=+4*S[6]*S[3];
S[10]=V[7]*V[7]*V[7]*V[7];
S[11]=V[6]*V[6]*V[6]*V[6];
S[12]=V[1]*V[1]*V[1]*V[1];
C[1]=+S[10]*S[11]*S[12];
S[13]=V[2]*V[2]*V[2]*V[2]*V[2]*V[2];
C[0]=+32*S[13];
}
REAL F119_a32a32_a10a10p1(double GG, REAL*V, REAL*DP,REAL*Q0,COMPLEX*Q1,REAL*
 Q2,REAL*cb_coeff,int gsw,int gtw)
{
REAL N,D,R; COMPLEX Prop;
static REAL C[22];REAL S[3];                                                
     
if(!DP){C119(V,C); return 0;} 
  REAL N_p1p2_=1/DP[0];
N=-C[0];
D=+C[1];
S[0]=DP[2]*DP[2];
S[1]=DP[3]*DP[3];
S[2]=DP[4]*DP[4];
R=+DP[0]*(C[2]*(DP[0]*(DP[2]+DP[3]+DP[4])+DP[1]*(-DP[2]-DP[3]-DP[4])-S[0]-
 S[1]-S[2])+DP[2]*(C[5]*(-DP[3]-DP[4])+C[15])+C[17]*DP[1]-C[20]-C[19]*DP[0]+
 C[12]*DP[3]+C[8]*DP[4])+DP[2]*(DP[3]*(C[2]*(DP[2]+DP[3]+DP[4])+C[5]*DP[1]-
 C[10])+DP[1]*(C[2]*DP[4]-C[14])+C[16]-C[14]*DP[2]-C[6]*DP[4])+DP[4]*(DP[1]*
 (C[2]*(DP[3]+DP[4])-C[7])+C[9]-C[4]*DP[3]-C[3]*DP[4])+DP[3]*(C[13]-C[11]*
 DP[1]-C[3]*DP[3])+C[18]*DP[1]-C[21];
R*=(N/D);
Prop=1*(gtw ? creal(Q1[8]):Q1[8])*(gsw ? creal(Q1[4]):Q1[4])*(gtw ? creal(
 Q1[10]):conj(Q1[10]))*(gsw ? creal(Q1[5]):conj(Q1[5]));
R*=creal(Prop);
 return R;
}
