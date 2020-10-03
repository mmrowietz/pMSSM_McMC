#include "intfile.hh"

dcmplx Pf5(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double x3=x[3];
dcmplx y[133];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=esx[0];
y[5]=2.*x0*y[1]*y[2];
y[6]=-(x0*y[1]*y[4]);
y[7]=2.*y[1]*y[2];
y[8]=2.*x1*y[1]*y[2];
y[9]=-(y[1]*y[4]);
y[10]=2.*x2*y[1]*y[2];
y[11]=2.*x3*y[1]*y[2];
y[12]=-(x2*y[1]*y[4]);
y[13]=-x1;
y[14]=1.+y[13];
y[15]=2.*x0*x2*y[1]*y[2];
y[16]=2.*x0*x3*y[1]*y[2];
y[17]=-x0;
y[18]=1.+y[17];
y[19]=x2*x2;
y[20]=x3*x3;
y[21]=lrs[0];
y[22]=2.*x1*x2*y[1]*y[2];
y[23]=y[1]*y[2]*y[19];
y[24]=2.*x1*x3*y[1]*y[2];
y[25]=2.*x2*x3*y[1]*y[2];
y[26]=y[1]*y[2]*y[20];
y[27]=-(x1*x2*y[1]*y[4]);
y[28]=-(x2*x3*y[1]*y[4]);
y[29]=2.*x0*y[1]*y[2]*y[19];
y[30]=x1*y[1]*y[2]*y[19];
y[31]=4.*x0*x2*x3*y[1]*y[2];
y[32]=2.*x1*x2*x3*y[1]*y[2];
y[33]=2.*x0*y[1]*y[2]*y[20];
y[34]=x1*y[1]*y[2]*y[20];
y[35]=-(x3*y[1]*y[4]);
y[36]=-2.*x0*x2*x3*y[1]*y[4];
y[37]=-(x1*x2*x3*y[1]*y[4]);
y[38]=y[10]+y[11]+y[22]+y[23]+y[24]+y[25]+y[26]+y[27]+y[28]+y[29]+y[30]+y[31\
]+y[32]+y[33]+y[34]+y[35]+y[36]+y[37];
y[39]=lrs[1];
y[40]=-x2;
y[41]=1.+y[40];
y[42]=4.*x0*x2*y[1]*y[2];
y[43]=4.*x0*x3*y[1]*y[2];
y[44]=lambda*lambda;
y[45]=-(x0*x2*y[1]*y[4]);
y[46]=y[5]+y[7]+y[8]+y[9]+y[15]+y[16]+y[45];
y[47]=y[10]+y[11]+y[12]+y[23]+y[25]+y[26]+y[28];
y[48]=-2.*x0*x2*y[1]*y[4];
y[49]=y[7]+y[8]+y[9]+y[10]+y[11]+y[12]+y[22]+y[24]+y[27]+y[42]+y[43]+y[48];
y[50]=x0*y[1]*y[2]*y[19];
y[51]=2.*x0*x2*x3*y[1]*y[2];
y[52]=x0*y[1]*y[2]*y[20];
y[53]=-(x0*x2*x3*y[1]*y[4]);
y[54]=y[7]+y[8]+y[9]+y[10]+y[11]+y[12]+y[15]+y[16]+y[22]+y[24]+y[35]+y[45]+y\
[50]+y[51]+y[52]+y[53];
y[55]=lrs[2];
y[56]=2.*y[1]*y[2]*y[3];
y[57]=2.*x0*x1*y[1]*y[2];
y[58]=-(y[1]*y[3]*y[4]);
y[59]=-(x0*x1*y[1]*y[4]);
y[60]=y[5]+y[6]+y[56]+y[57]+y[58]+y[59];
y[61]=2.*y[1]*y[2]*y[19];
y[62]=4.*x2*x3*y[1]*y[2];
y[63]=2.*y[1]*y[2]*y[20];
y[64]=-2.*x2*x3*y[1]*y[4];
y[65]=y[61]+y[62]+y[63]+y[64];
y[66]=-(lambda*MYI*x0*y[18]*y[21]*y[65]);
y[67]=-(lambda*MYI*y[18]*y[21]*y[38]);
y[68]=lambda*MYI*x0*y[21]*y[38];
y[69]=1.+y[66]+y[67]+y[68];
y[70]=y[7]+y[10]+y[11];
y[71]=-(lambda*MYI*x1*y[14]*y[39]*y[70]);
y[72]=-(lambda*MYI*y[14]*y[39]*y[54]);
y[73]=lambda*MYI*x1*y[39]*y[54];
y[74]=1.+y[71]+y[72]+y[73];
y[75]=-x3;
y[76]=1.+y[75];
y[77]=-(x1*y[1]*y[4]);
y[78]=-2.*x0*x3*y[1]*y[4];
y[79]=-(x1*x3*y[1]*y[4]);
y[80]=y[7]+y[8]+y[10]+y[11]+y[22]+y[24]+y[35]+y[42]+y[43]+y[77]+y[78]+y[79];
y[81]=-(x0*x3*y[1]*y[4]);
y[82]=y[5]+y[6]+y[7]+y[8]+y[9]+y[15]+y[16]+y[81];
y[83]=x0*x1*y[14]*y[18]*y[21]*y[39]*y[44]*y[47]*y[49];
y[84]=-(lambda*MYI*x1*y[14]*y[39]*y[46]*y[69]);
y[85]=y[83]+y[84];
y[86]=y[1]*y[2];
y[87]=x1*x1;
y[88]=y[1]*y[2]*y[87];
y[89]=2.*x2*y[1]*y[2]*y[3];
y[90]=2.*x0*x1*x2*y[1]*y[2];
y[91]=2.*x3*y[1]*y[2]*y[3];
y[92]=2.*x0*x1*x3*y[1]*y[2];
y[93]=-(x3*y[1]*y[3]*y[4]);
y[94]=-(x0*x1*x3*y[1]*y[4]);
y[95]=y[5]+y[8]+y[15]+y[16]+y[57]+y[59]+y[77]+y[81]+y[86]+y[88]+y[89]+y[90]+\
y[91]+y[92]+y[93]+y[94];
y[96]=lrs[3];
y[97]=x0*x1*y[14]*y[18]*y[21]*y[39]*y[44]*y[49]*y[82];
y[98]=-(x0*x1*y[14]*y[18]*y[21]*y[39]*y[44]*y[46]*y[80]);
y[99]=y[97]+y[98];
y[100]=-(x0*x1*y[14]*y[18]*y[21]*y[39]*y[44]*y[46]*y[47]);
y[101]=lambda*MYI*x0*y[18]*y[21]*y[49]*y[74];
y[102]=y[100]+y[101];
y[103]=y[5]+y[56]+y[57];
y[104]=-(lambda*MYI*x2*y[41]*y[55]*y[103]);
y[105]=-(lambda*MYI*y[41]*y[55]*y[95]);
y[106]=lambda*MYI*x2*y[55]*y[95];
y[107]=1.+y[104]+y[105]+y[106];
y[108]=x0*x1*y[14]*y[18]*y[21]*y[39]*y[44]*y[47]*y[80];
y[109]=-(lambda*MYI*x1*y[14]*y[39]*y[69]*y[82]);
y[110]=y[108]+y[109];
y[111]=-(x0*x1*y[14]*y[18]*y[21]*y[39]*y[44]*y[47]*y[82]);
y[112]=lambda*MYI*x0*y[18]*y[21]*y[74]*y[80];
y[113]=y[111]+y[112];
y[114]=pow(y[47],2);
y[115]=x0*x1*y[14]*y[18]*y[21]*y[39]*y[44]*y[114];
y[116]=y[69]*y[74];
y[117]=y[115]+y[116];
y[118]=-(x2*y[1]*y[3]*y[4]);
y[119]=-(x0*x1*x2*y[1]*y[4]);
y[120]=y[5]+y[6]+y[8]+y[15]+y[16]+y[45]+y[57]+y[77]+y[86]+y[88]+y[89]+y[90]+\
y[91]+y[92]+y[118]+y[119];
y[121]=-(lambda*MYI*x2*y[41]*y[55]*y[95]);
y[122]=x2+y[121];
y[123]=-(lambda*MYI*x1*y[14]*y[39]*y[54]);
y[124]=x1+y[123];
y[125]=-(lambda*MYI*x0*y[18]*y[21]*y[38]);
y[126]=x0+y[125];
y[127]=-(lambda*MYI*x3*y[76]*y[96]*y[120]);
y[128]=x3+y[127];
y[129]=pow(y[124],2);
y[130]=pow(y[122],2);
y[131]=pow(y[126],2);
y[132]=pow(y[128],2);
FOUT=(pow(bi,-2)*(-(lambda*MYI*x3*y[46]*y[76]*y[96]*(-(lambda*MYI*x2*y[41]*y\
[55]*y[80]*y[99])-y[85]*y[107]-lambda*MYI*x2*y[41]*y[55]*y[60]*y[110]))+lam\
bda*MYI*x3*y[49]*y[76]*y[96]*(-(lambda*MYI*x2*y[41]*y[55]*y[82]*y[99])-y[10\
2]*y[107]-lambda*MYI*x2*y[41]*y[55]*y[60]*y[113])+lambda*MYI*x3*y[60]*y[76]\
*y[96]*(lambda*MYI*x2*y[41]*y[55]*y[82]*y[85]-lambda*MYI*x2*y[41]*y[55]*y[8\
0]*y[102]-lambda*MYI*x2*y[41]*y[55]*y[60]*y[117])+(lambda*MYI*x2*y[41]*y[55\
]*y[82]*y[110]-lambda*MYI*x2*y[41]*y[55]*y[80]*y[113]+y[107]*y[117])*(1.-la\
mbda*MYI*x3*y[76]*y[96]*y[103]+lambda*MYI*x3*y[96]*y[120]-lambda*MYI*y[76]*\
y[96]*y[120])))/((y[1]+y[1]*y[122]+y[1]*y[124]+y[1]*y[122]*y[124]+y[1]*y[12\
2]*y[126]+y[1]*y[128]+y[1]*y[124]*y[128]+y[1]*y[126]*y[128])*(y[86]+y[1]*y[\
2]*y[122]+2.*y[1]*y[2]*y[124]-y[1]*y[4]*y[124]+2.*y[1]*y[2]*y[122]*y[124]-y\
[1]*y[4]*y[122]*y[124]+2.*y[1]*y[2]*y[122]*y[126]+2.*y[1]*y[2]*y[122]*y[124\
]*y[126]-y[1]*y[4]*y[122]*y[124]*y[126]+y[1]*y[2]*y[128]+2.*y[1]*y[2]*y[124\
]*y[128]-y[1]*y[4]*y[124]*y[128]+2.*y[1]*y[2]*y[126]*y[128]-y[1]*y[4]*y[126\
]*y[128]+2.*y[1]*y[2]*y[122]*y[126]*y[128]-y[1]*y[4]*y[122]*y[126]*y[128]+2\
.*y[1]*y[2]*y[124]*y[126]*y[128]+2.*y[1]*y[2]*y[122]*y[124]*y[126]*y[128]-y\
[1]*y[4]*y[122]*y[124]*y[126]*y[128]+y[1]*y[2]*y[129]+y[1]*y[2]*y[122]*y[12\
9]+y[1]*y[2]*y[128]*y[129]+y[1]*y[2]*y[126]*y[130]+y[1]*y[2]*y[124]*y[126]*\
y[130]+2.*y[1]*y[2]*y[122]*y[128]*y[131]-y[1]*y[4]*y[122]*y[128]*y[131]+y[1\
]*y[2]*y[130]*y[131]+y[1]*y[2]*y[126]*y[132]+y[1]*y[2]*y[124]*y[126]*y[132]\
+y[1]*y[2]*y[131]*y[132]));
return (FOUT);
}
