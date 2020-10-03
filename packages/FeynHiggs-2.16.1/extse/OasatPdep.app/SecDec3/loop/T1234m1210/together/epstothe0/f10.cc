#include "intfile.hh"

dcmplx Pf10(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[91];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=em[1];
y[5]=esx[0];
y[6]=-x1;
y[7]=1.+y[6];
y[8]=x0*y[1]*y[2];
y[9]=y[1]*y[2]*y[3];
y[10]=x0*y[1]*y[4];
y[11]=y[1]*y[3]*y[4];
y[12]=2.*x2*y[1]*y[3]*y[4];
y[13]=-(x0*y[1]*y[5]);
y[14]=-(y[1]*y[3]*y[5]);
y[15]=y[8]+y[9]+y[10]+y[11]+y[12]+y[13]+y[14];
y[16]=-x0;
y[17]=1.+y[16];
y[18]=x2*x2;
y[19]=lrs[0];
y[20]=y[1]*y[2];
y[21]=2.*x1*y[1]*y[2];
y[22]=2.*x0*x1*y[1]*y[2];
y[23]=x2*y[1]*y[2];
y[24]=x2*y[1]*y[4];
y[25]=-(x2*y[1]*y[5]);
y[26]=x1*x2*y[1]*y[2];
y[27]=2.*x0*x1*x2*y[1]*y[2];
y[28]=x1*x2*y[1]*y[4];
y[29]=2.*x0*x1*x2*y[1]*y[4];
y[30]=y[1]*y[4]*y[18];
y[31]=2.*x0*x1*y[1]*y[4]*y[18];
y[32]=-(x1*x2*y[1]*y[5]);
y[33]=-2.*x0*x1*x2*y[1]*y[5];
y[34]=y[20]+y[21]+y[22]+y[23]+y[24]+y[25]+y[26]+y[27]+y[28]+y[29]+y[30]+y[31\
]+y[32]+y[33];
y[35]=lrs[1];
y[36]=-x2;
y[37]=1.+y[36];
y[38]=x1*y[1]*y[2];
y[39]=y[1]*y[4];
y[40]=x1*y[1]*y[4];
y[41]=2.*x0*x1*y[1]*y[4];
y[42]=2.*x2*y[1]*y[4];
y[43]=4.*x0*x1*x2*y[1]*y[4];
y[44]=-(y[1]*y[5]);
y[45]=-(x1*y[1]*y[5]);
y[46]=-2.*x0*x1*y[1]*y[5];
y[47]=y[20]+y[22]+y[38]+y[39]+y[40]+y[41]+y[42]+y[43]+y[44]+y[45]+y[46];
y[48]=lambda*lambda;
y[49]=2.*y[1]*y[2];
y[50]=2.*x0*y[1]*y[2];
y[51]=2.*x0*x2*y[1]*y[2];
y[52]=2.*x0*x2*y[1]*y[4];
y[53]=2.*x0*y[1]*y[4]*y[18];
y[54]=-2.*x0*x2*y[1]*y[5];
y[55]=y[23]+y[24]+y[25]+y[49]+y[50]+y[51]+y[52]+y[53]+y[54];
y[56]=x0*x2*y[1]*y[2];
y[57]=x2*y[1]*y[2]*y[3];
y[58]=x0*x2*y[1]*y[4];
y[59]=x2*y[1]*y[3]*y[4];
y[60]=y[1]*y[3]*y[4]*y[18];
y[61]=-(x0*x2*y[1]*y[5]);
y[62]=-(x2*y[1]*y[3]*y[5]);
y[63]=y[9]+y[20]+y[50]+y[56]+y[57]+y[58]+y[59]+y[60]+y[61]+y[62];
y[64]=lrs[2];
y[65]=2.*x1*x2*y[1]*y[2];
y[66]=2.*x1*x2*y[1]*y[4];
y[67]=2.*x1*y[1]*y[4]*y[18];
y[68]=-2.*x1*x2*y[1]*y[5];
y[69]=y[21]+y[65]+y[66]+y[67]+y[68];
y[70]=-(lambda*MYI*x0*y[17]*y[19]*y[69]);
y[71]=-(lambda*MYI*y[17]*y[19]*y[34]);
y[72]=lambda*MYI*x0*y[19]*y[34];
y[73]=1.+y[70]+y[71]+y[72];
y[74]=-(lambda*MYI*y[7]*y[35]*y[63]);
y[75]=lambda*MYI*x1*y[35]*y[63];
y[76]=1.+y[74]+y[75];
y[77]=x0*x1*y[1]*y[2];
y[78]=x1*y[1]*y[2]*y[3];
y[79]=x0*x1*y[1]*y[4];
y[80]=x1*y[1]*y[3]*y[4];
y[81]=2.*x1*x2*y[1]*y[3]*y[4];
y[82]=-(x0*x1*y[1]*y[5]);
y[83]=-(x1*y[1]*y[3]*y[5]);
y[84]=y[8]+y[10]+y[13]+y[20]+y[52]+y[77]+y[78]+y[79]+y[80]+y[81]+y[82]+y[83]\
;
y[85]=-(lambda*MYI*x1*y[7]*y[35]*y[63]);
y[86]=x1+y[85];
y[87]=-(lambda*MYI*x0*y[17]*y[19]*y[34]);
y[88]=x0+y[87];
y[89]=-(lambda*MYI*x2*y[37]*y[64]*y[84]);
y[90]=x2+y[89];
FOUT=pow(bi,-2)*pow(y[1]+y[1]*y[86]+y[1]*y[86]*y[88]+y[1]*y[90]+y[1]*y[86]*y\
[88]*y[90],-2)*(lambda*MYI*x2*y[15]*y[37]*y[64]*(x0*x1*y[7]*y[17]*y[19]*y[3\
5]*y[47]*y[48]*y[55]-lambda*MYI*x1*y[7]*y[15]*y[35]*y[73])-lambda*MYI*x2*y[\
37]*y[47]*y[64]*(-(x0*x1*y[7]*y[15]*y[17]*y[19]*y[35]*y[48]*y[55])+lambda*M\
YI*x0*y[17]*y[19]*y[47]*y[76])+(x0*x1*pow(y[55],2)*y[7]*y[17]*y[19]*y[35]*y\
[48]+y[73]*y[76])*(1.-lambda*MYI*x2*(2.*x0*y[1]*y[4]+2.*x1*y[1]*y[3]*y[4])*\
y[37]*y[64]+lambda*MYI*x2*y[64]*y[84]-lambda*MYI*y[37]*y[64]*y[84]));
return (FOUT);
}
