#include "intfile.hh"

dcmplx Pf1(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[77];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[0];
y[3]=-x0;
y[4]=1.+y[3];
y[5]=2.*x1*y[1]*y[2];
y[6]=esx[0];
y[7]=lrs[0];
y[8]=y[1]*y[2];
y[9]=x1*x1;
y[10]=y[1]*y[2]*y[9];
y[11]=-(x1*y[1]*y[6]);
y[12]=y[5]+y[8]+y[10]+y[11];
y[13]=-x1;
y[14]=1.+y[13];
y[15]=2.*y[1]*y[2];
y[16]=lrs[1];
y[17]=2.*x0*y[1]*y[2];
y[18]=-(y[1]*y[6]);
y[19]=2.*x0*x1*y[1]*y[2];
y[20]=-(x0*y[1]*y[6]);
y[21]=y[5]+y[15]+y[17]+y[18]+y[19]+y[20];
y[22]=-(lambda*MYI*x0*y[4]*y[7]*y[12]);
y[23]=x0+y[22];
y[24]=-(lambda*MYI*x1*y[14]*y[16]*y[21]);
y[25]=x1+y[24];
y[26]=pow(bi,-2);
y[27]=lambda*lambda;
y[28]=y[5]+y[15]+y[18];
y[29]=pow(y[28],2);
y[30]=x0*x1*y[4]*y[7]*y[14]*y[16]*y[27]*y[29];
y[31]=-(lambda*MYI*y[4]*y[7]*y[12]);
y[32]=lambda*MYI*x0*y[7]*y[12];
y[33]=1.+y[31]+y[32];
y[34]=y[15]+y[17];
y[35]=-(lambda*MYI*x1*y[14]*y[16]*y[34]);
y[36]=-(lambda*MYI*y[14]*y[16]*y[21]);
y[37]=lambda*MYI*x1*y[16]*y[21];
y[38]=1.+y[35]+y[36]+y[37];
y[39]=y[33]*y[38];
y[40]=y[30]+y[39];
y[41]=y[1]*y[23];
y[42]=y[1]*y[25];
y[43]=y[1]*y[23]*y[25];
y[44]=y[1]+y[41]+y[42]+y[43];
y[45]=pow(y[44],-2);
y[46]=pow(y[25],2);
y[47]=1./x2;
y[48]=x1*y[1]*y[2];
y[49]=lrs[2];
y[50]=-x2;
y[51]=1.+y[50];
y[52]=y[8]+y[17];
y[53]=2.*x2*y[1]*y[2];
y[54]=2.*x1*x2*y[1]*y[2];
y[55]=x2*x2;
y[56]=y[1]*y[2]*y[55];
y[57]=-(x2*y[1]*y[6]);
y[58]=y[5]+y[8]+y[10]+y[11]+y[53]+y[54]+y[56]+y[57];
y[59]=y[5]+y[15]+y[18]+y[53];
y[60]=2.*x0*x2*y[1]*y[2];
y[61]=x2*y[1]*y[2];
y[62]=y[5]+y[15]+y[17]+y[18]+y[19]+y[20]+y[60]+y[61];
y[63]=pow(y[59],2);
y[64]=x0*x1*y[4]*y[7]*y[14]*y[16]*y[27]*y[63];
y[65]=-(lambda*MYI*y[4]*y[7]*y[58]);
y[66]=lambda*MYI*x0*y[7]*y[58];
y[67]=1.+y[65]+y[66];
y[68]=-(lambda*MYI*y[14]*y[16]*y[62]);
y[69]=lambda*MYI*x1*y[16]*y[62];
y[70]=1.+y[35]+y[68]+y[69];
y[71]=y[8]+y[17]+y[19]+y[20]+y[48]+y[60];
y[72]=-(lambda*MYI*y[49]*y[51]*y[71]);
y[73]=-(lambda*MYI*x0*y[4]*y[7]*y[58]);
y[74]=x0+y[73];
y[75]=-(lambda*MYI*x1*y[14]*y[16]*y[62]);
y[76]=x1+y[75];
FOUT=myLog(bi)*y[26]*y[40]*y[45]+3.*myLog(y[44])*y[26]*y[40]*y[45]-2.*myLog(\
y[8]+y[1]*y[2]*y[23]+2.*y[1]*y[2]*y[25]-y[1]*y[6]*y[25]+2.*y[1]*y[2]*y[23]*\
y[25]-y[1]*y[6]*y[23]*y[25]+y[1]*y[2]*y[46]+y[1]*y[2]*y[23]*y[46])*y[26]*y[\
40]*y[45]+myLog(1.-lambda*MYI*(y[8]+y[17]+y[19]+y[20]+y[48])*y[49])*y[26]*y\
[40]*y[45]-y[26]*y[40]*y[45]*y[47]+(pow(y[1]+y[1]*y[74]+y[1]*(x2-lambda*MYI\
*x2*y[49]*y[51]*y[71])*y[74]+y[1]*y[76]+y[1]*y[74]*y[76],-2)*y[26]*y[47]*(l\
ambda*MYI*x2*y[49]*y[51]*y[52]*(y[64]-lambda*MYI*x1*y[14]*y[16]*y[52]*y[67]\
)-lambda*MYI*x2*y[49]*y[51]*y[59]*(-(x0*x1*y[4]*y[7]*y[14]*y[16]*y[27]*y[52\
]*y[59])+lambda*MYI*x0*y[4]*y[7]*y[59]*y[70])+(y[64]+y[67]*y[70])*(1.-2.*la\
mbda*MYI*x0*x2*y[1]*y[2]*y[49]*y[51]+lambda*MYI*x2*y[49]*y[71]+y[72])))/(1.\
+y[72]);
return (FOUT);
}
