#include "intfile.hh"

dcmplx Pf4(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[123];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[0];
y[3]=-x0;
y[4]=1.+y[3];
y[5]=2.*y[1]*y[2];
y[6]=lrs[0];
y[7]=2.*x0*y[1]*y[2];
y[8]=2.*x1*y[1]*y[2];
y[9]=esx[0];
y[10]=-(y[1]*y[9]);
y[11]=2.*x0*x1*y[1]*y[2];
y[12]=-(x1*y[1]*y[9]);
y[13]=y[5]+y[7]+y[8]+y[10]+y[11]+y[12];
y[14]=-x1;
y[15]=1.+y[14];
y[16]=lrs[1];
y[17]=y[1]*y[2];
y[18]=x0*x0;
y[19]=y[1]*y[2]*y[18];
y[20]=-(x0*y[1]*y[9]);
y[21]=y[7]+y[17]+y[19]+y[20];
y[22]=-(lambda*MYI*x0*y[4]*y[6]*y[13]);
y[23]=x0+y[22];
y[24]=-(lambda*MYI*x1*y[15]*y[16]*y[21]);
y[25]=x1+y[24];
y[26]=1./x2;
y[27]=pow(bi,-2);
y[28]=lambda*lambda;
y[29]=y[5]+y[7]+y[10];
y[30]=pow(y[29],2);
y[31]=x0*x1*y[4]*y[6]*y[15]*y[16]*y[28]*y[30];
y[32]=y[5]+y[8];
y[33]=-(lambda*MYI*x0*y[4]*y[6]*y[32]);
y[34]=-(lambda*MYI*y[4]*y[6]*y[13]);
y[35]=lambda*MYI*x0*y[6]*y[13];
y[36]=1.+y[33]+y[34]+y[35];
y[37]=-(lambda*MYI*y[15]*y[16]*y[21]);
y[38]=lambda*MYI*x1*y[16]*y[21];
y[39]=1.+y[37]+y[38];
y[40]=y[36]*y[39];
y[41]=y[31]+y[40];
y[42]=y[1]*y[23];
y[43]=y[1]*y[25];
y[44]=y[1]*y[23]*y[25];
y[45]=y[1]+y[42]+y[43]+y[44];
y[46]=pow(y[45],-2);
y[47]=pow(y[23],2);
y[48]=myLog(bi);
y[49]=x0*y[1]*y[2];
y[50]=-(x0*x1*y[1]*y[9]);
y[51]=y[8]+y[11]+y[17]+y[49]+y[50];
y[52]=lrs[2];
y[53]=-(lambda*MYI*y[51]*y[52]);
y[54]=1.+y[53];
y[55]=myLog(y[54]);
y[56]=myLog(y[45]);
y[57]=2.*y[1]*y[2]*y[23];
y[58]=-(y[1]*y[9]*y[23]);
y[59]=y[1]*y[2]*y[47];
y[60]=y[1]*y[2]*y[25];
y[61]=2.*y[1]*y[2]*y[23]*y[25];
y[62]=-(y[1]*y[9]*y[23]*y[25]);
y[63]=y[1]*y[2]*y[25]*y[47];
y[64]=y[17]+y[57]+y[58]+y[59]+y[60]+y[61]+y[62]+y[63];
y[65]=myLog(y[64]);
y[66]=myLog(x2);
y[67]=-x2;
y[68]=1.+y[67];
y[69]=2.*x2*y[1]*y[2];
y[70]=y[5]+y[7]+y[20]+y[69];
y[71]=2.*x1*x2*y[1]*y[2];
y[72]=x2*y[1]*y[2];
y[73]=-(x1*x2*y[1]*y[9]);
y[74]=y[5]+y[7]+y[8]+y[10]+y[11]+y[12]+y[71]+y[72]+y[73];
y[75]=y[8]+y[12]+y[17];
y[76]=-(x2*y[1]*y[9]);
y[77]=y[5]+y[7]+y[10]+y[69]+y[76];
y[78]=2.*x0*x2*y[1]*y[2];
y[79]=x2*x2;
y[80]=y[1]*y[2]*y[79];
y[81]=-(x0*x2*y[1]*y[9]);
y[82]=y[7]+y[17]+y[19]+y[20]+y[69]+y[78]+y[80]+y[81];
y[83]=-(lambda*MYI*y[4]*y[6]*y[74]);
y[84]=lambda*MYI*x0*y[6]*y[74];
y[85]=1.+y[33]+y[83]+y[84];
y[86]=-(lambda*MYI*y[15]*y[16]*y[82]);
y[87]=lambda*MYI*x1*y[16]*y[82];
y[88]=1.+y[86]+y[87];
y[89]=y[8]+y[11]+y[17]+y[49]+y[50]+y[71];
y[90]=-(lambda*MYI*y[52]*y[68]*y[89]);
y[91]=-(lambda*MYI*x0*y[4]*y[6]*y[74]);
y[92]=x0+y[91];
y[93]=-(lambda*MYI*x1*y[15]*y[16]*y[82]);
y[94]=x1+y[93];
y[95]=1.+y[90];
y[96]=1./y[95];
y[97]=x0*x1*y[4]*y[6]*y[15]*y[16]*y[28]*y[75]*y[77];
y[98]=-(lambda*MYI*x1*y[15]*y[16]*y[70]*y[85]);
y[99]=y[97]+y[98];
y[100]=lambda*MYI*x2*y[52]*y[68]*y[70]*y[99];
y[101]=-(x0*x1*y[4]*y[6]*y[15]*y[16]*y[28]*y[70]*y[77]);
y[102]=lambda*MYI*x0*y[4]*y[6]*y[75]*y[88];
y[103]=y[101]+y[102];
y[104]=-(lambda*MYI*x2*y[52]*y[68]*y[75]*y[103]);
y[105]=pow(y[77],2);
y[106]=x0*x1*y[4]*y[6]*y[15]*y[16]*y[28]*y[105];
y[107]=y[85]*y[88];
y[108]=y[106]+y[107];
y[109]=-2.*lambda*MYI*x1*x2*y[1]*y[2]*y[52]*y[68];
y[110]=lambda*MYI*x2*y[52]*y[89];
y[111]=1.+y[90]+y[109]+y[110];
y[112]=y[108]*y[111];
y[113]=y[100]+y[104]+y[112];
y[114]=y[1]*y[92];
y[115]=y[1]*y[94];
y[116]=y[1]*y[92]*y[94];
y[117]=-(lambda*MYI*x2*y[52]*y[68]*y[89]);
y[118]=x2+y[117];
y[119]=y[1]*y[94]*y[118];
y[120]=y[1]+y[114]+y[115]+y[116]+y[119];
y[121]=pow(y[120],-2);
y[122]=pow(y[92],2);
FOUT=-(y[26]*(y[27]*y[41]*y[46]*y[48]+y[27]*y[41]*y[46]*y[55]+3.*y[27]*y[41]\
*y[46]*y[56]-2.*y[27]*y[41]*y[46]*y[65]))+0.5*(y[41]*y[46]*(pow(y[48],2)*y[\
27]+pow(y[55],2)*y[27]+2.*y[27]*y[48]*y[55])+2.*(y[27]*y[41]*y[48]+y[27]*y[\
41]*y[55])*(3.*y[46]*y[56]-2.*y[46]*y[65])+y[27]*y[41]*(9.*pow(y[56],2)*y[4\
6]+4.*pow(y[65],2)*y[46]-12.*y[46]*y[56]*y[65]))-y[26]*y[27]*y[41]*y[46]*y[\
66]+y[26]*y[27]*y[66]*y[96]*y[113]*y[121]+y[26]*(myLog(y[95])*y[27]*y[96]*y\
[113]*y[121]+3.*myLog(y[120])*y[27]*y[96]*y[113]*y[121]-2.*myLog(y[17]+2.*y\
[1]*y[2]*y[92]-y[1]*y[9]*y[92]+y[1]*y[2]*y[94]+pow(y[118],2)*y[1]*y[2]*y[94\
]+2.*y[1]*y[2]*y[92]*y[94]-y[1]*y[9]*y[92]*y[94]+y[1]*y[2]*y[118]+y[1]*y[2]\
*y[92]*y[118]+2.*y[1]*y[2]*y[94]*y[118]+2.*y[1]*y[2]*y[92]*y[94]*y[118]-y[1\
]*y[9]*y[92]*y[94]*y[118]+y[1]*y[2]*y[122]+y[1]*y[2]*y[94]*y[122])*y[27]*y[\
96]*y[113]*y[121]+y[27]*y[48]*y[96]*y[113]*y[121]);
return (FOUT);
}
