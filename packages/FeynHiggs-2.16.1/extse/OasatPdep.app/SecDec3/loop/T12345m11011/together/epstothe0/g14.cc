#include "intfile.hh"

double Pr14(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double x3=x[3];
double y[76];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x3*x3;
y[4]=esx[0];
y[5]=-x0;
y[6]=1.+y[5];
y[7]=y[1]*y[2];
y[8]=x1*y[1]*y[2];
y[9]=2.*x2*y[1]*y[2];
y[10]=2.*x0*x2*y[1]*y[2];
y[11]=2.*x1*x2*y[1]*y[2];
y[12]=2.*x3*y[1]*y[2];
y[13]=2.*x1*x3*y[1]*y[2];
y[14]=2.*x2*x3*y[1]*y[2];
y[15]=4.*x0*x2*x3*y[1]*y[2];
y[16]=2.*x1*x2*x3*y[1]*y[2];
y[17]=y[1]*y[2]*y[3];
y[18]=x1*y[1]*y[2]*y[3];
y[19]=2.*x0*x2*y[1]*y[2]*y[3];
y[20]=-(x2*y[1]*y[4]);
y[21]=-(x3*y[1]*y[4]);
y[22]=-(x1*x3*y[1]*y[4]);
y[23]=-2.*x0*x2*x3*y[1]*y[4];
y[24]=-(x1*x2*x3*y[1]*y[4]);
y[25]=y[7]+y[8]+y[9]+y[10]+y[11]+y[12]+y[13]+y[14]+y[15]+y[16]+y[17]+y[18]+y\
[19]+y[20]+y[21]+y[22]+y[23]+y[24];
y[26]=lrs[0];
y[27]=x0*x0;
y[28]=-x1;
y[29]=1.+y[28];
y[30]=2.*y[1]*y[2];
y[31]=x0*y[1]*y[2];
y[32]=2.*x1*y[1]*y[2];
y[33]=2.*x0*x3*y[1]*y[2];
y[34]=2.*x0*x2*x3*y[1]*y[2];
y[35]=x0*y[1]*y[2]*y[3];
y[36]=-(y[1]*y[4]);
y[37]=-(x0*x3*y[1]*y[4]);
y[38]=-(x0*x2*x3*y[1]*y[4]);
y[39]=y[9]+y[10]+y[11]+y[12]+y[13]+y[20]+y[21]+y[30]+y[31]+y[32]+y[33]+y[34]\
+y[35]+y[36]+y[37]+y[38];
y[40]=lrs[1];
y[41]=x1*x1;
y[42]=-x2;
y[43]=1.+y[42];
y[44]=2.*x0*y[1]*y[2];
y[45]=y[1]*y[2]*y[27];
y[46]=2.*x0*x1*y[1]*y[2];
y[47]=y[1]*y[2]*y[41];
y[48]=2.*x3*y[1]*y[2]*y[27];
y[49]=2.*x0*x1*x3*y[1]*y[2];
y[50]=y[1]*y[2]*y[3]*y[27];
y[51]=-(x0*y[1]*y[4]);
y[52]=-(x1*y[1]*y[4]);
y[53]=-(x3*y[1]*y[4]*y[27]);
y[54]=-(x0*x1*x3*y[1]*y[4]);
y[55]=y[7]+y[32]+y[33]+y[44]+y[45]+y[46]+y[47]+y[48]+y[49]+y[50]+y[51]+y[52]\
+y[53]+y[54];
y[56]=lrs[2];
y[57]=-x3;
y[58]=1.+y[57];
y[59]=2.*x2*y[1]*y[2]*y[27];
y[60]=2.*x0*x1*x2*y[1]*y[2];
y[61]=2.*x2*x3*y[1]*y[2]*y[27];
y[62]=-(x0*x1*y[1]*y[4]);
y[63]=-(x2*y[1]*y[4]*y[27]);
y[64]=-(x0*x1*x2*y[1]*y[4]);
y[65]=y[7]+y[10]+y[32]+y[33]+y[44]+y[46]+y[47]+y[49]+y[51]+y[52]+y[59]+y[60]\
+y[61]+y[62]+y[63]+y[64];
y[66]=lrs[3];
y[67]=pow(y[6],2);
y[68]=pow(y[25],2);
y[69]=pow(y[26],2);
y[70]=pow(y[29],2);
y[71]=pow(y[39],2);
y[72]=pow(y[40],2);
y[73]=pow(y[58],2);
y[74]=pow(y[65],2);
y[75]=pow(y[66],2);
FOUT=(2.*x0*x1*x2*y[1]*y[2]*y[6]*y[25]*y[26]*y[29]*y[39]*y[40]*y[43]*y[55]*y\
[56]+2.*x0*x1*x2*x3*y[1]*y[2]*y[6]*y[25]*y[26]*y[29]*y[39]*y[40]*y[43]*y[55\
]*y[56]-x0*x1*x2*x3*y[1]*y[4]*y[6]*y[25]*y[26]*y[29]*y[39]*y[40]*y[43]*y[55\
]*y[56]+2.*x0*x1*x3*y[1]*y[2]*y[6]*y[25]*y[26]*y[29]*y[39]*y[40]*y[58]*y[65\
]*y[66]+2.*x0*x1*x2*x3*y[1]*y[2]*y[6]*y[25]*y[26]*y[29]*y[39]*y[40]*y[58]*y\
[65]*y[66]+2.*x0*x1*y[1]*y[2]*y[3]*y[6]*y[25]*y[26]*y[29]*y[39]*y[40]*y[58]\
*y[65]*y[66]-x0*x1*x3*y[1]*y[4]*y[6]*y[25]*y[26]*y[29]*y[39]*y[40]*y[58]*y[\
65]*y[66]-x0*x1*x2*x3*y[1]*y[4]*y[6]*y[25]*y[26]*y[29]*y[39]*y[40]*y[58]*y[\
65]*y[66]+2.*x0*x2*x3*y[1]*y[2]*y[6]*y[25]*y[26]*y[43]*y[55]*y[56]*y[58]*y[\
65]*y[66]+2.*x0*x1*x2*x3*y[1]*y[2]*y[6]*y[25]*y[26]*y[43]*y[55]*y[56]*y[58]\
*y[65]*y[66]-x0*x1*x2*x3*y[1]*y[4]*y[6]*y[25]*y[26]*y[43]*y[55]*y[56]*y[58]\
*y[65]*y[66]+4.*x2*x3*y[1]*y[2]*y[6]*y[25]*y[26]*y[27]*y[43]*y[55]*y[56]*y[\
58]*y[65]*y[66]+4.*x2*y[1]*y[2]*y[3]*y[6]*y[25]*y[26]*y[27]*y[43]*y[55]*y[5\
6]*y[58]*y[65]*y[66]-2.*x2*x3*y[1]*y[4]*y[6]*y[25]*y[26]*y[27]*y[43]*y[55]*\
y[56]*y[58]*y[65]*y[66]+2.*x0*x1*x2*x3*y[1]*y[2]*y[29]*y[39]*y[40]*y[43]*y[\
55]*y[56]*y[58]*y[65]*y[66]-x0*x1*x2*x3*y[1]*y[4]*y[29]*y[39]*y[40]*y[43]*y\
[55]*y[56]*y[58]*y[65]*y[66]+x2*y[1]*y[2]*y[27]*y[43]*y[55]*y[56]*y[67]*y[6\
8]*y[69]+2.*x2*x3*y[1]*y[2]*y[27]*y[43]*y[55]*y[56]*y[67]*y[68]*y[69]+x2*y[\
1]*y[2]*y[3]*y[27]*y[43]*y[55]*y[56]*y[67]*y[68]*y[69]-x2*x3*y[1]*y[4]*y[27\
]*y[43]*y[55]*y[56]*y[67]*y[68]*y[69]+2.*x2*x3*y[1]*y[2]*y[27]*y[58]*y[65]*\
y[66]*y[67]*y[68]*y[69]+2.*x2*y[1]*y[2]*y[3]*y[27]*y[58]*y[65]*y[66]*y[67]*\
y[68]*y[69]-x2*x3*y[1]*y[4]*y[27]*y[58]*y[65]*y[66]*y[67]*y[68]*y[69]+x2*y[\
1]*y[2]*y[41]*y[43]*y[55]*y[56]*y[70]*y[71]*y[72]+x3*y[1]*y[2]*y[41]*y[58]*\
y[65]*y[66]*y[70]*y[71]*y[72]+x0*y[1]*y[2]*y[3]*y[6]*y[25]*y[26]*y[73]*y[74\
]*y[75]+x0*x1*y[1]*y[2]*y[3]*y[6]*y[25]*y[26]*y[73]*y[74]*y[75]+2.*x2*y[1]*\
y[2]*y[3]*y[6]*y[25]*y[26]*y[27]*y[73]*y[74]*y[75]+x0*x1*y[1]*y[2]*y[3]*y[2\
9]*y[39]*y[40]*y[73]*y[74]*y[75]+x2*y[1]*y[2]*y[3]*y[27]*y[43]*y[55]*y[56]*\
y[73]*y[74]*y[75])/(-(x0*y[1]*y[2]*y[6]*y[25]*y[26])-x0*x1*y[1]*y[2]*y[6]*y\
[25]*y[26]-2.*x0*x2*y[1]*y[2]*y[6]*y[25]*y[26]-2.*x0*x1*x2*y[1]*y[2]*y[6]*y\
[25]*y[26]-2.*x0*x3*y[1]*y[2]*y[6]*y[25]*y[26]-2.*x0*x1*x3*y[1]*y[2]*y[6]*y\
[25]*y[26]-2.*x0*x2*x3*y[1]*y[2]*y[6]*y[25]*y[26]-2.*x0*x1*x2*x3*y[1]*y[2]*\
y[6]*y[25]*y[26]-x0*y[1]*y[2]*y[3]*y[6]*y[25]*y[26]-x0*x1*y[1]*y[2]*y[3]*y[\
6]*y[25]*y[26]+x0*x2*y[1]*y[4]*y[6]*y[25]*y[26]+x0*x3*y[1]*y[4]*y[6]*y[25]*\
y[26]+x0*x1*x3*y[1]*y[4]*y[6]*y[25]*y[26]+x0*x1*x2*x3*y[1]*y[4]*y[6]*y[25]*\
y[26]-2.*x2*y[1]*y[2]*y[6]*y[25]*y[26]*y[27]-4.*x2*x3*y[1]*y[2]*y[6]*y[25]*\
y[26]*y[27]-2.*x2*y[1]*y[2]*y[3]*y[6]*y[25]*y[26]*y[27]+2.*x2*x3*y[1]*y[4]*\
y[6]*y[25]*y[26]*y[27]-2.*x1*y[1]*y[2]*y[29]*y[39]*y[40]-x0*x1*y[1]*y[2]*y[\
29]*y[39]*y[40]-2.*x1*x2*y[1]*y[2]*y[29]*y[39]*y[40]-2.*x0*x1*x2*y[1]*y[2]*\
y[29]*y[39]*y[40]-2.*x1*x3*y[1]*y[2]*y[29]*y[39]*y[40]-2.*x0*x1*x3*y[1]*y[2\
]*y[29]*y[39]*y[40]-2.*x0*x1*x2*x3*y[1]*y[2]*y[29]*y[39]*y[40]-x0*x1*y[1]*y\
[2]*y[3]*y[29]*y[39]*y[40]+x1*y[1]*y[4]*y[29]*y[39]*y[40]+x1*x2*y[1]*y[4]*y\
[29]*y[39]*y[40]+x1*x3*y[1]*y[4]*y[29]*y[39]*y[40]+x0*x1*x3*y[1]*y[4]*y[29]\
*y[39]*y[40]+x0*x1*x2*x3*y[1]*y[4]*y[29]*y[39]*y[40]-2.*y[1]*y[2]*y[29]*y[3\
9]*y[40]*y[41]-2.*x2*y[1]*y[2]*y[29]*y[39]*y[40]*y[41]-2.*x3*y[1]*y[2]*y[29\
]*y[39]*y[40]*y[41]-x2*y[1]*y[2]*y[43]*y[55]*y[56]-2.*x0*x2*y[1]*y[2]*y[43]\
*y[55]*y[56]-2.*x1*x2*y[1]*y[2]*y[43]*y[55]*y[56]-2.*x0*x1*x2*y[1]*y[2]*y[4\
3]*y[55]*y[56]-2.*x0*x2*x3*y[1]*y[2]*y[43]*y[55]*y[56]-2.*x0*x1*x2*x3*y[1]*\
y[2]*y[43]*y[55]*y[56]+x0*x2*y[1]*y[4]*y[43]*y[55]*y[56]+x1*x2*y[1]*y[4]*y[\
43]*y[55]*y[56]+x0*x1*x2*x3*y[1]*y[4]*y[43]*y[55]*y[56]-x2*y[1]*y[2]*y[27]*\
y[43]*y[55]*y[56]-2.*x2*x3*y[1]*y[2]*y[27]*y[43]*y[55]*y[56]-x2*y[1]*y[2]*y\
[3]*y[27]*y[43]*y[55]*y[56]+x2*x3*y[1]*y[4]*y[27]*y[43]*y[55]*y[56]-x2*y[1]\
*y[2]*y[41]*y[43]*y[55]*y[56]-x3*y[1]*y[2]*y[58]*y[65]*y[66]-2.*x0*x3*y[1]*\
y[2]*y[58]*y[65]*y[66]-2.*x1*x3*y[1]*y[2]*y[58]*y[65]*y[66]-2.*x0*x1*x3*y[1\
]*y[2]*y[58]*y[65]*y[66]-2.*x0*x2*x3*y[1]*y[2]*y[58]*y[65]*y[66]-2.*x0*x1*x\
2*x3*y[1]*y[2]*y[58]*y[65]*y[66]-2.*x0*y[1]*y[2]*y[3]*y[58]*y[65]*y[66]-2.*\
x0*x1*y[1]*y[2]*y[3]*y[58]*y[65]*y[66]+x0*x3*y[1]*y[4]*y[58]*y[65]*y[66]+x1\
*x3*y[1]*y[4]*y[58]*y[65]*y[66]+x0*x1*x3*y[1]*y[4]*y[58]*y[65]*y[66]+x0*x1*\
x2*x3*y[1]*y[4]*y[58]*y[65]*y[66]-2.*x2*x3*y[1]*y[2]*y[27]*y[58]*y[65]*y[66\
]-2.*x2*y[1]*y[2]*y[3]*y[27]*y[58]*y[65]*y[66]+x2*x3*y[1]*y[4]*y[27]*y[58]*\
y[65]*y[66]-x3*y[1]*y[2]*y[41]*y[58]*y[65]*y[66]);
return (FOUT);
}
double Pm14(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double x3=x[3];
double y[76];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=esx[0];
y[5]=y[1]*y[2];
y[6]=2.*x0*y[1]*y[2];
y[7]=2.*x1*y[1]*y[2];
y[8]=2.*x0*x1*y[1]*y[2];
y[9]=x1*x1;
y[10]=y[1]*y[2]*y[9];
y[11]=2.*x0*x3*y[1]*y[2];
y[12]=2.*x0*x1*x3*y[1]*y[2];
y[13]=x3*x3;
y[14]=-(x0*y[1]*y[4]);
y[15]=-(x1*y[1]*y[4]);
y[16]=2.*x0*x2*y[1]*y[2];
y[17]=-x0;
y[18]=1.+y[17];
y[19]=x1*y[1]*y[2];
y[20]=2.*x2*y[1]*y[2];
y[21]=2.*x1*x2*y[1]*y[2];
y[22]=2.*x3*y[1]*y[2];
y[23]=2.*x1*x3*y[1]*y[2];
y[24]=2.*x2*x3*y[1]*y[2];
y[25]=4.*x0*x2*x3*y[1]*y[2];
y[26]=2.*x1*x2*x3*y[1]*y[2];
y[27]=y[1]*y[2]*y[13];
y[28]=x1*y[1]*y[2]*y[13];
y[29]=2.*x0*x2*y[1]*y[2]*y[13];
y[30]=-(x2*y[1]*y[4]);
y[31]=-(x3*y[1]*y[4]);
y[32]=-(x1*x3*y[1]*y[4]);
y[33]=-2.*x0*x2*x3*y[1]*y[4];
y[34]=-(x1*x2*x3*y[1]*y[4]);
y[35]=y[5]+y[16]+y[19]+y[20]+y[21]+y[22]+y[23]+y[24]+y[25]+y[26]+y[27]+y[28]\
+y[29]+y[30]+y[31]+y[32]+y[33]+y[34];
y[36]=lrs[0];
y[37]=-x1;
y[38]=1.+y[37];
y[39]=2.*y[1]*y[2];
y[40]=x0*y[1]*y[2];
y[41]=2.*x0*x2*x3*y[1]*y[2];
y[42]=x0*y[1]*y[2]*y[13];
y[43]=-(y[1]*y[4]);
y[44]=-(x0*x3*y[1]*y[4]);
y[45]=-(x0*x2*x3*y[1]*y[4]);
y[46]=y[7]+y[11]+y[16]+y[20]+y[21]+y[22]+y[23]+y[30]+y[31]+y[39]+y[40]+y[41]\
+y[42]+y[43]+y[44]+y[45];
y[47]=lrs[1];
y[48]=-x2;
y[49]=1.+y[48];
y[50]=y[1]*y[2]*y[3];
y[51]=2.*x3*y[1]*y[2]*y[3];
y[52]=y[1]*y[2]*y[3]*y[13];
y[53]=-(x3*y[1]*y[3]*y[4]);
y[54]=-(x0*x1*x3*y[1]*y[4]);
y[55]=y[5]+y[6]+y[7]+y[8]+y[10]+y[11]+y[12]+y[14]+y[15]+y[50]+y[51]+y[52]+y[\
53]+y[54];
y[56]=lrs[2];
y[57]=-x3;
y[58]=1.+y[57];
y[59]=2.*x2*y[1]*y[2]*y[3];
y[60]=2.*x0*x1*x2*y[1]*y[2];
y[61]=2.*x2*x3*y[1]*y[2]*y[3];
y[62]=-(x0*x1*y[1]*y[4]);
y[63]=-(x2*y[1]*y[3]*y[4]);
y[64]=-(x0*x1*x2*y[1]*y[4]);
y[65]=y[5]+y[6]+y[7]+y[8]+y[10]+y[11]+y[12]+y[14]+y[15]+y[16]+y[59]+y[60]+y[\
61]+y[62]+y[63]+y[64];
y[66]=lrs[3];
y[67]=pow(y[18],2);
y[68]=pow(y[35],2);
y[69]=pow(y[36],2);
y[70]=pow(y[38],2);
y[71]=pow(y[46],2);
y[72]=pow(y[47],2);
y[73]=pow(y[58],2);
y[74]=pow(y[65],2);
y[75]=pow(y[66],2);
FOUT=pow(lambda*(-(x0*y[1]*y[2]*y[18]*y[35]*y[36])-x0*x1*y[1]*y[2]*y[18]*y[3\
5]*y[36]-2.*x0*x2*y[1]*y[2]*y[18]*y[35]*y[36]-2.*x0*x1*x2*y[1]*y[2]*y[18]*y\
[35]*y[36]-2.*x0*x3*y[1]*y[2]*y[18]*y[35]*y[36]-2.*x0*x1*x3*y[1]*y[2]*y[18]\
*y[35]*y[36]-2.*x0*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]-2.*x0*x1*x2*x3*y[1]*y[\
2]*y[18]*y[35]*y[36]-2.*x2*y[1]*y[2]*y[3]*y[18]*y[35]*y[36]-4.*x2*x3*y[1]*y\
[2]*y[3]*y[18]*y[35]*y[36]+x0*x2*y[1]*y[4]*y[18]*y[35]*y[36]+x0*x3*y[1]*y[4\
]*y[18]*y[35]*y[36]+x0*x1*x3*y[1]*y[4]*y[18]*y[35]*y[36]+x0*x1*x2*x3*y[1]*y\
[4]*y[18]*y[35]*y[36]+2.*x2*x3*y[1]*y[3]*y[4]*y[18]*y[35]*y[36]-x0*y[1]*y[2\
]*y[13]*y[18]*y[35]*y[36]-x0*x1*y[1]*y[2]*y[13]*y[18]*y[35]*y[36]-2.*x2*y[1\
]*y[2]*y[3]*y[13]*y[18]*y[35]*y[36]-2.*x1*y[1]*y[2]*y[38]*y[46]*y[47]-x0*x1\
*y[1]*y[2]*y[38]*y[46]*y[47]-2.*x1*x2*y[1]*y[2]*y[38]*y[46]*y[47]-2.*x0*x1*\
x2*y[1]*y[2]*y[38]*y[46]*y[47]-2.*x1*x3*y[1]*y[2]*y[38]*y[46]*y[47]-2.*x0*x\
1*x3*y[1]*y[2]*y[38]*y[46]*y[47]-2.*x0*x1*x2*x3*y[1]*y[2]*y[38]*y[46]*y[47]\
+x1*y[1]*y[4]*y[38]*y[46]*y[47]+x1*x2*y[1]*y[4]*y[38]*y[46]*y[47]+x1*x3*y[1\
]*y[4]*y[38]*y[46]*y[47]+x0*x1*x3*y[1]*y[4]*y[38]*y[46]*y[47]+x0*x1*x2*x3*y\
[1]*y[4]*y[38]*y[46]*y[47]-2.*y[1]*y[2]*y[9]*y[38]*y[46]*y[47]-2.*x2*y[1]*y\
[2]*y[9]*y[38]*y[46]*y[47]-2.*x3*y[1]*y[2]*y[9]*y[38]*y[46]*y[47]-x0*x1*y[1\
]*y[2]*y[13]*y[38]*y[46]*y[47]-x2*y[1]*y[2]*y[49]*y[55]*y[56]-2.*x0*x2*y[1]\
*y[2]*y[49]*y[55]*y[56]-2.*x1*x2*y[1]*y[2]*y[49]*y[55]*y[56]-2.*x0*x1*x2*y[\
1]*y[2]*y[49]*y[55]*y[56]-2.*x0*x2*x3*y[1]*y[2]*y[49]*y[55]*y[56]-2.*x0*x1*\
x2*x3*y[1]*y[2]*y[49]*y[55]*y[56]-x2*y[1]*y[2]*y[3]*y[49]*y[55]*y[56]-2.*x2\
*x3*y[1]*y[2]*y[3]*y[49]*y[55]*y[56]+x0*x2*y[1]*y[4]*y[49]*y[55]*y[56]+x1*x\
2*y[1]*y[4]*y[49]*y[55]*y[56]+x0*x1*x2*x3*y[1]*y[4]*y[49]*y[55]*y[56]+x2*x3\
*y[1]*y[3]*y[4]*y[49]*y[55]*y[56]-x2*y[1]*y[2]*y[9]*y[49]*y[55]*y[56]-x2*y[\
1]*y[2]*y[3]*y[13]*y[49]*y[55]*y[56]-x3*y[1]*y[2]*y[58]*y[65]*y[66]-2.*x0*x\
3*y[1]*y[2]*y[58]*y[65]*y[66]-2.*x1*x3*y[1]*y[2]*y[58]*y[65]*y[66]-2.*x0*x1\
*x3*y[1]*y[2]*y[58]*y[65]*y[66]-2.*x0*x2*x3*y[1]*y[2]*y[58]*y[65]*y[66]-2.*\
x0*x1*x2*x3*y[1]*y[2]*y[58]*y[65]*y[66]-2.*x2*x3*y[1]*y[2]*y[3]*y[58]*y[65]\
*y[66]+x0*x3*y[1]*y[4]*y[58]*y[65]*y[66]+x1*x3*y[1]*y[4]*y[58]*y[65]*y[66]+\
x0*x1*x3*y[1]*y[4]*y[58]*y[65]*y[66]+x0*x1*x2*x3*y[1]*y[4]*y[58]*y[65]*y[66\
]+x2*x3*y[1]*y[3]*y[4]*y[58]*y[65]*y[66]-x3*y[1]*y[2]*y[9]*y[58]*y[65]*y[66\
]-2.*x0*y[1]*y[2]*y[13]*y[58]*y[65]*y[66]-2.*x0*x1*y[1]*y[2]*y[13]*y[58]*y[\
65]*y[66]-2.*x2*y[1]*y[2]*y[3]*y[13]*y[58]*y[65]*y[66])-x2*pow(lambda,5)*y[\
1]*y[2]*y[3]*y[13]*y[49]*y[55]*y[56]*y[67]*y[68]*y[69]*y[73]*y[74]*y[75]+po\
w(lambda,3)*(2.*x0*x1*x2*y[1]*y[2]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]*y[49\
]*y[55]*y[56]+2.*x0*x1*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]*\
y[49]*y[55]*y[56]-x0*x1*x2*x3*y[1]*y[4]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]\
*y[49]*y[55]*y[56]+2.*x0*x1*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47\
]*y[58]*y[65]*y[66]+2.*x0*x1*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[38]*y[46]*\
y[47]*y[58]*y[65]*y[66]-x0*x1*x3*y[1]*y[4]*y[18]*y[35]*y[36]*y[38]*y[46]*y[\
47]*y[58]*y[65]*y[66]-x0*x1*x2*x3*y[1]*y[4]*y[18]*y[35]*y[36]*y[38]*y[46]*y\
[47]*y[58]*y[65]*y[66]+2.*x0*x1*y[1]*y[2]*y[13]*y[18]*y[35]*y[36]*y[38]*y[4\
6]*y[47]*y[58]*y[65]*y[66]+2.*x0*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[49]*y[\
55]*y[56]*y[58]*y[65]*y[66]+2.*x0*x1*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[49\
]*y[55]*y[56]*y[58]*y[65]*y[66]+4.*x2*x3*y[1]*y[2]*y[3]*y[18]*y[35]*y[36]*y\
[49]*y[55]*y[56]*y[58]*y[65]*y[66]-x0*x1*x2*x3*y[1]*y[4]*y[18]*y[35]*y[36]*\
y[49]*y[55]*y[56]*y[58]*y[65]*y[66]-2.*x2*x3*y[1]*y[3]*y[4]*y[18]*y[35]*y[3\
6]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]+4.*x2*y[1]*y[2]*y[3]*y[13]*y[18]*y[3\
5]*y[36]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]+2.*x0*x1*x2*x3*y[1]*y[2]*y[38]\
*y[46]*y[47]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]-x0*x1*x2*x3*y[1]*y[4]*y[38\
]*y[46]*y[47]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]+x2*y[1]*y[2]*y[3]*y[49]*y\
[55]*y[56]*y[67]*y[68]*y[69]+2.*x2*x3*y[1]*y[2]*y[3]*y[49]*y[55]*y[56]*y[67\
]*y[68]*y[69]-x2*x3*y[1]*y[3]*y[4]*y[49]*y[55]*y[56]*y[67]*y[68]*y[69]+x2*y\
[1]*y[2]*y[3]*y[13]*y[49]*y[55]*y[56]*y[67]*y[68]*y[69]+2.*x2*x3*y[1]*y[2]*\
y[3]*y[58]*y[65]*y[66]*y[67]*y[68]*y[69]-x2*x3*y[1]*y[3]*y[4]*y[58]*y[65]*y\
[66]*y[67]*y[68]*y[69]+2.*x2*y[1]*y[2]*y[3]*y[13]*y[58]*y[65]*y[66]*y[67]*y\
[68]*y[69]+x2*y[1]*y[2]*y[9]*y[49]*y[55]*y[56]*y[70]*y[71]*y[72]+x3*y[1]*y[\
2]*y[9]*y[58]*y[65]*y[66]*y[70]*y[71]*y[72]+x0*y[1]*y[2]*y[13]*y[18]*y[35]*\
y[36]*y[73]*y[74]*y[75]+x0*x1*y[1]*y[2]*y[13]*y[18]*y[35]*y[36]*y[73]*y[74]\
*y[75]+2.*x2*y[1]*y[2]*y[3]*y[13]*y[18]*y[35]*y[36]*y[73]*y[74]*y[75]+x0*x1\
*y[1]*y[2]*y[13]*y[38]*y[46]*y[47]*y[73]*y[74]*y[75]+x2*y[1]*y[2]*y[3]*y[13\
]*y[49]*y[55]*y[56]*y[73]*y[74]*y[75]),2)+pow(x0*x1*y[1]*y[2]+x2*y[1]*y[2]+\
x3*y[1]*y[2]+2.*x0*x1*x2*x3*y[1]*y[2]+x2*y[1]*y[2]*y[3]-x0*x2*y[1]*y[4]-x1*\
x2*y[1]*y[4]-x0*x1*x2*x3*y[1]*y[4]-x2*x3*y[1]*y[3]*y[4]+y[5]+y[7]+x2*y[1]*y\
[2]*y[9]+x3*y[1]*y[2]*y[9]+y[10]+y[11]+y[12]+x0*x1*y[1]*y[2]*y[13]+x2*y[1]*\
y[2]*y[3]*y[13]+y[15]+y[16]+y[21]+y[23]+y[32]+y[40]+y[41]+y[42]+y[44]+y[54]\
+y[60]+y[61]+lambda*lambda*(-(x0*x1*y[1]*y[2]*y[18]*y[35]*y[36]*y[38]*y[46]\
*y[47])-2.*x0*x1*x2*y[1]*y[2]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]-2.*x0*x1*\
x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]-2.*x0*x1*x2*x3*y[1]*y[2]*y\
[18]*y[35]*y[36]*y[38]*y[46]*y[47]+x0*x1*x3*y[1]*y[4]*y[18]*y[35]*y[36]*y[3\
8]*y[46]*y[47]+x0*x1*x2*x3*y[1]*y[4]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]-x0\
*x1*y[1]*y[2]*y[13]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]-2.*x0*x2*y[1]*y[2]*\
y[18]*y[35]*y[36]*y[49]*y[55]*y[56]-2.*x0*x1*x2*y[1]*y[2]*y[18]*y[35]*y[36]\
*y[49]*y[55]*y[56]-2.*x0*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[49]*y[55]*y[56\
]-2.*x0*x1*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[49]*y[55]*y[56]-2.*x2*y[1]*y\
[2]*y[3]*y[18]*y[35]*y[36]*y[49]*y[55]*y[56]-4.*x2*x3*y[1]*y[2]*y[3]*y[18]*\
y[35]*y[36]*y[49]*y[55]*y[56]+x0*x2*y[1]*y[4]*y[18]*y[35]*y[36]*y[49]*y[55]\
*y[56]+x0*x1*x2*x3*y[1]*y[4]*y[18]*y[35]*y[36]*y[49]*y[55]*y[56]+2.*x2*x3*y\
[1]*y[3]*y[4]*y[18]*y[35]*y[36]*y[49]*y[55]*y[56]-2.*x2*y[1]*y[2]*y[3]*y[13\
]*y[18]*y[35]*y[36]*y[49]*y[55]*y[56]-2.*x1*x2*y[1]*y[2]*y[38]*y[46]*y[47]*\
y[49]*y[55]*y[56]-2.*x0*x1*x2*y[1]*y[2]*y[38]*y[46]*y[47]*y[49]*y[55]*y[56]\
-2.*x0*x1*x2*x3*y[1]*y[2]*y[38]*y[46]*y[47]*y[49]*y[55]*y[56]+x1*x2*y[1]*y[\
4]*y[38]*y[46]*y[47]*y[49]*y[55]*y[56]+x0*x1*x2*x3*y[1]*y[4]*y[38]*y[46]*y[\
47]*y[49]*y[55]*y[56]-2.*x2*y[1]*y[2]*y[9]*y[38]*y[46]*y[47]*y[49]*y[55]*y[\
56]-2.*x0*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[58]*y[65]*y[66]-2.*x0*x1*x3*y[1]\
*y[2]*y[18]*y[35]*y[36]*y[58]*y[65]*y[66]-2.*x0*x2*x3*y[1]*y[2]*y[18]*y[35]\
*y[36]*y[58]*y[65]*y[66]-2.*x0*x1*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[58]*y\
[65]*y[66]-4.*x2*x3*y[1]*y[2]*y[3]*y[18]*y[35]*y[36]*y[58]*y[65]*y[66]+x0*x\
3*y[1]*y[4]*y[18]*y[35]*y[36]*y[58]*y[65]*y[66]+x0*x1*x3*y[1]*y[4]*y[18]*y[\
35]*y[36]*y[58]*y[65]*y[66]+x0*x1*x2*x3*y[1]*y[4]*y[18]*y[35]*y[36]*y[58]*y\
[65]*y[66]+2.*x2*x3*y[1]*y[3]*y[4]*y[18]*y[35]*y[36]*y[58]*y[65]*y[66]-2.*x\
0*y[1]*y[2]*y[13]*y[18]*y[35]*y[36]*y[58]*y[65]*y[66]-2.*x0*x1*y[1]*y[2]*y[\
13]*y[18]*y[35]*y[36]*y[58]*y[65]*y[66]-4.*x2*y[1]*y[2]*y[3]*y[13]*y[18]*y[\
35]*y[36]*y[58]*y[65]*y[66]-2.*x1*x3*y[1]*y[2]*y[38]*y[46]*y[47]*y[58]*y[65\
]*y[66]-2.*x0*x1*x3*y[1]*y[2]*y[38]*y[46]*y[47]*y[58]*y[65]*y[66]-2.*x0*x1*\
x2*x3*y[1]*y[2]*y[38]*y[46]*y[47]*y[58]*y[65]*y[66]+x1*x3*y[1]*y[4]*y[38]*y\
[46]*y[47]*y[58]*y[65]*y[66]+x0*x1*x3*y[1]*y[4]*y[38]*y[46]*y[47]*y[58]*y[6\
5]*y[66]+x0*x1*x2*x3*y[1]*y[4]*y[38]*y[46]*y[47]*y[58]*y[65]*y[66]-2.*x3*y[\
1]*y[2]*y[9]*y[38]*y[46]*y[47]*y[58]*y[65]*y[66]-2.*x0*x1*y[1]*y[2]*y[13]*y\
[38]*y[46]*y[47]*y[58]*y[65]*y[66]-2.*x0*x2*x3*y[1]*y[2]*y[49]*y[55]*y[56]*\
y[58]*y[65]*y[66]-2.*x0*x1*x2*x3*y[1]*y[2]*y[49]*y[55]*y[56]*y[58]*y[65]*y[\
66]-2.*x2*x3*y[1]*y[2]*y[3]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]+x0*x1*x2*x3\
*y[1]*y[4]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]+x2*x3*y[1]*y[3]*y[4]*y[49]*y\
[55]*y[56]*y[58]*y[65]*y[66]-2.*x2*y[1]*y[2]*y[3]*y[13]*y[49]*y[55]*y[56]*y\
[58]*y[65]*y[66]-x2*y[1]*y[2]*y[3]*y[67]*y[68]*y[69]-2.*x2*x3*y[1]*y[2]*y[3\
]*y[67]*y[68]*y[69]+x2*x3*y[1]*y[3]*y[4]*y[67]*y[68]*y[69]-x2*y[1]*y[2]*y[3\
]*y[13]*y[67]*y[68]*y[69]-y[1]*y[2]*y[9]*y[70]*y[71]*y[72]-x2*y[1]*y[2]*y[9\
]*y[70]*y[71]*y[72]-x3*y[1]*y[2]*y[9]*y[70]*y[71]*y[72]-x0*y[1]*y[2]*y[13]*\
y[73]*y[74]*y[75]-x0*x1*y[1]*y[2]*y[13]*y[73]*y[74]*y[75]-x2*y[1]*y[2]*y[3]\
*y[13]*y[73]*y[74]*y[75])+pow(lambda,4)*(2.*x0*x1*x2*x3*y[1]*y[2]*y[18]*y[3\
5]*y[36]*y[38]*y[46]*y[47]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]-x0*x1*x2*x3*\
y[1]*y[4]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]*y[49]*y[55]*y[56]*y[58]*y[65]\
*y[66]+2.*x2*x3*y[1]*y[2]*y[3]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]*y[67]*y[\
68]*y[69]-x2*x3*y[1]*y[3]*y[4]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]*y[67]*y[\
68]*y[69]+2.*x2*y[1]*y[2]*y[3]*y[13]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]*y[\
67]*y[68]*y[69]+x0*x1*y[1]*y[2]*y[13]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]*y\
[73]*y[74]*y[75]+2.*x2*y[1]*y[2]*y[3]*y[13]*y[18]*y[35]*y[36]*y[49]*y[55]*y\
[56]*y[73]*y[74]*y[75]+x2*y[1]*y[2]*y[3]*y[13]*y[67]*y[68]*y[69]*y[73]*y[74\
]*y[75]),2);
return (FOUT);
}
double Ps14(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double x3=x[3];
double y[76];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=esx[0];
y[5]=y[1]*y[2];
y[6]=2.*x0*y[1]*y[2];
y[7]=2.*x1*y[1]*y[2];
y[8]=2.*x0*x1*y[1]*y[2];
y[9]=x1*x1;
y[10]=y[1]*y[2]*y[9];
y[11]=2.*x0*x3*y[1]*y[2];
y[12]=2.*x0*x1*x3*y[1]*y[2];
y[13]=x3*x3;
y[14]=-(x0*y[1]*y[4]);
y[15]=-(x1*y[1]*y[4]);
y[16]=2.*x0*x2*y[1]*y[2];
y[17]=-x0;
y[18]=1.+y[17];
y[19]=x1*y[1]*y[2];
y[20]=2.*x2*y[1]*y[2];
y[21]=2.*x1*x2*y[1]*y[2];
y[22]=2.*x3*y[1]*y[2];
y[23]=2.*x1*x3*y[1]*y[2];
y[24]=2.*x2*x3*y[1]*y[2];
y[25]=4.*x0*x2*x3*y[1]*y[2];
y[26]=2.*x1*x2*x3*y[1]*y[2];
y[27]=y[1]*y[2]*y[13];
y[28]=x1*y[1]*y[2]*y[13];
y[29]=2.*x0*x2*y[1]*y[2]*y[13];
y[30]=-(x2*y[1]*y[4]);
y[31]=-(x3*y[1]*y[4]);
y[32]=-(x1*x3*y[1]*y[4]);
y[33]=-2.*x0*x2*x3*y[1]*y[4];
y[34]=-(x1*x2*x3*y[1]*y[4]);
y[35]=y[5]+y[16]+y[19]+y[20]+y[21]+y[22]+y[23]+y[24]+y[25]+y[26]+y[27]+y[28]\
+y[29]+y[30]+y[31]+y[32]+y[33]+y[34];
y[36]=lrs[0];
y[37]=-x1;
y[38]=1.+y[37];
y[39]=2.*y[1]*y[2];
y[40]=x0*y[1]*y[2];
y[41]=2.*x0*x2*x3*y[1]*y[2];
y[42]=x0*y[1]*y[2]*y[13];
y[43]=-(y[1]*y[4]);
y[44]=-(x0*x3*y[1]*y[4]);
y[45]=-(x0*x2*x3*y[1]*y[4]);
y[46]=y[7]+y[11]+y[16]+y[20]+y[21]+y[22]+y[23]+y[30]+y[31]+y[39]+y[40]+y[41]\
+y[42]+y[43]+y[44]+y[45];
y[47]=lrs[1];
y[48]=-x2;
y[49]=1.+y[48];
y[50]=y[1]*y[2]*y[3];
y[51]=2.*x3*y[1]*y[2]*y[3];
y[52]=y[1]*y[2]*y[3]*y[13];
y[53]=-(x3*y[1]*y[3]*y[4]);
y[54]=-(x0*x1*x3*y[1]*y[4]);
y[55]=y[5]+y[6]+y[7]+y[8]+y[10]+y[11]+y[12]+y[14]+y[15]+y[50]+y[51]+y[52]+y[\
53]+y[54];
y[56]=lrs[2];
y[57]=-x3;
y[58]=1.+y[57];
y[59]=2.*x2*y[1]*y[2]*y[3];
y[60]=2.*x0*x1*x2*y[1]*y[2];
y[61]=2.*x2*x3*y[1]*y[2]*y[3];
y[62]=-(x0*x1*y[1]*y[4]);
y[63]=-(x2*y[1]*y[3]*y[4]);
y[64]=-(x0*x1*x2*y[1]*y[4]);
y[65]=y[5]+y[6]+y[7]+y[8]+y[10]+y[11]+y[12]+y[14]+y[15]+y[16]+y[59]+y[60]+y[\
61]+y[62]+y[63]+y[64];
y[66]=lrs[3];
y[67]=pow(y[18],2);
y[68]=pow(y[35],2);
y[69]=pow(y[36],2);
y[70]=pow(y[38],2);
y[71]=pow(y[46],2);
y[72]=pow(y[47],2);
y[73]=pow(y[58],2);
y[74]=pow(y[65],2);
y[75]=pow(y[66],2);
FOUT=lambda*(-(x0*y[1]*y[2]*y[18]*y[35]*y[36])-x0*x1*y[1]*y[2]*y[18]*y[35]*y\
[36]-2.*x0*x2*y[1]*y[2]*y[18]*y[35]*y[36]-2.*x0*x1*x2*y[1]*y[2]*y[18]*y[35]\
*y[36]-2.*x0*x3*y[1]*y[2]*y[18]*y[35]*y[36]-2.*x0*x1*x3*y[1]*y[2]*y[18]*y[3\
5]*y[36]-2.*x0*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]-2.*x0*x1*x2*x3*y[1]*y[2]*y\
[18]*y[35]*y[36]-2.*x2*y[1]*y[2]*y[3]*y[18]*y[35]*y[36]-4.*x2*x3*y[1]*y[2]*\
y[3]*y[18]*y[35]*y[36]+x0*x2*y[1]*y[4]*y[18]*y[35]*y[36]+x0*x3*y[1]*y[4]*y[\
18]*y[35]*y[36]+x0*x1*x3*y[1]*y[4]*y[18]*y[35]*y[36]+x0*x1*x2*x3*y[1]*y[4]*\
y[18]*y[35]*y[36]+2.*x2*x3*y[1]*y[3]*y[4]*y[18]*y[35]*y[36]-x0*y[1]*y[2]*y[\
13]*y[18]*y[35]*y[36]-x0*x1*y[1]*y[2]*y[13]*y[18]*y[35]*y[36]-2.*x2*y[1]*y[\
2]*y[3]*y[13]*y[18]*y[35]*y[36]-2.*x1*y[1]*y[2]*y[38]*y[46]*y[47]-x0*x1*y[1\
]*y[2]*y[38]*y[46]*y[47]-2.*x1*x2*y[1]*y[2]*y[38]*y[46]*y[47]-2.*x0*x1*x2*y\
[1]*y[2]*y[38]*y[46]*y[47]-2.*x1*x3*y[1]*y[2]*y[38]*y[46]*y[47]-2.*x0*x1*x3\
*y[1]*y[2]*y[38]*y[46]*y[47]-2.*x0*x1*x2*x3*y[1]*y[2]*y[38]*y[46]*y[47]+x1*\
y[1]*y[4]*y[38]*y[46]*y[47]+x1*x2*y[1]*y[4]*y[38]*y[46]*y[47]+x1*x3*y[1]*y[\
4]*y[38]*y[46]*y[47]+x0*x1*x3*y[1]*y[4]*y[38]*y[46]*y[47]+x0*x1*x2*x3*y[1]*\
y[4]*y[38]*y[46]*y[47]-2.*y[1]*y[2]*y[9]*y[38]*y[46]*y[47]-2.*x2*y[1]*y[2]*\
y[9]*y[38]*y[46]*y[47]-2.*x3*y[1]*y[2]*y[9]*y[38]*y[46]*y[47]-x0*x1*y[1]*y[\
2]*y[13]*y[38]*y[46]*y[47]-x2*y[1]*y[2]*y[49]*y[55]*y[56]-2.*x0*x2*y[1]*y[2\
]*y[49]*y[55]*y[56]-2.*x1*x2*y[1]*y[2]*y[49]*y[55]*y[56]-2.*x0*x1*x2*y[1]*y\
[2]*y[49]*y[55]*y[56]-2.*x0*x2*x3*y[1]*y[2]*y[49]*y[55]*y[56]-2.*x0*x1*x2*x\
3*y[1]*y[2]*y[49]*y[55]*y[56]-x2*y[1]*y[2]*y[3]*y[49]*y[55]*y[56]-2.*x2*x3*\
y[1]*y[2]*y[3]*y[49]*y[55]*y[56]+x0*x2*y[1]*y[4]*y[49]*y[55]*y[56]+x1*x2*y[\
1]*y[4]*y[49]*y[55]*y[56]+x0*x1*x2*x3*y[1]*y[4]*y[49]*y[55]*y[56]+x2*x3*y[1\
]*y[3]*y[4]*y[49]*y[55]*y[56]-x2*y[1]*y[2]*y[9]*y[49]*y[55]*y[56]-x2*y[1]*y\
[2]*y[3]*y[13]*y[49]*y[55]*y[56]-x3*y[1]*y[2]*y[58]*y[65]*y[66]-2.*x0*x3*y[\
1]*y[2]*y[58]*y[65]*y[66]-2.*x1*x3*y[1]*y[2]*y[58]*y[65]*y[66]-2.*x0*x1*x3*\
y[1]*y[2]*y[58]*y[65]*y[66]-2.*x0*x2*x3*y[1]*y[2]*y[58]*y[65]*y[66]-2.*x0*x\
1*x2*x3*y[1]*y[2]*y[58]*y[65]*y[66]-2.*x2*x3*y[1]*y[2]*y[3]*y[58]*y[65]*y[6\
6]+x0*x3*y[1]*y[4]*y[58]*y[65]*y[66]+x1*x3*y[1]*y[4]*y[58]*y[65]*y[66]+x0*x\
1*x3*y[1]*y[4]*y[58]*y[65]*y[66]+x0*x1*x2*x3*y[1]*y[4]*y[58]*y[65]*y[66]+x2\
*x3*y[1]*y[3]*y[4]*y[58]*y[65]*y[66]-x3*y[1]*y[2]*y[9]*y[58]*y[65]*y[66]-2.\
*x0*y[1]*y[2]*y[13]*y[58]*y[65]*y[66]-2.*x0*x1*y[1]*y[2]*y[13]*y[58]*y[65]*\
y[66]-2.*x2*y[1]*y[2]*y[3]*y[13]*y[58]*y[65]*y[66])-x2*pow(lambda,5)*y[1]*y\
[2]*y[3]*y[13]*y[49]*y[55]*y[56]*y[67]*y[68]*y[69]*y[73]*y[74]*y[75]+pow(la\
mbda,3)*(2.*x0*x1*x2*y[1]*y[2]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]*y[49]*y[\
55]*y[56]+2.*x0*x1*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]*y[49\
]*y[55]*y[56]-x0*x1*x2*x3*y[1]*y[4]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]*y[4\
9]*y[55]*y[56]+2.*x0*x1*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]*y[\
58]*y[65]*y[66]+2.*x0*x1*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47\
]*y[58]*y[65]*y[66]-x0*x1*x3*y[1]*y[4]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]*\
y[58]*y[65]*y[66]-x0*x1*x2*x3*y[1]*y[4]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]\
*y[58]*y[65]*y[66]+2.*x0*x1*y[1]*y[2]*y[13]*y[18]*y[35]*y[36]*y[38]*y[46]*y\
[47]*y[58]*y[65]*y[66]+2.*x0*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[49]*y[55]*\
y[56]*y[58]*y[65]*y[66]+2.*x0*x1*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[49]*y[\
55]*y[56]*y[58]*y[65]*y[66]+4.*x2*x3*y[1]*y[2]*y[3]*y[18]*y[35]*y[36]*y[49]\
*y[55]*y[56]*y[58]*y[65]*y[66]-x0*x1*x2*x3*y[1]*y[4]*y[18]*y[35]*y[36]*y[49\
]*y[55]*y[56]*y[58]*y[65]*y[66]-2.*x2*x3*y[1]*y[3]*y[4]*y[18]*y[35]*y[36]*y\
[49]*y[55]*y[56]*y[58]*y[65]*y[66]+4.*x2*y[1]*y[2]*y[3]*y[13]*y[18]*y[35]*y\
[36]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]+2.*x0*x1*x2*x3*y[1]*y[2]*y[38]*y[4\
6]*y[47]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]-x0*x1*x2*x3*y[1]*y[4]*y[38]*y[\
46]*y[47]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]+x2*y[1]*y[2]*y[3]*y[49]*y[55]\
*y[56]*y[67]*y[68]*y[69]+2.*x2*x3*y[1]*y[2]*y[3]*y[49]*y[55]*y[56]*y[67]*y[\
68]*y[69]-x2*x3*y[1]*y[3]*y[4]*y[49]*y[55]*y[56]*y[67]*y[68]*y[69]+x2*y[1]*\
y[2]*y[3]*y[13]*y[49]*y[55]*y[56]*y[67]*y[68]*y[69]+2.*x2*x3*y[1]*y[2]*y[3]\
*y[58]*y[65]*y[66]*y[67]*y[68]*y[69]-x2*x3*y[1]*y[3]*y[4]*y[58]*y[65]*y[66]\
*y[67]*y[68]*y[69]+2.*x2*y[1]*y[2]*y[3]*y[13]*y[58]*y[65]*y[66]*y[67]*y[68]\
*y[69]+x2*y[1]*y[2]*y[9]*y[49]*y[55]*y[56]*y[70]*y[71]*y[72]+x3*y[1]*y[2]*y\
[9]*y[58]*y[65]*y[66]*y[70]*y[71]*y[72]+x0*y[1]*y[2]*y[13]*y[18]*y[35]*y[36\
]*y[73]*y[74]*y[75]+x0*x1*y[1]*y[2]*y[13]*y[18]*y[35]*y[36]*y[73]*y[74]*y[7\
5]+2.*x2*y[1]*y[2]*y[3]*y[13]*y[18]*y[35]*y[36]*y[73]*y[74]*y[75]+x0*x1*y[1\
]*y[2]*y[13]*y[38]*y[46]*y[47]*y[73]*y[74]*y[75]+x2*y[1]*y[2]*y[3]*y[13]*y[\
49]*y[55]*y[56]*y[73]*y[74]*y[75]);
return (FOUT);
}
double Pa14(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double x3=x[3];
double y[76];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=esx[0];
y[5]=y[1]*y[2];
y[6]=2.*x0*y[1]*y[2];
y[7]=2.*x1*y[1]*y[2];
y[8]=2.*x0*x1*y[1]*y[2];
y[9]=x1*x1;
y[10]=y[1]*y[2]*y[9];
y[11]=2.*x0*x3*y[1]*y[2];
y[12]=2.*x0*x1*x3*y[1]*y[2];
y[13]=x3*x3;
y[14]=-(x0*y[1]*y[4]);
y[15]=-(x1*y[1]*y[4]);
y[16]=2.*x0*x2*y[1]*y[2];
y[17]=-x0;
y[18]=1.+y[17];
y[19]=x1*y[1]*y[2];
y[20]=2.*x2*y[1]*y[2];
y[21]=2.*x1*x2*y[1]*y[2];
y[22]=2.*x3*y[1]*y[2];
y[23]=2.*x1*x3*y[1]*y[2];
y[24]=2.*x2*x3*y[1]*y[2];
y[25]=4.*x0*x2*x3*y[1]*y[2];
y[26]=2.*x1*x2*x3*y[1]*y[2];
y[27]=y[1]*y[2]*y[13];
y[28]=x1*y[1]*y[2]*y[13];
y[29]=2.*x0*x2*y[1]*y[2]*y[13];
y[30]=-(x2*y[1]*y[4]);
y[31]=-(x3*y[1]*y[4]);
y[32]=-(x1*x3*y[1]*y[4]);
y[33]=-2.*x0*x2*x3*y[1]*y[4];
y[34]=-(x1*x2*x3*y[1]*y[4]);
y[35]=y[5]+y[16]+y[19]+y[20]+y[21]+y[22]+y[23]+y[24]+y[25]+y[26]+y[27]+y[28]\
+y[29]+y[30]+y[31]+y[32]+y[33]+y[34];
y[36]=lrs[0];
y[37]=-x1;
y[38]=1.+y[37];
y[39]=2.*y[1]*y[2];
y[40]=x0*y[1]*y[2];
y[41]=2.*x0*x2*x3*y[1]*y[2];
y[42]=x0*y[1]*y[2]*y[13];
y[43]=-(y[1]*y[4]);
y[44]=-(x0*x3*y[1]*y[4]);
y[45]=-(x0*x2*x3*y[1]*y[4]);
y[46]=y[7]+y[11]+y[16]+y[20]+y[21]+y[22]+y[23]+y[30]+y[31]+y[39]+y[40]+y[41]\
+y[42]+y[43]+y[44]+y[45];
y[47]=lrs[1];
y[48]=-x2;
y[49]=1.+y[48];
y[50]=y[1]*y[2]*y[3];
y[51]=2.*x3*y[1]*y[2]*y[3];
y[52]=y[1]*y[2]*y[3]*y[13];
y[53]=-(x3*y[1]*y[3]*y[4]);
y[54]=-(x0*x1*x3*y[1]*y[4]);
y[55]=y[5]+y[6]+y[7]+y[8]+y[10]+y[11]+y[12]+y[14]+y[15]+y[50]+y[51]+y[52]+y[\
53]+y[54];
y[56]=lrs[2];
y[57]=-x3;
y[58]=1.+y[57];
y[59]=2.*x2*y[1]*y[2]*y[3];
y[60]=2.*x0*x1*x2*y[1]*y[2];
y[61]=2.*x2*x3*y[1]*y[2]*y[3];
y[62]=-(x0*x1*y[1]*y[4]);
y[63]=-(x2*y[1]*y[3]*y[4]);
y[64]=-(x0*x1*x2*y[1]*y[4]);
y[65]=y[5]+y[6]+y[7]+y[8]+y[10]+y[11]+y[12]+y[14]+y[15]+y[16]+y[59]+y[60]+y[\
61]+y[62]+y[63]+y[64];
y[66]=lrs[3];
y[67]=pow(y[18],2);
y[68]=pow(y[35],2);
y[69]=pow(y[36],2);
y[70]=pow(y[38],2);
y[71]=pow(y[46],2);
y[72]=pow(y[47],2);
y[73]=pow(y[58],2);
y[74]=pow(y[65],2);
y[75]=pow(y[66],2);
FOUT=(lambda*(-(x0*y[1]*y[2]*y[18]*y[35]*y[36])-x0*x1*y[1]*y[2]*y[18]*y[35]*\
y[36]-2.*x0*x2*y[1]*y[2]*y[18]*y[35]*y[36]-2.*x0*x1*x2*y[1]*y[2]*y[18]*y[35\
]*y[36]-2.*x0*x3*y[1]*y[2]*y[18]*y[35]*y[36]-2.*x0*x1*x3*y[1]*y[2]*y[18]*y[\
35]*y[36]-2.*x0*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]-2.*x0*x1*x2*x3*y[1]*y[2]*\
y[18]*y[35]*y[36]-2.*x2*y[1]*y[2]*y[3]*y[18]*y[35]*y[36]-4.*x2*x3*y[1]*y[2]\
*y[3]*y[18]*y[35]*y[36]+x0*x2*y[1]*y[4]*y[18]*y[35]*y[36]+x0*x3*y[1]*y[4]*y\
[18]*y[35]*y[36]+x0*x1*x3*y[1]*y[4]*y[18]*y[35]*y[36]+x0*x1*x2*x3*y[1]*y[4]\
*y[18]*y[35]*y[36]+2.*x2*x3*y[1]*y[3]*y[4]*y[18]*y[35]*y[36]-x0*y[1]*y[2]*y\
[13]*y[18]*y[35]*y[36]-x0*x1*y[1]*y[2]*y[13]*y[18]*y[35]*y[36]-2.*x2*y[1]*y\
[2]*y[3]*y[13]*y[18]*y[35]*y[36]-2.*x1*y[1]*y[2]*y[38]*y[46]*y[47]-x0*x1*y[\
1]*y[2]*y[38]*y[46]*y[47]-2.*x1*x2*y[1]*y[2]*y[38]*y[46]*y[47]-2.*x0*x1*x2*\
y[1]*y[2]*y[38]*y[46]*y[47]-2.*x1*x3*y[1]*y[2]*y[38]*y[46]*y[47]-2.*x0*x1*x\
3*y[1]*y[2]*y[38]*y[46]*y[47]-2.*x0*x1*x2*x3*y[1]*y[2]*y[38]*y[46]*y[47]+x1\
*y[1]*y[4]*y[38]*y[46]*y[47]+x1*x2*y[1]*y[4]*y[38]*y[46]*y[47]+x1*x3*y[1]*y\
[4]*y[38]*y[46]*y[47]+x0*x1*x3*y[1]*y[4]*y[38]*y[46]*y[47]+x0*x1*x2*x3*y[1]\
*y[4]*y[38]*y[46]*y[47]-2.*y[1]*y[2]*y[9]*y[38]*y[46]*y[47]-2.*x2*y[1]*y[2]\
*y[9]*y[38]*y[46]*y[47]-2.*x3*y[1]*y[2]*y[9]*y[38]*y[46]*y[47]-x0*x1*y[1]*y\
[2]*y[13]*y[38]*y[46]*y[47]-x2*y[1]*y[2]*y[49]*y[55]*y[56]-2.*x0*x2*y[1]*y[\
2]*y[49]*y[55]*y[56]-2.*x1*x2*y[1]*y[2]*y[49]*y[55]*y[56]-2.*x0*x1*x2*y[1]*\
y[2]*y[49]*y[55]*y[56]-2.*x0*x2*x3*y[1]*y[2]*y[49]*y[55]*y[56]-2.*x0*x1*x2*\
x3*y[1]*y[2]*y[49]*y[55]*y[56]-x2*y[1]*y[2]*y[3]*y[49]*y[55]*y[56]-2.*x2*x3\
*y[1]*y[2]*y[3]*y[49]*y[55]*y[56]+x0*x2*y[1]*y[4]*y[49]*y[55]*y[56]+x1*x2*y\
[1]*y[4]*y[49]*y[55]*y[56]+x0*x1*x2*x3*y[1]*y[4]*y[49]*y[55]*y[56]+x2*x3*y[\
1]*y[3]*y[4]*y[49]*y[55]*y[56]-x2*y[1]*y[2]*y[9]*y[49]*y[55]*y[56]-x2*y[1]*\
y[2]*y[3]*y[13]*y[49]*y[55]*y[56]-x3*y[1]*y[2]*y[58]*y[65]*y[66]-2.*x0*x3*y\
[1]*y[2]*y[58]*y[65]*y[66]-2.*x1*x3*y[1]*y[2]*y[58]*y[65]*y[66]-2.*x0*x1*x3\
*y[1]*y[2]*y[58]*y[65]*y[66]-2.*x0*x2*x3*y[1]*y[2]*y[58]*y[65]*y[66]-2.*x0*\
x1*x2*x3*y[1]*y[2]*y[58]*y[65]*y[66]-2.*x2*x3*y[1]*y[2]*y[3]*y[58]*y[65]*y[\
66]+x0*x3*y[1]*y[4]*y[58]*y[65]*y[66]+x1*x3*y[1]*y[4]*y[58]*y[65]*y[66]+x0*\
x1*x3*y[1]*y[4]*y[58]*y[65]*y[66]+x0*x1*x2*x3*y[1]*y[4]*y[58]*y[65]*y[66]+x\
2*x3*y[1]*y[3]*y[4]*y[58]*y[65]*y[66]-x3*y[1]*y[2]*y[9]*y[58]*y[65]*y[66]-2\
.*x0*y[1]*y[2]*y[13]*y[58]*y[65]*y[66]-2.*x0*x1*y[1]*y[2]*y[13]*y[58]*y[65]\
*y[66]-2.*x2*y[1]*y[2]*y[3]*y[13]*y[58]*y[65]*y[66])-x2*pow(lambda,5)*y[1]*\
y[2]*y[3]*y[13]*y[49]*y[55]*y[56]*y[67]*y[68]*y[69]*y[73]*y[74]*y[75]+pow(l\
ambda,3)*(2.*x0*x1*x2*y[1]*y[2]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]*y[49]*y\
[55]*y[56]+2.*x0*x1*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]*y[4\
9]*y[55]*y[56]-x0*x1*x2*x3*y[1]*y[4]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]*y[\
49]*y[55]*y[56]+2.*x0*x1*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]*y\
[58]*y[65]*y[66]+2.*x0*x1*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[38]*y[46]*y[4\
7]*y[58]*y[65]*y[66]-x0*x1*x3*y[1]*y[4]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]\
*y[58]*y[65]*y[66]-x0*x1*x2*x3*y[1]*y[4]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47\
]*y[58]*y[65]*y[66]+2.*x0*x1*y[1]*y[2]*y[13]*y[18]*y[35]*y[36]*y[38]*y[46]*\
y[47]*y[58]*y[65]*y[66]+2.*x0*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[49]*y[55]\
*y[56]*y[58]*y[65]*y[66]+2.*x0*x1*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[49]*y\
[55]*y[56]*y[58]*y[65]*y[66]+4.*x2*x3*y[1]*y[2]*y[3]*y[18]*y[35]*y[36]*y[49\
]*y[55]*y[56]*y[58]*y[65]*y[66]-x0*x1*x2*x3*y[1]*y[4]*y[18]*y[35]*y[36]*y[4\
9]*y[55]*y[56]*y[58]*y[65]*y[66]-2.*x2*x3*y[1]*y[3]*y[4]*y[18]*y[35]*y[36]*\
y[49]*y[55]*y[56]*y[58]*y[65]*y[66]+4.*x2*y[1]*y[2]*y[3]*y[13]*y[18]*y[35]*\
y[36]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]+2.*x0*x1*x2*x3*y[1]*y[2]*y[38]*y[\
46]*y[47]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]-x0*x1*x2*x3*y[1]*y[4]*y[38]*y\
[46]*y[47]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]+x2*y[1]*y[2]*y[3]*y[49]*y[55\
]*y[56]*y[67]*y[68]*y[69]+2.*x2*x3*y[1]*y[2]*y[3]*y[49]*y[55]*y[56]*y[67]*y\
[68]*y[69]-x2*x3*y[1]*y[3]*y[4]*y[49]*y[55]*y[56]*y[67]*y[68]*y[69]+x2*y[1]\
*y[2]*y[3]*y[13]*y[49]*y[55]*y[56]*y[67]*y[68]*y[69]+2.*x2*x3*y[1]*y[2]*y[3\
]*y[58]*y[65]*y[66]*y[67]*y[68]*y[69]-x2*x3*y[1]*y[3]*y[4]*y[58]*y[65]*y[66\
]*y[67]*y[68]*y[69]+2.*x2*y[1]*y[2]*y[3]*y[13]*y[58]*y[65]*y[66]*y[67]*y[68\
]*y[69]+x2*y[1]*y[2]*y[9]*y[49]*y[55]*y[56]*y[70]*y[71]*y[72]+x3*y[1]*y[2]*\
y[9]*y[58]*y[65]*y[66]*y[70]*y[71]*y[72]+x0*y[1]*y[2]*y[13]*y[18]*y[35]*y[3\
6]*y[73]*y[74]*y[75]+x0*x1*y[1]*y[2]*y[13]*y[18]*y[35]*y[36]*y[73]*y[74]*y[\
75]+2.*x2*y[1]*y[2]*y[3]*y[13]*y[18]*y[35]*y[36]*y[73]*y[74]*y[75]+x0*x1*y[\
1]*y[2]*y[13]*y[38]*y[46]*y[47]*y[73]*y[74]*y[75]+x2*y[1]*y[2]*y[3]*y[13]*y\
[49]*y[55]*y[56]*y[73]*y[74]*y[75]))/(lambda*(x0*x1*y[1]*y[2]+x2*y[1]*y[2]+\
x3*y[1]*y[2]+2.*x0*x1*x2*x3*y[1]*y[2]+x2*y[1]*y[2]*y[3]-x0*x2*y[1]*y[4]-x1*\
x2*y[1]*y[4]-x0*x1*x2*x3*y[1]*y[4]-x2*x3*y[1]*y[3]*y[4]+y[5]+y[7]+x2*y[1]*y\
[2]*y[9]+x3*y[1]*y[2]*y[9]+y[10]+y[11]+y[12]+x0*x1*y[1]*y[2]*y[13]+x2*y[1]*\
y[2]*y[3]*y[13]+y[15]+y[16]+y[21]+y[23]+y[32]+y[40]+y[41]+y[42]+y[44]+y[54]\
+y[60]+y[61]+lambda*lambda*(-(x0*x1*y[1]*y[2]*y[18]*y[35]*y[36]*y[38]*y[46]\
*y[47])-2.*x0*x1*x2*y[1]*y[2]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]-2.*x0*x1*\
x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]-2.*x0*x1*x2*x3*y[1]*y[2]*y\
[18]*y[35]*y[36]*y[38]*y[46]*y[47]+x0*x1*x3*y[1]*y[4]*y[18]*y[35]*y[36]*y[3\
8]*y[46]*y[47]+x0*x1*x2*x3*y[1]*y[4]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]-x0\
*x1*y[1]*y[2]*y[13]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]-2.*x0*x2*y[1]*y[2]*\
y[18]*y[35]*y[36]*y[49]*y[55]*y[56]-2.*x0*x1*x2*y[1]*y[2]*y[18]*y[35]*y[36]\
*y[49]*y[55]*y[56]-2.*x0*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[49]*y[55]*y[56\
]-2.*x0*x1*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[49]*y[55]*y[56]-2.*x2*y[1]*y\
[2]*y[3]*y[18]*y[35]*y[36]*y[49]*y[55]*y[56]-4.*x2*x3*y[1]*y[2]*y[3]*y[18]*\
y[35]*y[36]*y[49]*y[55]*y[56]+x0*x2*y[1]*y[4]*y[18]*y[35]*y[36]*y[49]*y[55]\
*y[56]+x0*x1*x2*x3*y[1]*y[4]*y[18]*y[35]*y[36]*y[49]*y[55]*y[56]+2.*x2*x3*y\
[1]*y[3]*y[4]*y[18]*y[35]*y[36]*y[49]*y[55]*y[56]-2.*x2*y[1]*y[2]*y[3]*y[13\
]*y[18]*y[35]*y[36]*y[49]*y[55]*y[56]-2.*x1*x2*y[1]*y[2]*y[38]*y[46]*y[47]*\
y[49]*y[55]*y[56]-2.*x0*x1*x2*y[1]*y[2]*y[38]*y[46]*y[47]*y[49]*y[55]*y[56]\
-2.*x0*x1*x2*x3*y[1]*y[2]*y[38]*y[46]*y[47]*y[49]*y[55]*y[56]+x1*x2*y[1]*y[\
4]*y[38]*y[46]*y[47]*y[49]*y[55]*y[56]+x0*x1*x2*x3*y[1]*y[4]*y[38]*y[46]*y[\
47]*y[49]*y[55]*y[56]-2.*x2*y[1]*y[2]*y[9]*y[38]*y[46]*y[47]*y[49]*y[55]*y[\
56]-2.*x0*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[58]*y[65]*y[66]-2.*x0*x1*x3*y[1]\
*y[2]*y[18]*y[35]*y[36]*y[58]*y[65]*y[66]-2.*x0*x2*x3*y[1]*y[2]*y[18]*y[35]\
*y[36]*y[58]*y[65]*y[66]-2.*x0*x1*x2*x3*y[1]*y[2]*y[18]*y[35]*y[36]*y[58]*y\
[65]*y[66]-4.*x2*x3*y[1]*y[2]*y[3]*y[18]*y[35]*y[36]*y[58]*y[65]*y[66]+x0*x\
3*y[1]*y[4]*y[18]*y[35]*y[36]*y[58]*y[65]*y[66]+x0*x1*x3*y[1]*y[4]*y[18]*y[\
35]*y[36]*y[58]*y[65]*y[66]+x0*x1*x2*x3*y[1]*y[4]*y[18]*y[35]*y[36]*y[58]*y\
[65]*y[66]+2.*x2*x3*y[1]*y[3]*y[4]*y[18]*y[35]*y[36]*y[58]*y[65]*y[66]-2.*x\
0*y[1]*y[2]*y[13]*y[18]*y[35]*y[36]*y[58]*y[65]*y[66]-2.*x0*x1*y[1]*y[2]*y[\
13]*y[18]*y[35]*y[36]*y[58]*y[65]*y[66]-4.*x2*y[1]*y[2]*y[3]*y[13]*y[18]*y[\
35]*y[36]*y[58]*y[65]*y[66]-2.*x1*x3*y[1]*y[2]*y[38]*y[46]*y[47]*y[58]*y[65\
]*y[66]-2.*x0*x1*x3*y[1]*y[2]*y[38]*y[46]*y[47]*y[58]*y[65]*y[66]-2.*x0*x1*\
x2*x3*y[1]*y[2]*y[38]*y[46]*y[47]*y[58]*y[65]*y[66]+x1*x3*y[1]*y[4]*y[38]*y\
[46]*y[47]*y[58]*y[65]*y[66]+x0*x1*x3*y[1]*y[4]*y[38]*y[46]*y[47]*y[58]*y[6\
5]*y[66]+x0*x1*x2*x3*y[1]*y[4]*y[38]*y[46]*y[47]*y[58]*y[65]*y[66]-2.*x3*y[\
1]*y[2]*y[9]*y[38]*y[46]*y[47]*y[58]*y[65]*y[66]-2.*x0*x1*y[1]*y[2]*y[13]*y\
[38]*y[46]*y[47]*y[58]*y[65]*y[66]-2.*x0*x2*x3*y[1]*y[2]*y[49]*y[55]*y[56]*\
y[58]*y[65]*y[66]-2.*x0*x1*x2*x3*y[1]*y[2]*y[49]*y[55]*y[56]*y[58]*y[65]*y[\
66]-2.*x2*x3*y[1]*y[2]*y[3]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]+x0*x1*x2*x3\
*y[1]*y[4]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]+x2*x3*y[1]*y[3]*y[4]*y[49]*y\
[55]*y[56]*y[58]*y[65]*y[66]-2.*x2*y[1]*y[2]*y[3]*y[13]*y[49]*y[55]*y[56]*y\
[58]*y[65]*y[66]-x2*y[1]*y[2]*y[3]*y[67]*y[68]*y[69]-2.*x2*x3*y[1]*y[2]*y[3\
]*y[67]*y[68]*y[69]+x2*x3*y[1]*y[3]*y[4]*y[67]*y[68]*y[69]-x2*y[1]*y[2]*y[3\
]*y[13]*y[67]*y[68]*y[69]-y[1]*y[2]*y[9]*y[70]*y[71]*y[72]-x2*y[1]*y[2]*y[9\
]*y[70]*y[71]*y[72]-x3*y[1]*y[2]*y[9]*y[70]*y[71]*y[72]-x0*y[1]*y[2]*y[13]*\
y[73]*y[74]*y[75]-x0*x1*y[1]*y[2]*y[13]*y[73]*y[74]*y[75]-x2*y[1]*y[2]*y[3]\
*y[13]*y[73]*y[74]*y[75])+pow(lambda,4)*(2.*x0*x1*x2*x3*y[1]*y[2]*y[18]*y[3\
5]*y[36]*y[38]*y[46]*y[47]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]-x0*x1*x2*x3*\
y[1]*y[4]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]*y[49]*y[55]*y[56]*y[58]*y[65]\
*y[66]+2.*x2*x3*y[1]*y[2]*y[3]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]*y[67]*y[\
68]*y[69]-x2*x3*y[1]*y[3]*y[4]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]*y[67]*y[\
68]*y[69]+2.*x2*y[1]*y[2]*y[3]*y[13]*y[49]*y[55]*y[56]*y[58]*y[65]*y[66]*y[\
67]*y[68]*y[69]+x0*x1*y[1]*y[2]*y[13]*y[18]*y[35]*y[36]*y[38]*y[46]*y[47]*y\
[73]*y[74]*y[75]+2.*x2*y[1]*y[2]*y[3]*y[13]*y[18]*y[35]*y[36]*y[49]*y[55]*y\
[56]*y[73]*y[74]*y[75]+x2*y[1]*y[2]*y[3]*y[13]*y[67]*y[68]*y[69]*y[73]*y[74\
]*y[75])));
return (FOUT);
}
double Pt14t1(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double x3=x[3];
double y[5];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x3*x3;
y[4]=esx[0];
FOUT=(1.-x0)*x0*(y[1]*y[2]+x1*y[1]*y[2]+2.*x2*y[1]*y[2]+2.*x0*x2*y[1]*y[2]+2\
.*x1*x2*y[1]*y[2]+2.*x3*y[1]*y[2]+2.*x1*x3*y[1]*y[2]+2.*x2*x3*y[1]*y[2]+4.*\
x0*x2*x3*y[1]*y[2]+2.*x1*x2*x3*y[1]*y[2]+y[1]*y[2]*y[3]+x1*y[1]*y[2]*y[3]+2\
.*x0*x2*y[1]*y[2]*y[3]-x2*y[1]*y[4]-x3*y[1]*y[4]-x1*x3*y[1]*y[4]-2.*x0*x2*x\
3*y[1]*y[4]-x1*x2*x3*y[1]*y[4]);
return (FOUT);
}
double Pt14t2(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double x3=x[3];
double y[4];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=esx[0];
FOUT=(1.-x1)*x1*(2.*y[1]*y[2]+x0*y[1]*y[2]+2.*x1*y[1]*y[2]+2.*x2*y[1]*y[2]+2\
.*x0*x2*y[1]*y[2]+2.*x1*x2*y[1]*y[2]+2.*x3*y[1]*y[2]+2.*x0*x3*y[1]*y[2]+2.*\
x1*x3*y[1]*y[2]+2.*x0*x2*x3*y[1]*y[2]+x0*(x3*x3)*y[1]*y[2]-y[1]*y[3]-x2*y[1\
]*y[3]-x3*y[1]*y[3]-x0*x3*y[1]*y[3]-x0*x2*x3*y[1]*y[3]);
return (FOUT);
}
double Pt14t3(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double x3=x[3];
double y[5];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=esx[0];
FOUT=(1.-x2)*x2*(y[1]*y[2]+2.*x0*y[1]*y[2]+2.*x1*y[1]*y[2]+2.*x0*x1*y[1]*y[2\
]+2.*x0*x3*y[1]*y[2]+2.*x0*x1*x3*y[1]*y[2]+x1*x1*y[1]*y[2]+y[1]*y[2]*y[3]+2\
.*x3*y[1]*y[2]*y[3]+x3*x3*y[1]*y[2]*y[3]-x0*y[1]*y[4]-x1*y[1]*y[4]-x0*x1*x3\
*y[1]*y[4]-x3*y[1]*y[3]*y[4]);
return (FOUT);
}
double Pt14t4(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double x3=x[3];
double y[5];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=esx[0];
FOUT=(1.-x3)*x3*(y[1]*y[2]+2.*x0*y[1]*y[2]+2.*x1*y[1]*y[2]+2.*x0*x1*y[1]*y[2\
]+2.*x0*x2*y[1]*y[2]+2.*x0*x1*x2*y[1]*y[2]+2.*x0*x3*y[1]*y[2]+2.*x0*x1*x3*y\
[1]*y[2]+x1*x1*y[1]*y[2]+2.*x2*y[1]*y[2]*y[3]+2.*x2*x3*y[1]*y[2]*y[3]-x0*y[\
1]*y[4]-x1*y[1]*y[4]-x0*x1*y[1]*y[4]-x0*x1*x2*y[1]*y[4]-x2*y[1]*y[3]*y[4]);
return (FOUT);
}