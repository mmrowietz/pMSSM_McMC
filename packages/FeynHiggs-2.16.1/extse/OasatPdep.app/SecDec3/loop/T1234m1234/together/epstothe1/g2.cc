#include "intfile.hh"

double Pr2(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[68];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=em[1];
y[4]=em[2];
y[5]=em[3];
y[6]=x2*x2;
y[7]=esx[0];
y[8]=-x0;
y[9]=1.+y[8];
y[10]=y[1]*y[2];
y[11]=x1*y[1]*y[2];
y[12]=x2*y[1]*y[2];
y[13]=x1*y[1]*y[3];
y[14]=x1*x1;
y[15]=y[1]*y[3]*y[14];
y[16]=x1*x2*y[1]*y[3];
y[17]=x2*y[1]*y[4];
y[18]=2.*x0*x2*y[1]*y[4];
y[19]=x1*x2*y[1]*y[4];
y[20]=2.*x0*x1*x2*y[1]*y[4];
y[21]=2.*x0*y[1]*y[4]*y[6];
y[22]=x2*y[1]*y[5];
y[23]=x1*x2*y[1]*y[5];
y[24]=y[1]*y[5]*y[6];
y[25]=-(x1*y[1]*y[7]);
y[26]=-(x1*x2*y[1]*y[7]);
y[27]=y[10]+y[11]+y[12]+y[13]+y[15]+y[16]+y[17]+y[18]+y[19]+y[20]+y[21]+y[22\
]+y[23]+y[24]+y[25]+y[26];
y[28]=lrs[0];
y[29]=x0*x0;
y[30]=-x1;
y[31]=1.+y[30];
y[32]=x0*y[1]*y[2];
y[33]=y[1]*y[3];
y[34]=x0*y[1]*y[3];
y[35]=2.*x1*y[1]*y[3];
y[36]=2.*x0*x1*y[1]*y[3];
y[37]=x0*x2*y[1]*y[3];
y[38]=x0*x2*y[1]*y[4];
y[39]=x2*y[1]*y[4]*y[29];
y[40]=x0*x2*y[1]*y[5];
y[41]=-(y[1]*y[7]);
y[42]=-(x0*y[1]*y[7]);
y[43]=-(x0*x2*y[1]*y[7]);
y[44]=y[10]+y[22]+y[32]+y[33]+y[34]+y[35]+y[36]+y[37]+y[38]+y[39]+y[40]+y[41\
]+y[42]+y[43];
y[45]=lrs[1];
y[46]=-x2;
y[47]=1.+y[46];
y[48]=x0*x1*y[1]*y[3];
y[49]=x0*y[1]*y[4];
y[50]=y[1]*y[4]*y[29];
y[51]=x0*x1*y[1]*y[4];
y[52]=x1*y[1]*y[4]*y[29];
y[53]=2.*x2*y[1]*y[4]*y[29];
y[54]=y[1]*y[5];
y[55]=x0*y[1]*y[5];
y[56]=x1*y[1]*y[5];
y[57]=x0*x1*y[1]*y[5];
y[58]=2.*x0*x2*y[1]*y[5];
y[59]=-(x0*x1*y[1]*y[7]);
y[60]=y[32]+y[48]+y[49]+y[50]+y[51]+y[52]+y[53]+y[54]+y[55]+y[56]+y[57]+y[58\
]+y[59];
y[61]=lrs[2];
y[62]=pow(y[9],2);
y[63]=pow(y[27],2);
y[64]=pow(y[28],2);
y[65]=pow(y[47],2);
y[66]=pow(y[60],2);
y[67]=pow(y[61],2);
FOUT=(x0*pow(y[31],2)*pow(y[44],2)*pow(y[45],2)*y[1]*y[3]*y[9]*y[14]*y[27]*y\
[28]+x0*x1*x2*y[1]*y[3]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61\
]+x0*x1*x2*y[1]*y[4]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]+x\
0*x1*x2*y[1]*y[5]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]-x0*x\
1*x2*y[1]*y[7]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]+2.*x1*x\
2*y[1]*y[4]*y[9]*y[27]*y[28]*y[29]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]+x1*x\
2*y[1]*y[4]*y[29]*y[31]*y[44]*y[45]*y[62]*y[63]*y[64]+x2*y[1]*y[4]*y[29]*y[\
47]*y[60]*y[61]*y[62]*y[63]*y[64]+x1*x2*y[1]*y[4]*y[29]*y[47]*y[60]*y[61]*y\
[62]*y[63]*y[64]+2.*y[1]*y[4]*y[6]*y[29]*y[47]*y[60]*y[61]*y[62]*y[63]*y[64\
]+x0*y[1]*y[5]*y[6]*y[9]*y[27]*y[28]*y[65]*y[66]*y[67]+2.*y[1]*y[4]*y[6]*y[\
9]*y[27]*y[28]*y[29]*y[65]*y[66]*y[67])/(-(x0*y[1]*y[2]*y[9]*y[27]*y[28])-x\
0*x1*y[1]*y[2]*y[9]*y[27]*y[28]-x0*x2*y[1]*y[2]*y[9]*y[27]*y[28]-x0*x1*y[1]\
*y[3]*y[9]*y[27]*y[28]-x0*x1*x2*y[1]*y[3]*y[9]*y[27]*y[28]-x0*x2*y[1]*y[4]*\
y[9]*y[27]*y[28]-x0*x1*x2*y[1]*y[4]*y[9]*y[27]*y[28]-x0*x2*y[1]*y[5]*y[9]*y\
[27]*y[28]-x0*x1*x2*y[1]*y[5]*y[9]*y[27]*y[28]-x0*y[1]*y[5]*y[6]*y[9]*y[27]\
*y[28]+x0*x1*y[1]*y[7]*y[9]*y[27]*y[28]+x0*x1*x2*y[1]*y[7]*y[9]*y[27]*y[28]\
-x0*y[1]*y[3]*y[9]*y[14]*y[27]*y[28]-2.*x2*y[1]*y[4]*y[9]*y[27]*y[28]*y[29]\
-2.*x1*x2*y[1]*y[4]*y[9]*y[27]*y[28]*y[29]-2.*y[1]*y[4]*y[6]*y[9]*y[27]*y[2\
8]*y[29]-x1*y[1]*y[2]*y[31]*y[44]*y[45]-x0*x1*y[1]*y[2]*y[31]*y[44]*y[45]-x\
1*y[1]*y[3]*y[31]*y[44]*y[45]-x0*x1*y[1]*y[3]*y[31]*y[44]*y[45]-x0*x1*x2*y[\
1]*y[3]*y[31]*y[44]*y[45]-x0*x1*x2*y[1]*y[4]*y[31]*y[44]*y[45]-x1*x2*y[1]*y\
[5]*y[31]*y[44]*y[45]-x0*x1*x2*y[1]*y[5]*y[31]*y[44]*y[45]+x1*y[1]*y[7]*y[3\
1]*y[44]*y[45]+x0*x1*y[1]*y[7]*y[31]*y[44]*y[45]+x0*x1*x2*y[1]*y[7]*y[31]*y\
[44]*y[45]-2.*y[1]*y[3]*y[14]*y[31]*y[44]*y[45]-2.*x0*y[1]*y[3]*y[14]*y[31]\
*y[44]*y[45]-x1*x2*y[1]*y[4]*y[29]*y[31]*y[44]*y[45]-x0*x2*y[1]*y[2]*y[47]*\
y[60]*y[61]-x0*x1*x2*y[1]*y[3]*y[47]*y[60]*y[61]-x0*x2*y[1]*y[4]*y[47]*y[60\
]*y[61]-x0*x1*x2*y[1]*y[4]*y[47]*y[60]*y[61]-x2*y[1]*y[5]*y[47]*y[60]*y[61]\
-x0*x2*y[1]*y[5]*y[47]*y[60]*y[61]-x1*x2*y[1]*y[5]*y[47]*y[60]*y[61]-x0*x1*\
x2*y[1]*y[5]*y[47]*y[60]*y[61]-2.*x0*y[1]*y[5]*y[6]*y[47]*y[60]*y[61]+x0*x1\
*x2*y[1]*y[7]*y[47]*y[60]*y[61]-x2*y[1]*y[4]*y[29]*y[47]*y[60]*y[61]-x1*x2*\
y[1]*y[4]*y[29]*y[47]*y[60]*y[61]-2.*y[1]*y[4]*y[6]*y[29]*y[47]*y[60]*y[61]\
);
return (FOUT);
}
double Pm2(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[71];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=em[1];
y[4]=em[2];
y[5]=em[3];
y[6]=x2*x2;
y[7]=esx[0];
y[8]=-x0;
y[9]=1.+y[8];
y[10]=y[1]*y[2];
y[11]=x1*y[1]*y[2];
y[12]=x2*y[1]*y[2];
y[13]=x1*y[1]*y[3];
y[14]=x1*x1;
y[15]=y[1]*y[3]*y[14];
y[16]=x1*x2*y[1]*y[3];
y[17]=x2*y[1]*y[4];
y[18]=2.*x0*x2*y[1]*y[4];
y[19]=x1*x2*y[1]*y[4];
y[20]=2.*x0*x1*x2*y[1]*y[4];
y[21]=2.*x0*y[1]*y[4]*y[6];
y[22]=x2*y[1]*y[5];
y[23]=x1*x2*y[1]*y[5];
y[24]=y[1]*y[5]*y[6];
y[25]=-(x1*y[1]*y[7]);
y[26]=-(x1*x2*y[1]*y[7]);
y[27]=y[10]+y[11]+y[12]+y[13]+y[15]+y[16]+y[17]+y[18]+y[19]+y[20]+y[21]+y[22\
]+y[23]+y[24]+y[25]+y[26];
y[28]=lrs[0];
y[29]=x0*x0;
y[30]=-x1;
y[31]=1.+y[30];
y[32]=x0*y[1]*y[2];
y[33]=y[1]*y[3];
y[34]=x0*y[1]*y[3];
y[35]=2.*x1*y[1]*y[3];
y[36]=2.*x0*x1*y[1]*y[3];
y[37]=x0*x2*y[1]*y[3];
y[38]=x0*x2*y[1]*y[4];
y[39]=x2*y[1]*y[4]*y[29];
y[40]=x0*x2*y[1]*y[5];
y[41]=-(y[1]*y[7]);
y[42]=-(x0*y[1]*y[7]);
y[43]=-(x0*x2*y[1]*y[7]);
y[44]=y[10]+y[22]+y[32]+y[33]+y[34]+y[35]+y[36]+y[37]+y[38]+y[39]+y[40]+y[41\
]+y[42]+y[43];
y[45]=lrs[1];
y[46]=-x2;
y[47]=1.+y[46];
y[48]=x0*x1*y[1]*y[3];
y[49]=x0*y[1]*y[4];
y[50]=y[1]*y[4]*y[29];
y[51]=x0*x1*y[1]*y[4];
y[52]=x1*y[1]*y[4]*y[29];
y[53]=2.*x2*y[1]*y[4]*y[29];
y[54]=y[1]*y[5];
y[55]=x0*y[1]*y[5];
y[56]=x1*y[1]*y[5];
y[57]=x0*x1*y[1]*y[5];
y[58]=2.*x0*x2*y[1]*y[5];
y[59]=-(x0*x1*y[1]*y[7]);
y[60]=y[32]+y[48]+y[49]+y[50]+y[51]+y[52]+y[53]+y[54]+y[55]+y[56]+y[57]+y[58\
]+y[59];
y[61]=lrs[2];
y[62]=pow(y[9],2);
y[63]=pow(y[27],2);
y[64]=pow(y[28],2);
y[65]=pow(y[47],2);
y[66]=pow(y[60],2);
y[67]=pow(y[61],2);
y[68]=pow(y[31],2);
y[69]=pow(y[44],2);
y[70]=pow(y[45],2);
FOUT=pow(x0*x1*y[1]*y[2]+x0*x2*y[1]*y[2]+x0*x1*x2*y[1]*y[3]+x0*x1*x2*y[1]*y[\
4]+x0*x1*x2*y[1]*y[5]+x0*y[1]*y[5]*y[6]-x0*x1*x2*y[1]*y[7]+y[10]+y[11]+y[13\
]+x0*y[1]*y[3]*y[14]+y[15]+y[22]+y[23]+y[25]+x1*x2*y[1]*y[4]*y[29]+y[1]*y[4\
]*y[6]*y[29]+y[32]+y[38]+y[39]+y[40]+y[48]+y[59]+pow(lambda,4)*(x1*x2*y[1]*\
y[4]*y[29]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]*y[62]*y[63]*y[64]+y[1]*y[4]*\
y[6]*y[29]*y[62]*y[63]*y[64]*y[65]*y[66]*y[67])+lambda*lambda*(-(x0*x1*y[1]\
*y[2]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45])-x0*x1*y[1]*y[3]*y[9]*y[27]*y[28]*\
y[31]*y[44]*y[45]-x0*x1*x2*y[1]*y[3]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45]-x0*\
x1*x2*y[1]*y[4]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45]-x0*x1*x2*y[1]*y[5]*y[9]*\
y[27]*y[28]*y[31]*y[44]*y[45]+x0*x1*y[1]*y[7]*y[9]*y[27]*y[28]*y[31]*y[44]*\
y[45]+x0*x1*x2*y[1]*y[7]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45]-2.*x0*y[1]*y[3]\
*y[9]*y[14]*y[27]*y[28]*y[31]*y[44]*y[45]-2.*x1*x2*y[1]*y[4]*y[9]*y[27]*y[2\
8]*y[29]*y[31]*y[44]*y[45]-x0*x2*y[1]*y[2]*y[9]*y[27]*y[28]*y[47]*y[60]*y[6\
1]-x0*x1*x2*y[1]*y[3]*y[9]*y[27]*y[28]*y[47]*y[60]*y[61]-x0*x2*y[1]*y[4]*y[\
9]*y[27]*y[28]*y[47]*y[60]*y[61]-x0*x1*x2*y[1]*y[4]*y[9]*y[27]*y[28]*y[47]*\
y[60]*y[61]-x0*x2*y[1]*y[5]*y[9]*y[27]*y[28]*y[47]*y[60]*y[61]-x0*x1*x2*y[1\
]*y[5]*y[9]*y[27]*y[28]*y[47]*y[60]*y[61]-2.*x0*y[1]*y[5]*y[6]*y[9]*y[27]*y\
[28]*y[47]*y[60]*y[61]+x0*x1*x2*y[1]*y[7]*y[9]*y[27]*y[28]*y[47]*y[60]*y[61\
]-2.*x2*y[1]*y[4]*y[9]*y[27]*y[28]*y[29]*y[47]*y[60]*y[61]-2.*x1*x2*y[1]*y[\
4]*y[9]*y[27]*y[28]*y[29]*y[47]*y[60]*y[61]-4.*y[1]*y[4]*y[6]*y[9]*y[27]*y[\
28]*y[29]*y[47]*y[60]*y[61]-x0*x1*x2*y[1]*y[3]*y[31]*y[44]*y[45]*y[47]*y[60\
]*y[61]-x0*x1*x2*y[1]*y[4]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]-x1*x2*y[1]*y\
[5]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]-x0*x1*x2*y[1]*y[5]*y[31]*y[44]*y[45\
]*y[47]*y[60]*y[61]+x0*x1*x2*y[1]*y[7]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]-\
x1*x2*y[1]*y[4]*y[29]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]-x2*y[1]*y[4]*y[29\
]*y[62]*y[63]*y[64]-x1*x2*y[1]*y[4]*y[29]*y[62]*y[63]*y[64]-y[1]*y[4]*y[6]*\
y[29]*y[62]*y[63]*y[64]-x0*y[1]*y[5]*y[6]*y[65]*y[66]*y[67]-y[1]*y[4]*y[6]*\
y[29]*y[65]*y[66]*y[67]-y[1]*y[3]*y[14]*y[68]*y[69]*y[70]-x0*y[1]*y[3]*y[14\
]*y[68]*y[69]*y[70]),2)+pow(lambda*(-(x0*y[1]*y[2]*y[9]*y[27]*y[28])-x0*x1*\
y[1]*y[2]*y[9]*y[27]*y[28]-x0*x2*y[1]*y[2]*y[9]*y[27]*y[28]-x0*x1*y[1]*y[3]\
*y[9]*y[27]*y[28]-x0*x1*x2*y[1]*y[3]*y[9]*y[27]*y[28]-x0*x2*y[1]*y[4]*y[9]*\
y[27]*y[28]-x0*x1*x2*y[1]*y[4]*y[9]*y[27]*y[28]-x0*x2*y[1]*y[5]*y[9]*y[27]*\
y[28]-x0*x1*x2*y[1]*y[5]*y[9]*y[27]*y[28]-x0*y[1]*y[5]*y[6]*y[9]*y[27]*y[28\
]+x0*x1*y[1]*y[7]*y[9]*y[27]*y[28]+x0*x1*x2*y[1]*y[7]*y[9]*y[27]*y[28]-x0*y\
[1]*y[3]*y[9]*y[14]*y[27]*y[28]-2.*x2*y[1]*y[4]*y[9]*y[27]*y[28]*y[29]-2.*x\
1*x2*y[1]*y[4]*y[9]*y[27]*y[28]*y[29]-2.*y[1]*y[4]*y[6]*y[9]*y[27]*y[28]*y[\
29]-x1*y[1]*y[2]*y[31]*y[44]*y[45]-x0*x1*y[1]*y[2]*y[31]*y[44]*y[45]-x1*y[1\
]*y[3]*y[31]*y[44]*y[45]-x0*x1*y[1]*y[3]*y[31]*y[44]*y[45]-x0*x1*x2*y[1]*y[\
3]*y[31]*y[44]*y[45]-x0*x1*x2*y[1]*y[4]*y[31]*y[44]*y[45]-x1*x2*y[1]*y[5]*y\
[31]*y[44]*y[45]-x0*x1*x2*y[1]*y[5]*y[31]*y[44]*y[45]+x1*y[1]*y[7]*y[31]*y[\
44]*y[45]+x0*x1*y[1]*y[7]*y[31]*y[44]*y[45]+x0*x1*x2*y[1]*y[7]*y[31]*y[44]*\
y[45]-2.*y[1]*y[3]*y[14]*y[31]*y[44]*y[45]-2.*x0*y[1]*y[3]*y[14]*y[31]*y[44\
]*y[45]-x1*x2*y[1]*y[4]*y[29]*y[31]*y[44]*y[45]-x0*x2*y[1]*y[2]*y[47]*y[60]\
*y[61]-x0*x1*x2*y[1]*y[3]*y[47]*y[60]*y[61]-x0*x2*y[1]*y[4]*y[47]*y[60]*y[6\
1]-x0*x1*x2*y[1]*y[4]*y[47]*y[60]*y[61]-x2*y[1]*y[5]*y[47]*y[60]*y[61]-x0*x\
2*y[1]*y[5]*y[47]*y[60]*y[61]-x1*x2*y[1]*y[5]*y[47]*y[60]*y[61]-x0*x1*x2*y[\
1]*y[5]*y[47]*y[60]*y[61]-2.*x0*y[1]*y[5]*y[6]*y[47]*y[60]*y[61]+x0*x1*x2*y\
[1]*y[7]*y[47]*y[60]*y[61]-x2*y[1]*y[4]*y[29]*y[47]*y[60]*y[61]-x1*x2*y[1]*\
y[4]*y[29]*y[47]*y[60]*y[61]-2.*y[1]*y[4]*y[6]*y[29]*y[47]*y[60]*y[61])+pow\
(lambda,3)*(x0*x1*x2*y[1]*y[3]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45]*y[47]*y[6\
0]*y[61]+x0*x1*x2*y[1]*y[4]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45]*y[47]*y[60]*\
y[61]+x0*x1*x2*y[1]*y[5]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45]*y[47]*y[60]*y[6\
1]-x0*x1*x2*y[1]*y[7]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]+\
2.*x1*x2*y[1]*y[4]*y[9]*y[27]*y[28]*y[29]*y[31]*y[44]*y[45]*y[47]*y[60]*y[6\
1]+x1*x2*y[1]*y[4]*y[29]*y[31]*y[44]*y[45]*y[62]*y[63]*y[64]+x2*y[1]*y[4]*y\
[29]*y[47]*y[60]*y[61]*y[62]*y[63]*y[64]+x1*x2*y[1]*y[4]*y[29]*y[47]*y[60]*\
y[61]*y[62]*y[63]*y[64]+2.*y[1]*y[4]*y[6]*y[29]*y[47]*y[60]*y[61]*y[62]*y[6\
3]*y[64]+x0*y[1]*y[5]*y[6]*y[9]*y[27]*y[28]*y[65]*y[66]*y[67]+2.*y[1]*y[4]*\
y[6]*y[9]*y[27]*y[28]*y[29]*y[65]*y[66]*y[67]+x0*y[1]*y[3]*y[9]*y[14]*y[27]\
*y[28]*y[68]*y[69]*y[70]),2);
return (FOUT);
}
double Ps2(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[68];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=em[1];
y[4]=em[2];
y[5]=em[3];
y[6]=x2*x2;
y[7]=esx[0];
y[8]=-x0;
y[9]=1.+y[8];
y[10]=y[1]*y[2];
y[11]=x1*y[1]*y[2];
y[12]=x2*y[1]*y[2];
y[13]=x1*y[1]*y[3];
y[14]=x1*x1;
y[15]=y[1]*y[3]*y[14];
y[16]=x1*x2*y[1]*y[3];
y[17]=x2*y[1]*y[4];
y[18]=2.*x0*x2*y[1]*y[4];
y[19]=x1*x2*y[1]*y[4];
y[20]=2.*x0*x1*x2*y[1]*y[4];
y[21]=2.*x0*y[1]*y[4]*y[6];
y[22]=x2*y[1]*y[5];
y[23]=x1*x2*y[1]*y[5];
y[24]=y[1]*y[5]*y[6];
y[25]=-(x1*y[1]*y[7]);
y[26]=-(x1*x2*y[1]*y[7]);
y[27]=y[10]+y[11]+y[12]+y[13]+y[15]+y[16]+y[17]+y[18]+y[19]+y[20]+y[21]+y[22\
]+y[23]+y[24]+y[25]+y[26];
y[28]=lrs[0];
y[29]=x0*x0;
y[30]=-x1;
y[31]=1.+y[30];
y[32]=x0*y[1]*y[2];
y[33]=y[1]*y[3];
y[34]=x0*y[1]*y[3];
y[35]=2.*x1*y[1]*y[3];
y[36]=2.*x0*x1*y[1]*y[3];
y[37]=x0*x2*y[1]*y[3];
y[38]=x0*x2*y[1]*y[4];
y[39]=x2*y[1]*y[4]*y[29];
y[40]=x0*x2*y[1]*y[5];
y[41]=-(y[1]*y[7]);
y[42]=-(x0*y[1]*y[7]);
y[43]=-(x0*x2*y[1]*y[7]);
y[44]=y[10]+y[22]+y[32]+y[33]+y[34]+y[35]+y[36]+y[37]+y[38]+y[39]+y[40]+y[41\
]+y[42]+y[43];
y[45]=lrs[1];
y[46]=-x2;
y[47]=1.+y[46];
y[48]=x0*x1*y[1]*y[3];
y[49]=x0*y[1]*y[4];
y[50]=y[1]*y[4]*y[29];
y[51]=x0*x1*y[1]*y[4];
y[52]=x1*y[1]*y[4]*y[29];
y[53]=2.*x2*y[1]*y[4]*y[29];
y[54]=y[1]*y[5];
y[55]=x0*y[1]*y[5];
y[56]=x1*y[1]*y[5];
y[57]=x0*x1*y[1]*y[5];
y[58]=2.*x0*x2*y[1]*y[5];
y[59]=-(x0*x1*y[1]*y[7]);
y[60]=y[32]+y[48]+y[49]+y[50]+y[51]+y[52]+y[53]+y[54]+y[55]+y[56]+y[57]+y[58\
]+y[59];
y[61]=lrs[2];
y[62]=pow(y[9],2);
y[63]=pow(y[27],2);
y[64]=pow(y[28],2);
y[65]=pow(y[47],2);
y[66]=pow(y[60],2);
y[67]=pow(y[61],2);
FOUT=lambda*(-(x0*y[1]*y[2]*y[9]*y[27]*y[28])-x0*x1*y[1]*y[2]*y[9]*y[27]*y[2\
8]-x0*x2*y[1]*y[2]*y[9]*y[27]*y[28]-x0*x1*y[1]*y[3]*y[9]*y[27]*y[28]-x0*x1*\
x2*y[1]*y[3]*y[9]*y[27]*y[28]-x0*x2*y[1]*y[4]*y[9]*y[27]*y[28]-x0*x1*x2*y[1\
]*y[4]*y[9]*y[27]*y[28]-x0*x2*y[1]*y[5]*y[9]*y[27]*y[28]-x0*x1*x2*y[1]*y[5]\
*y[9]*y[27]*y[28]-x0*y[1]*y[5]*y[6]*y[9]*y[27]*y[28]+x0*x1*y[1]*y[7]*y[9]*y\
[27]*y[28]+x0*x1*x2*y[1]*y[7]*y[9]*y[27]*y[28]-x0*y[1]*y[3]*y[9]*y[14]*y[27\
]*y[28]-2.*x2*y[1]*y[4]*y[9]*y[27]*y[28]*y[29]-2.*x1*x2*y[1]*y[4]*y[9]*y[27\
]*y[28]*y[29]-2.*y[1]*y[4]*y[6]*y[9]*y[27]*y[28]*y[29]-x1*y[1]*y[2]*y[31]*y\
[44]*y[45]-x0*x1*y[1]*y[2]*y[31]*y[44]*y[45]-x1*y[1]*y[3]*y[31]*y[44]*y[45]\
-x0*x1*y[1]*y[3]*y[31]*y[44]*y[45]-x0*x1*x2*y[1]*y[3]*y[31]*y[44]*y[45]-x0*\
x1*x2*y[1]*y[4]*y[31]*y[44]*y[45]-x1*x2*y[1]*y[5]*y[31]*y[44]*y[45]-x0*x1*x\
2*y[1]*y[5]*y[31]*y[44]*y[45]+x1*y[1]*y[7]*y[31]*y[44]*y[45]+x0*x1*y[1]*y[7\
]*y[31]*y[44]*y[45]+x0*x1*x2*y[1]*y[7]*y[31]*y[44]*y[45]-2.*y[1]*y[3]*y[14]\
*y[31]*y[44]*y[45]-2.*x0*y[1]*y[3]*y[14]*y[31]*y[44]*y[45]-x1*x2*y[1]*y[4]*\
y[29]*y[31]*y[44]*y[45]-x0*x2*y[1]*y[2]*y[47]*y[60]*y[61]-x0*x1*x2*y[1]*y[3\
]*y[47]*y[60]*y[61]-x0*x2*y[1]*y[4]*y[47]*y[60]*y[61]-x0*x1*x2*y[1]*y[4]*y[\
47]*y[60]*y[61]-x2*y[1]*y[5]*y[47]*y[60]*y[61]-x0*x2*y[1]*y[5]*y[47]*y[60]*\
y[61]-x1*x2*y[1]*y[5]*y[47]*y[60]*y[61]-x0*x1*x2*y[1]*y[5]*y[47]*y[60]*y[61\
]-2.*x0*y[1]*y[5]*y[6]*y[47]*y[60]*y[61]+x0*x1*x2*y[1]*y[7]*y[47]*y[60]*y[6\
1]-x2*y[1]*y[4]*y[29]*y[47]*y[60]*y[61]-x1*x2*y[1]*y[4]*y[29]*y[47]*y[60]*y\
[61]-2.*y[1]*y[4]*y[6]*y[29]*y[47]*y[60]*y[61])+pow(lambda,3)*(x0*pow(y[31]\
,2)*pow(y[44],2)*pow(y[45],2)*y[1]*y[3]*y[9]*y[14]*y[27]*y[28]+x0*x1*x2*y[1\
]*y[3]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]+x0*x1*x2*y[1]*y\
[4]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]+x0*x1*x2*y[1]*y[5]\
*y[9]*y[27]*y[28]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]-x0*x1*x2*y[1]*y[7]*y[\
9]*y[27]*y[28]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]+2.*x1*x2*y[1]*y[4]*y[9]*\
y[27]*y[28]*y[29]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]+x1*x2*y[1]*y[4]*y[29]\
*y[31]*y[44]*y[45]*y[62]*y[63]*y[64]+x2*y[1]*y[4]*y[29]*y[47]*y[60]*y[61]*y\
[62]*y[63]*y[64]+x1*x2*y[1]*y[4]*y[29]*y[47]*y[60]*y[61]*y[62]*y[63]*y[64]+\
2.*y[1]*y[4]*y[6]*y[29]*y[47]*y[60]*y[61]*y[62]*y[63]*y[64]+x0*y[1]*y[5]*y[\
6]*y[9]*y[27]*y[28]*y[65]*y[66]*y[67]+2.*y[1]*y[4]*y[6]*y[9]*y[27]*y[28]*y[\
29]*y[65]*y[66]*y[67]);
return (FOUT);
}
double Pa2(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[71];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=em[1];
y[4]=em[2];
y[5]=em[3];
y[6]=x2*x2;
y[7]=esx[0];
y[8]=-x0;
y[9]=1.+y[8];
y[10]=y[1]*y[2];
y[11]=x1*y[1]*y[2];
y[12]=x2*y[1]*y[2];
y[13]=x1*y[1]*y[3];
y[14]=x1*x1;
y[15]=y[1]*y[3]*y[14];
y[16]=x1*x2*y[1]*y[3];
y[17]=x2*y[1]*y[4];
y[18]=2.*x0*x2*y[1]*y[4];
y[19]=x1*x2*y[1]*y[4];
y[20]=2.*x0*x1*x2*y[1]*y[4];
y[21]=2.*x0*y[1]*y[4]*y[6];
y[22]=x2*y[1]*y[5];
y[23]=x1*x2*y[1]*y[5];
y[24]=y[1]*y[5]*y[6];
y[25]=-(x1*y[1]*y[7]);
y[26]=-(x1*x2*y[1]*y[7]);
y[27]=y[10]+y[11]+y[12]+y[13]+y[15]+y[16]+y[17]+y[18]+y[19]+y[20]+y[21]+y[22\
]+y[23]+y[24]+y[25]+y[26];
y[28]=lrs[0];
y[29]=x0*x0;
y[30]=-x1;
y[31]=1.+y[30];
y[32]=x0*y[1]*y[2];
y[33]=y[1]*y[3];
y[34]=x0*y[1]*y[3];
y[35]=2.*x1*y[1]*y[3];
y[36]=2.*x0*x1*y[1]*y[3];
y[37]=x0*x2*y[1]*y[3];
y[38]=x0*x2*y[1]*y[4];
y[39]=x2*y[1]*y[4]*y[29];
y[40]=x0*x2*y[1]*y[5];
y[41]=-(y[1]*y[7]);
y[42]=-(x0*y[1]*y[7]);
y[43]=-(x0*x2*y[1]*y[7]);
y[44]=y[10]+y[22]+y[32]+y[33]+y[34]+y[35]+y[36]+y[37]+y[38]+y[39]+y[40]+y[41\
]+y[42]+y[43];
y[45]=lrs[1];
y[46]=-x2;
y[47]=1.+y[46];
y[48]=x0*x1*y[1]*y[3];
y[49]=x0*y[1]*y[4];
y[50]=y[1]*y[4]*y[29];
y[51]=x0*x1*y[1]*y[4];
y[52]=x1*y[1]*y[4]*y[29];
y[53]=2.*x2*y[1]*y[4]*y[29];
y[54]=y[1]*y[5];
y[55]=x0*y[1]*y[5];
y[56]=x1*y[1]*y[5];
y[57]=x0*x1*y[1]*y[5];
y[58]=2.*x0*x2*y[1]*y[5];
y[59]=-(x0*x1*y[1]*y[7]);
y[60]=y[32]+y[48]+y[49]+y[50]+y[51]+y[52]+y[53]+y[54]+y[55]+y[56]+y[57]+y[58\
]+y[59];
y[61]=lrs[2];
y[62]=pow(y[9],2);
y[63]=pow(y[27],2);
y[64]=pow(y[28],2);
y[65]=pow(y[47],2);
y[66]=pow(y[60],2);
y[67]=pow(y[61],2);
y[68]=pow(y[31],2);
y[69]=pow(y[44],2);
y[70]=pow(y[45],2);
FOUT=(lambda*(-(x0*y[1]*y[2]*y[9]*y[27]*y[28])-x0*x1*y[1]*y[2]*y[9]*y[27]*y[\
28]-x0*x2*y[1]*y[2]*y[9]*y[27]*y[28]-x0*x1*y[1]*y[3]*y[9]*y[27]*y[28]-x0*x1\
*x2*y[1]*y[3]*y[9]*y[27]*y[28]-x0*x2*y[1]*y[4]*y[9]*y[27]*y[28]-x0*x1*x2*y[\
1]*y[4]*y[9]*y[27]*y[28]-x0*x2*y[1]*y[5]*y[9]*y[27]*y[28]-x0*x1*x2*y[1]*y[5\
]*y[9]*y[27]*y[28]-x0*y[1]*y[5]*y[6]*y[9]*y[27]*y[28]+x0*x1*y[1]*y[7]*y[9]*\
y[27]*y[28]+x0*x1*x2*y[1]*y[7]*y[9]*y[27]*y[28]-x0*y[1]*y[3]*y[9]*y[14]*y[2\
7]*y[28]-2.*x2*y[1]*y[4]*y[9]*y[27]*y[28]*y[29]-2.*x1*x2*y[1]*y[4]*y[9]*y[2\
7]*y[28]*y[29]-2.*y[1]*y[4]*y[6]*y[9]*y[27]*y[28]*y[29]-x1*y[1]*y[2]*y[31]*\
y[44]*y[45]-x0*x1*y[1]*y[2]*y[31]*y[44]*y[45]-x1*y[1]*y[3]*y[31]*y[44]*y[45\
]-x0*x1*y[1]*y[3]*y[31]*y[44]*y[45]-x0*x1*x2*y[1]*y[3]*y[31]*y[44]*y[45]-x0\
*x1*x2*y[1]*y[4]*y[31]*y[44]*y[45]-x1*x2*y[1]*y[5]*y[31]*y[44]*y[45]-x0*x1*\
x2*y[1]*y[5]*y[31]*y[44]*y[45]+x1*y[1]*y[7]*y[31]*y[44]*y[45]+x0*x1*y[1]*y[\
7]*y[31]*y[44]*y[45]+x0*x1*x2*y[1]*y[7]*y[31]*y[44]*y[45]-2.*y[1]*y[3]*y[14\
]*y[31]*y[44]*y[45]-2.*x0*y[1]*y[3]*y[14]*y[31]*y[44]*y[45]-x1*x2*y[1]*y[4]\
*y[29]*y[31]*y[44]*y[45]-x0*x2*y[1]*y[2]*y[47]*y[60]*y[61]-x0*x1*x2*y[1]*y[\
3]*y[47]*y[60]*y[61]-x0*x2*y[1]*y[4]*y[47]*y[60]*y[61]-x0*x1*x2*y[1]*y[4]*y\
[47]*y[60]*y[61]-x2*y[1]*y[5]*y[47]*y[60]*y[61]-x0*x2*y[1]*y[5]*y[47]*y[60]\
*y[61]-x1*x2*y[1]*y[5]*y[47]*y[60]*y[61]-x0*x1*x2*y[1]*y[5]*y[47]*y[60]*y[6\
1]-2.*x0*y[1]*y[5]*y[6]*y[47]*y[60]*y[61]+x0*x1*x2*y[1]*y[7]*y[47]*y[60]*y[\
61]-x2*y[1]*y[4]*y[29]*y[47]*y[60]*y[61]-x1*x2*y[1]*y[4]*y[29]*y[47]*y[60]*\
y[61]-2.*y[1]*y[4]*y[6]*y[29]*y[47]*y[60]*y[61])+pow(lambda,3)*(x0*x1*x2*y[\
1]*y[3]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]+x0*x1*x2*y[1]*\
y[4]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]+x0*x1*x2*y[1]*y[5\
]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]-x0*x1*x2*y[1]*y[7]*y\
[9]*y[27]*y[28]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]+2.*x1*x2*y[1]*y[4]*y[9]\
*y[27]*y[28]*y[29]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]+x1*x2*y[1]*y[4]*y[29\
]*y[31]*y[44]*y[45]*y[62]*y[63]*y[64]+x2*y[1]*y[4]*y[29]*y[47]*y[60]*y[61]*\
y[62]*y[63]*y[64]+x1*x2*y[1]*y[4]*y[29]*y[47]*y[60]*y[61]*y[62]*y[63]*y[64]\
+2.*y[1]*y[4]*y[6]*y[29]*y[47]*y[60]*y[61]*y[62]*y[63]*y[64]+x0*y[1]*y[5]*y\
[6]*y[9]*y[27]*y[28]*y[65]*y[66]*y[67]+2.*y[1]*y[4]*y[6]*y[9]*y[27]*y[28]*y\
[29]*y[65]*y[66]*y[67]+x0*y[1]*y[3]*y[9]*y[14]*y[27]*y[28]*y[68]*y[69]*y[70\
]))/(lambda*(x0*x1*y[1]*y[2]+x0*x2*y[1]*y[2]+x0*x1*x2*y[1]*y[3]+x0*x1*x2*y[\
1]*y[4]+x0*x1*x2*y[1]*y[5]+x0*y[1]*y[5]*y[6]-x0*x1*x2*y[1]*y[7]+y[10]+y[11]\
+y[13]+x0*y[1]*y[3]*y[14]+y[15]+y[22]+y[23]+y[25]+x1*x2*y[1]*y[4]*y[29]+y[1\
]*y[4]*y[6]*y[29]+y[32]+y[38]+y[39]+y[40]+y[48]+y[59]+pow(lambda,4)*(x1*x2*\
y[1]*y[4]*y[29]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]*y[62]*y[63]*y[64]+y[1]*\
y[4]*y[6]*y[29]*y[62]*y[63]*y[64]*y[65]*y[66]*y[67])+lambda*lambda*(-(x0*x1\
*y[1]*y[2]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45])-x0*x1*y[1]*y[3]*y[9]*y[27]*y\
[28]*y[31]*y[44]*y[45]-x0*x1*x2*y[1]*y[3]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45\
]-x0*x1*x2*y[1]*y[4]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45]-x0*x1*x2*y[1]*y[5]*\
y[9]*y[27]*y[28]*y[31]*y[44]*y[45]+x0*x1*y[1]*y[7]*y[9]*y[27]*y[28]*y[31]*y\
[44]*y[45]+x0*x1*x2*y[1]*y[7]*y[9]*y[27]*y[28]*y[31]*y[44]*y[45]-2.*x0*y[1]\
*y[3]*y[9]*y[14]*y[27]*y[28]*y[31]*y[44]*y[45]-2.*x1*x2*y[1]*y[4]*y[9]*y[27\
]*y[28]*y[29]*y[31]*y[44]*y[45]-x0*x2*y[1]*y[2]*y[9]*y[27]*y[28]*y[47]*y[60\
]*y[61]-x0*x1*x2*y[1]*y[3]*y[9]*y[27]*y[28]*y[47]*y[60]*y[61]-x0*x2*y[1]*y[\
4]*y[9]*y[27]*y[28]*y[47]*y[60]*y[61]-x0*x1*x2*y[1]*y[4]*y[9]*y[27]*y[28]*y\
[47]*y[60]*y[61]-x0*x2*y[1]*y[5]*y[9]*y[27]*y[28]*y[47]*y[60]*y[61]-x0*x1*x\
2*y[1]*y[5]*y[9]*y[27]*y[28]*y[47]*y[60]*y[61]-2.*x0*y[1]*y[5]*y[6]*y[9]*y[\
27]*y[28]*y[47]*y[60]*y[61]+x0*x1*x2*y[1]*y[7]*y[9]*y[27]*y[28]*y[47]*y[60]\
*y[61]-2.*x2*y[1]*y[4]*y[9]*y[27]*y[28]*y[29]*y[47]*y[60]*y[61]-2.*x1*x2*y[\
1]*y[4]*y[9]*y[27]*y[28]*y[29]*y[47]*y[60]*y[61]-4.*y[1]*y[4]*y[6]*y[9]*y[2\
7]*y[28]*y[29]*y[47]*y[60]*y[61]-x0*x1*x2*y[1]*y[3]*y[31]*y[44]*y[45]*y[47]\
*y[60]*y[61]-x0*x1*x2*y[1]*y[4]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]-x1*x2*y\
[1]*y[5]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]-x0*x1*x2*y[1]*y[5]*y[31]*y[44]\
*y[45]*y[47]*y[60]*y[61]+x0*x1*x2*y[1]*y[7]*y[31]*y[44]*y[45]*y[47]*y[60]*y\
[61]-x1*x2*y[1]*y[4]*y[29]*y[31]*y[44]*y[45]*y[47]*y[60]*y[61]-x2*y[1]*y[4]\
*y[29]*y[62]*y[63]*y[64]-x1*x2*y[1]*y[4]*y[29]*y[62]*y[63]*y[64]-y[1]*y[4]*\
y[6]*y[29]*y[62]*y[63]*y[64]-x0*y[1]*y[5]*y[6]*y[65]*y[66]*y[67]-y[1]*y[4]*\
y[6]*y[29]*y[65]*y[66]*y[67]-y[1]*y[3]*y[14]*y[68]*y[69]*y[70]-x0*y[1]*y[3]\
*y[14]*y[68]*y[69]*y[70])));
return (FOUT);
}
double Pt2t1(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[8];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=em[1];
y[4]=em[2];
y[5]=em[3];
y[6]=x2*x2;
y[7]=esx[0];
FOUT=(1.-x0)*x0*(y[1]*y[2]+x1*y[1]*y[2]+x2*y[1]*y[2]+x1*y[1]*y[3]+x1*x2*y[1]\
*y[3]+x1*x1*y[1]*y[3]+x2*y[1]*y[4]+2.*x0*x2*y[1]*y[4]+x1*x2*y[1]*y[4]+2.*x0\
*x1*x2*y[1]*y[4]+x2*y[1]*y[5]+x1*x2*y[1]*y[5]+2.*x0*y[1]*y[4]*y[6]+y[1]*y[5\
]*y[6]-x1*y[1]*y[7]-x1*x2*y[1]*y[7]);
return (FOUT);
}
double Pt2t2(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[7];
double FOUT;
y[1]=1./bi;
y[2]=em[0];
y[3]=em[1];
y[4]=em[2];
y[5]=em[3];
y[6]=esx[0];
FOUT=(1.-x1)*x1*(y[1]*y[2]+x0*y[1]*y[2]+y[1]*y[3]+x0*y[1]*y[3]+2.*x1*y[1]*y[\
3]+2.*x0*x1*y[1]*y[3]+x0*x2*y[1]*y[3]+x0*x2*y[1]*y[4]+x2*(x0*x0)*y[1]*y[4]+\
x2*y[1]*y[5]+x0*x2*y[1]*y[5]-y[1]*y[6]-x0*y[1]*y[6]-x0*x2*y[1]*y[6]);
return (FOUT);
}
double Pt2t3(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double y[5];
double FOUT;
y[1]=1./bi;
y[2]=em[2];
y[3]=x0*x0;
y[4]=em[3];
FOUT=(1.-x2)*x2*(x0*em[0]*y[1]+x0*x1*em[1]*y[1]-x0*x1*esx[0]*y[1]+x0*y[1]*y[\
2]+x0*x1*y[1]*y[2]+y[1]*y[2]*y[3]+x1*y[1]*y[2]*y[3]+2.*x2*y[1]*y[2]*y[3]+y[\
1]*y[4]+x0*y[1]*y[4]+x1*y[1]*y[4]+x0*x1*y[1]*y[4]+2.*x0*x2*y[1]*y[4]);
return (FOUT);
}
