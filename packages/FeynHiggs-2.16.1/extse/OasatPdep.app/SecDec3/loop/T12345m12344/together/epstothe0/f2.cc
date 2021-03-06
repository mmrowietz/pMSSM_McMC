#include "intfile.hh"

dcmplx Pf2(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double x3=x[3];
dcmplx y[202];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[0];
y[3]=x0*x0;
y[4]=em[1];
y[5]=esx[0];
y[6]=x0*y[1]*y[4];
y[7]=em[3];
y[8]=2.*x0*x2*y[1]*y[4];
y[9]=em[2];
y[10]=x2*y[1]*y[4];
y[11]=y[1]*y[9];
y[12]=-x1;
y[13]=1.+y[12];
y[14]=x0*y[1]*y[2];
y[15]=x0*y[1]*y[9];
y[16]=2.*y[1]*y[7];
y[17]=x0*y[1]*y[7];
y[18]=2.*x1*y[1]*y[7];
y[19]=-(y[1]*y[5]);
y[20]=-(x0*y[1]*y[5]);
y[21]=-x0;
y[22]=1.+y[21];
y[23]=x3*x3;
y[24]=x2*x2;
y[25]=lrs[0];
y[26]=x3*y[1]*y[2];
y[27]=x2*x3*y[1]*y[2];
y[28]=y[1]*y[2]*y[23];
y[29]=x1*x2*y[1]*y[4];
y[30]=y[1]*y[4]*y[24];
y[31]=x2*x3*y[1]*y[4];
y[32]=x1*y[1]*y[9];
y[33]=x2*y[1]*y[9];
y[34]=x3*y[1]*y[9];
y[35]=x2*y[1]*y[7];
y[36]=x3*y[1]*y[7];
y[37]=-(x2*y[1]*y[5]);
y[38]=-(x2*x3*y[1]*y[5]);
y[39]=x1*x3*y[1]*y[2];
y[40]=2.*x0*x2*x3*y[1]*y[2];
y[41]=x1*x2*x3*y[1]*y[2];
y[42]=2.*x0*y[1]*y[2]*y[23];
y[43]=x1*y[1]*y[2]*y[23];
y[44]=2.*x0*y[1]*y[4]*y[24];
y[45]=x1*y[1]*y[4]*y[24];
y[46]=2.*x0*x2*x3*y[1]*y[4];
y[47]=x1*x2*x3*y[1]*y[4];
y[48]=2.*x0*x2*y[1]*y[9];
y[49]=x1*x2*y[1]*y[9];
y[50]=2.*x0*x3*y[1]*y[9];
y[51]=x1*x3*y[1]*y[9];
y[52]=x1*x2*y[1]*y[7];
y[53]=x1*x3*y[1]*y[7];
y[54]=-(x1*x3*y[1]*y[5]);
y[55]=-2.*x0*x2*x3*y[1]*y[5];
y[56]=-(x1*x2*x3*y[1]*y[5]);
y[57]=y[10]+y[11]+y[26]+y[27]+y[28]+y[29]+y[30]+y[31]+y[32]+y[33]+y[34]+y[35\
]+y[36]+y[37]+y[38]+y[39]+y[40]+y[41]+y[42]+y[43]+y[44]+y[45]+y[46]+y[47]+y\
[48]+y[49]+y[50]+y[51]+y[52]+y[53]+y[54]+y[55]+y[56];
y[58]=lrs[1];
y[59]=-x2;
y[60]=1.+y[59];
y[61]=2.*x0*x3*y[1]*y[2];
y[62]=2.*x0*y[1]*y[9];
y[63]=y[1]*y[7];
y[64]=x1*y[1]*y[7];
y[65]=-(x3*y[1]*y[5]);
y[66]=lambda*lambda;
y[67]=x0*x2*y[1]*y[2];
y[68]=x0*x2*y[1]*y[4];
y[69]=-(x0*x2*y[1]*y[5]);
y[70]=y[14]+y[15]+y[16]+y[17]+y[18]+y[19]+y[20]+y[61]+y[67]+y[68]+y[69];
y[71]=y[10]+y[11]+y[26]+y[27]+y[28]+y[30]+y[31]+y[33]+y[34]+y[35]+y[36]+y[38\
]+y[65];
y[72]=y[1]*y[2];
y[73]=x1*y[1]*y[2];
y[74]=x2*y[1]*y[2];
y[75]=2.*x0*x2*y[1]*y[2];
y[76]=x1*x2*y[1]*y[2];
y[77]=2.*x3*y[1]*y[2];
y[78]=4.*x0*x3*y[1]*y[2];
y[79]=2.*x1*x3*y[1]*y[2];
y[80]=-(x1*y[1]*y[5]);
y[81]=-2.*x0*x2*y[1]*y[5];
y[82]=-(x1*x2*y[1]*y[5]);
y[83]=y[8]+y[10]+y[11]+y[29]+y[32]+y[37]+y[62]+y[63]+y[64]+y[72]+y[73]+y[74]\
+y[75]+y[76]+y[77]+y[78]+y[79]+y[80]+y[81]+y[82];
y[84]=x0*x3*y[1]*y[2];
y[85]=2.*x2*y[1]*y[7];
y[86]=2.*x3*y[1]*y[7];
y[87]=-(x0*x3*y[1]*y[5]);
y[88]=x0*x2*x3*y[1]*y[2];
y[89]=x0*y[1]*y[2]*y[23];
y[90]=x0*y[1]*y[4]*y[24];
y[91]=x0*x2*x3*y[1]*y[4];
y[92]=x0*x2*y[1]*y[9];
y[93]=x0*x3*y[1]*y[9];
y[94]=x0*x2*y[1]*y[7];
y[95]=2.*x1*x2*y[1]*y[7];
y[96]=x0*x3*y[1]*y[7];
y[97]=2.*x1*x3*y[1]*y[7];
y[98]=-(x0*x2*x3*y[1]*y[5]);
y[99]=y[15]+y[16]+y[18]+y[19]+y[37]+y[65]+y[68]+y[84]+y[85]+y[86]+y[87]+y[88\
]+y[89]+y[90]+y[91]+y[92]+y[93]+y[94]+y[95]+y[96]+y[97]+y[98];
y[100]=lrs[2];
y[101]=y[1]*y[2]*y[3];
y[102]=x0*x1*y[1]*y[2];
y[103]=y[1]*y[3]*y[4];
y[104]=x0*x1*y[1]*y[4];
y[105]=-(y[1]*y[3]*y[5]);
y[106]=-(x0*x1*y[1]*y[5]);
y[107]=y[6]+y[14]+y[20]+y[101]+y[102]+y[103]+y[104]+y[105]+y[106];
y[108]=2.*x2*x3*y[1]*y[2];
y[109]=2.*y[1]*y[2]*y[23];
y[110]=2.*y[1]*y[4]*y[24];
y[111]=2.*x2*x3*y[1]*y[4];
y[112]=2.*x2*y[1]*y[9];
y[113]=2.*x3*y[1]*y[9];
y[114]=-2.*x2*x3*y[1]*y[5];
y[115]=y[108]+y[109]+y[110]+y[111]+y[112]+y[113]+y[114];
y[116]=-(lambda*MYI*x0*y[22]*y[25]*y[115]);
y[117]=-(lambda*MYI*y[22]*y[25]*y[57]);
y[118]=lambda*MYI*x0*y[25]*y[57];
y[119]=1.+y[116]+y[117]+y[118];
y[120]=y[16]+y[85]+y[86];
y[121]=-(lambda*MYI*x1*y[13]*y[58]*y[120]);
y[122]=-(lambda*MYI*y[13]*y[58]*y[99]);
y[123]=lambda*MYI*x1*y[58]*y[99];
y[124]=1.+y[121]+y[122]+y[123];
y[125]=-x3;
y[126]=1.+y[125];
y[127]=y[1]*y[4];
y[128]=x1*y[1]*y[4];
y[129]=2.*x2*y[1]*y[4];
y[130]=4.*x0*x2*y[1]*y[4];
y[131]=2.*x1*x2*y[1]*y[4];
y[132]=x3*y[1]*y[4];
y[133]=2.*x0*x3*y[1]*y[4];
y[134]=x1*x3*y[1]*y[4];
y[135]=-2.*x0*x3*y[1]*y[5];
y[136]=y[11]+y[19]+y[26]+y[32]+y[39]+y[54]+y[61]+y[62]+y[63]+y[64]+y[65]+y[1\
27]+y[128]+y[129]+y[130]+y[131]+y[132]+y[133]+y[134]+y[135];
y[137]=x0*x3*y[1]*y[4];
y[138]=y[6]+y[8]+y[15]+y[16]+y[17]+y[18]+y[19]+y[84]+y[87]+y[137];
y[139]=x0*x1*y[13]*y[22]*y[25]*y[58]*y[66]*y[71]*y[83];
y[140]=-(lambda*MYI*x1*y[13]*y[58]*y[70]*y[119]);
y[141]=y[139]+y[140];
y[142]=x3*y[1]*y[2]*y[3];
y[143]=x0*x1*x3*y[1]*y[2];
y[144]=2.*x2*y[1]*y[3]*y[4];
y[145]=2.*x0*x1*x2*y[1]*y[4];
y[146]=x3*y[1]*y[3]*y[4];
y[147]=x0*x1*x3*y[1]*y[4];
y[148]=y[1]*y[3]*y[9];
y[149]=x0*x1*y[1]*y[9];
y[150]=x0*x1*y[1]*y[7];
y[151]=x1*x1;
y[152]=y[1]*y[7]*y[151];
y[153]=-(x3*y[1]*y[3]*y[5]);
y[154]=-(x0*x1*x3*y[1]*y[5]);
y[155]=y[6]+y[8]+y[15]+y[17]+y[18]+y[20]+y[63]+y[80]+y[84]+y[87]+y[104]+y[13\
7]+y[142]+y[143]+y[144]+y[145]+y[146]+y[147]+y[148]+y[149]+y[150]+y[152]+y[\
153]+y[154];
y[156]=lrs[3];
y[157]=x0*x1*y[13]*y[22]*y[25]*y[58]*y[66]*y[83]*y[138];
y[158]=-(x0*x1*y[13]*y[22]*y[25]*y[58]*y[66]*y[70]*y[136]);
y[159]=y[157]+y[158];
y[160]=-(x0*x1*y[13]*y[22]*y[25]*y[58]*y[66]*y[70]*y[71]);
y[161]=lambda*MYI*x0*y[22]*y[25]*y[83]*y[124];
y[162]=y[160]+y[161];
y[163]=2.*x0*y[1]*y[4];
y[164]=2.*y[1]*y[3]*y[4];
y[165]=2.*x0*x1*y[1]*y[4];
y[166]=y[163]+y[164]+y[165];
y[167]=-(lambda*MYI*x2*y[60]*y[100]*y[166]);
y[168]=-(lambda*MYI*y[60]*y[100]*y[155]);
y[169]=lambda*MYI*x2*y[100]*y[155];
y[170]=1.+y[167]+y[168]+y[169];
y[171]=x0*x1*y[13]*y[22]*y[25]*y[58]*y[66]*y[71]*y[136];
y[172]=-(lambda*MYI*x1*y[13]*y[58]*y[119]*y[138]);
y[173]=y[171]+y[172];
y[174]=-(x0*x1*y[13]*y[22]*y[25]*y[58]*y[66]*y[71]*y[138]);
y[175]=lambda*MYI*x0*y[22]*y[25]*y[124]*y[136];
y[176]=y[174]+y[175];
y[177]=pow(y[71],2);
y[178]=x0*x1*y[13]*y[22]*y[25]*y[58]*y[66]*y[177];
y[179]=y[119]*y[124];
y[180]=y[178]+y[179];
y[181]=x2*y[1]*y[2]*y[3];
y[182]=x0*x1*x2*y[1]*y[2];
y[183]=2.*x3*y[1]*y[2]*y[3];
y[184]=2.*x0*x1*x3*y[1]*y[2];
y[185]=x2*y[1]*y[3]*y[4];
y[186]=x0*x1*x2*y[1]*y[4];
y[187]=-(x2*y[1]*y[3]*y[5]);
y[188]=-(x0*x1*x2*y[1]*y[5]);
y[189]=y[14]+y[15]+y[17]+y[18]+y[61]+y[63]+y[67]+y[68]+y[69]+y[80]+y[102]+y[\
106]+y[148]+y[149]+y[150]+y[152]+y[181]+y[182]+y[183]+y[184]+y[185]+y[186]+\
y[187]+y[188];
y[190]=-(lambda*MYI*x2*y[60]*y[100]*y[155]);
y[191]=x2+y[190];
y[192]=-(lambda*MYI*x1*y[13]*y[58]*y[99]);
y[193]=x1+y[192];
y[194]=-(lambda*MYI*x0*y[22]*y[25]*y[57]);
y[195]=x0+y[194];
y[196]=-(lambda*MYI*x3*y[126]*y[156]*y[189]);
y[197]=x3+y[196];
y[198]=pow(y[193],2);
y[199]=pow(y[195],2);
y[200]=pow(y[191],2);
y[201]=pow(y[197],2);
FOUT=(pow(bi,-2)*(-(lambda*MYI*x3*y[70]*y[126]*y[156]*(-(lambda*MYI*x2*y[60]\
*y[100]*y[136]*y[159])-y[141]*y[170]-lambda*MYI*x2*y[60]*y[100]*y[107]*y[17\
3]))+lambda*MYI*x3*y[83]*y[126]*y[156]*(-(lambda*MYI*x2*y[60]*y[100]*y[138]\
*y[159])-y[162]*y[170]-lambda*MYI*x2*y[60]*y[100]*y[107]*y[176])+lambda*MYI\
*x3*y[107]*y[126]*y[156]*(lambda*MYI*x2*y[60]*y[100]*y[138]*y[141]-lambda*M\
YI*x2*y[60]*y[100]*y[136]*y[162]-lambda*MYI*x2*y[60]*y[100]*y[107]*y[180])+\
(lambda*MYI*x2*y[60]*y[100]*y[138]*y[173]-lambda*MYI*x2*y[60]*y[100]*y[136]\
*y[176]+y[170]*y[180])*(1.-lambda*MYI*x3*(2.*x0*y[1]*y[2]+2.*x0*x1*y[1]*y[2\
]+2.*y[1]*y[2]*y[3])*y[126]*y[156]+lambda*MYI*x3*y[156]*y[189]-lambda*MYI*y\
[126]*y[156]*y[189])))/((y[1]+y[1]*y[191]+y[1]*y[193]+y[1]*y[191]*y[193]+y[\
1]*y[191]*y[195]+y[1]*y[197]+y[1]*y[193]*y[197]+y[1]*y[195]*y[197])*(y[63]+\
y[1]*y[7]*y[191]-y[1]*y[5]*y[193]+2.*y[1]*y[7]*y[193]-y[1]*y[5]*y[191]*y[19\
3]+2.*y[1]*y[7]*y[191]*y[193]+y[1]*y[9]*y[195]+y[1]*y[4]*y[191]*y[195]-y[1]\
*y[5]*y[191]*y[195]+y[1]*y[7]*y[191]*y[195]+y[1]*y[9]*y[191]*y[195]+y[1]*y[\
9]*y[193]*y[195]+y[1]*y[4]*y[191]*y[193]*y[195]+y[1]*y[7]*y[191]*y[193]*y[1\
95]+y[1]*y[9]*y[191]*y[193]*y[195]+y[1]*y[7]*y[197]-y[1]*y[5]*y[193]*y[197]\
+2.*y[1]*y[7]*y[193]*y[197]+y[1]*y[2]*y[195]*y[197]+y[1]*y[7]*y[195]*y[197]\
+y[1]*y[9]*y[195]*y[197]+y[1]*y[2]*y[191]*y[195]*y[197]+y[1]*y[4]*y[191]*y[\
195]*y[197]-y[1]*y[5]*y[191]*y[195]*y[197]+y[1]*y[2]*y[193]*y[195]*y[197]-y\
[1]*y[5]*y[193]*y[195]*y[197]+y[1]*y[7]*y[193]*y[195]*y[197]+y[1]*y[9]*y[19\
3]*y[195]*y[197]+y[1]*y[2]*y[191]*y[193]*y[195]*y[197]+y[1]*y[4]*y[191]*y[1\
93]*y[195]*y[197]-y[1]*y[5]*y[191]*y[193]*y[195]*y[197]+y[1]*y[7]*y[198]+y[\
1]*y[7]*y[191]*y[198]+y[1]*y[7]*y[197]*y[198]+y[1]*y[9]*y[191]*y[199]+y[1]*\
y[9]*y[197]*y[199]+y[1]*y[2]*y[191]*y[197]*y[199]+y[1]*y[4]*y[191]*y[197]*y\
[199]-y[1]*y[5]*y[191]*y[197]*y[199]+y[1]*y[4]*y[195]*y[200]+y[1]*y[4]*y[19\
3]*y[195]*y[200]+y[1]*y[4]*y[199]*y[200]+y[1]*y[2]*y[195]*y[201]+y[1]*y[2]*\
y[193]*y[195]*y[201]+y[1]*y[2]*y[199]*y[201]));
return (FOUT);
}
