#include "intfile.hh"

dcmplx Pf8(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
double x3=x[3];
dcmplx y[204];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[1];
y[3]=em[3];
y[4]=esx[0];
y[5]=em[0];
y[6]=y[1]*y[5];
y[7]=y[1]*y[2];
y[8]=2.*x2*y[1]*y[2];
y[9]=em[2];
y[10]=-(y[1]*y[4]);
y[11]=-(x0*y[1]*y[4]);
y[12]=2.*x0*y[1]*y[3];
y[13]=-(x2*y[1]*y[4]);
y[14]=-x1;
y[15]=1.+y[14];
y[16]=x0*y[1]*y[9];
y[17]=x0*y[1]*y[3];
y[18]=x0*x0;
y[19]=-x0;
y[20]=1.+y[19];
y[21]=x1*x1;
y[22]=x3*x3;
y[23]=lrs[0];
y[24]=x1*y[1]*y[5];
y[25]=x1*x2*y[1]*y[2];
y[26]=x1*y[1]*y[9];
y[27]=2.*x0*y[1]*y[9]*y[21];
y[28]=x1*x2*y[1]*y[9];
y[29]=y[1]*y[3];
y[30]=x1*y[1]*y[3];
y[31]=x2*y[1]*y[3];
y[32]=x1*x2*y[1]*y[3];
y[33]=2.*x3*y[1]*y[3];
y[34]=4.*x0*x1*x3*y[1]*y[3];
y[35]=2.*x2*x3*y[1]*y[3];
y[36]=-(x3*y[1]*y[4]);
y[37]=x1*x3*y[1]*y[5];
y[38]=x1*x2*x3*y[1]*y[2];
y[39]=y[1]*y[9]*y[21];
y[40]=x2*y[1]*y[9]*y[21];
y[41]=x1*x3*y[1]*y[9];
y[42]=2.*x0*x3*y[1]*y[9]*y[21];
y[43]=x1*x2*x3*y[1]*y[9];
y[44]=2.*x0*x1*y[1]*y[3];
y[45]=x1*x3*y[1]*y[3];
y[46]=x1*x2*x3*y[1]*y[3];
y[47]=y[1]*y[3]*y[22];
y[48]=2.*x0*x1*y[1]*y[3]*y[22];
y[49]=x2*y[1]*y[3]*y[22];
y[50]=-(x1*x2*y[1]*y[4]);
y[51]=-(x1*x3*y[1]*y[4]);
y[52]=-2.*x0*x1*x3*y[1]*y[4];
y[53]=-(x2*x3*y[1]*y[4]);
y[54]=y[24]+y[25]+y[26]+y[27]+y[28]+y[29]+y[30]+y[31]+y[32]+y[33]+y[34]+y[35\
]+y[36]+y[37]+y[38]+y[39]+y[40]+y[41]+y[42]+y[43]+y[44]+y[45]+y[46]+y[47]+y\
[48]+y[49]+y[50]+y[51]+y[52]+y[53];
y[55]=lrs[1];
y[56]=-x2;
y[57]=1.+y[56];
y[58]=-(x1*y[1]*y[4]);
y[59]=lambda*lambda;
y[60]=x0*y[1]*y[5];
y[61]=x0*x2*y[1]*y[2];
y[62]=2.*x1*y[1]*y[9]*y[18];
y[63]=x0*x2*y[1]*y[9];
y[64]=2.*y[1]*y[3]*y[18];
y[65]=x0*x2*y[1]*y[3];
y[66]=2.*x3*y[1]*y[3]*y[18];
y[67]=-(y[1]*y[4]*y[18]);
y[68]=y[11]+y[16]+y[17]+y[60]+y[61]+y[62]+y[63]+y[64]+y[65]+y[66]+y[67];
y[69]=x3*y[1]*y[5];
y[70]=x2*y[1]*y[2];
y[71]=x2*x3*y[1]*y[2];
y[72]=y[1]*y[9];
y[73]=2.*x1*y[1]*y[9];
y[74]=4.*x0*x1*y[1]*y[9];
y[75]=x2*y[1]*y[9];
y[76]=2.*x1*x2*y[1]*y[9];
y[77]=x3*y[1]*y[9];
y[78]=4.*x0*x1*x3*y[1]*y[9];
y[79]=x2*x3*y[1]*y[9];
y[80]=x3*y[1]*y[3];
y[81]=4.*x0*x3*y[1]*y[3];
y[82]=x2*x3*y[1]*y[3];
y[83]=2.*x0*y[1]*y[3]*y[22];
y[84]=-2.*x0*x3*y[1]*y[4];
y[85]=y[6]+y[12]+y[13]+y[29]+y[31]+y[36]+y[69]+y[70]+y[71]+y[72]+y[73]+y[74]\
+y[75]+y[76]+y[77]+y[78]+y[79]+y[80]+y[81]+y[82]+y[83]+y[84];
y[86]=2.*y[1]*y[3];
y[87]=4.*x0*x1*y[1]*y[3];
y[88]=2.*x2*y[1]*y[3];
y[89]=-2.*x0*x1*y[1]*y[4];
y[90]=y[10]+y[13]+y[24]+y[25]+y[26]+y[27]+y[28]+y[30]+y[32]+y[33]+y[34]+y[35\
]+y[58]+y[86]+y[87]+y[88]+y[89];
y[91]=2.*x0*x1*y[1]*y[9];
y[92]=x0*x3*y[1]*y[9];
y[93]=x0*x3*y[1]*y[3];
y[94]=x2*y[1]*y[5];
y[95]=x0*x3*y[1]*y[5];
y[96]=x2*x2;
y[97]=y[1]*y[2]*y[96];
y[98]=x0*x2*x3*y[1]*y[2];
y[99]=2.*x0*x1*x2*y[1]*y[9];
y[100]=2.*x1*x3*y[1]*y[9]*y[18];
y[101]=x0*x2*x3*y[1]*y[9];
y[102]=y[1]*y[3]*y[18];
y[103]=x0*x2*x3*y[1]*y[3];
y[104]=y[1]*y[3]*y[18]*y[22];
y[105]=-(x0*x2*y[1]*y[4]);
y[106]=-(x0*x3*y[1]*y[4]);
y[107]=-(x3*y[1]*y[4]*y[18]);
y[108]=y[6]+y[13]+y[16]+y[17]+y[60]+y[61]+y[62]+y[63]+y[65]+y[66]+y[70]+y[91\
]+y[92]+y[93]+y[94]+y[95]+y[97]+y[98]+y[99]+y[100]+y[101]+y[102]+y[103]+y[1\
04]+y[105]+y[106]+y[107];
y[109]=lrs[2];
y[110]=x0*x1*y[1]*y[2];
y[111]=x0*x1*y[1]*y[9];
y[112]=x0*x1*y[1]*y[3];
y[113]=2.*x0*x3*y[1]*y[3];
y[114]=y[6]+y[7]+y[8]+y[10]+y[11]+y[12]+y[110]+y[111]+y[112]+y[113];
y[115]=2.*y[1]*y[9]*y[21];
y[116]=2.*x3*y[1]*y[9]*y[21];
y[117]=2.*x1*y[1]*y[3];
y[118]=4.*x1*x3*y[1]*y[3];
y[119]=2.*x1*y[1]*y[3]*y[22];
y[120]=-2.*x1*x3*y[1]*y[4];
y[121]=y[115]+y[116]+y[117]+y[118]+y[119]+y[120];
y[122]=-(lambda*MYI*x0*y[20]*y[23]*y[121]);
y[123]=-(lambda*MYI*y[20]*y[23]*y[54]);
y[124]=lambda*MYI*x0*y[23]*y[54];
y[125]=1.+y[122]+y[123]+y[124];
y[126]=2.*x0*y[1]*y[9];
y[127]=2.*y[1]*y[9]*y[18];
y[128]=2.*x0*x2*y[1]*y[9];
y[129]=2.*x3*y[1]*y[9]*y[18];
y[130]=y[126]+y[127]+y[128]+y[129];
y[131]=-(lambda*MYI*x1*y[15]*y[55]*y[130]);
y[132]=-(lambda*MYI*y[15]*y[55]*y[108]);
y[133]=lambda*MYI*x1*y[55]*y[108];
y[134]=1.+y[131]+y[132]+y[133];
y[135]=-x3;
y[136]=1.+y[135];
y[137]=x1*y[1]*y[2];
y[138]=x1*x3*y[1]*y[2];
y[139]=y[26]+y[29]+y[30]+y[33]+y[36]+y[39]+y[41]+y[45]+y[47]+y[58]+y[137]+y[\
138];
y[140]=x0*y[1]*y[2];
y[141]=x0*x3*y[1]*y[2];
y[142]=y[6]+y[7]+y[8]+y[10]+y[11]+y[16]+y[17]+y[91]+y[92]+y[93]+y[140]+y[141\
];
y[143]=x0*x1*y[15]*y[20]*y[23]*y[55]*y[59]*y[85]*y[90];
y[144]=-(lambda*MYI*x1*y[15]*y[55]*y[68]*y[125]);
y[145]=y[143]+y[144];
y[146]=2.*x1*x2*y[1]*y[2];
y[147]=x3*y[1]*y[2];
y[148]=x0*x1*x3*y[1]*y[2];
y[149]=2.*x2*x3*y[1]*y[2];
y[150]=x0*y[1]*y[9]*y[21];
y[151]=x0*x1*x3*y[1]*y[9];
y[152]=x0*x1*x3*y[1]*y[3];
y[153]=x0*y[1]*y[3]*y[22];
y[154]=-(x0*x1*y[1]*y[4]);
y[155]=y[6]+y[7]+y[8]+y[10]+y[17]+y[24]+y[36]+y[58]+y[69]+y[106]+y[110]+y[11\
1]+y[112]+y[113]+y[137]+y[146]+y[147]+y[148]+y[149]+y[150]+y[151]+y[152]+y[\
153]+y[154];
y[156]=lrs[3];
y[157]=x0*x1*y[15]*y[20]*y[23]*y[55]*y[59]*y[90]*y[142];
y[158]=-(x0*x1*y[15]*y[20]*y[23]*y[55]*y[59]*y[68]*y[139]);
y[159]=y[157]+y[158];
y[160]=-(x0*x1*y[15]*y[20]*y[23]*y[55]*y[59]*y[68]*y[85]);
y[161]=lambda*MYI*x0*y[20]*y[23]*y[90]*y[134];
y[162]=y[160]+y[161];
y[163]=2.*y[1]*y[2];
y[164]=2.*x1*y[1]*y[2];
y[165]=2.*x3*y[1]*y[2];
y[166]=y[163]+y[164]+y[165];
y[167]=-(lambda*MYI*x2*y[57]*y[109]*y[166]);
y[168]=-(lambda*MYI*y[57]*y[109]*y[155]);
y[169]=lambda*MYI*x2*y[109]*y[155];
y[170]=1.+y[167]+y[168]+y[169];
y[171]=x0*x1*y[15]*y[20]*y[23]*y[55]*y[59]*y[85]*y[139];
y[172]=-(lambda*MYI*x1*y[15]*y[55]*y[125]*y[142]);
y[173]=y[171]+y[172];
y[174]=-(x0*x1*y[15]*y[20]*y[23]*y[55]*y[59]*y[85]*y[142]);
y[175]=lambda*MYI*x0*y[20]*y[23]*y[134]*y[139];
y[176]=y[174]+y[175];
y[177]=pow(y[85],2);
y[178]=x0*x1*y[15]*y[20]*y[23]*y[55]*y[59]*y[177];
y[179]=y[125]*y[134];
y[180]=y[178]+y[179];
y[181]=2.*x1*y[1]*y[3]*y[18];
y[182]=2.*x0*x2*y[1]*y[3];
y[183]=x0*x1*y[1]*y[5];
y[184]=x0*x1*x2*y[1]*y[2];
y[185]=y[1]*y[9]*y[18]*y[21];
y[186]=x0*x1*x2*y[1]*y[9];
y[187]=x0*x1*x2*y[1]*y[3];
y[188]=2.*x1*x3*y[1]*y[3]*y[18];
y[189]=2.*x0*x2*x3*y[1]*y[3];
y[190]=-(x1*y[1]*y[4]*y[18]);
y[191]=y[6]+y[11]+y[12]+y[13]+y[70]+y[94]+y[97]+y[105]+y[111]+y[112]+y[113]+\
y[154]+y[181]+y[182]+y[183]+y[184]+y[185]+y[186]+y[187]+y[188]+y[189]+y[190\
];
y[192]=-(lambda*MYI*x1*y[15]*y[55]*y[108]);
y[193]=x1+y[192];
y[194]=-(lambda*MYI*x2*y[57]*y[109]*y[155]);
y[195]=x2+y[194];
y[196]=-(lambda*MYI*x0*y[20]*y[23]*y[54]);
y[197]=x0+y[196];
y[198]=-(lambda*MYI*x3*y[136]*y[156]*y[191]);
y[199]=x3+y[198];
y[200]=pow(y[197],2);
y[201]=pow(y[193],2);
y[202]=pow(y[195],2);
y[203]=pow(y[199],2);
FOUT=(pow(bi,-2)*(-(lambda*MYI*x3*y[68]*y[136]*y[156]*(-(lambda*MYI*x2*y[57]\
*y[109]*y[139]*y[159])-y[145]*y[170]-lambda*MYI*x2*y[57]*y[109]*y[114]*y[17\
3]))+lambda*MYI*x3*y[90]*y[136]*y[156]*(-(lambda*MYI*x2*y[57]*y[109]*y[142]\
*y[159])-y[162]*y[170]-lambda*MYI*x2*y[57]*y[109]*y[114]*y[176])+lambda*MYI\
*x3*y[114]*y[136]*y[156]*(lambda*MYI*x2*y[57]*y[109]*y[142]*y[145]-lambda*M\
YI*x2*y[57]*y[109]*y[139]*y[162]-lambda*MYI*x2*y[57]*y[109]*y[114]*y[180])+\
(lambda*MYI*x2*y[57]*y[109]*y[142]*y[173]-lambda*MYI*x2*y[57]*y[109]*y[139]\
*y[176]+y[170]*y[180])*(1.-lambda*MYI*x3*y[136]*y[156]*(y[12]+y[181]+y[182]\
)+lambda*MYI*x3*y[156]*y[191]-lambda*MYI*y[136]*y[156]*y[191])))/((y[1]+y[1\
]*y[193]+y[1]*y[195]+y[1]*y[193]*y[195]+y[1]*y[193]*y[197]+y[1]*y[199]+y[1]\
*y[195]*y[199]+y[1]*y[193]*y[197]*y[199])*(y[6]+y[1]*y[5]*y[193]+y[1]*y[2]*\
y[195]-y[1]*y[4]*y[195]+y[1]*y[5]*y[195]+y[1]*y[2]*y[193]*y[195]-y[1]*y[4]*\
y[193]*y[195]+y[1]*y[5]*y[193]*y[195]+y[1]*y[3]*y[197]+y[1]*y[3]*y[193]*y[1\
97]+y[1]*y[5]*y[193]*y[197]+y[1]*y[9]*y[193]*y[197]+y[1]*y[3]*y[195]*y[197]\
+y[1]*y[2]*y[193]*y[195]*y[197]+y[1]*y[3]*y[193]*y[195]*y[197]-y[1]*y[4]*y[\
193]*y[195]*y[197]+y[1]*y[9]*y[193]*y[195]*y[197]+y[1]*y[5]*y[199]+y[1]*y[2\
]*y[195]*y[199]-y[1]*y[4]*y[195]*y[199]+y[1]*y[5]*y[195]*y[199]+2.*y[1]*y[3\
]*y[197]*y[199]-y[1]*y[4]*y[197]*y[199]+y[1]*y[3]*y[193]*y[197]*y[199]-y[1]\
*y[4]*y[193]*y[197]*y[199]+y[1]*y[5]*y[193]*y[197]*y[199]+y[1]*y[9]*y[193]*\
y[197]*y[199]+2.*y[1]*y[3]*y[195]*y[197]*y[199]-y[1]*y[4]*y[195]*y[197]*y[1\
99]+y[1]*y[2]*y[193]*y[195]*y[197]*y[199]+y[1]*y[3]*y[193]*y[195]*y[197]*y[\
199]+y[1]*y[9]*y[193]*y[195]*y[197]*y[199]+y[1]*y[3]*y[193]*y[200]+2.*y[1]*\
y[3]*y[193]*y[199]*y[200]-y[1]*y[4]*y[193]*y[199]*y[200]+y[1]*y[9]*y[197]*y\
[201]+y[1]*y[9]*y[195]*y[197]*y[201]+y[1]*y[9]*y[200]*y[201]+y[1]*y[9]*y[19\
9]*y[200]*y[201]+y[1]*y[2]*y[202]+y[1]*y[2]*y[193]*y[202]+y[1]*y[2]*y[199]*\
y[202]+y[1]*y[3]*y[197]*y[203]+y[1]*y[3]*y[195]*y[197]*y[203]+y[1]*y[3]*y[1\
93]*y[200]*y[203]));
return (FOUT);
}
