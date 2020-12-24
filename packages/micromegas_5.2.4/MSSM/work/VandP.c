#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../CalcHEP_src/include/extern.h"
#include "../../CalcHEP_src/include/VandP.h"
#include "autoprot.h"
extern int  FError;
/*  Special model functions  */

int nModelParticles=48;
static ModelPrtclsStr ModelPrtcls_[48]=
{
  {"A","A", 22, "0","0",2,1,0}
, {"Z","Z", 23, "MZ","wZ",2,1,0}
, {"W+","W-", 24, "MW","wW",2,1,3}
, {"G","G", 21, "0","0",2,8,0}
, {"ne","Ne", 12, "0","0",1,1,0}
, {"e","E", 11, "0","0",1,1,-3}
, {"nm","Nm", 14, "0","0",1,1,0}
, {"m","M", 13, "0","0",1,1,-3}
, {"nl","Nl", 16, "0","0",1,1,0}
, {"l","L", 15, "Ml","0",1,1,-3}
, {"u","U", 2, "Mq","0",1,3,2}
, {"d","D", 1, "Mq","0",1,3,-1}
, {"s","S", 3, "Mq","0",1,3,-1}
, {"c","C", 4, "Mc","0",1,3,2}
, {"t","T", 6, "Mt","wt",1,3,2}
, {"b","B", 5, "Mb","0",1,3,-1}
, {"h","h", 25, "Mh","wh",0,1,0}
, {"H","H", 35, "MH","wH",0,1,0}
, {"H3","H3", 36, "MH3","wH3",0,1,0}
, {"H+","H-", 37, "MHc","wHc",0,1,3}
, {"~1+","~1-", 1000024, "MC1","wC1",1,1,3}
, {"~2+","~2-", 1000037, "MC2","wC2",1,1,3}
, {"~o1","~o1", 1000022, "MNE1","0",1,1,0}
, {"~o2","~o2", 1000023, "MNE2","wNE2",1,1,0}
, {"~o3","~o3", 1000025, "MNE3","wNE3",1,1,0}
, {"~o4","~o4", 1000035, "MNE4","wNE4",1,1,0}
, {"~g","~g", 1000021, "MSG","wSG",1,8,0}
, {"~eL","~EL", 1000011, "MSeL","wSeL",0,1,-3}
, {"~eR","~ER", 2000011, "MSeR","wSeR",0,1,-3}
, {"~mL","~ML", 1000013, "MSmL","wSmL",0,1,-3}
, {"~mR","~MR", 2000013, "MSmR","wSmR",0,1,-3}
, {"~l1","~L1", 1000015, "MSl1","wSl1",0,1,-3}
, {"~l2","~L2", 2000015, "MSl2","wSl2",0,1,-3}
, {"~ne","~Ne", 1000012, "MSne","wSne",0,1,0}
, {"~nm","~Nm", 1000014, "MSnm","wSnm",0,1,0}
, {"~nl","~Nl", 1000016, "MSnl","wSnl",0,1,0}
, {"~uL","~UL", 1000002, "MSuL","wSuL",0,3,2}
, {"~uR","~UR", 2000002, "MSuR","wSuR",0,3,2}
, {"~dL","~DL", 1000001, "MSdL","wSdL",0,3,-1}
, {"~dR","~DR", 2000001, "MSdR","wSdR",0,3,-1}
, {"~cL","~CL", 1000004, "MScL","wScL",0,3,2}
, {"~cR","~CR", 2000004, "MScR","wScR",0,3,2}
, {"~sL","~SL", 1000003, "MSsL","wSsL",0,3,-1}
, {"~sR","~SR", 2000003, "MSsR","wSsR",0,3,-1}
, {"~t1","~T1", 1000006, "MSt1","wSt1",0,3,2}
, {"~t2","~T2", 2000006, "MSt2","wSt2",0,3,2}
, {"~b1","~B1", 1000005, "MSb1","wSb1",0,3,-1}
, {"~b2","~B2", 2000005, "MSb2","wSb2",0,3,-1}
};
ModelPrtclsStr *ModelPrtcls=ModelPrtcls_; 
int nModelVars=38;
int nModelFunc=233;
static int nCurrentVars=37;
int*currentVarPtr=&nCurrentVars;
static char*varNames_[271]={
 "Am","alfSMZ","MZ","MW","Mtp","MbMb","McMc","Q","EE","Ml"
,"Mq","Au","Ad","tb","MH3","mu","MG1","MG2","MG3","Ml1"
,"Ml2","Ml3","Mr1","Mr2","Mr3","Mq1","Mq2","Mq3","Mu1","Mu2"
,"Mu3","Md1","Md2","Md3","At","Ab","Al","Maux","Ms","CW"
,"SW","S2W","C2W","VEV","GF","Lqcd","Mbp","Mcp","tB","sb"
,"cb","c2b","MSne","MSnm","MSeL","MSeR","MSmL","MSmR","MSdL","MSdR"
,"MSuL","MSuR","MSsL","MSsR","MScL","MScR","MSnl","MSl1","MSl2","MSb1"
,"MSb2","MSt1","MSt2","QSUSY","Zt11","Zt12","Zt21","Zt22","Zb11","Zb12"
,"Zb21","Zb22","Zl11","Zl12","Zl21","Zl22","MbMM","MtMM","MSG","MNE1"
,"MNE2","MNE3","MNE4","MC1","MC2","Zn11","Zn12","Zn13","Zn14","Zn21"
,"Zn22","Zn23","Zn24","Zn31","Zn32","Zn33","Zn34","Zn41","Zn42","Zn43"
,"Zn44","Zu11","Zu12","Zu21","Zu22","Zv11","Zv12","Zv21","Zv22","Mh"
,"MH","MHc","alpha","ca","sa","calcL","la1","la3","la6","Zh11"
,"Zh12","Zh21","Zh22","Mhh11","Mhh12","Mhh22","la5","la4","la7","la2"
,"dMb","dMl","Td3","Tl3","Mb","Mt","Mc","dMd","Td2","fiuu"
,"fidd","ficc","fiss","Zuu11","Zuu12","Zuu21","Zuu22","Zdd11","Zdd12","Zdd21"
,"Zdd22","Zcc11","Zcc12","Zcc21","Zcc22","Zss11","Zss12","Zss21","Zss22","nmm11"
,"nmm12","nmm13","nmm14","nmm22","nmm23","nmm24","nmm33","nmm34","nmm44","NMM34"
,"NMM23","NMM13","NMM14","NMM24","ahF_c","ahF_b","ahF_t","ahF_l","aHF_c","aHF_b"
,"aHF_t","aHF_l","ahS_eL","ahS_eR","ahS_mL","ahS_mR","ahS_uL","ahS_uR","ahS_cL","ahS_cR"
,"ahS_dL","ahS_dR","ahS_sL","ahS_sR","ahS_l1","ahS_l2","ahS_t1","ahS_t2","ahS_b1","ahS_b2"
,"aHS_eL","aHS_eR","aHS_mL","aHS_mR","aHS_uL","aHS_uR","aHS_cL","aHS_cR","aHS_dL","aHS_dR"
,"aHS_sL","aHS_sR","aHS_l1","aHS_l2","aHS_t1","aHS_t2","aHS_b1","aHS_b2","ahV_W","B00000"
,"B00001","B00002","ahS_Hc","aHV_W","B00003","B00004","B00005","aHS_Hc","ahF_c1","ahF_c2"
,"aHF_c1","aHF_c2","aAF_c","aAF_b","aAF_t","aAF_l","aAF_c1","aAF_c2","ahF_b0","aHF_b0"
,"aAF_b0","ahF_l0","aHF_l0","aAF_l0","alphaE0","Quq","Qdq","aSMhF_f","aSMhV_W","aQCDh"
,"Rqcdh","aQCDH","RqcdH","aQCDH3","RqcdH3","LGGh","LGGH","LAAh","LAAH","LGGH3"
,"LAAH3"};
char**varNames=varNames_;
static REAL varValues_[271]={
   0.000000E+00,  1.184000E-01,  9.118700E+01,  8.020000E+01,  1.730700E+02,  4.230000E+00,  1.270000E+00,  1.000000E+02,  3.123000E-01,  1.777000E+00
,  5.000000E-02,  0.000000E+00,  0.000000E+00,  1.000000E+01,  1.000000E+03,  3.500000E+02,  2.000000E+02,  4.000000E+02,  8.000000E+02,  5.000000E+02
,  5.000000E+02,  5.000000E+02,  2.000000E+02,  2.000000E+02,  2.000000E+02,  1.000000E+03,  1.000000E+03,  1.000000E+03,  3.000000E+02,  3.000000E+02
,  3.000000E+02,  3.000000E+02,  3.000000E+02,  3.000000E+02, -1.000000E+03,  0.000000E+00,  0.000000E+00,  1.000000E+00};
REAL*varValues=varValues_;
int calcMainFunc(void)
{
   int i;
   static REAL * VV=NULL;
   static int iQ=-1;
   static int cErr=1;
   REAL *V=varValues;
   FError=0;
   if(VV && cErr==0)
   { for(i=0;i<nModelVars;i++) if(i!=iQ && VV[i]!=V[i]) break;
     if(i==nModelVars)      {if(iQ>=0 && VV[iQ]!=V[iQ]) goto FirstQ; else return 0;} 
   }
  cErr=1;
   nCurrentVars=38;
   V[38]=V[10];

   nCurrentVars=39;
   V[39]=V[3]/(V[2]);
   if(!isfinite(V[39]) || FError) return 39;
   nCurrentVars=40;
   V[40]=Sqrt(1-Pow(V[39],2));
   if(!isfinite(V[40]) || FError) return 40;
   nCurrentVars=41;
   V[41]=2*V[40]*V[39];
   if(!isfinite(V[41]) || FError) return 41;
   nCurrentVars=42;
   V[42]=Pow(V[39],2)-Pow(V[40],2);
   if(!isfinite(V[42]) || FError) return 42;
   nCurrentVars=43;
   V[43]=2*V[3]*V[40]/(V[8]);
   if(!isfinite(V[43]) || FError) return 43;
   nCurrentVars=44;
   V[44]=Pow(V[8],2)/(Pow(2*V[40]*V[3],2))/(M_SQRT2);
   if(!isfinite(V[44]) || FError) return 44;
   nCurrentVars=45;
   V[45]=initQCD5(V[1],V[6],V[5],V[4]);
   if(!isfinite(V[45]) || FError) return 45;
   nCurrentVars=46;
   V[46]=bPoleMass();
   if(!isfinite(V[46]) || FError) return 46;
   nCurrentVars=47;
   V[47]=V[6]*(1+4/(double)((3))*alphaQCD(V[6])/(M_PI));
   if(!isfinite(V[47]) || FError) return 47;
   nCurrentVars=48;
   V[48]=slhaVal("HMIX",0,1,2);
   if(!isfinite(V[48]) || FError) return 48;
   nCurrentVars=49;
   V[49]=V[48]/(Sqrt(1+Pow(V[48],2)));
   if(!isfinite(V[49]) || FError) return 49;
   nCurrentVars=50;
   V[50]=Sqrt(1-Pow(V[49],2));
   if(!isfinite(V[50]) || FError) return 50;
   nCurrentVars=51;
   V[51]=Pow(V[50],2)-Pow(V[49],2);
   if(!isfinite(V[51]) || FError) return 51;
   nCurrentVars=52;
   V[52]=slhaVal("MASS",V[2],1,1000012);
   if(!isfinite(V[52]) || FError) return 52;
   nCurrentVars=53;
   V[53]=slhaVal("MASS",V[2],1,1000014);
   if(!isfinite(V[53]) || FError) return 53;
   nCurrentVars=54;
   V[54]=slhaVal("MASS",V[2],1,1000011);
   if(!isfinite(V[54]) || FError) return 54;
   nCurrentVars=55;
   V[55]=slhaVal("MASS",V[2],1,2000011);
   if(!isfinite(V[55]) || FError) return 55;
   nCurrentVars=56;
   V[56]=slhaVal("MASS",V[2],1,1000013);
   if(!isfinite(V[56]) || FError) return 56;
   nCurrentVars=57;
   V[57]=slhaVal("MASS",V[2],1,2000013);
   if(!isfinite(V[57]) || FError) return 57;
   nCurrentVars=58;
   V[58]=slhaVal("MASS",V[2],1,1000001);
   if(!isfinite(V[58]) || FError) return 58;
   nCurrentVars=59;
   V[59]=slhaVal("MASS",V[2],1,2000001);
   if(!isfinite(V[59]) || FError) return 59;
   nCurrentVars=60;
   V[60]=slhaVal("MASS",V[2],1,1000002);
   if(!isfinite(V[60]) || FError) return 60;
   nCurrentVars=61;
   V[61]=slhaVal("MASS",V[2],1,2000002);
   if(!isfinite(V[61]) || FError) return 61;
   nCurrentVars=62;
   V[62]=slhaVal("MASS",V[2],1,1000003);
   if(!isfinite(V[62]) || FError) return 62;
   nCurrentVars=63;
   V[63]=slhaVal("MASS",V[2],1,2000003);
   if(!isfinite(V[63]) || FError) return 63;
   nCurrentVars=64;
   V[64]=slhaVal("MASS",V[2],1,1000004);
   if(!isfinite(V[64]) || FError) return 64;
   nCurrentVars=65;
   V[65]=slhaVal("MASS",V[2],1,2000004);
   if(!isfinite(V[65]) || FError) return 65;
   nCurrentVars=66;
   V[66]=slhaVal("MASS",V[2],1,1000016);
   if(!isfinite(V[66]) || FError) return 66;
   nCurrentVars=67;
   V[67]=slhaVal("MASS",V[2],1,1000015);
   if(!isfinite(V[67]) || FError) return 67;
   nCurrentVars=68;
   V[68]=slhaVal("MASS",V[2],1,2000015);
   if(!isfinite(V[68]) || FError) return 68;
   nCurrentVars=69;
   V[69]=slhaVal("MASS",V[2],1,1000005);
   if(!isfinite(V[69]) || FError) return 69;
   nCurrentVars=70;
   V[70]=slhaVal("MASS",V[2],1,2000005);
   if(!isfinite(V[70]) || FError) return 70;
   nCurrentVars=71;
   V[71]=slhaVal("MASS",V[2],1,1000006);
   if(!isfinite(V[71]) || FError) return 71;
   nCurrentVars=72;
   V[72]=slhaVal("MASS",V[2],1,2000006);
   if(!isfinite(V[72]) || FError) return 72;
   nCurrentVars=73;
   V[73]=Sqrt(V[71]*V[72]);
   if(!isfinite(V[73]) || FError) return 73;
   nCurrentVars=74;
   V[74]=slhaVal("STOPMIX",V[73],2,1,1);
   if(!isfinite(V[74]) || FError) return 74;
   nCurrentVars=75;
   V[75]=slhaVal("STOPMIX",V[73],2,1,2);
   if(!isfinite(V[75]) || FError) return 75;
   nCurrentVars=76;
   V[76]=slhaVal("STOPMIX",V[73],2,2,1);
   if(!isfinite(V[76]) || FError) return 76;
   nCurrentVars=77;
   V[77]=slhaVal("STOPMIX",V[73],2,2,2);
   if(!isfinite(V[77]) || FError) return 77;
   nCurrentVars=78;
   V[78]=slhaVal("SBOTMIX",V[73],2,1,1);
   if(!isfinite(V[78]) || FError) return 78;
   nCurrentVars=79;
   V[79]=slhaVal("SBOTMIX",V[73],2,1,2);
   if(!isfinite(V[79]) || FError) return 79;
   nCurrentVars=80;
   V[80]=slhaVal("SBOTMIX",V[73],2,2,1);
   if(!isfinite(V[80]) || FError) return 80;
   nCurrentVars=81;
   V[81]=slhaVal("SBOTMIX",V[73],2,2,2);
   if(!isfinite(V[81]) || FError) return 81;
   nCurrentVars=82;
   V[82]=slhaVal("STAUMIX",V[73],2,1,1);
   if(!isfinite(V[82]) || FError) return 82;
   nCurrentVars=83;
   V[83]=slhaVal("STAUMIX",V[73],2,1,2);
   if(!isfinite(V[83]) || FError) return 83;
   nCurrentVars=84;
   V[84]=slhaVal("STAUMIX",V[73],2,2,1);
   if(!isfinite(V[84]) || FError) return 84;
   nCurrentVars=85;
   V[85]=slhaVal("STAUMIX",V[73],2,2,2);
   if(!isfinite(V[85]) || FError) return 85;
   nCurrentVars=86;
   V[86]=MbRun(V[73]);
   if(!isfinite(V[86]) || FError) return 86;
   nCurrentVars=87;
   V[87]=MtRun(V[73]);
   if(!isfinite(V[87]) || FError) return 87;
   nCurrentVars=88;
   V[88]=slhaVal("MASS",V[2],1,1000021);
   if(!isfinite(V[88]) || FError) return 88;
   nCurrentVars=89;
   V[89]=slhaVal("MASS",V[2],1,1000022);
   if(!isfinite(V[89]) || FError) return 89;
   nCurrentVars=90;
   V[90]=slhaVal("MASS",V[2],1,1000023);
   if(!isfinite(V[90]) || FError) return 90;
   nCurrentVars=91;
   V[91]=slhaVal("MASS",V[2],1,1000025);
   if(!isfinite(V[91]) || FError) return 91;
   nCurrentVars=92;
   V[92]=slhaVal("MASS",V[2],1,1000035);
   if(!isfinite(V[92]) || FError) return 92;
   nCurrentVars=93;
   V[93]=slhaVal("MASS",V[2],1,1000024);
   if(!isfinite(V[93]) || FError) return 93;
   nCurrentVars=94;
   V[94]=slhaVal("MASS",V[2],1,1000037);
   if(!isfinite(V[94]) || FError) return 94;
   nCurrentVars=95;
   V[95]=slhaVal("NMIX",V[73],2,1,1);
   if(!isfinite(V[95]) || FError) return 95;
   nCurrentVars=96;
   V[96]=slhaVal("NMIX",V[73],2,1,2);
   if(!isfinite(V[96]) || FError) return 96;
   nCurrentVars=97;
   V[97]=slhaVal("NMIX",V[73],2,1,3);
   if(!isfinite(V[97]) || FError) return 97;
   nCurrentVars=98;
   V[98]=slhaVal("NMIX",V[73],2,1,4);
   if(!isfinite(V[98]) || FError) return 98;
   nCurrentVars=99;
   V[99]=slhaVal("NMIX",V[73],2,2,1);
   if(!isfinite(V[99]) || FError) return 99;
   nCurrentVars=100;
   V[100]=slhaVal("NMIX",V[73],2,2,2);
   if(!isfinite(V[100]) || FError) return 100;
   nCurrentVars=101;
   V[101]=slhaVal("NMIX",V[73],2,2,3);
   if(!isfinite(V[101]) || FError) return 101;
   nCurrentVars=102;
   V[102]=slhaVal("NMIX",V[73],2,2,4);
   if(!isfinite(V[102]) || FError) return 102;
   nCurrentVars=103;
   V[103]=slhaVal("NMIX",V[73],2,3,1);
   if(!isfinite(V[103]) || FError) return 103;
   nCurrentVars=104;
   V[104]=slhaVal("NMIX",V[73],2,3,2);
   if(!isfinite(V[104]) || FError) return 104;
   nCurrentVars=105;
   V[105]=slhaVal("NMIX",V[73],2,3,3);
   if(!isfinite(V[105]) || FError) return 105;
   nCurrentVars=106;
   V[106]=slhaVal("NMIX",V[73],2,3,4);
   if(!isfinite(V[106]) || FError) return 106;
   nCurrentVars=107;
   V[107]=slhaVal("NMIX",V[73],2,4,1);
   if(!isfinite(V[107]) || FError) return 107;
   nCurrentVars=108;
   V[108]=slhaVal("NMIX",V[73],2,4,2);
   if(!isfinite(V[108]) || FError) return 108;
   nCurrentVars=109;
   V[109]=slhaVal("NMIX",V[73],2,4,3);
   if(!isfinite(V[109]) || FError) return 109;
   nCurrentVars=110;
   V[110]=slhaVal("NMIX",V[73],2,4,4);
   if(!isfinite(V[110]) || FError) return 110;
   nCurrentVars=111;
   V[111]=slhaVal("UMIX",V[73],2,1,1);
   if(!isfinite(V[111]) || FError) return 111;
   nCurrentVars=112;
   V[112]=slhaVal("UMIX",V[73],2,1,2);
   if(!isfinite(V[112]) || FError) return 112;
   nCurrentVars=113;
   V[113]=slhaVal("UMIX",V[73],2,2,1);
   if(!isfinite(V[113]) || FError) return 113;
   nCurrentVars=114;
   V[114]=slhaVal("UMIX",V[73],2,2,2);
   if(!isfinite(V[114]) || FError) return 114;
   nCurrentVars=115;
   V[115]=slhaVal("VMIX",V[73],2,1,1);
   if(!isfinite(V[115]) || FError) return 115;
   nCurrentVars=116;
   V[116]=slhaVal("VMIX",V[73],2,1,2);
   if(!isfinite(V[116]) || FError) return 116;
   nCurrentVars=117;
   V[117]=slhaVal("VMIX",V[73],2,2,1);
   if(!isfinite(V[117]) || FError) return 117;
   nCurrentVars=118;
   V[118]=slhaVal("VMIX",V[73],2,2,2);
   if(!isfinite(V[118]) || FError) return 118;
   nCurrentVars=119;
   V[119]=slhaVal("MASS",V[2],1,25);
   if(!isfinite(V[119]) || FError) return 119;
   nCurrentVars=120;
   V[120]=slhaVal("MASS",V[2],1,35);
   if(!isfinite(V[120]) || FError) return 120;
   nCurrentVars=121;
   V[121]=slhaVal("MASS",V[2],1,37);
   if(!isfinite(V[121]) || FError) return 121;
   nCurrentVars=122;
   V[122]=slhaVal("ALPHA",V[73],0);
   if(!isfinite(V[122]) || FError) return 122;
   nCurrentVars=123;
   V[123]=Cos(V[122]);
   if(!isfinite(V[123]) || FError) return 123;
   nCurrentVars=124;
   V[124]=Sin(V[122]);
   if(!isfinite(V[124]) || FError) return 124;
   nCurrentVars=125;
   V[125]=calcLambdas(V[48],V[14],V[15],V[73],V[4],V[88],V[94],V[27],V[30],V[33],V[35],V[34]);
   if(!isfinite(V[125]) || FError) return 125;
   nCurrentVars=126;
   V[126]=Lambda1();
   if(!isfinite(V[126]) || FError) return 126;
   nCurrentVars=127;
   V[127]=Lambda3();
   if(!isfinite(V[127]) || FError) return 127;
   nCurrentVars=128;
   V[128]=-Lambda6();
   if(!isfinite(V[128]) || FError) return 128;
   nCurrentVars=129;
   V[129]=V[123];

   nCurrentVars=130;
   V[130]=V[124];

   nCurrentVars=131;
   V[131]=-V[124];
   if(!isfinite(V[131]) || FError) return 131;
   nCurrentVars=132;
   V[132]=V[123];

   nCurrentVars=133;
   V[133]=Pow(V[119],2)*Pow(V[130],2)+Pow(V[120],2)*Pow(V[129],2);
   if(!isfinite(V[133]) || FError) return 133;
   nCurrentVars=134;
   V[134]=Pow(V[120],2)*V[129]*V[130]+Pow(V[119],2)*V[132]*V[131];
   if(!isfinite(V[134]) || FError) return 134;
   nCurrentVars=135;
   V[135]=Pow(V[119],2)*Pow(V[132],2)+Pow(V[120],2)*Pow(V[131],2);
   if(!isfinite(V[135]) || FError) return 135;
   nCurrentVars=136;
   V[136]=2*V[128]*V[50]/(V[49])-V[126]*Pow(V[50],2)/(Pow(V[49],2))-Pow(V[8],2)*(Pow(V[14]*V[49],2)-V[133])/(Pow(2*V[3]*V[40]*V[49],2));
   if(!isfinite(V[136]) || FError) return 136;
   nCurrentVars=137;
   V[137]=V[136]-1/(double)((2))*Pow(V[8],2)*(Pow(V[121],2)-Pow(V[14],2))/(Pow(V[3]*V[40],2));
   if(!isfinite(V[137]) || FError) return 137;
   nCurrentVars=138;
   V[138]=V[50]/(V[49])*(V[127]+V[137]-V[50]/(V[49])*V[128]-1/(double)((4))*Pow(V[8],2)*(V[134]/(V[50])/(V[49])+Pow(V[14],2))/(Pow(V[3]*V[40],2)));
   if(!isfinite(V[138]) || FError) return 138;
   nCurrentVars=139;
   V[139]=2*V[138]*V[50]/(V[49])-V[136]*Pow(V[50],2)/(Pow(V[49],2))-Pow(V[8],2)*(Pow(V[14]*V[50],2)-V[135])/(Pow(2*V[3]*V[40]*V[49],2));
   if(!isfinite(V[139]) || FError) return 139;
   nCurrentVars=140;
   V[140]=deltaMb();
   if(!isfinite(V[140]) || FError) return 140;
   nCurrentVars=141;
   V[141]=deltaMl();
   if(!isfinite(V[141]) || FError) return 141;
   nCurrentVars=142;
   V[142]=V[140]/(1+V[140]);
   if(!isfinite(V[142]) || FError) return 142;
   nCurrentVars=143;
   V[143]=V[141]/(1+V[141]);
   if(!isfinite(V[143]) || FError) return 143;
 FirstQ:
 cErr=1;
   nCurrentVars=144;
   V[144]=MbEff(V[7]);
   if(!isfinite(V[144]) || FError) return 144;
   nCurrentVars=145;
   V[145]=MtEff(V[7]);
   if(!isfinite(V[145]) || FError) return 145;
   nCurrentVars=146;
   V[146]=McEff(V[7]);
   if(!isfinite(V[146]) || FError) return 146;
   nCurrentVars=147;
   V[147]=deltaMd();
   if(!isfinite(V[147]) || FError) return 147;
   nCurrentVars=148;
   V[148]=V[147]/(1+V[147]);
   if(!isfinite(V[148]) || FError) return 148;
   nCurrentVars=149;
   V[149]=Atan(-2*V[10]*(V[11]-V[15]/(V[48]))/(Pow(V[61],2)-Pow(V[60],2)))/(2);
   if(!isfinite(V[149]) || FError) return 149;
   nCurrentVars=150;
   V[150]=Atan(-2*V[10]*(1-V[148])*(V[12]-V[15]*V[48])/(Pow(V[59],2)-Pow(V[58],2)))/(2);
   if(!isfinite(V[150]) || FError) return 150;
   nCurrentVars=151;
   V[151]=Atan(-2*V[146]*(V[11]-V[15]/(V[48]))/(Pow(V[65],2)-Pow(V[64],2)))/(2);
   if(!isfinite(V[151]) || FError) return 151;
   nCurrentVars=152;
   V[152]=Atan(-2*V[10]*(1-V[148])*(V[12]-V[15]*V[48])/(Pow(V[63],2)-Pow(V[62],2)))/(2);
   if(!isfinite(V[152]) || FError) return 152;
   nCurrentVars=153;
   V[153]=Cos(V[149]);
   if(!isfinite(V[153]) || FError) return 153;
   nCurrentVars=154;
   V[154]=Sin(V[149]);
   if(!isfinite(V[154]) || FError) return 154;
   nCurrentVars=155;
   V[155]=-Sin(V[149]);
   if(!isfinite(V[155]) || FError) return 155;
   nCurrentVars=156;
   V[156]=Cos(V[149]);
   if(!isfinite(V[156]) || FError) return 156;
   nCurrentVars=157;
   V[157]=Cos(V[150]);
   if(!isfinite(V[157]) || FError) return 157;
   nCurrentVars=158;
   V[158]=Sin(V[150]);
   if(!isfinite(V[158]) || FError) return 158;
   nCurrentVars=159;
   V[159]=-Sin(V[150]);
   if(!isfinite(V[159]) || FError) return 159;
   nCurrentVars=160;
   V[160]=Cos(V[150]);
   if(!isfinite(V[160]) || FError) return 160;
   nCurrentVars=161;
   V[161]=Cos(V[151]);
   if(!isfinite(V[161]) || FError) return 161;
   nCurrentVars=162;
   V[162]=Sin(V[151]);
   if(!isfinite(V[162]) || FError) return 162;
   nCurrentVars=163;
   V[163]=-Sin(V[151]);
   if(!isfinite(V[163]) || FError) return 163;
   nCurrentVars=164;
   V[164]=Cos(V[151]);
   if(!isfinite(V[164]) || FError) return 164;
   nCurrentVars=165;
   V[165]=Cos(V[152]);
   if(!isfinite(V[165]) || FError) return 165;
   nCurrentVars=166;
   V[166]=Sin(V[152]);
   if(!isfinite(V[166]) || FError) return 166;
   nCurrentVars=167;
   V[167]=-Sin(V[152]);
   if(!isfinite(V[167]) || FError) return 167;
   nCurrentVars=168;
   V[168]=Cos(V[152]);
   if(!isfinite(V[168]) || FError) return 168;
   nCurrentVars=169;
   V[169]=V[95]*V[89]*V[95]+V[99]*V[90]*V[99]+V[103]*V[91]*V[103]+V[107]*V[92]*V[107];
   if(!isfinite(V[169]) || FError) return 169;
   nCurrentVars=170;
   V[170]=V[95]*V[89]*V[96]+V[99]*V[90]*V[100]+V[103]*V[91]*V[104]+V[107]*V[92]*V[108];
   if(!isfinite(V[170]) || FError) return 170;
   nCurrentVars=171;
   V[171]=V[95]*V[89]*V[97]+V[99]*V[90]*V[101]+V[103]*V[91]*V[105]+V[107]*V[92]*V[109];
   if(!isfinite(V[171]) || FError) return 171;
   nCurrentVars=172;
   V[172]=V[95]*V[89]*V[98]+V[99]*V[90]*V[102]+V[103]*V[91]*V[106]+V[107]*V[92]*V[110];
   if(!isfinite(V[172]) || FError) return 172;
   nCurrentVars=173;
   V[173]=V[96]*V[89]*V[96]+V[100]*V[90]*V[100]+V[104]*V[91]*V[104]+V[108]*V[92]*V[108];
   if(!isfinite(V[173]) || FError) return 173;
   nCurrentVars=174;
   V[174]=V[96]*V[89]*V[97]+V[100]*V[90]*V[101]+V[104]*V[91]*V[105]+V[108]*V[92]*V[109];
   if(!isfinite(V[174]) || FError) return 174;
   nCurrentVars=175;
   V[175]=V[96]*V[89]*V[98]+V[100]*V[90]*V[102]+V[104]*V[91]*V[106]+V[108]*V[92]*V[110];
   if(!isfinite(V[175]) || FError) return 175;
   nCurrentVars=176;
   V[176]=V[97]*V[89]*V[97]+V[101]*V[90]*V[101]+V[105]*V[91]*V[105]+V[109]*V[92]*V[109];
   if(!isfinite(V[176]) || FError) return 176;
   nCurrentVars=177;
   V[177]=V[97]*V[89]*V[98]+V[101]*V[90]*V[102]+V[105]*V[91]*V[106]+V[109]*V[92]*V[110];
   if(!isfinite(V[177]) || FError) return 177;
   nCurrentVars=178;
   V[178]=V[98]*V[89]*V[98]+V[102]*V[90]*V[102]+V[106]*V[91]*V[106]+V[110]*V[92]*V[110];
   if(!isfinite(V[178]) || FError) return 178;
   nCurrentVars=179;
   V[179]=-V[15];
   if(!isfinite(V[179]) || FError) return 179;
   nCurrentVars=180;
   V[180]=V[2]*V[50]*V[39];
   if(!isfinite(V[180]) || FError) return 180;
   nCurrentVars=181;
   V[181]=-V[2]*V[50]*V[40];
   if(!isfinite(V[181]) || FError) return 181;
   nCurrentVars=182;
   V[182]=V[2]*V[49]*V[40];
   if(!isfinite(V[182]) || FError) return 182;
   nCurrentVars=183;
   V[183]=-V[2]*V[49]*V[39];
   if(!isfinite(V[183]) || FError) return 183;
   nCurrentVars=184;
   V[184]=-V[8]/(V[3])*V[146]/(V[40])*V[132]/(V[49])/(2)/(V[146]);
   if(!isfinite(V[184]) || FError) return 184;
   nCurrentVars=185;
   V[185]=-V[8]/(V[3])*V[144]/(V[40])/(V[50])/(V[49])*(V[49]*V[131]-V[49]*V[142]*V[131]+V[142]*V[132]*V[50])/(2)/(V[144]);
   if(!isfinite(V[185]) || FError) return 185;
   nCurrentVars=186;
   V[186]=-V[8]/(V[3])*V[145]/(V[40])*V[132]/(V[49])/(2)/(V[145]);
   if(!isfinite(V[186]) || FError) return 186;
   nCurrentVars=187;
   V[187]=-V[8]/(V[3])*V[9]/(V[40])/(V[50])/(V[49])*(V[49]*V[131]-V[49]*V[143]*V[131]+V[143]*V[132]*V[50])/(2)/(V[9]);
   if(!isfinite(V[187]) || FError) return 187;
   nCurrentVars=188;
   V[188]=-V[8]/(V[3])*V[146]/(V[40])*V[130]/(V[49])/(2)/(V[146]);
   if(!isfinite(V[188]) || FError) return 188;
   nCurrentVars=189;
   V[189]=-V[8]/(V[3])*V[144]/(V[40])/(V[50])/(V[49])*(V[49]*V[129]-V[49]*V[142]*V[129]+V[142]*V[130]*V[50])/(2)/(V[144]);
   if(!isfinite(V[189]) || FError) return 189;
   nCurrentVars=190;
   V[190]=-V[8]/(V[3])*V[145]/(V[40])*V[130]/(V[49])/(2)/(V[145]);
   if(!isfinite(V[190]) || FError) return 190;
   nCurrentVars=191;
   V[191]=-V[8]/(V[3])*V[9]/(V[40])/(V[50])/(V[49])*(V[49]*V[129]-V[49]*V[143]*V[129]+V[143]*V[130]*V[50])/(2)/(V[9]);
   if(!isfinite(V[191]) || FError) return 191;
   nCurrentVars=192;
   V[192]=1/(Pow(V[39],2))*V[8]*V[3]/(V[40])*(V[42]*V[131]*V[50]-V[42]*V[49]*V[132])/(2)/(Pow(V[54],2));
   if(!isfinite(V[192]) || FError) return 192;
   nCurrentVars=193;
   V[193]=1/(Pow(V[39],2))*V[8]*V[3]*V[40]*(V[131]*V[50]-V[49]*V[132])/(Pow(V[55],2));
   if(!isfinite(V[193]) || FError) return 193;
   nCurrentVars=194;
   V[194]=1/(Pow(V[39],2))*V[8]*V[3]/(V[40])*(V[42]*V[131]*V[50]-V[42]*V[49]*V[132])/(2)/(Pow(V[56],2));
   if(!isfinite(V[194]) || FError) return 194;
   nCurrentVars=195;
   V[195]=1/(Pow(V[39],2))*V[8]*V[3]*V[40]*(V[131]*V[50]-V[49]*V[132])/(Pow(V[57],2));
   if(!isfinite(V[195]) || FError) return 195;
   nCurrentVars=196;
   V[196]=1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[49])*(4*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[131]*Pow(V[153],2)*V[50]-3*V[49]*Pow(V[3],2)*V[131]*Pow(V[153],2)*V[50]-4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*V[132]*Pow(V[153],2)+3*Pow(V[49],2)*Pow(V[3],2)*V[132]*Pow(V[153],2)-4*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[131]*Pow(V[154],2)*V[50]+4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*V[132]*Pow(V[154],2)+6*Pow(V[39],2)*V[10]*V[131]*V[153]*V[154]*V[15]-6*Pow(V[39],2)*V[11]*V[10]*V[132]*V[153]*V[154]-6*Pow(V[39],2)*Pow(V[10],2)*V[132])/(6)/(Pow(V[60],2));
   if(!isfinite(V[196]) || FError) return 196;
   nCurrentVars=197;
   V[197]=1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[49])*(4*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[131]*Pow(V[155],2)*V[50]-3*V[49]*Pow(V[3],2)*V[131]*Pow(V[155],2)*V[50]-4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*V[132]*Pow(V[155],2)+3*Pow(V[49],2)*Pow(V[3],2)*V[132]*Pow(V[155],2)-4*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[131]*Pow(V[156],2)*V[50]+4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*V[132]*Pow(V[156],2)+6*Pow(V[39],2)*V[10]*V[131]*V[155]*V[156]*V[15]-6*Pow(V[39],2)*V[11]*V[10]*V[132]*V[155]*V[156]-6*Pow(V[39],2)*Pow(V[10],2)*V[132])/(6)/(Pow(V[61],2));
   if(!isfinite(V[197]) || FError) return 197;
   nCurrentVars=198;
   V[198]=1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[49])*(4*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[161],2)*V[131]*V[50]-3*V[49]*Pow(V[3],2)*Pow(V[161],2)*V[131]*V[50]-4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*Pow(V[161],2)*V[132]+3*Pow(V[49],2)*Pow(V[3],2)*Pow(V[161],2)*V[132]-4*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[162],2)*V[131]*V[50]+4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*Pow(V[162],2)*V[132]+6*Pow(V[39],2)*V[146]*V[161]*V[162]*V[131]*V[15]-6*Pow(V[39],2)*V[11]*V[146]*V[161]*V[162]*V[132]-6*Pow(V[39],2)*Pow(V[146],2)*V[132])/(6)/(Pow(V[64],2));
   if(!isfinite(V[198]) || FError) return 198;
   nCurrentVars=199;
   V[199]=1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[49])*(4*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[163],2)*V[131]*V[50]-3*V[49]*Pow(V[3],2)*Pow(V[163],2)*V[131]*V[50]-4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*Pow(V[163],2)*V[132]+3*Pow(V[49],2)*Pow(V[3],2)*Pow(V[163],2)*V[132]-4*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[164],2)*V[131]*V[50]+4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*Pow(V[164],2)*V[132]+6*Pow(V[39],2)*V[146]*V[163]*V[164]*V[131]*V[15]-6*Pow(V[39],2)*V[11]*V[146]*V[163]*V[164]*V[132]-6*Pow(V[39],2)*Pow(V[146],2)*V[132])/(6)/(Pow(V[65],2));
   if(!isfinite(V[199]) || FError) return 199;
   nCurrentVars=200;
   V[200]=-1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[50])*(2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*Pow(V[157],2)*V[131]-3*Pow(V[50],2)*Pow(V[3],2)*Pow(V[157],2)*V[131]-2*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[157],2)*V[132]*V[50]+3*V[49]*Pow(V[3],2)*Pow(V[157],2)*V[132]*V[50]-2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*Pow(V[158],2)*V[131]+2*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[158],2)*V[132]*V[50]-6*Pow(V[39],2)*V[10]*V[157]*V[158]*V[132]*V[15]+6*Pow(V[39],2)*V[12]*V[10]*V[157]*V[158]*V[131]+6*Pow(V[39],2)*Pow(V[10],2)*V[131])/(6)/(Pow(V[58],2));
   if(!isfinite(V[200]) || FError) return 200;
   nCurrentVars=201;
   V[201]=-1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[50])*(2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*Pow(V[159],2)*V[131]-3*Pow(V[50],2)*Pow(V[3],2)*Pow(V[159],2)*V[131]-2*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[159],2)*V[132]*V[50]+3*V[49]*Pow(V[3],2)*Pow(V[159],2)*V[132]*V[50]-2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*Pow(V[160],2)*V[131]+2*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[160],2)*V[132]*V[50]-6*Pow(V[39],2)*V[10]*V[159]*V[160]*V[132]*V[15]+6*Pow(V[39],2)*V[12]*V[10]*V[159]*V[160]*V[131]+6*Pow(V[39],2)*Pow(V[10],2)*V[131])/(6)/(Pow(V[59],2));
   if(!isfinite(V[201]) || FError) return 201;
   nCurrentVars=202;
   V[202]=-1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[50])*(2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*V[131]*Pow(V[165],2)-3*Pow(V[50],2)*Pow(V[3],2)*V[131]*Pow(V[165],2)-2*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[132]*Pow(V[165],2)*V[50]+3*V[49]*Pow(V[3],2)*V[132]*Pow(V[165],2)*V[50]-2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*V[131]*Pow(V[166],2)+2*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[132]*Pow(V[166],2)*V[50]-6*Pow(V[39],2)*V[10]*V[132]*V[165]*V[166]*V[15]+6*Pow(V[39],2)*V[12]*V[10]*V[131]*V[165]*V[166]+6*Pow(V[39],2)*Pow(V[10],2)*V[131])/(6)/(Pow(V[62],2));
   if(!isfinite(V[202]) || FError) return 202;
   nCurrentVars=203;
   V[203]=-1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[50])*(2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*V[131]*Pow(V[167],2)-3*Pow(V[50],2)*Pow(V[3],2)*V[131]*Pow(V[167],2)-2*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[132]*Pow(V[167],2)*V[50]+3*V[49]*Pow(V[3],2)*V[132]*Pow(V[167],2)*V[50]-2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*V[131]*Pow(V[168],2)+2*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[132]*Pow(V[168],2)*V[50]-6*Pow(V[39],2)*V[10]*V[132]*V[167]*V[168]*V[15]+6*Pow(V[39],2)*V[12]*V[10]*V[131]*V[167]*V[168]+6*Pow(V[39],2)*Pow(V[10],2)*V[131])/(6)/(Pow(V[63],2));
   if(!isfinite(V[203]) || FError) return 203;
   nCurrentVars=204;
   V[204]=1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[50])*(V[42]*Pow(V[50],2)*Pow(V[3],2)*V[131]*Pow(V[82],2)-V[42]*V[49]*Pow(V[3],2)*V[132]*Pow(V[82],2)*V[50]+2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*V[131]*Pow(V[83],2)-2*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[132]*Pow(V[83],2)*V[50]+2*Pow(V[39],2)*V[9]*V[132]*V[82]*V[83]*V[15]-2*Pow(V[39],2)*V[36]*V[9]*V[131]*V[82]*V[83]-2*Pow(V[39],2)*Pow(V[9],2)*V[131])/(2)/(Pow(V[67],2));
   if(!isfinite(V[204]) || FError) return 204;
   nCurrentVars=205;
   V[205]=1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[50])*(V[42]*Pow(V[50],2)*Pow(V[3],2)*V[131]*Pow(V[84],2)-V[42]*V[49]*Pow(V[3],2)*V[132]*Pow(V[84],2)*V[50]+2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*V[131]*Pow(V[85],2)-2*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[132]*Pow(V[85],2)*V[50]+2*Pow(V[39],2)*V[9]*V[132]*V[84]*V[85]*V[15]-2*Pow(V[39],2)*V[36]*V[9]*V[131]*V[84]*V[85]-2*Pow(V[39],2)*Pow(V[9],2)*V[131])/(2)/(Pow(V[68],2));
   if(!isfinite(V[205]) || FError) return 205;
   nCurrentVars=206;
   V[206]=1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[49])*(4*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[131]*Pow(V[74],2)*V[50]-3*V[49]*Pow(V[3],2)*V[131]*Pow(V[74],2)*V[50]-4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*V[132]*Pow(V[74],2)+3*Pow(V[49],2)*Pow(V[3],2)*V[132]*Pow(V[74],2)-4*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[131]*Pow(V[75],2)*V[50]+4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*V[132]*Pow(V[75],2)+6*Pow(V[39],2)*V[87]*V[131]*V[74]*V[75]*V[15]-6*Pow(V[39],2)*V[34]*V[87]*V[132]*V[74]*V[75]-6*Pow(V[39],2)*Pow(V[87],2)*V[132])/(6)/(Pow(V[71],2));
   if(!isfinite(V[206]) || FError) return 206;
   nCurrentVars=207;
   V[207]=1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[49])*(4*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[131]*Pow(V[76],2)*V[50]-3*V[49]*Pow(V[3],2)*V[131]*Pow(V[76],2)*V[50]-4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*V[132]*Pow(V[76],2)+3*Pow(V[49],2)*Pow(V[3],2)*V[132]*Pow(V[76],2)-4*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[131]*Pow(V[77],2)*V[50]+4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*V[132]*Pow(V[77],2)+6*Pow(V[39],2)*V[87]*V[131]*V[76]*V[77]*V[15]-6*Pow(V[39],2)*V[34]*V[87]*V[132]*V[76]*V[77]-6*Pow(V[39],2)*Pow(V[87],2)*V[132])/(6)/(Pow(V[72],2));
   if(!isfinite(V[207]) || FError) return 207;
   nCurrentVars=208;
   V[208]=-1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[50])*(2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*Pow(V[78],2)*V[131]-3*Pow(V[50],2)*Pow(V[3],2)*Pow(V[78],2)*V[131]-2*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[78],2)*V[132]*V[50]+3*V[49]*Pow(V[3],2)*Pow(V[78],2)*V[132]*V[50]-2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*Pow(V[79],2)*V[131]+2*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[79],2)*V[132]*V[50]-6*Pow(V[39],2)*V[86]*V[78]*V[79]*V[132]*V[15]+6*Pow(V[39],2)*V[35]*V[86]*V[78]*V[79]*V[131]+6*Pow(V[39],2)*Pow(V[86],2)*V[131])/(6)/(Pow(V[69],2));
   if(!isfinite(V[208]) || FError) return 208;
   nCurrentVars=209;
   V[209]=-1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[50])*(2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*Pow(V[80],2)*V[131]-3*Pow(V[50],2)*Pow(V[3],2)*Pow(V[80],2)*V[131]-2*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[80],2)*V[132]*V[50]+3*V[49]*Pow(V[3],2)*Pow(V[80],2)*V[132]*V[50]-2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*Pow(V[81],2)*V[131]+2*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[81],2)*V[132]*V[50]-6*Pow(V[39],2)*V[86]*V[80]*V[81]*V[132]*V[15]+6*Pow(V[39],2)*V[35]*V[86]*V[80]*V[81]*V[131]+6*Pow(V[39],2)*Pow(V[86],2)*V[131])/(6)/(Pow(V[70],2));
   if(!isfinite(V[209]) || FError) return 209;
   nCurrentVars=210;
   V[210]=1/(Pow(V[39],2))*V[8]*V[3]/(V[40])*(V[42]*V[129]*V[50]-V[42]*V[49]*V[130])/(2)/(Pow(V[54],2));
   if(!isfinite(V[210]) || FError) return 210;
   nCurrentVars=211;
   V[211]=1/(Pow(V[39],2))*V[8]*V[3]*V[40]*(V[129]*V[50]-V[49]*V[130])/(Pow(V[55],2));
   if(!isfinite(V[211]) || FError) return 211;
   nCurrentVars=212;
   V[212]=1/(Pow(V[39],2))*V[8]*V[3]/(V[40])*(V[42]*V[129]*V[50]-V[42]*V[49]*V[130])/(2)/(Pow(V[56],2));
   if(!isfinite(V[212]) || FError) return 212;
   nCurrentVars=213;
   V[213]=1/(Pow(V[39],2))*V[8]*V[3]*V[40]*(V[129]*V[50]-V[49]*V[130])/(Pow(V[57],2));
   if(!isfinite(V[213]) || FError) return 213;
   nCurrentVars=214;
   V[214]=1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[49])*(4*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[129]*Pow(V[153],2)*V[50]-3*V[49]*Pow(V[3],2)*V[129]*Pow(V[153],2)*V[50]-4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*V[130]*Pow(V[153],2)+3*Pow(V[49],2)*Pow(V[3],2)*V[130]*Pow(V[153],2)-4*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[129]*Pow(V[154],2)*V[50]+4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*V[130]*Pow(V[154],2)+6*Pow(V[39],2)*V[10]*V[129]*V[153]*V[154]*V[15]-6*Pow(V[39],2)*V[11]*V[10]*V[130]*V[153]*V[154]-6*Pow(V[39],2)*Pow(V[10],2)*V[130])/(6)/(Pow(V[60],2));
   if(!isfinite(V[214]) || FError) return 214;
   nCurrentVars=215;
   V[215]=1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[49])*(4*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[129]*Pow(V[155],2)*V[50]-3*V[49]*Pow(V[3],2)*V[129]*Pow(V[155],2)*V[50]-4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*V[130]*Pow(V[155],2)+3*Pow(V[49],2)*Pow(V[3],2)*V[130]*Pow(V[155],2)-4*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[129]*Pow(V[156],2)*V[50]+4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*V[130]*Pow(V[156],2)+6*Pow(V[39],2)*V[10]*V[129]*V[155]*V[156]*V[15]-6*Pow(V[39],2)*V[11]*V[10]*V[130]*V[155]*V[156]-6*Pow(V[39],2)*Pow(V[10],2)*V[130])/(6)/(Pow(V[61],2));
   if(!isfinite(V[215]) || FError) return 215;
   nCurrentVars=216;
   V[216]=1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[49])*(4*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[161],2)*V[129]*V[50]-3*V[49]*Pow(V[3],2)*Pow(V[161],2)*V[129]*V[50]-4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*Pow(V[161],2)*V[130]+3*Pow(V[49],2)*Pow(V[3],2)*Pow(V[161],2)*V[130]-4*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[162],2)*V[129]*V[50]+4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*Pow(V[162],2)*V[130]+6*Pow(V[39],2)*V[146]*V[161]*V[162]*V[129]*V[15]-6*Pow(V[39],2)*V[11]*V[146]*V[161]*V[162]*V[130]-6*Pow(V[39],2)*Pow(V[146],2)*V[130])/(6)/(Pow(V[64],2));
   if(!isfinite(V[216]) || FError) return 216;
   nCurrentVars=217;
   V[217]=1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[49])*(4*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[163],2)*V[129]*V[50]-3*V[49]*Pow(V[3],2)*Pow(V[163],2)*V[129]*V[50]-4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*Pow(V[163],2)*V[130]+3*Pow(V[49],2)*Pow(V[3],2)*Pow(V[163],2)*V[130]-4*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[164],2)*V[129]*V[50]+4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*Pow(V[164],2)*V[130]+6*Pow(V[39],2)*V[146]*V[163]*V[164]*V[129]*V[15]-6*Pow(V[39],2)*V[11]*V[146]*V[163]*V[164]*V[130]-6*Pow(V[39],2)*Pow(V[146],2)*V[130])/(6)/(Pow(V[65],2));
   if(!isfinite(V[217]) || FError) return 217;
   nCurrentVars=218;
   V[218]=-1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[50])*(2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*Pow(V[157],2)*V[129]-3*Pow(V[50],2)*Pow(V[3],2)*Pow(V[157],2)*V[129]-2*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[157],2)*V[130]*V[50]+3*V[49]*Pow(V[3],2)*Pow(V[157],2)*V[130]*V[50]-2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*Pow(V[158],2)*V[129]+2*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[158],2)*V[130]*V[50]-6*Pow(V[39],2)*V[10]*V[157]*V[158]*V[130]*V[15]+6*Pow(V[39],2)*V[12]*V[10]*V[157]*V[158]*V[129]+6*Pow(V[39],2)*Pow(V[10],2)*V[129])/(6)/(Pow(V[58],2));
   if(!isfinite(V[218]) || FError) return 218;
   nCurrentVars=219;
   V[219]=-1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[50])*(2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*Pow(V[159],2)*V[129]-3*Pow(V[50],2)*Pow(V[3],2)*Pow(V[159],2)*V[129]-2*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[159],2)*V[130]*V[50]+3*V[49]*Pow(V[3],2)*Pow(V[159],2)*V[130]*V[50]-2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*Pow(V[160],2)*V[129]+2*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[160],2)*V[130]*V[50]-6*Pow(V[39],2)*V[10]*V[159]*V[160]*V[130]*V[15]+6*Pow(V[39],2)*V[12]*V[10]*V[159]*V[160]*V[129]+6*Pow(V[39],2)*Pow(V[10],2)*V[129])/(6)/(Pow(V[59],2));
   if(!isfinite(V[219]) || FError) return 219;
   nCurrentVars=220;
   V[220]=-1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[50])*(2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*V[129]*Pow(V[165],2)-3*Pow(V[50],2)*Pow(V[3],2)*V[129]*Pow(V[165],2)-2*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[130]*Pow(V[165],2)*V[50]+3*V[49]*Pow(V[3],2)*V[130]*Pow(V[165],2)*V[50]-2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*V[129]*Pow(V[166],2)+2*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[130]*Pow(V[166],2)*V[50]-6*Pow(V[39],2)*V[10]*V[130]*V[165]*V[166]*V[15]+6*Pow(V[39],2)*V[12]*V[10]*V[129]*V[165]*V[166]+6*Pow(V[39],2)*Pow(V[10],2)*V[129])/(6)/(Pow(V[62],2));
   if(!isfinite(V[220]) || FError) return 220;
   nCurrentVars=221;
   V[221]=-1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[50])*(2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*V[129]*Pow(V[167],2)-3*Pow(V[50],2)*Pow(V[3],2)*V[129]*Pow(V[167],2)-2*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[130]*Pow(V[167],2)*V[50]+3*V[49]*Pow(V[3],2)*V[130]*Pow(V[167],2)*V[50]-2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*V[129]*Pow(V[168],2)+2*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[130]*Pow(V[168],2)*V[50]-6*Pow(V[39],2)*V[10]*V[130]*V[167]*V[168]*V[15]+6*Pow(V[39],2)*V[12]*V[10]*V[129]*V[167]*V[168]+6*Pow(V[39],2)*Pow(V[10],2)*V[129])/(6)/(Pow(V[63],2));
   if(!isfinite(V[221]) || FError) return 221;
   nCurrentVars=222;
   V[222]=1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[50])*(V[42]*Pow(V[50],2)*Pow(V[3],2)*V[129]*Pow(V[82],2)-V[42]*V[49]*Pow(V[3],2)*V[130]*Pow(V[82],2)*V[50]+2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*V[129]*Pow(V[83],2)-2*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[130]*Pow(V[83],2)*V[50]+2*Pow(V[39],2)*V[9]*V[130]*V[82]*V[83]*V[15]-2*Pow(V[39],2)*V[36]*V[9]*V[129]*V[82]*V[83]-2*Pow(V[39],2)*Pow(V[9],2)*V[129])/(2)/(Pow(V[67],2));
   if(!isfinite(V[222]) || FError) return 222;
   nCurrentVars=223;
   V[223]=1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[50])*(V[42]*Pow(V[50],2)*Pow(V[3],2)*V[129]*Pow(V[84],2)-V[42]*V[49]*Pow(V[3],2)*V[130]*Pow(V[84],2)*V[50]+2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*V[129]*Pow(V[85],2)-2*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[130]*Pow(V[85],2)*V[50]+2*Pow(V[39],2)*V[9]*V[130]*V[84]*V[85]*V[15]-2*Pow(V[39],2)*V[36]*V[9]*V[129]*V[84]*V[85]-2*Pow(V[39],2)*Pow(V[9],2)*V[129])/(2)/(Pow(V[68],2));
   if(!isfinite(V[223]) || FError) return 223;
   nCurrentVars=224;
   V[224]=1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[49])*(4*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[129]*Pow(V[74],2)*V[50]-3*V[49]*Pow(V[3],2)*V[129]*Pow(V[74],2)*V[50]-4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*V[130]*Pow(V[74],2)+3*Pow(V[49],2)*Pow(V[3],2)*V[130]*Pow(V[74],2)-4*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[129]*Pow(V[75],2)*V[50]+4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*V[130]*Pow(V[75],2)+6*Pow(V[39],2)*V[87]*V[129]*V[74]*V[75]*V[15]-6*Pow(V[39],2)*V[34]*V[87]*V[130]*V[74]*V[75]-6*Pow(V[39],2)*Pow(V[87],2)*V[130])/(6)/(Pow(V[71],2));
   if(!isfinite(V[224]) || FError) return 224;
   nCurrentVars=225;
   V[225]=1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[49])*(4*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[129]*Pow(V[76],2)*V[50]-3*V[49]*Pow(V[3],2)*V[129]*Pow(V[76],2)*V[50]-4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*V[130]*Pow(V[76],2)+3*Pow(V[49],2)*Pow(V[3],2)*V[130]*Pow(V[76],2)-4*Pow(V[40],2)*V[49]*Pow(V[3],2)*V[129]*Pow(V[77],2)*V[50]+4*Pow(V[40],2)*Pow(V[49],2)*Pow(V[3],2)*V[130]*Pow(V[77],2)+6*Pow(V[39],2)*V[87]*V[129]*V[76]*V[77]*V[15]-6*Pow(V[39],2)*V[34]*V[87]*V[130]*V[76]*V[77]-6*Pow(V[39],2)*Pow(V[87],2)*V[130])/(6)/(Pow(V[72],2));
   if(!isfinite(V[225]) || FError) return 225;
   nCurrentVars=226;
   V[226]=-1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[50])*(2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*Pow(V[78],2)*V[129]-3*Pow(V[50],2)*Pow(V[3],2)*Pow(V[78],2)*V[129]-2*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[78],2)*V[130]*V[50]+3*V[49]*Pow(V[3],2)*Pow(V[78],2)*V[130]*V[50]-2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*Pow(V[79],2)*V[129]+2*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[79],2)*V[130]*V[50]-6*Pow(V[39],2)*V[86]*V[78]*V[79]*V[130]*V[15]+6*Pow(V[39],2)*V[35]*V[86]*V[78]*V[79]*V[129]+6*Pow(V[39],2)*Pow(V[86],2)*V[129])/(6)/(Pow(V[69],2));
   if(!isfinite(V[226]) || FError) return 226;
   nCurrentVars=227;
   V[227]=-1/(Pow(V[39],2))*V[8]/(V[3])/(V[40])/(V[50])*(2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*Pow(V[80],2)*V[129]-3*Pow(V[50],2)*Pow(V[3],2)*Pow(V[80],2)*V[129]-2*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[80],2)*V[130]*V[50]+3*V[49]*Pow(V[3],2)*Pow(V[80],2)*V[130]*V[50]-2*Pow(V[40],2)*Pow(V[50],2)*Pow(V[3],2)*Pow(V[81],2)*V[129]+2*Pow(V[40],2)*V[49]*Pow(V[3],2)*Pow(V[81],2)*V[130]*V[50]-6*Pow(V[39],2)*V[86]*V[80]*V[81]*V[130]*V[15]+6*Pow(V[39],2)*V[35]*V[86]*V[80]*V[81]*V[129]+6*Pow(V[39],2)*Pow(V[86],2)*V[129])/(6)/(Pow(V[70],2));
   if(!isfinite(V[227]) || FError) return 227;
   nCurrentVars=228;
   V[228]=V[8]*V[3]/(V[40])*(V[131]*V[50]+V[49]*V[132])/(Pow(V[3],2));
   if(!isfinite(V[228]) || FError) return 228;
   nCurrentVars=229;
   V[229]=Pow(V[49],2)*V[131]*V[50]*V[126]+Pow(V[50],2)*V[49]*V[132]*V[139]+Pow(V[50],2)*V[131]*V[50]*V[127]+Pow(V[49],3)*V[132]*V[127]-Pow(V[50],2)*V[49]*V[132]*V[137];
   if(!isfinite(V[229]) || FError) return 229;
   nCurrentVars=230;
   V[230]=V[229]-Pow(V[49],2)*V[131]*V[50]*V[137]-Pow(V[50],2)*V[49]*V[132]*V[136]-Pow(V[49],2)*V[131]*V[50]*V[136]-3*Pow(V[49],3)*V[131]*V[128]+2*V[49]*V[131]*V[128];
   if(!isfinite(V[230]) || FError) return 230;
   nCurrentVars=231;
   V[231]=V[230]-Pow(V[49],2)*V[132]*V[50]*V[128]+3*Pow(V[49],2)*V[132]*V[50]*V[138]-V[132]*V[50]*V[138]-Pow(V[50],2)*V[49]*V[131]*V[138];
   if(!isfinite(V[231]) || FError) return 231;
   nCurrentVars=232;
   V[232]=-2/(V[8])*V[3]*V[40]*V[231]/(Pow(V[121],2));
   if(!isfinite(V[232]) || FError) return 232;
   nCurrentVars=233;
   V[233]=V[8]*V[3]/(V[40])*(V[129]*V[50]+V[49]*V[130])/(Pow(V[3],2));
   if(!isfinite(V[233]) || FError) return 233;
   nCurrentVars=234;
   V[234]=Pow(V[49],2)*V[129]*V[50]*V[126]+Pow(V[50],2)*V[49]*V[130]*V[139]+Pow(V[50],2)*V[129]*V[50]*V[127]+Pow(V[49],3)*V[130]*V[127]-Pow(V[50],2)*V[49]*V[130]*V[137];
   if(!isfinite(V[234]) || FError) return 234;
   nCurrentVars=235;
   V[235]=V[234]-Pow(V[49],2)*V[129]*V[50]*V[137]-Pow(V[50],2)*V[49]*V[130]*V[136]-Pow(V[49],2)*V[129]*V[50]*V[136]-3*Pow(V[49],3)*V[129]*V[128]+2*V[49]*V[129]*V[128];
   if(!isfinite(V[235]) || FError) return 235;
   nCurrentVars=236;
   V[236]=V[235]-Pow(V[49],2)*V[130]*V[50]*V[128]+3*Pow(V[49],2)*V[130]*V[50]*V[138]-V[130]*V[50]*V[138]-Pow(V[50],2)*V[49]*V[129]*V[138];
   if(!isfinite(V[236]) || FError) return 236;
   nCurrentVars=237;
   V[237]=-2/(V[8])*V[3]*V[40]*V[236]/(Pow(V[121],2));
   if(!isfinite(V[237]) || FError) return 237;
   nCurrentVars=238;
   V[238]=-V[8]/(V[3])/(V[40])*M_SQRT2/(V[50])/(V[49])*(V[49]*V[180]*V[131]*V[112]*V[115]-V[183]*V[132]*V[111]*V[116]*V[50])/(2)/(V[93]);
   if(!isfinite(V[238]) || FError) return 238;
   nCurrentVars=239;
   V[239]=-V[8]/(V[3])/(V[40])*M_SQRT2/(V[50])/(V[49])*(V[49]*V[180]*V[131]*V[114]*V[117]-V[183]*V[132]*V[113]*V[118]*V[50])/(2)/(V[94]);
   if(!isfinite(V[239]) || FError) return 239;
   nCurrentVars=240;
   V[240]=-V[8]/(V[3])/(V[40])*M_SQRT2/(V[50])/(V[49])*(V[49]*V[180]*V[129]*V[112]*V[115]-V[183]*V[130]*V[111]*V[116]*V[50])/(2)/(V[93]);
   if(!isfinite(V[240]) || FError) return 240;
   nCurrentVars=241;
   V[241]=-V[8]/(V[3])/(V[40])*M_SQRT2/(V[50])/(V[49])*(V[49]*V[180]*V[129]*V[114]*V[117]-V[183]*V[130]*V[113]*V[118]*V[50])/(2)/(V[94]);
   if(!isfinite(V[241]) || FError) return 241;
   nCurrentVars=242;
   V[242]=V[8]/(V[3])*V[146]/(V[40])/(V[48])/(2)/(V[146])/(2);
   if(!isfinite(V[242]) || FError) return 242;
   nCurrentVars=243;
   V[243]=-V[8]/(V[3])*V[144]/(V[40])/(V[50])/(V[49])*(V[142]-Pow(V[49],2))/(2)/(V[144])/(2);
   if(!isfinite(V[243]) || FError) return 243;
   nCurrentVars=244;
   V[244]=V[8]/(V[3])*V[145]/(V[40])/(V[48])/(2)/(V[145])/(2);
   if(!isfinite(V[244]) || FError) return 244;
   nCurrentVars=245;
   V[245]=-V[8]/(V[3])*V[9]/(V[40])/(V[50])/(V[49])*(V[143]-Pow(V[49],2))/(2)/(V[9])/(2);
   if(!isfinite(V[245]) || FError) return 245;
   nCurrentVars=246;
   V[246]=V[8]/(V[3])/(V[40])*M_SQRT2/(V[50])/(V[49])*(Pow(V[50],2)*V[183]*V[111]*V[116]-Pow(V[49],2)*V[180]*V[112]*V[115])/(2)/(V[93])/(2);
   if(!isfinite(V[246]) || FError) return 246;
   nCurrentVars=247;
   V[247]=V[8]/(V[3])/(V[40])*M_SQRT2/(V[50])/(V[49])*(Pow(V[50],2)*V[183]*V[113]*V[118]-Pow(V[49],2)*V[180]*V[114]*V[117])/(2)/(V[94])/(2);
   if(!isfinite(V[247]) || FError) return 247;
   nCurrentVars=248;
   V[248]=V[185]/(1-V[142]-V[142]*V[123]/(V[124])/(V[48]));
   if(!isfinite(V[248]) || FError) return 248;
   nCurrentVars=249;
   V[249]=V[189]/(1-V[142]+V[142]*V[124]/(V[123])/(V[48]));
   if(!isfinite(V[249]) || FError) return 249;
   nCurrentVars=250;
   V[250]=V[243]/(1-V[142]-V[142]/(Pow(V[48],2)));
   if(!isfinite(V[250]) || FError) return 250;
   nCurrentVars=251;
   V[251]=V[187]/(1-V[143]-V[143]*V[123]/(V[124])/(V[48]));
   if(!isfinite(V[251]) || FError) return 251;
   nCurrentVars=252;
   V[252]=V[191]/(1-V[143]+V[143]*V[124]/(V[123])/(V[48]));
   if(!isfinite(V[252]) || FError) return 252;
   nCurrentVars=253;
   V[253]=V[245]/(1-V[143]-V[143]/(Pow(V[48],2)));
   if(!isfinite(V[253]) || FError) return 253;
   nCurrentVars=254;
   V[254]=1/(137.036);
   if(!isfinite(V[254]) || FError) return 254;
   nCurrentVars=255;
   V[255]=4/(double)((9));
   if(!isfinite(V[255]) || FError) return 255;
   nCurrentVars=256;
   V[256]=1/(double)((9));
   if(!isfinite(V[256]) || FError) return 256;
   nCurrentVars=257;
   V[257]=-V[8]/(V[3])/(V[40])/(2);
   if(!isfinite(V[257]) || FError) return 257;
   nCurrentVars=258;
   V[258]=V[8]/(V[40])/(V[3])/(2);
   if(!isfinite(V[258]) || FError) return 258;
   nCurrentVars=259;
   V[259]=alphaQCD(V[119])/(M_PI);
   if(!isfinite(V[259]) || FError) return 259;
   nCurrentVars=260;
   V[260]=Sqrt(1+V[259]*(149/(double)((12))+V[259]*(68.6482-V[259]*212.447)));
   if(!isfinite(V[260]) || FError) return 260;
   nCurrentVars=261;
   V[261]=alphaQCD(V[120])/(M_PI);
   if(!isfinite(V[261]) || FError) return 261;
   nCurrentVars=262;
   V[262]=Sqrt(1+V[261]*(149/(double)((12))+V[261]*(68.6482-V[261]*212.447)));
   if(!isfinite(V[262]) || FError) return 262;
   nCurrentVars=263;
   V[263]=alphaQCD(V[14])/(M_PI);
   if(!isfinite(V[263]) || FError) return 263;
   nCurrentVars=264;
   V[264]=Sqrt(1+V[263]*(149/(double)((12))+V[263]*(68.6482-V[263]*212.447)));
   if(!isfinite(V[264]) || FError) return 264;
   nCurrentVars=265;
   V[265]=-Cabs(hGGeven(V[119],V[259],15,0,3,V[60],V[196],0,3,V[61],V[197],0,3,V[58],V[200],0,3,V[59],V[201],0,3,V[62],V[202],0,3,V[63],V[203],0,3,V[64],V[198],0,3,V[65],V[199],0,3,V[69],V[208],0,3,V[70],V[209],0,3,V[71],V[206],0,3,V[72],V[207],1,3,V[46],V[248],1,3,V[47],V[184],1,3,V[4],V[186]));
   if(!isfinite(V[265]) || FError) return 265;
   nCurrentVars=266;
   V[266]=-Cabs(hGGeven(V[120],V[261],15,0,3,V[60],V[214],0,3,V[61],V[215],0,3,V[58],V[218],0,3,V[59],V[219],0,3,V[62],V[220],0,3,V[63],V[221],0,3,V[64],V[216],0,3,V[65],V[217],0,3,V[69],V[226],0,3,V[70],V[227],0,3,V[71],V[224],0,3,V[72],V[225],1,3,V[46],V[249],1,3,V[47],V[188],1,3,V[4],V[190]));
   if(!isfinite(V[266]) || FError) return 266;
   nCurrentVars=267;
   V[267]=-Cabs(V[255]*hAAeven(V[119],V[259],8,1,3,V[4],V[186],1,3,V[47],V[184],0,3,V[60],V[196],0,3,V[61],V[197],0,3,V[64],V[198],0,3,V[65],V[199],0,3,V[71],V[206],0,3,V[72],V[207])+V[256]*hAAeven(V[119],V[259],7,1,3,V[46],V[248],0,3,V[58],V[200],0,3,V[59],V[201],0,3,V[62],V[202],0,3,V[63],V[203],0,3,V[69],V[208],0,3,V[70],V[209])+hAAeven(V[119],V[259],11,2,1,V[3],V[228],1,1,V[93],V[238],1,1,V[94],V[239],0,1,V[54],V[192],0,1,V[55],V[193],0,1,V[56],V[194],0,1,V[57],V[195],0,1,V[67],V[204],0,1,V[68],V[205],1,1,V[9],V[251],0,1,V[121],V[232]));
   if(!isfinite(V[267]) || FError) return 267;
   nCurrentVars=268;
   V[268]=-Cabs(V[255]*hAAeven(V[120],V[261],8,1,3,V[4],V[190],1,3,V[47],V[188],0,3,V[60],V[214],0,3,V[61],V[215],0,3,V[64],V[216],0,3,V[65],V[217],0,3,V[71],V[224],0,3,V[72],V[225])+V[256]*hAAeven(V[120],V[261],7,1,3,V[46],V[249],0,3,V[58],V[218],0,3,V[59],V[219],0,3,V[62],V[220],0,3,V[63],V[221],0,3,V[69],V[226],0,3,V[70],V[227])+hAAeven(V[120],V[261],11,2,1,V[3],V[233],1,1,V[93],V[240],1,1,V[94],V[241],0,1,V[54],V[210],0,1,V[55],V[211],0,1,V[56],V[212],0,1,V[57],V[213],0,1,V[67],V[222],0,1,V[68],V[223],1,1,V[9],V[252],0,1,V[121],V[237]));
   if(!isfinite(V[268]) || FError) return 268;
   nCurrentVars=269;
   V[269]=-Cabs(hGGodd(V[14],V[263],3,1,3,V[47],V[242],1,3,V[46],V[250],1,3,V[4],V[244]));
   if(!isfinite(V[269]) || FError) return 269;
   nCurrentVars=270;
   V[270]=-Cabs(V[256]*hAAodd(V[14],V[263],1,1,3,V[46],V[250])+V[255]*hAAodd(V[14],V[263],2,1,3,V[4],V[244],1,3,V[47],V[242])+hAAodd(V[14],V[263],3,1,1,V[9],V[253],1,1,V[93],V[246],1,1,V[94],V[247]));
   if(!isfinite(V[270]) || FError) return 270;
   if(VV==NULL) 
   {  VV=malloc(sizeof(REAL)*nModelVars);
      for(i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) iQ=i;
   }
   for(i=0;i<nModelVars;i++) VV[i]=V[i];
   cErr=0;
   nCurrentVars++;
   return 0;
}
