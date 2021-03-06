#ifndef RECORDINDICES_H
#define RECORDINDICES_H

#define iVar 1
#define iLower 2
#define iUpper 3
#define iStep 4
#define iAdmin 1
#define FHRecordR 2
#define iinvAlfa0 2
#define iinvAlfaMZ 3
#define iAlfasMZ 4
#define iGF 5
#define iME 6
#define iMU 7
#define iMD 8
#define iMM 9
#define iMC 10
#define iMS 11
#define iML 12
#define iMT 13
#define iMB 14
#define iMW 15
#define iMZ 16
#define iGammaW 17
#define iGammaZ 18
#define iCKMlambda 19
#define iCKMA 20
#define iCKMrhobar 21
#define iCKMetabar 22
#define iTB 23
#define iMA0 24
#define iMHp 25
#define iMSusy 26
#define iM1SL 27
#define iM1SE 28
#define iM1SQ 29
#define iM1SU 30
#define iM1SD 31
#define iM2SL 32
#define iM2SE 33
#define iM2SQ 34
#define iM2SU 35
#define iM2SD 36
#define iM3SL 37
#define iM3SE 38
#define iM3SQ 39
#define iM3SU 40
#define iM3SD 41
#define iQtau 42
#define iQt 43
#define iQb 44
#define iscalefactor 45
#define iprodSqrts 46
#define FHRecordC 47
#define iAe 47
#define iAu 51
#define iAd 55
#define iAmu 59
#define iAc 63
#define iAs 67
#define iAtau 71
#define iAt 75
#define iAb 79
#define iXtau 83
#define iXt 87
#define iXb 91
#define iMUE 95
#define iM1 99
#define iM2 103
#define iM3 107
#define ideltaLLL12 111
#define ideltaLLL23 115
#define ideltaLLL13 119
#define ideltaELR12 123
#define ideltaELR23 127
#define ideltaELR13 131
#define ideltaERL12 135
#define ideltaERL23 139
#define ideltaERL13 143
#define ideltaERR12 147
#define ideltaERR23 151
#define ideltaERR13 155
#define ideltaQLL12 159
#define ideltaQLL23 163
#define ideltaQLL13 167
#define ideltaULR12 171
#define ideltaULR23 175
#define ideltaULR13 179
#define ideltaURL12 183
#define ideltaURL23 187
#define ideltaURL13 191
#define ideltaURR12 195
#define ideltaURR23 199
#define ideltaURR13 203
#define ideltaDLR12 207
#define ideltaDLR23 211
#define ideltaDLR13 215
#define ideltaDRL12 219
#define ideltaDRL23 223
#define ideltaDRL13 227
#define ideltaDRR12 231
#define ideltaDRR23 235
#define ideltaDRR13 239
#define FHRecordE 243
#define FHRecordN 242

#endif
* FHRecord.h.in
* the data structures for a FH record
* this file is part of FeynHiggs
* last modified 26 Jan 16 th


#ifndef FHRangeR
#define FHRangeR FHRecordR, FHRecordC - 1
#define FHRangeC FHRecordC, FHRecordN, 4
#define FHRangeA FHRecordR, FHRecordN

#define iRe(i) i
#define iIm(i) i+1
#define iAbs(i) i+2
#define iArg(i) i+3

#define iMSS(n,g) iM1SL+(n-1)*(iM1SE-iM1SL)+(g-1)*(iM2SL-iM1SL)
#define iMf(t,g) iMU+(t-3)*(iMD-iMU)+(g-1)*(iMC-iMU)
#define iAf(t,g) iAu+(t-3)*(iAd-iAu)+(g-1)*(iAc-iAu)
#define iQSf(t) iQt+(t-3)*(iQb-iQt)

#define FHNameR(i) FHName(i)(1:len_trim(FHName(i)))
#define FHNameC(i) FHName(i)(4:index(FHName(i),")")-1)

#define RecordDecl(rec) RealType rec(FHRecordN,4)
#endif

	RealType unset, default, bytable
	parameter (unset = -999)
	parameter (default = -888)
	parameter (bytable = 777)

	character*16 FHName(FHRecordR:FHRecordN)
	common /fhrecnames/ FHName

	integer maxcols, maxrows
	parameter (maxcols = FHRecordN, maxrows = 10000)

	RealType tabledata(maxcols,maxrows)
	integer tableflag(0:maxcols), tablerows
	common /fhtable/ tabledata, tableflag, tablerows

