#if 0
* hgaga-v1_All4_Ind1_Ind2.F
* form factors for {S[10, {hno}]} -> {V[1], V[1]}
* generated by FormCalc 9.8 (22 Apr 2019) on 20-Dec-2019 10:10
* this file is part of FeynHiggs
#endif

#include "hgaga-vars.h"
#include "inline.h"

	subroutine hgaga_v1_All4_Ind1_Ind2
	implicit none

#include "hgaga-vars.h"
#include "inline.h"

	Cloop(1) = Cloop(1) + 
     &    1/Pi*(Alfa0*EL0*(1/4.D0*
     &           (Sub22(All4,Ind1,Ind2)*Sub3(All4)) + 
     &          1/36.D0*(Sub25(All4)*Sub44(All4,Ind1,Ind2)) - 
     &          1/9.D0*(Sub47(All4)*Sub66(All4,Ind1,Ind2))))/(MW*SW)

#ifdef DETAILED_DEBUG
	DCOUP "hgaga-v1_All4_Ind1_Ind2: Cloop(1) =", Cloop(1) ENDL
#endif

	end

