#if 0
* hee-v1_Neu4_Neu5.F
* form factors for {S[10, {hno}]} -> {F[2, {gno1}], -F[2, {gno2}]}
* generated by FormCalc 9.9 (29 Jan 2020) on 21-Apr-2020 10:38
* this file is part of FeynHiggs
#endif

#include "hee-vars.h"
#include "inline.h"

	subroutine hee_v1_Neu4_Neu5
	implicit none

#include "hee-vars.h"
#include "inline.h"

	Cloop(1) = Cloop(1) + 
     &    cI/(16.D0*Pi)*(AlfaMZ*ELMZ*Delta(gno1,gno2)*
     &        (A0i8(Neu5)*Mf(2,gno1)*Sub233(Neu4,Neu5) + 
     &          MW*(B1i8(Neu4,Neu5)*Sub238(Neu4,Neu5) + 
     &             B0i9(Neu4,Neu5)*Sub240(Neu4,Neu5))))/
     &      (CW2*MW*(m12 - MZ2)*SW*SW2)

#ifdef DETAILED_DEBUG
	DCOUP "hee-v1_Neu4_Neu5: Cloop(1) =", Cloop(1) ENDL
#endif

	end

