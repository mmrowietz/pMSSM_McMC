#if 0
* hdd-v1_Cha4_Cha5.F
* form factors for {S[10, {hno}]} -> {F[4, {gno1, Col2}], -F[4, {gno2, Col3}]}
* generated by FormCalc 9.9 (29 Jan 2020) on 21-Apr-2020 10:42
* this file is part of FeynHiggs
#endif

#include "hdd-vars.h"
#include "inline.h"

	subroutine hdd_v1_Cha4_Cha5
	implicit none

#include "hdd-vars.h"
#include "inline.h"

	Cloop(1) = Cloop(1) - 
     &    cI/(24.D0*Pi)*(AlfaMZ*ELMZ*Delta(gno1,gno2)*
     &        (3*A0i4(Cha5)*Mf(bH,gno1)*Sub147(Cha4,Cha5) + 
     &          MW*(B1i4(Cha4,Cha5)*Sub152(Cha4,Cha5) + 
     &             B0i4(Cha4,Cha5)*Sub154(Cha4,Cha5))))/
     &      (MW*(m12 - MZ2)*SW*SW2)

#ifdef DETAILED_DEBUG
	DCOUP "hdd-v1_Cha4_Cha5: Cloop(1) =", Cloop(1) ENDL
#endif

	end

