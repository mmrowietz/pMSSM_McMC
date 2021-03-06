#if 0
* huu-v0.F
* form factors for {S[10, {hno}]} -> {F[3, {gno1, Col2}], -F[3, {gno2, Col3}]}
* generated by FormCalc 9.9 (29 Jan 2020) on 21-Apr-2020 10:40
* this file is part of FeynHiggs
#endif

#include "huu-vars.h"
#include "inline.h"

	subroutine huu_v0
	implicit none

#include "huu-vars.h"
#include "inline.h"

	Ctree(1) = Ctree(1) + 
     &    1/2.D0*(ELMZ*Delta(gno1,gno2)*Mf(tH,gno1)*
     &        ((F1 + F2)*Sub1*TB + SB*Sub2*XHiggs(hno,3,uzext1)))/
     &      (MW*SB*SW*TB)

#ifdef DETAILED_DEBUG
	DCOUP "huu-v0: Ctree(1) =", Ctree(1) ENDL
#endif

	end

