#if 0
* hgaZSM-v1.F
* form factors for {S[1]} -> {V[1], V[2]}
* generated by FormCalc 9.8 (22 Apr 2019) on 20-Dec-2019 10:10
* this file is part of FeynHiggs
#endif

#include "hgaZSM-vars.h"
#include "inline.h"

	subroutine hgaZSM_v1
	implicit none

#include "hgaZSM-vars.h"
#include "inline.h"

	Cloop(1) = Cloop(1) + 
     &    EL0*(1/4.D0*(MW*Sub13) + 1/(8.D0*Pi)*(Alfa0*Sub23))

#ifdef DETAILED_DEBUG
	DCOUP "hgaZSM-v1: Cloop(1) =", Cloop(1) ENDL
#endif

	end

