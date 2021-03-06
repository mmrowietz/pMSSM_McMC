* dZ_HHHH_exp.F
* generated 6-Aug-2018 12:46
* this file is part of FeynHiggs
* please do not edit directly

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine dZ_HHHH_exp(se,k2)
	implicit none
	ComplexType se
	RealType k2

#include "FH.h"
#include "looptools.h"

	RealType MSusy, xt

#ifdef DETAILED_DEBUG
	DHIGGS "dZ_HHHH_exp k2 =", k2 ENDL
#endif

	MSusy = sqrt(Re(MSS(tQ(3),3)*MSS(tU(3),3)))

	xt = Re(Xtgl)/MSusy

	se = 0

	if( mssmpart .le. 1 ) return

	if( mssmpart .le. 2 ) return

	if( mssmpart .le. 3 ) return

	end
