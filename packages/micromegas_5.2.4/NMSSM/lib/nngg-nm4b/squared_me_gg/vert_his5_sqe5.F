#if 0
* vert_his5_sqe5.F
* form factors for {F[6, {1}], F[6, {1}]} -> {V[1], V[1]}
* generated by FormCalc 8.0 on 6-Apr-2015 1:04
#endif

#include "vars.h"


	subroutine vert_his5_sqe5
	implicit none

#include "vars.h"


#include "inline.h"

        Cloop(1) = Cloop(1) + 
     &    1/Pi**2*(AAABR(373)*
     &        (Pair3*(-(1/288.D0*pave15(bb0,sqe5)) + 
     &             1/72.D0*pave16(cc00,sqe5)) - 
     &          1/72.D0*(Pair1*Pair2*pave16(cc12,sqe5)))*
     &        MTR026(his5,sqe5)*Sub7(his5))/(S - hisMass(his5)**2)

#ifdef DEBUG
	DEB("vert_his5_sqe5: Cloop(1) =", Cloop(1))
#endif

	end
