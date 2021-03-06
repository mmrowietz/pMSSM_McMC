#if 0
* hgg-abbr1h.F
* abbreviations for {S[10, {hno}]} -> {V[5, {Glu2}], V[5, {Glu3}]}
* generated by FormCalc 9.8 (22 Apr 2019) on 20-Dec-2019 10:11
* this file is part of FeynHiggs
#endif

#include "hgg-vars.h"
#include "inline.h"

	subroutine hgg_abbr1h
	implicit none

#include "hgg-vars.h"
#include "inline.h"

	Eps1 = Eps(ec(2),ec(3),k(1),k(2))

	Pair1 = Pair(ec(2),ec(3))

	Pair2 = Pair(ec(2),k(1))

	Pair3 = Pair(ec(3),k(1))

	Abb1 = -(m12*Pair1) + 4*Pair2*Pair3

	Abb2 = -(m12*Pair1) + 2*(Eps1 + Pair2*Pair3)

	Abb3 = -2*Eps1 - m12*Pair1 + 2*Pair2*Pair3

	LOOP(All4, 1,6,1)

	Sub3(All4) = Pair1*Sub1(All4) + 4*Pair2*Pair3*Sub2(All4)

	Sub25(All4) = Pair1*Sub23(All4) + 4*Pair2*Pair3*Sub24(All4)

	ENDLOOP(All4)

	LOOP(Gen4, 1,3,1)

	Sub47(Gen4) = Abb1*C0i3(Epsi(cc2),Gen4) - 
     &    Pair1*Sub45(Gen4) + 4*Pair2*Pair3*Sub46(Gen4)

	Sub51(Gen4) = cI*
     &     (-(Abb1*Im(C0i4(Epsi(cc2),Gen4))) - 
     &       4*Pair2*Pair3*Sub50(Gen4)) + Pair1*Sub49(Gen4)

	Sub59(Gen4) = cI*
     &     (-(Abb1*Re(C0i4(Epsi(cc2),Gen4))) - 
     &       4*Pair2*Pair3*Sub58(Gen4)) + Pair1*Sub57(Gen4)

	Sub72(Gen4) = C0i3(Epsi(cc0),Gen4)*Mf2(tH,Gen4)*
     &     (Abb2*Sub62(Gen4) + Abb3*Sub63(Gen4)) + 
     &    Mf(4,Gen4)*Mf(bBR,Gen4)*
     &     (cI*(Im(C0i4(Epsi(cc0),Gen4))*
     &          (Abb2*Sub65(Gen4) + Abb3*Sub67(Gen4))) + 
     &       Re(C0i4(Epsi(cc0),Gen4))*
     &        (Abb2*Sub69(Gen4) + Abb3*Sub71(Gen4)))

	ENDLOOP(Gen4)

	end

