#if 0
* hgaZSM-abbr1h.F
* abbreviations for {S[1]} -> {V[1], V[2]}
* generated by FormCalc 9.8 (22 Apr 2019) on 20-Dec-2019 10:10
* this file is part of FeynHiggs
#endif

#include "hgaZSM-vars.h"
#include "inline.h"

	subroutine hgaZSM_abbr1h
	implicit none

#include "hgaZSM-vars.h"
#include "inline.h"

	Eps1 = Eps(ec(2),ec(3),k(1),k(2))

	Pair1 = Pair(ec(2),ec(3))

	Pair2 = Pair(ec(2),k(1))

	Pair3 = Pair(ec(3),k(1))

	Abb1 = -2*Eps1 - (m12 - m32)*Pair1 + 4*Pair2*Pair3

	Abb2 = 2*Eps1 + (-m12 + m32)*Pair1 + 4*Pair2*Pair3

	Abb3 = -2*Eps1 - (m12 - m32)*Pair1 + 2*Pair2*Pair3

	Abb4 = (-m12 + m32)*Pair1 + 2*(Eps1 + Pair2*Pair3)

	Abb5 = -((m12 - m32)*Pair1)

	Abb6 = (m12 + m32)*Pair1

	Abb7 = (3*m12 + m32)*Pair1 + 2*Pair2*Pair3

	Abb8 = (3*m12 - m32)*Pair1 + 2*Pair2*Pair3

	Abb9 = -((m12 - m32)*Pair1) - 2*Pair2*Pair3

	Abb10 = (m12 + m32)*Pair1 + 2*Pair2*Pair3

	Sub12 = (Abb9 + Abb5*CW2 + 2*Pair2*Pair3)/CW

	Sub13 = -(1/Pi*
     &      (-(Alfa0*Sub12*C0i1(Epsi(cc1))) + 
     &         Pi*(Pair1*Sub11*SW2))/SW2)

	Sub14 = ((2*Abb10 - 5*Abb5 + Abb7)*CW2 + 2*Abb8*SW2)/(CW*SW2)

	Sub16 = (2*Pair2*Pair3*Sub15 + CW*MW2*Sub14*SW2)/(CW*MW*SW2)

	Sub17 = Abb5*(7*CW2 - 4*SW2) + 
     &    8*Pair2*Pair3*(3*CW2 - SW2) + 
     &    2*Pair1*(3*MW2 + SW2*MHiggs2(hno))

	Sub22 = -(((4*Pair2*Pair3*Sub18*Sub19 + Pair1*Sub21)/
     &         MW - MW*(Abb6*CW2 - Sub17)*C0i1(Epsi(cc0)))/CW)

	Sub23 = (Sub22 - Sub16*SW2*C0i1(Epsi(cc2)))/SW2

	LOOP(Gen4, 1,3,1)

	Sub3(Gen4) = -((Abb3 + Abb4)*C0i2(Epsi(cc0),Gen4)) - 
     &    2*(Abb1 + Abb2)*C0i2(Epsi(cc2),Gen4) + 
     &    4*Pair1*Sub1(Gen4) - 16*Pair2*Pair3*Sub2(Gen4)

	Sub6(Gen4) = -((Abb3 + Abb4)*C0i3(Epsi(cc0),Gen4)) - 
     &    2*(Abb1 + Abb2)*C0i3(Epsi(cc2),Gen4) + 
     &    4*Pair1*Sub4(Gen4) - 16*Pair2*Pair3*Sub5(Gen4)

	Sub9(Gen4) = -((Abb3 + Abb4)*C0i4(Epsi(cc0),Gen4)) - 
     &    2*(Abb1 + Abb2)*C0i4(Epsi(cc2),Gen4) + 
     &    4*Pair1*Sub7(Gen4) - 16*Pair2*Pair3*Sub8(Gen4)

	ENDLOOP(Gen4)

	end

