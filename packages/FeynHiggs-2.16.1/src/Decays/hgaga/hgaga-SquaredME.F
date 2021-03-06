*#define CHECK

#if 0
* hgaga-SquaredME.F
* assembly of squared matrix element
* generated by FormCalc 9.8 (22 Apr 2019) on 20-Dec-2019 10:10
* this file is part of FeynHiggs
#endif

#include "hgaga-vars.h"
#include "inline.h"

************************************************************************

	subroutine hgaga_SquaredMEHel(res)
	implicit none
	RealType res(*)

#include "hgaga-vars.h"

	ComplexType amp
	RealType ampsq

* hgaga-BEGIN ABBR_HEL
	call hgaga_abbr1h
* hgaga-END ABBR_HEL

	ampsq = 0
	res(1) = ampsq

	ampsq = 0
* hgaga-BEGIN FF_LOOP
	Cloop(1) = 0

	call hgaga_v1

	LOOP(Cha4, 1,2,1)
	call hgaga_v1_Cha4
	ENDLOOP(Cha4)

	LOOP(Gen4, 1,3,1)
	call hgaga_v1_Gen4
	ENDLOOP(Gen4)

	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)
	LOOP(All4, 1,6,1)
	call hgaga_v1_All4_Ind1_Ind2
	ENDLOOP(All4)
	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
* hgaga-END FF_LOOP

* hgaga-BEGIN M2_LOOP
	amp = 0
	amp = amp + Cloop(1)
	ampsq = ampsq + Re(amp*Conjugate(Cloop(1)))
* hgaga-END M2_LOOP
	res(2) = ampsq

	end

************************************************************************

	subroutine hgaga_SquaredME(result)
	implicit none
	RealType result(*)

#include "hgaga-vars.h"

* hgaga-BEGIN VAR_DECL
	integer i, h, hmax
	parameter (hmax = 2)
	RealType res(2,hmax)
	RealType rtree, rloop
	external hgaga_SquaredMEHel
* hgaga-END VAR_DECL

#include "inline.h"

* hgaga-BEGIN INVARIANTS
* hgaga-END INVARIANTS

* hgaga-BEGIN ABBR_S
	call hgaga_abbr1s
* hgaga-END ABBR_S

* hgaga-BEGIN ABBR_ANGLE
* hgaga-END ABBR_ANGLE

* hgaga-BEGIN HEL_LOOP

	h = 0
	do i = 0, hmax - 1

#define hval(s,m,a,b) a*mod(i/s,m)-b

	  Hel0(1) = 0
	  Hel0(2) = hval(1,2,2,1)
	  Hel0(3) = Hel0(2)
	  h = h + 1
	  call hgaga_SquaredMEHel(res(1,h))
	enddo

* hgaga-END HEL_LOOP

* hgaga-BEGIN RESULT
	rtree = 0
	rloop = 0
	do i = 1, h
	  rtree = rtree + res(1,i)
	  rloop = rloop + res(2,i)
	enddo
	result(1) = rtree
	result(2) = rloop

* hgaga-END RESULT

#ifdef CHECK
	print *, 'tree =', rtree
	print *, 'loop =', rloop
	stop
#endif

	end

