* Pi12t.F
* the stop1-stop2 self-energy
* this file is part of FeynHiggs
* last modified 1 Jan 17 th

#include "externals.h"
#include "types.h"


	RealType function Pi12t(p2, Q2)
	implicit none
	RealType p2, Q2

#include "TLps.h"

	Pi12t = 1/2.D0*((Atr*c2tt*
     &         (2*ctb2*ht2*pMT + hb*ht*MUEr*s2tb) + 
     &        Atr**2*ctb2*ht2*s2tt - 
     &        s2tt*(ctb2*ht2*pMT2 + hb*ht*pMT*MUEr*s2tb + 
     &           hb2*MUEr2*stb2))*myB0(p2,MSb12,0D0,Q2) + 
     &     (ctb2*ht2*MUEr2*s2tt + 
     &        hb*ht*MUEr*s2tb*(Abr*c2tt + pMT*s2tt) - 
     &        hb2*(-2*Abr*c2tt*pMT + Abr**2*s2tt - pMT2*s2tt)*stb2)*
     &      myB0(p2,MSb12,MA02,Q2) + 
     &     (-(ctb2*hb2*MUEr2*s2tt) + 
     &        hb*ht*MUEr*s2tb*(-(Atr*c2tt) + pMT*s2tt) + 
     &        ht2*(2*Atr*c2tt*pMT + Atr**2*s2tt - pMT2*s2tt)*stb2)*
     &      myB0(p2,MSb22,0D0,Q2) + 
     &     (-(hb*ht*MUEr*s2tb*(Abr*c2tt + pMT*s2tt)) + 
     &        ctb2*hb2*(2*Abr*c2tt*pMT - Abr**2*s2tt + pMT2*s2tt) + 
     &        ht2*MUEr2*s2tt*stb2)*myB0(p2,MSb22,MA02,Q2) + 
     &     Atr*c2tt*ht2*((2*pMT + Atr*s2tt)*myB0(p2,MSt12,0D0,Q2) + 
     &        (2*pMT - Atr*s2tt)*myB0(p2,MSt22,0D0,Q2)) + 
     &     s2tt*(-(hb2*myAA(MA02,Q2)) + 
     &        (ctb2*ht2 - hb2*stb2)*myAA(MSb12,Q2) + 
     &        (-(ctb2*hb2) + ht2*stb2)*myAA(MSb22,Q2) + 
     &        c2tt*ht2*(1 + Nc)*myAA(MSt12,Q2) - 
     &        c2tt*ht2*(1 + Nc)*myAA(MSt22,Q2) + 
     &        (hb2 - ht2)*(myAA(0D0,Q2) + myAA(MUEr2,Q2)) + 
     &        (hb2 - ht2)*(MUEr2 - p2)*myB0(p2,0D0,MUEr2,Q2) + 
     &        c2tt*ht2*MUEr2*myB0(p2,MSt12,MA02,Q2) - 
     &        c2tt*ht2*MUEr2*myB0(p2,MSt22,MA02,Q2)))
	end

