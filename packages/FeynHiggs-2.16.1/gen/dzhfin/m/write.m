CodeExpr[{dAt, dCB1, dMHinsq1, dMHinsq2, dSB1, dTA01, dTA02, dTB1, dTB1div, 
  dTB2div, dTG01, dTh01, dTh02, dTHH1, dTHH2, dZ11H1gl, dZ22H1gl, 
  dMHiggs1gl[5, 6], dMHiggs2gl[1, 3], dMHiggsZ2gl[5, 5], dZHiggs1gl[6, 6], 
  dZHiggs2gl[2, 3], se["HHA0"], seshift["HHA0"]}, 
 {A0q1, A0q10, A0q11, A0q12, A0q13, A0q14, A0q15, A0q16, A0q17, A0q2, A0q3, 
  A0q4, A0q5, A0q6, A0q7, A0q8, A0q9, B0q1, B0q10, B0q11, B0q12, B0q13, 
  B0q14, B0q15, B0q16, B0q17, B0q18, B0q19, B0q2, B0q20, B0q21, B0q22, B0q23, 
  B0q24, B0q25, B0q26, B0q3, B0q4, B0q5, B0q6, B0q7, B0q8, B0q9, dup100, 
  dup101, dup102, dup103, dup104, dup105, dup106, dup107, dup108, dup109, 
  dup110, dup111, dup112, dup113, dup114, dup115, dup116, dup117, dup17, 
  dup18, dup19, dup20, dup21, dup22, dup23, dup24, dup25, dup26, dup28, 
  dup29, dup30, dup31, dup32, dup33, dup34, dup35, dup36, dup37, dup39, dup4, 
  dup40, dup41, dup42, dup43, dup44, dup45, dup46, dup47, dup48, dup49, dup5, 
  dup50, dup51, dup52, dup53, dup54, dup55, dup56, dup57, dup58, dup59, dup6, 
  dup60, dup61, dup62, dup63, dup64, dup65, dup66, dup67, dup68, dup69, dup7, 
  dup70, dup71, dup72, dup73, dup74, dup75, dup76, dup77, dup78, dup79, dup8, 
  dup80, dup81, dup82, dup83, dup84, dup85, dup86, dup87, dup88, dup89, 
  dup90, dup91, dup92, dup93, dup94, dup95, dup96, dup97, dup98, dup99, tmp1, 
  tmp2, U2c11, U2c12, U2c13, U2c14, U2c21, U2c22, U2c23, U2c24, U2c25, U2c26, 
  U2c27, U2c28, U2s11, U2s12, U2s13, U2s14, U2s21, U2s22, U2s23, U2s24, 
  U2s25, U2s26, U2s27, U2s28}, {B0q12 -> B0q[0, MStgl2[2], MStgl2[1], Q2], 
  B0q13 -> B0q[0, MSbgl2[2], MSbgl2[1], Q2], U2s25 -> U2s2[UCStgl, Ytgl], 
  U2s26 -> U2s2[UCSbgl, Ybgl], 
  dTA01 -> (3*EL1L*(MTgl*U2s25*MStgl2[3]*Re[B0q12] + MBgl*TB2*U2s26*MSbgl2[3]*
       Re[B0q13]))/(16*MW*Pi^2*SW*TB), DebugLine[1, dTA01], 
  U2s27 -> U2s2[UCStgl, Xtgl], U2s28 -> U2s2[UCSbgl, Xbgl], 
  dTG01 -> (3*EL1L*(MTgl*U2s27*MStgl2[3]*Re[B0q12] - 
      MBgl*U2s28*MSbgl2[3]*Re[B0q13]))/(16*MW*Pi^2*SW), DebugLine[1, dTG01], 
  A0q12 -> A0q[MBgl2, Q2], A0q13 -> A0q[MTgl2, Q2], 
  A0q14 -> A0q[MStgl2[2], Q2], A0q15 -> A0q[MStgl2[1], Q2], 
  A0q16 -> A0q[MSbgl2[2], Q2], A0q17 -> A0q[MSbgl2[1], Q2], 
  U2s11 -> U2s1[UCStgl, Xtgl], U2s12 -> U2s1[UCSbgl, Xbgl], 
  dTh01 -> (EL1L*((3*(MBgl2*Re[A0q12] + MTgl2*Re[A0q13]))/8 - 
      (3*(MTgl*((MTgl - U2s11)*Re[A0q14] + (MTgl + U2s11)*Re[A0q15]) + 
         MBgl*((MBgl - U2s12)*Re[A0q16] + (MBgl + U2s12)*Re[A0q17])))/16))/
    (MW*Pi^2*SW), DebugLine[1, dTh01], 
  B0q14 -> B0q[0, MSbgl2[1], MSbgl2[2], Q2], U2s13 -> U2s1[UCStgl, Ytgl], 
  U2s14 -> U2s1[UCSbgl, Ybgl], dup74 -> 2*Re[A0q12] - Re[A0q16] - Re[A0q17], 
  dTHH1 -> (3*EL1L*(CB2*(MTgl2*(-2*Re[A0q13] + Re[A0q14] + Re[A0q15]) - 
        MTgl*U2s13*MStgl2[3]*Re[B0q12]) + 
      SB2*(dup74*MBgl2 + MBgl*U2s14*MSbgl2[3]*Re[B0q14])))/
    (8*MW*Pi^2*S2B*SW), DebugLine[1, dTHH1], 
  dZ11H1gl -> (-3*Alfa1L*Divergence*((CB2*MTgl2)/SB2 + MBgl2*TB2))/
     (8*MW2*Pi*SW2) + dZ11H1fingl[zM], DebugLine[1, dZ11H1gl], 
  dZ22H1gl -> (-3*Alfa1L*Divergence*(MBgl2 + MTgl2))/(8*MW2*Pi*SW2) + 
    dZ22H1fingl[zM], DebugLine[1, dZ22H1gl], dMHiggs1gl[1, 1] -> 
   -(dTh01*HoldCode[EL1L/(2*MW*SW)]), DebugLine[1, dMHiggs1gl[1, 1]], 
  dMHiggs1gl[1, 3] -> -(dTA01*HoldCode[EL1L/(2*MW*SW)]), 
  DebugLine[1, dMHiggs1gl[1, 3]], dMHiggs1gl[2, 4] -> 
   -(dTA01*HoldCode[EL1L/(2*MW*SW)]), DebugLine[1, dMHiggs1gl[2, 4]], 
  dZHiggs1gl[1, 1] -> CB2*Re[dZ11H1gl] + SB2*Re[dZ22H1gl] + 
    S2B*Re[dZ12H1fingl[zM]], DebugLine[1, dZHiggs1gl[1, 1]], 
  dZHiggs1gl[1, 2] -> (S2B*(Re[dZ11H1gl] - Re[dZ22H1gl]))/2 - 
    C2B*Re[dZ12H1fingl[zM]], DebugLine[1, dZHiggs1gl[1, 2]], 
  dZHiggs1gl[1, 3] -> -Im[dZ12H1fingl[zM]], DebugLine[1, dZHiggs1gl[1, 3]], 
  dZHiggs1gl[2, 2] -> SB2*Re[dZ11H1gl] + CB2*Re[dZ22H1gl] - 
    S2B*Re[dZ12H1fingl[zM]], DebugLine[1, dZHiggs1gl[2, 2]], 
  dZHiggs1gl[2, 4] -> -Im[dZ12H1fingl[zM]], DebugLine[1, dZHiggs1gl[2, 4]], 
  dZHiggs1gl[3, 4] -> -dZHiggs1gl[1, 2], DebugLine[1, dZHiggs1gl[3, 4]], 
  dZHiggs1gl[5, 6] -> dZHiggs1gl[3, 4] + I*Im[dZ12H1fingl[zM]], 
  DebugLine[1, dZHiggs1gl[5, 6]], dZHiggs1gl[6, 5] -> 
   dZHiggs1gl[3, 4] - I*Im[dZ12H1fingl[zM]], DebugLine[1, dZHiggs1gl[6, 5]], 
  dZHiggs2gl[1, 2] -> (S2B*(-Re[dZ11H1gl^2] + Re[dZ22H1gl^2]) + 
     2*C2B*Re[dZ12H1fingl[zM]^2])/8, DebugLine[1, dZHiggs2gl[1, 2]], 
  dZHiggs2gl[1, 3] -> Im[dZ12H1fingl[zM]^2]/4, 
  DebugLine[1, dZHiggs2gl[1, 3]], dZHiggs2gl[2, 2] -> 
   (-(SB2*Re[dZ11H1gl^2]) - CB2*Re[dZ22H1gl^2] + S2B*Re[dZ12H1fingl[zM]^2])/
    4, DebugLine[1, dZHiggs2gl[2, 2]], 
  dTB1div -> (SB*(-Re[dZ11H1gl] + Re[dZ22H1gl]) + 
     CB*(1 - TB2)*Re[dZ12H1fingl[zM]])/(2*CB), DebugLine[1, dTB1div], 
  dTB2div -> (-Conjugate[dZ12H1fingl[zM]]^2 - dZ12H1fingl[zM]^2 - 
     2*dZ11H1gl*(Conjugate[dZ12H1fingl[zM]] + dZ12H1fingl[zM]) + 
     TB2*(Conjugate[dZ12H1fingl[zM]]^2 + dZ12H1fingl[zM]^2 + 
       (4*dZ11H1gl - 2*dZ22H1gl)*(Conjugate[dZ12H1fingl[zM]] + 
         dZ12H1fingl[zM])) + (SB*(6*dZ11H1gl^2 - 4*dZ11H1gl*dZ22H1gl - 
        Conjugate[dZ12H1fingl[zM]]^2 - dZ12H1fingl[zM]^2 + 
        TB2*(Conjugate[dZ12H1fingl[zM]] + dZ12H1fingl[zM])^2 - 
        2*(dZ22H1gl^2 + Sq[dZ12H1fingl[zM]])))/CB)/16, DebugLine[1, dTB2div], 
  dMHiggs1gl[1, 2] -> CB2*dTB1div*MHin2 - dTHH1*HoldCode[EL1L/(2*MW*SW)], 
  DebugLine[1, dMHiggs1gl[1, 2]], dMHiggs1gl[3, 4] -> -dMHiggs1gl[1, 2], 
  DebugLine[1, dMHiggs1gl[3, 4]], dMHiggs1gl[5, 6] -> 
   -(CB2*dTB1div*MHin2) + (I*dTA01 + dTHH1)*HoldCode[EL1L/(2*MW*SW)], 
  DebugLine[1, dMHiggs1gl[5, 6]], dTB1 -> dTB1div + dTB1fingl[zM], 
  DebugLine[1, dTB1], dAt -> -((dTB1*MUEC)/TB2), DebugLine[1, dAt], 
  dCB1 -> -(CB2*dTB1*SB), DebugLine[1, dCB1], dSB1 -> CB*CB2*dTB1, 
  DebugLine[1, dSB1], dup60 -> dZHiggs1gl[1, 3] + I*dZHiggs1gl[3, 4], 
  dup61 -> I*dZHiggs1gl[1, 3] + dZHiggs1gl[3, 4], 
  dup64 -> -4*dCB1*SB*Ybgl + S2B*(2*dTB1*MUEC + Ybgl*dZHiggs1gl[2, 2]) + 
    CB2*Xbgl*((-2*I)*dZHiggs1gl[1, 3] - 2*dZHiggs1gl[3, 4]), 
  dup65 -> -2*MUE*S2B*Conjugate[dTB1] + 
    YbglC*(4*dCB1*SB - S2B*dZHiggs1gl[2, 2]) + 
    CB2*XbglC*((-2*I)*dZHiggs1gl[1, 3] + 2*dZHiggs1gl[3, 4]), 
  dup89 -> -2*(dup61*SB2*XtglC + S2B*Conjugate[dAt]) + 
    YtglC*(4*CB*dSB1 - S2B*dZHiggs1gl[2, 2]), 
  dup94 -> (2*I)*dup60*SB2*Xtgl + 4*CB*dSB1*Ytgl + 
    S2B*(-2*dAt - Ytgl*dZHiggs1gl[2, 2]), 
  dup116 -> (-4*I)*MTgl2*dZHiggs1gl[1, 3] + 
    MTgl*((8*I)*U2s27*dBn1gl[zM] + (dup89*UCStgl[1, 3] - dup94*UCStglC[1, 3])/
       SB2), dup117 -> (4*I)*MTgl2*dZHiggs1gl[1, 3] + 
    MTgl*((8*I)*U2s27*dBn1gl[zM] + (dup89*UCStgl[1, 3] - dup94*UCStglC[1, 3])/
       SB2), dTA02 -> (-(dTh01*dZHiggs1gl[1, 3]) - dTA01*dZHiggs1gl[2, 2] - 
      dTG01*dZHiggs1gl[3, 4])/2 + 
    (3*EL1L*(8*dZHiggs1gl[1, 3]*(MBgl2*Re[A0q12] + MTgl2*Re[A0q13]) - 
       I*(-(dup117*Re[A0q14]) + dup116*Re[A0q15] - (4*I)*MBgl2*
          dZHiggs1gl[1, 3]*(Re[A0q16] + Re[A0q17]) + 
         (MBgl*MSbgl2[3]*Re[B0q14]*((8*I)*CB2*U2s28*dBn1gl[zM] - 
            dup65*UCSbgl[1, 3] - dup64*UCSbglC[1, 3]))/CB2)))/
     (128*MW*Pi^2*SW), DebugLine[1, dTA02], 
  dup41 -> 8*dSB1*SB - 4*SB2*dZHiggs1gl[1, 1], 
  dup42 -> 2*MUE*Conjugate[dTB1] + XbglC*dZHiggs1gl[1, 1], 
  dup47 -> 2*dA1gl[zM] - dZHiggs1gl[1, 2] - I*dZHiggs1gl[1, 3], 
  dup48 -> dZHiggs1gl[1, 2] - I*dZHiggs1gl[1, 3], 
  dup49 -> 2*dA1gl[zM] - dZHiggs1gl[1, 2] + I*dZHiggs1gl[1, 3], 
  dup50 -> dZHiggs1gl[1, 2] + I*dZHiggs1gl[1, 3], 
  dup83 -> -(dup47*S2B*Ybgl) + CB2*(4*dTB1*MUEC + 2*Xbgl*dZHiggs1gl[1, 1]), 
  dup84 -> 4*dSB1*SB*Xtgl + SB2*(-4*dAt - 2*Xtgl*dZHiggs1gl[1, 1]), 
  dup86 -> -2*SB2*Conjugate[dAt] + XtglC*(2*dSB1*SB - SB2*dZHiggs1gl[1, 1]), 
  dup87 -> 2*dSB1 - SB*dZHiggs1gl[1, 1] + 
    CB*(-2*dA1gl[zM] + dZHiggs1gl[1, 2]), 
  dup88 -> 4*CB*dCB1 - 2*CB2*dZHiggs1gl[1, 1] + 
    S2B*(2*dA1gl[zM] - dZHiggs1gl[1, 2]), 
  dup102 -> 2*MTgl2*dZHiggs1gl[1, 2] - MTgl*(dup50*YtglC*UCStgl[1, 3] + 
      dup48*Ytgl*UCStglC[1, 3]), dup103 -> 2*MTgl2*dZHiggs1gl[1, 2] + 
    MTgl*(dup50*YtglC*UCStgl[1, 3] + dup48*Ytgl*UCStglC[1, 3]), 
  dup113 -> (8*CB*MTgl*(MTgl - U2s13)*dA1gl[zM])/SB + 
    (-(dup41*MTgl2) - dup102*S2B + 2*dup86*MTgl*UCStgl[1, 3] + 
      dup84*MTgl*UCStglC[1, 3])/SB2, 
  dup115 -> (8*CB*MTgl*(MTgl + U2s13)*dA1gl[zM])/SB - 
    (dup41*MTgl2 + dup103*S2B + MTgl*(2*dup86*UCStgl[1, 3] + 
        dup84*UCStglC[1, 3]))/SB2, 
  dTh02 -> (-(dTh01*dZHiggs1gl[1, 1]) - dTHH1*dZHiggs1gl[1, 2] - 
     dTA01*dZHiggs1gl[1, 3] - (EL1L*((24*dup87*MTgl2*Re[A0q13])/SB + 
        3*(dup113*Re[A0q14] + dup115*Re[A0q15]) + 
        (6*dup74*dup88*MBgl2 + 3*MBgl*MSbgl2[3]*Re[B0q14]*(8*CB*dCB1*U2s12 - 
            2*CB2*dup42*UCSbgl[1, 3] + dup49*S2B*YbglC*UCSbgl[1, 3] - 
            dup83*UCSbglC[1, 3]))/CB2))/(64*MW*Pi^2*SW))/2, 
  DebugLine[1, dTh02], B0q15 -> B0q[0, MStgl2[1], MStgl2[2], Q2], 
  dup45 -> dTB1*MUEC*S2B + CB2*Xbgl*dZHiggs1gl[1, 2], 
  dup46 -> MUE*S2B*Conjugate[dTB1] + CB2*XbglC*dZHiggs1gl[1, 2], 
  dup51 -> 4*CB*dSB1 - S2B*dZHiggs1gl[2, 2], 
  dup52 -> 4*dCB1*SB - S2B*dZHiggs1gl[2, 2], 
  dup56 -> dZHiggs1gl[1, 2] + I*dZHiggs1gl[2, 4], 
  dup57 -> 2*(dZHiggs1gl[1, 2] - I*dZHiggs1gl[2, 4]), 
  dup78 -> 4*(CB*dSB1 + SB2*dA1gl[zM]) + 2*SB2*dZHiggs1gl[1, 2] - 
    S2B*dZHiggs1gl[2, 2], 
  dup79 -> 4*dCB1*SB - 2*CB2*(2*dA1gl[zM] + dZHiggs1gl[1, 2]) - 
    S2B*dZHiggs1gl[2, 2], dup93 -> 2*dup56*SB2*Xtgl + 4*CB*dSB1*Ytgl + 
    S2B*(-2*dAt - Ytgl*dZHiggs1gl[2, 2]), 
  dup98 -> 8*CB2*U2s12*dA1gl[zM] + (2*dup46 - dup52*YbglC + 
      (2*I)*CB2*XbglC*dZHiggs1gl[2, 4])*UCSbgl[1, 3] + 
    (2*dup45 - dup52*Ybgl - (2*I)*CB2*Xbgl*dZHiggs1gl[2, 4])*UCSbglC[1, 3], 
  dTHH2 -> (-(dTh01*dZHiggs1gl[1, 2]) - dTHH1*dZHiggs1gl[2, 2] - 
      dTG01*dZHiggs1gl[2, 4])/2 + 
    (3*EL1L*((dup79*MBgl2*(-4*Re[A0q12] + 2*(Re[A0q16] + Re[A0q17])) + 
         dup98*MBgl*MSbgl2[3]*Re[B0q14])/CB2 + 
       (dup78*MTgl2*(4*Re[A0q13] - 2*(Re[A0q14] + Re[A0q15])) + 
         MTgl*MStgl2[3]*Re[B0q15]*((dup51*YtglC - 2*S2B*Conjugate[dAt])*
            UCStgl[1, 3] + SB2*(8*U2s11*dA1gl[zM] + dup57*XtglC*
              UCStgl[1, 3]) + dup93*UCStglC[1, 3]))/SB2))/(128*MW*Pi^2*SW), 
  DebugLine[1, dTHH2], dMHiggs2gl[1, 1] -> CB2^2*dTB1div^2*MHin2 - 
    (dTh02 + CB2*dTB1div*dTHH1)*HoldCode[EL1L/(2*MW*SW)], 
  DebugLine[1, dMHiggs2gl[1, 1]], dMHiggs2gl[1, 3] -> 
   -(dTA02*HoldCode[EL1L/(2*MW*SW)]), DebugLine[1, dMHiggs2gl[1, 3]], 
  dMHiggsZ2gl[1, 1] -> dMHiggs2gl[1, 1] + dMHiggs1gl[1, 1]*dZHiggs1gl[1, 1] + 
    dMHiggs1gl[1, 2]*dZHiggs1gl[1, 2] + dMHiggs1gl[1, 3]*dZHiggs1gl[1, 3] + 
    (MHin2*(dZHiggs1gl[1, 2]^2 + dZHiggs1gl[1, 3]^2))/4, 
  DebugLine[1, dMHiggsZ2gl[1, 1]], 
  dup39 -> 8*dCB1*SB*YbglC - 4*MUE*S2B*Conjugate[dTB1], 
  dup54 -> 2*MUE*Conjugate[dTB1] + YbglC*dZHiggs1gl[2, 2], 
  dup59 -> (-2*I)*dZHiggs1gl[1, 3] - 2*dZHiggs1gl[3, 4], 
  dup62 -> 2*((-I)*dZHiggs1gl[1, 3] + dZHiggs1gl[3, 4]), 
  dup63 -> 2*SB2*dZHiggs1gl[2, 2] - S2B*dZHiggs1gl[3, 4], 
  dup67 -> 2*SB2*dZHiggs1gl[2, 2] - S2B*dZHiggs1gl[5, 6], 
  dup68 -> 2*CB2*dZHiggs1gl[2, 2] + S2B*dZHiggs1gl[5, 6], 
  dup69 -> 2*SB2*YbglC*dZHiggs1gl[2, 2] - S2B*XbglC*dZHiggs1gl[5, 6], 
  dup70 -> 2*(dBc1gl[zM] + dBn1gl[zM]) - dZHiggs1gl[5, 6] - dZHiggs1gl[6, 5], 
  dup71 -> 2*CB2*dZHiggs1gl[2, 2] + S2B*dZHiggs1gl[6, 5], 
  dup72 -> 8*dCB1*SB*SB2 + S2B^2*dZHiggs1gl[6, 5], 
  dup73 -> 8*dCB1*SB*SB2*Ybgl + Abgl*S2B^2*dZHiggs1gl[6, 5], 
  dup77 -> 4*(CB*CB2*dSB1 + dCB1*SB*SB2) - S2B*dZHiggs1gl[2, 2], 
  dup80 -> S2B*dZHiggs1gl[2, 2] + SB2*(dZHiggs1gl[5, 6] + dZHiggs1gl[6, 5]), 
  dup95 -> CB2*(4*dAt + 2*Ytgl*dZHiggs1gl[2, 2]) + S2B*Xtgl*dZHiggs1gl[5, 6], 
  IndexIf[inputmass == A0A0, {A0q1 -> A0q[MTgl2, Q2], 
    A0q2 -> A0q[MStgl2[1], Q2], A0q3 -> A0q[MStgl2[2], Q2], 
    A0q4 -> A0q[MBgl2, Q2], A0q5 -> A0q[MSbgl2[1], Q2], 
    A0q6 -> A0q[MSbgl2[2], Q2], B0q1 -> B0q[0, MStgl2[1], MStgl2[2], Q2], 
    B0q2 -> B0q[0, MStgl2[1], MStgl2[1], Q2], 
    B0q3 -> B0q[0, MStgl2[2], MStgl2[2], Q2], 
    B0q4 -> B0q[0, MSbgl2[1], MSbgl2[2], Q2], 
    B0q5 -> B0q[0, MSbgl2[1], MSbgl2[1], Q2], 
    B0q6 -> B0q[0, MSbgl2[2], MSbgl2[2], Q2], U2c21 -> U2c2[UCStgl, Ytgl], 
    U2c22 -> U2c2[UCSbgl, Ybgl], U2s21 -> U2s2[UCStgl, Ytgl], 
    U2s22 -> U2s2[UCSbgl, Ybgl], dMHinsq1 -> -(MHin2*dZHiggs1gl[2, 2]) - 
      (3*Alfa1L*(MTgl2*(S2B^2*Re[A0q1] - 2*CB2^2*TB2*(Re[A0q2] + Re[A0q3] + 
             2*U2s21^2*(Re[B0q2] + Re[B0q3]) + 4*Re[B0q1]*Sq[U2c21])) + 
         MBgl2*TB2*(S2B^2*TB2*Re[A0q4] - 2*SB2^2*(Re[A0q5] + Re[A0q6] + 
             2*U2s22^2*(Re[B0q5] + Re[B0q6]) + 4*Re[B0q4]*Sq[U2c22]))))/
       (16*MW2*Pi*SB2^2*SW2), DebugLine[1, dMHinsq1], 
    dMHiggs1gl[3, 3] -> dMHinsq1, DebugLine[1, dMHiggs1gl[3, 3]], 
    U2c23 -> U2c2[UCSbgl, Xbgl], U2c24 -> U2c2[UCStgl, Xtgl], 
    U2s23 -> U2s2[UCSbgl, Xbgl], U2s24 -> U2s2[UCStgl, Xtgl], 
    dup4 -> CB2*dup59*Xbgl - 4*dCB1*SB*Ybgl + 
      S2B*(2*dTB1*MUEC + Ybgl*dZHiggs1gl[2, 2]), 
    dup5 -> CB2*dup62*XbglC + dup52*YbglC - 2*MUE*S2B*Conjugate[dTB1], 
    dup6 -> dup51*YtglC - 2*(dup61*SB2*XtglC + S2B*Conjugate[dAt]), 
    dup7 -> (2*I)*dup60*SB2*Xtgl + 4*CB*dSB1*Ytgl + 
      S2B*(-2*dAt - Ytgl*dZHiggs1gl[2, 2]), 
    dup8 -> (-(dup54*S2B) + 4*dCB1*SB*YbglC)*UCSbgl[1, 3] - 
      (2*I)*CB2*(4*U2s23*dBn1gl[zM] + dup60*XbglC*UCSbgl[1, 3]) + 
      dup4*UCSbglC[1, 3], dMHinsq2 -> -(dMHiggs1gl[1, 3]*dZHiggs1gl[1, 3]) - 
      dMHiggs1gl[3, 3]*dZHiggs1gl[2, 2] - dMHiggs1gl[3, 4]*dZHiggs1gl[3, 4] - 
      (MHin2*(dZHiggs1gl[2, 2]^2 + 4*dZHiggs2gl[2, 2]))/4 - 
      (3*Alfa1L*(-(MTgl2*S2B*(8*CB*CB2*dSB1 + S2B*(2*S2B*dBn1gl[zM] - 
              2*CB2*dZHiggs1gl[2, 2] - S2B*dZHiggs1gl[3, 4]))*Re[A0q1]) + 
         2*CB2*MTgl2*SB*(4*CB2*dSB1 + CB*(4*SB2*dBn1gl[zM] - 
             S2B*dZHiggs1gl[2, 2] - 2*SB2*dZHiggs1gl[3, 4]))*
          (Re[A0q2] + Re[A0q3]) + MBgl2*SB*TB*(dup63*S2B*SB - 8*dCB1*SB2^2 + 
           2*S2B^2*SB*dBn1gl[zM])*(2*Re[A0q4] - Re[A0q5] - Re[A0q6]) + 
         2*MBgl*SB2^2*TB*U2s22*(((-I)*dup8*MBgl - 4*CB2*MBgl2*
              dZHiggs1gl[1, 3])*Re[B0q5] + ((-I)*dup8*MBgl + 
             4*CB2*MBgl2*dZHiggs1gl[1, 3])*Re[B0q6]) - 2*MBgl^2*SB2^2*TB*
          Re[B0q4]*(Conjugate[U2c22]*(8*CB2*U2c23*dBn1gl[zM] - 
             dup5*UCSbgl[3, 3] - dup4*UCSbgl[3, 4]) + 
           U2c22*(8*CB2*Conjugate[U2c23]*dBn1gl[zM] + dup4*UCSbglC[3, 3] + 
             dup5*UCSbglC[3, 4])) + CB2*MTgl*S2B*
          (I*U2s21*(Re[B0q3]*((-8*I)*MTgl*SB2*U2s24*dBn1gl[zM] - (4*I)*MTgl2*
                SB2*dZHiggs1gl[1, 3] - dup6*MTgl*UCStgl[1, 3] + dup7*MTgl*
                UCStglC[1, 3]) + Re[B0q2]*((-8*I)*MTgl*SB2*U2s24*dBn1gl[zM] + 
               (4*I)*MTgl2*SB2*dZHiggs1gl[1, 3] - dup6*MTgl*UCStgl[1, 3] + 
               dup7*MTgl*UCStglC[1, 3])) + MTgl*Re[B0q1]*
            (Conjugate[U2c21]*(8*SB2*U2c24*dBn1gl[zM] + dup6*UCStgl[3, 3] - 
               dup7*UCStgl[3, 4]) + U2c21*(8*SB2*Conjugate[U2c24]*
                dBn1gl[zM] + dup7*UCStglC[3, 3] - (4*CB*dSB1*YtglC - 
                 2*S2B*Conjugate[dAt] - (2*I)*SB2*XtglC*dZHiggs1gl[1, 3] - 
                 S2B*YtglC*dZHiggs1gl[2, 2] - 2*SB2*XtglC*dZHiggs1gl[3, 4])*
                UCStglC[3, 4])))))/(8*MW2*Pi*S2B^2*SB2*SW2), 
    DebugLine[1, dMHinsq2]}, True, {A0q7 -> A0q[MBgl2, Q2], 
    A0q8 -> A0q[MSbgl2[2], Q2], A0q9 -> A0q[MSbgl2[1], Q2], 
    A0q10 -> A0q[MStgl2[1], Q2], A0q11 -> A0q[MStgl2[2], Q2], 
    B0q7 -> B0q[0, MBgl2, MTgl2, Q2], B0q8 -> B0q[0, MStgl2[1], MSbgl2[1], 
      Q2], B0q9 -> B0q[0, MStgl2[1], MSbgl2[2], Q2], 
    B0q10 -> B0q[0, MSbgl2[1], MStgl2[2], Q2], 
    B0q11 -> B0q[0, MStgl2[2], MSbgl2[2], Q2], 
    dup17 -> (CB2*Ytgl*USbglC[1, 1] + MBgl*USbglC[1, 2])*UStgl[1, 2], 
    dup18 -> (CB2*Ytgl*USbglC[2, 1] + MBgl*USbglC[2, 2])*UStgl[1, 2], 
    dup19 -> (CB2*Ytgl*USbglC[1, 1] + MBgl*USbglC[1, 2])*UStgl[2, 2], 
    dup20 -> (CB2*Ytgl*USbglC[2, 1] + MBgl*USbglC[2, 2])*UStgl[2, 2], 
    dup21 -> (CB2*YtglC*USbgl[1, 1] + MBgl*USbgl[1, 2])*UStglC[1, 2], 
    dup22 -> (CB2*YtglC*USbgl[2, 1] + MBgl*USbgl[2, 2])*UStglC[1, 2], 
    dup23 -> (CB2*YtglC*USbgl[1, 1] + MBgl*USbgl[1, 2])*UStglC[2, 2], 
    dup24 -> (CB2*YtglC*USbgl[2, 1] + MBgl*USbgl[2, 2])*UStglC[2, 2], 
    dup30 -> dup17*MTgl + ((CB2*MTgl2 + MBgl2*SB2)*USbglC[1, 1] + 
        MBgl*SB2*YbglC*USbglC[1, 2])*UStgl[1, 1], 
    dup31 -> dup18*MTgl + ((CB2*MTgl2 + MBgl2*SB2)*USbglC[2, 1] + 
        MBgl*SB2*YbglC*USbglC[2, 2])*UStgl[1, 1], 
    dup32 -> dup19*MTgl + ((CB2*MTgl2 + MBgl2*SB2)*USbglC[1, 1] + 
        MBgl*SB2*YbglC*USbglC[1, 2])*UStgl[2, 1], 
    dup33 -> dup20*MTgl + ((CB2*MTgl2 + MBgl2*SB2)*USbglC[2, 1] + 
        MBgl*SB2*YbglC*USbglC[2, 2])*UStgl[2, 1], 
    dup34 -> dup21*MTgl + ((CB2*MTgl2 + MBgl2*SB2)*USbgl[1, 1] + 
        MBgl*SB2*Ybgl*USbgl[1, 2])*UStglC[1, 1], 
    dup35 -> dup22*MTgl + ((CB2*MTgl2 + MBgl2*SB2)*USbgl[2, 1] + 
        MBgl*SB2*Ybgl*USbgl[2, 2])*UStglC[1, 1], 
    dup36 -> dup23*MTgl + ((CB2*MTgl2 + MBgl2*SB2)*USbgl[1, 1] + 
        MBgl*SB2*Ybgl*USbgl[1, 2])*UStglC[2, 1], 
    dup37 -> dup24*MTgl + ((CB2*MTgl2 + MBgl2*SB2)*USbgl[2, 1] + 
        MBgl*SB2*Ybgl*USbgl[2, 2])*UStglC[2, 1], 
    dMHinsq1 -> -(MHin2*dZHiggs1gl[2, 2]) - 
      (3*Alfa1L*(((MTgl2 + MBgl2*TB2^2)*Re[A0q7] + 
           MTgl2*(MTgl2 + MBgl2*TB2*(2 + TB2))*Re[B0q7])/TB2 - 
         (2*(dup32*dup36*Re[B0q10] + dup33*dup37*Re[B0q11] + 
            dup30*dup34*Re[B0q8] + dup31*dup35*Re[B0q9] + 
            (MBgl2*SB2^2*Re[A0q8] + CB2^2*MTgl2*Re[A0q9])*USbgl2[1, 1] + 
            (CB2^2*MTgl2*Re[A0q8] + MBgl2*SB2^2*Re[A0q9])*USbgl2[1, 2] + 
            (MBgl2*SB2^2*Re[A0q10] + CB2^2*MTgl2*Re[A0q11])*UStgl2[1, 1] + 
            (CB2^2*MTgl2*Re[A0q10] + MBgl2*SB2^2*Re[A0q11])*UStgl2[1, 2]))/
          S2B^2))/(4*MW2*Pi*SW2), DebugLine[1, dMHinsq1], 
    dup25 -> dup72*MBgl2 + MTgl2*(8*CB*CB2*dSB1 - dup71*S2B), 
    dup26 -> MTgl2*(8*CB*CB2*dSB1 - dup68*S2B) + 
      MBgl2*(-(dup67*S2B) + 8*dCB1*SB*SB2), 
    dup28 -> dup70*S2B + 4*SB2*dZHiggs1gl[2, 2], 
    dup29 -> 4*CB*dSB1 + dup70*SB2 - S2B*dZHiggs1gl[2, 2], 
    tmp1 -> Re[B0q8]*(dup34*(MBgl*USbglC[1, 2]*((-(dup69*S2B) + dup39*SB2)*
             UStgl[1, 1] + 2*dup77*MTgl*UStgl[1, 2]) + USbglC[1, 1]*
           (dup26*UStgl[1, 1] + MTgl*(-(dup95*S2B) + 8*CB*CB2*dSB1*Ytgl)*
             UStgl[1, 2]) + 2*S2B^2*dBc1gl[zM]*
           ((MTgl2*USbglC[1, 1] - MBgl*(MBgl*USbglC[1, 1] + XbglC*
                 USbglC[1, 2]))*UStgl[1, 1] + MTgl*Xtgl*USbglC[1, 1]*
             UStgl[1, 2])) + dup30*((dup25*USbgl[1, 1] + 
            dup73*MBgl*USbgl[1, 2])*UStglC[1, 1] + 
          MTgl*((8*CB*CB2*dSB1*YtglC - S2B^2*XtglC*dZHiggs1gl[6, 5])*
             USbgl[1, 1] + 8*dCB1*MBgl*SB*SB2*USbgl[1, 2])*UStglC[1, 2] - 
          2*(MBgl*USbgl[1, 2]*(MUEC*S2B*SB2*dZHiggs1gl[6, 5]*UStglC[1, 1] - 
              4*CB*CB2*dSB1*MTgl*UStglC[1, 2]) + 
            S2B*(2*CB2*MTgl*Conjugate[dAt]*USbgl[1, 1]*UStglC[1, 2] + 
              MBgl*USbgl[1, 2]*(2*dTB1*MUEC*SB2*UStglC[1, 1] + 
                MTgl*dZHiggs1gl[2, 2]*UStglC[1, 2]) + S2B*dBn1gl[zM]*(
                ((MBgl2 - MTgl2)*USbgl[1, 1] + MBgl*Xbgl*USbgl[1, 2])*
                 UStglC[1, 1] - MTgl*XtglC*USbgl[1, 1]*UStglC[1, 2]) + 
              dZHiggs1gl[2, 2]*(SB2*(MBgl2*USbgl[1, 1] + MBgl*Ybgl*
                   USbgl[1, 2])*UStglC[1, 1] + CB2*MTgl*YtglC*USbgl[1, 1]*
                 UStglC[1, 2]))))) + Re[B0q9]*
       (dup35*(MBgl*USbglC[2, 2]*((-(dup69*S2B) + dup39*SB2)*UStgl[1, 1] + 
            2*dup77*MTgl*UStgl[1, 2]) + USbglC[2, 1]*(dup26*UStgl[1, 1] + 
            MTgl*(-(dup95*S2B) + 8*CB*CB2*dSB1*Ytgl)*UStgl[1, 2]) + 
          2*S2B^2*dBc1gl[zM]*((MTgl2*USbglC[2, 1] - MBgl*(MBgl*USbglC[2, 1] + 
                XbglC*USbglC[2, 2]))*UStgl[1, 1] + MTgl*Xtgl*USbglC[2, 1]*
             UStgl[1, 2])) + dup31*((dup25*USbgl[2, 1] + 
            dup73*MBgl*USbgl[2, 2])*UStglC[1, 1] + 
          MTgl*((8*CB*CB2*dSB1*YtglC - S2B^2*XtglC*dZHiggs1gl[6, 5])*
             USbgl[2, 1] + 8*dCB1*MBgl*SB*SB2*USbgl[2, 2])*UStglC[1, 2] - 
          2*(MBgl*USbgl[2, 2]*(MUEC*S2B*SB2*dZHiggs1gl[6, 5]*UStglC[1, 1] - 
              4*CB*CB2*dSB1*MTgl*UStglC[1, 2]) + 
            S2B*(2*CB2*MTgl*Conjugate[dAt]*USbgl[2, 1]*UStglC[1, 2] + 
              MBgl*USbgl[2, 2]*(2*dTB1*MUEC*SB2*UStglC[1, 1] + 
                MTgl*dZHiggs1gl[2, 2]*UStglC[1, 2]) + S2B*dBn1gl[zM]*(
                ((MBgl2 - MTgl2)*USbgl[2, 1] + MBgl*Xbgl*USbgl[2, 2])*
                 UStglC[1, 1] - MTgl*XtglC*USbgl[2, 1]*UStglC[1, 2]) + 
              dZHiggs1gl[2, 2]*(SB2*(MBgl2*USbgl[2, 1] + MBgl*Ybgl*
                   USbgl[2, 2])*UStglC[1, 1] + CB2*MTgl*YtglC*USbgl[2, 1]*
                 UStglC[1, 2]))))), RuleAdd[tmp1, 
     Re[B0q10]*(dup36*(MBgl*USbglC[1, 2]*((-(dup69*S2B) + dup39*SB2)*
             UStgl[2, 1] + 2*dup77*MTgl*UStgl[2, 2]) + USbglC[1, 1]*
           (dup26*UStgl[2, 1] + MTgl*(-(dup95*S2B) + 8*CB*CB2*dSB1*Ytgl)*
             UStgl[2, 2]) + 2*S2B^2*dBc1gl[zM]*
           ((MTgl2*USbglC[1, 1] - MBgl*(MBgl*USbglC[1, 1] + XbglC*
                 USbglC[1, 2]))*UStgl[2, 1] + MTgl*Xtgl*USbglC[1, 1]*
             UStgl[2, 2])) + dup32*((dup25*USbgl[1, 1] + 
            dup73*MBgl*USbgl[1, 2])*UStglC[2, 1] + 
          MTgl*((8*CB*CB2*dSB1*YtglC - S2B^2*XtglC*dZHiggs1gl[6, 5])*
             USbgl[1, 1] + 8*dCB1*MBgl*SB*SB2*USbgl[1, 2])*UStglC[2, 2] - 
          2*(MBgl*USbgl[1, 2]*(MUEC*S2B*SB2*dZHiggs1gl[6, 5]*UStglC[2, 1] - 
              4*CB*CB2*dSB1*MTgl*UStglC[2, 2]) + 
            S2B*(2*CB2*MTgl*Conjugate[dAt]*USbgl[1, 1]*UStglC[2, 2] + 
              MBgl*USbgl[1, 2]*(2*dTB1*MUEC*SB2*UStglC[2, 1] + 
                MTgl*dZHiggs1gl[2, 2]*UStglC[2, 2]) + S2B*dBn1gl[zM]*(
                ((MBgl2 - MTgl2)*USbgl[1, 1] + MBgl*Xbgl*USbgl[1, 2])*
                 UStglC[2, 1] - MTgl*XtglC*USbgl[1, 1]*UStglC[2, 2]) + 
              dZHiggs1gl[2, 2]*(SB2*(MBgl2*USbgl[1, 1] + MBgl*Ybgl*
                   USbgl[1, 2])*UStglC[2, 1] + CB2*MTgl*YtglC*USbgl[1, 1]*
                 UStglC[2, 2]))))) + Re[B0q11]*
       (dup37*(MBgl*USbglC[2, 2]*((-(dup69*S2B) + dup39*SB2)*UStgl[2, 1] + 
            2*dup77*MTgl*UStgl[2, 2]) + USbglC[2, 1]*(dup26*UStgl[2, 1] + 
            MTgl*(-(dup95*S2B) + 8*CB*CB2*dSB1*Ytgl)*UStgl[2, 2]) + 
          2*S2B^2*dBc1gl[zM]*((MTgl2*USbglC[2, 1] - MBgl*(MBgl*USbglC[2, 1] + 
                XbglC*USbglC[2, 2]))*UStgl[2, 1] + MTgl*Xtgl*USbglC[2, 1]*
             UStgl[2, 2])) + dup33*((dup25*USbgl[2, 1] + 
            dup73*MBgl*USbgl[2, 2])*UStglC[2, 1] + 
          MTgl*((8*CB*CB2*dSB1*YtglC - S2B^2*XtglC*dZHiggs1gl[6, 5])*
             USbgl[2, 1] + 8*dCB1*MBgl*SB*SB2*USbgl[2, 2])*UStglC[2, 2] - 
          2*(MBgl*USbgl[2, 2]*(MUEC*S2B*SB2*dZHiggs1gl[6, 5]*UStglC[2, 1] - 
              4*CB*CB2*dSB1*MTgl*UStglC[2, 2]) + 
            S2B*(2*CB2*MTgl*Conjugate[dAt]*USbgl[2, 1]*UStglC[2, 2] + 
              MBgl*USbgl[2, 2]*(2*dTB1*MUEC*SB2*UStglC[2, 1] + 
                MTgl*dZHiggs1gl[2, 2]*UStglC[2, 2]) + S2B*dBn1gl[zM]*(
                ((MBgl2 - MTgl2)*USbgl[2, 1] + MBgl*Xbgl*USbgl[2, 2])*
                 UStglC[2, 1] - MTgl*XtglC*USbgl[2, 1]*UStglC[2, 2]) + 
              dZHiggs1gl[2, 2]*(SB2*(MBgl2*USbgl[2, 1] + MBgl*Ybgl*
                   USbgl[2, 2])*UStglC[2, 1] + CB2*MTgl*YtglC*USbgl[2, 1]*
                 UStglC[2, 2])))))], dMHinsq2 -> 
     -(dMHinsq1*dZHiggs1gl[2, 2]) + 
      (-(Conjugate[dMHiggs1gl[5, 6]]*dZHiggs1gl[5, 6]) - 
        dMHiggs1gl[5, 6]*dZHiggs1gl[6, 5])/2 - 
      (MHin2*(dZHiggs1gl[2, 2]^2 + 4*dZHiggs2gl[2, 2]))/4 + 
      (3*Alfa1L*(2*((CB*MTgl2*(8*CB*dSB1 - 2*dup80 + 4*SB2*(dBc1gl[zM] + 
                dBn1gl[zM])))/(SB*SB2) - (2*MBgl2*SB*(CB2*dup70 - 4*dCB1*SB + 
              S2B*dZHiggs1gl[2, 2]))/(CB*CB2))*Re[A0q7] - 
         (2*MTgl2*(4*MBgl2*SB2*(CB2*dup70 - 4*(CB*dSB1 + dCB1*SB*(1 + TB2)) + 
               S2B*(2 + TB2)*dZHiggs1gl[2, 2]) - MTgl2*(16*CB*CB2*dSB1 + S2B*
                (dup70*S2B - 4*CB2*dZHiggs1gl[2, 2])))*Re[B0q7] + 
           (tmp1 + 4*(Re[A0q8]*(-(MBgl2*SB2^2*(CB2*dup70 - 4*dCB1*SB + 
                    S2B*dZHiggs1gl[2, 2])*USbgl2[1, 1]) + CB2^2*dup29*MTgl2*
                  USbgl2[1, 2]) + Re[A0q9]*(CB2^2*dup29*MTgl2*USbgl2[1, 1] - 
                 MBgl2*SB2^2*(CB2*dup70 - 4*dCB1*SB + S2B*dZHiggs1gl[2, 2])*
                  USbgl2[1, 2])) + Re[A0q10]*(MBgl2*SB2*(-(dup28*S2B) + 
                 16*dCB1*SB*SB2)*UStgl2[1, 1] + 4*CB2^2*dup29*MTgl2*
                UStgl2[1, 2]) + Re[A0q11]*(4*CB2^2*dup29*MTgl2*UStgl2[1, 1] + 
               MBgl2*SB2*(-(dup28*S2B) + 16*dCB1*SB*SB2)*UStgl2[1, 2]))/CB2)/
          (S2B*SB2)))/(32*MW2*Pi*SW2), DebugLine[1, dMHinsq2], 
    dMHiggs1gl[3, 3] -> dMHinsq1, DebugLine[1, dMHiggs1gl[3, 3]]}], 
  dMHiggs2gl[1, 2] -> 
   CB2*(dMHinsq1*dTB1div + MHin2*(dTB2div - (dTB1div^2*S2B)/2)) - 
    dTHH2*HoldCode[EL1L/(2*MW*SW)], DebugLine[1, dMHiggs2gl[1, 2]], 
  dMHiggsZ2gl[1, 2] -> dMHiggs2gl[1, 2] + 
    ((dMHinsq1 + dMHiggs1gl[1, 1])*dZHiggs1gl[1, 2] + 
      dMHiggs1gl[1, 2]*(dZHiggs1gl[1, 1] + dZHiggs1gl[2, 2]))/2 + 
    (MHin2*(dZHiggs1gl[1, 2]*dZHiggs1gl[2, 2] + 2*dZHiggs2gl[1, 2]))/4, 
  DebugLine[1, dMHiggsZ2gl[1, 2]], dMHiggsZ2gl[1, 3] -> 
   dMHiggs2gl[1, 3] + ((dMHiggs1gl[1, 1] + dMHiggs1gl[3, 3])*
       dZHiggs1gl[1, 3] + dMHiggs1gl[1, 3]*(dZHiggs1gl[1, 1] + 
        dZHiggs1gl[2, 2]))/2 + (MHin2*(dZHiggs1gl[1, 3]*dZHiggs1gl[2, 2] + 
       2*dZHiggs2gl[1, 3]))/4, DebugLine[1, dMHiggsZ2gl[1, 3]], 
  dup81 -> dZHiggs1gl[2, 2]^2/4 + dZHiggs2gl[2, 2], 
  dMHiggsZ2gl[2, 2] -> dMHinsq2 + dup81*MHin2 + 
    dMHiggs1gl[1, 2]*dZHiggs1gl[1, 2] + dMHinsq1*dZHiggs1gl[2, 2] + 
    dMHiggs1gl[2, 4]*dZHiggs1gl[2, 4], DebugLine[1, dMHiggsZ2gl[2, 2]], 
  dMHiggsZ2gl[2, 3] -> (dMHiggs1gl[1, 3]*dZHiggs1gl[1, 2] + 
     dMHiggs1gl[1, 2]*dZHiggs1gl[1, 3] + dMHiggs1gl[3, 4]*dZHiggs1gl[2, 4] + 
     dMHiggs1gl[2, 4]*dZHiggs1gl[3, 4])/2, DebugLine[1, dMHiggsZ2gl[2, 3]], 
  dMHiggsZ2gl[3, 3] -> dMHinsq2 + dup81*MHin2 + 
    dMHiggs1gl[1, 3]*dZHiggs1gl[1, 3] + dMHiggs1gl[3, 3]*dZHiggs1gl[2, 2] + 
    dMHiggs1gl[3, 4]*dZHiggs1gl[3, 4], DebugLine[1, dMHiggsZ2gl[3, 3]], 
  dMHiggsZ2gl[5, 5] -> dMHinsq2 + dup81*MHin2 + dMHinsq1*dZHiggs1gl[2, 2] + 
    (Conjugate[dMHiggs1gl[5, 6]]*dZHiggs1gl[5, 6] + 
      dMHiggs1gl[5, 6]*dZHiggs1gl[6, 5])/2, DebugLine[1, dMHiggsZ2gl[5, 5]], 
  B0q16 -> B0q[0, MBgl2, MBgl2, Q2], 
  B0q17 -> B0q[0, MSbgl2[1], MSbgl2[1], Q2], 
  B0q18 -> B0q[0, MSbgl2[2], MSbgl2[2], Q2], 
  B0q19 -> B0q[0, MTgl2, MTgl2, Q2], 
  B0q20 -> B0q[0, MStgl2[2], MStgl2[2], Q2], 
  B0q21 -> B0q[0, MStgl2[1], MStgl2[1], Q2], U2c11 -> U2c1[UCSbgl, Xbgl], 
  U2c12 -> U2c1[UCSbgl, Ybgl], U2c13 -> U2c1[UCStgl, Xtgl], 
  U2c14 -> U2c1[UCStgl, Ytgl], dup40 -> 2*CB*dCB1 - CB2*dZHiggs1gl[1, 1], 
  dup43 -> 4*Conjugate[dAt] + 2*XtglC*dZHiggs1gl[1, 1], 
  dup85 -> 4*dSB1*SB*Xtgl + dup48*S2B*Ytgl + 
    SB2*(-4*dAt - 2*Xtgl*dZHiggs1gl[1, 1]), 
  dup97 -> 4*(S2B*U2c12*dA1gl[zM] + CB*dCB1*XbglC*UCSbgl[3, 3]) + 
    2*dup40*Xbgl*UCSbgl[3, 4] - 2*CB2*(dup42*UCSbgl[3, 3] + 
      2*dTB1*MUEC*UCSbgl[3, 4]) - S2B*(dup48*YbglC*UCSbgl[3, 3] + 
      dup50*Ybgl*UCSbgl[3, 4]), 
  dup101 -> (4*CB2*dTB1*MUEC - 2*dup40*Xbgl)*UCSbglC[3, 3] + 
    (2*CB2*dup42 - 4*CB*dCB1*XbglC)*UCSbglC[3, 4] + 
    S2B*(-4*Conjugate[U2c12]*dA1gl[zM] + dup50*Ybgl*UCSbglC[3, 3] + 
      dup48*YbglC*UCSbglC[3, 4]), dup106 -> 8*CB*SB2*U2c14*dA1gl[zM] - 
    SB*((-(dup43*SB2) + 4*dSB1*SB*XtglC + dup50*S2B*YtglC)*UCStgl[3, 3] + 
      dup85*UCStgl[3, 4]), dup108 -> 8*CB*SB2*Conjugate[U2c14]*dA1gl[zM] - 
    SB*(dup85*UCStglC[3, 3] + (-(dup43*SB2) + 4*dSB1*SB*XtglC + 
        dup50*S2B*YtglC)*UCStglC[3, 4]), 
  se["h0h0"] -> 
   (3*Alfa1L*(CB2*(MTgl*(SB*SB2*(B0q20*dup113*(MTgl - U2s11) + 
            B0q21*dup115*(MTgl + U2s11)) + B0q15*MTgl*(dup108*U2c13 + 
            dup106*Conjugate[U2c13])) - MTgl2*(A0q14 + A0q15 - 4*B0q19*MTgl2)*
         SB*(4*dSB1*SB - 2*(S2B*dA1gl[zM] + SB2*dZHiggs1gl[1, 1]) + 
          S2B*dZHiggs1gl[1, 2])) + 
      SB2*(2*(2*A0q13*CB2*dup87*MTgl2 + MBgl2*(2*A0q12 - A0q16 - A0q17 + 
            4*B0q16*MBgl2)*(dCB1*S2B - CB2*SB*dZHiggs1gl[1, 1] + 
            CB*SB2*(2*dA1gl[zM] - dZHiggs1gl[1, 2]))) - 
        MBgl*SB*(B0q14*MBgl*(-(dup101*U2c11) + dup97*Conjugate[U2c11]) + 
          B0q17*(MBgl + U2s12)*(2*dup88*MBgl2 + MBgl*(8*CB*dCB1*U2s12 + 
              (-2*CB2*dup42 + dup49*S2B*YbglC)*UCSbgl[1, 3] - 
              dup83*UCSbglC[1, 3])) + B0q18*(MBgl - U2s12)*
           (2*dup88*MBgl2 + MBgl*(-8*CB*dCB1*U2s12 + (-(dup49*S2B*YbglC) + 
                CB2*(4*MUE*Conjugate[dTB1] + 2*XbglC*dZHiggs1gl[1, 1]))*
               UCSbgl[1, 3] + dup83*UCSbglC[1, 3]))))))/
    (4*MW2*Pi*S2B^2*SB*SW2), DebugLine[1, se["h0h0"]], 
  dup53 -> 4*CB*dSB1 + 2*SB2*dZHiggs1gl[1, 2] - S2B*dZHiggs1gl[2, 2], 
  dup55 -> dZHiggs1gl[1, 2] - I*dZHiggs1gl[2, 4], 
  dup91 -> 4*dCB1*Ybgl - 2*CB*(2*dTB1*MUEC + Ybgl*dZHiggs1gl[2, 2]), 
  dup96 -> 8*CB2*U2c11*dA1gl[zM] + (2*dup46 - dup52*YbglC + 
      (2*I)*CB2*XbglC*dZHiggs1gl[2, 4])*UCSbgl[3, 3] + 
    (2*dup45 - dup52*Ybgl - (2*I)*CB2*Xbgl*dZHiggs1gl[2, 4])*UCSbgl[3, 4], 
  dup100 -> 8*CB2*Conjugate[U2c11]*dA1gl[zM] + 
    (2*dup45 - dup52*Ybgl - (2*I)*CB2*Xbgl*dZHiggs1gl[2, 4])*UCSbglC[3, 3] + 
    (2*dup46 - dup52*YbglC + (2*I)*CB2*XbglC*dZHiggs1gl[2, 4])*UCSbglC[3, 4], 
  dup104 -> (-(dup54*S2B) + 4*dCB1*SB*YbglC)*UCSbgl[1, 3] + 
    dup91*SB*UCSbglC[1, 3] - 2*CB2*(4*U2s12*dA1gl[zM] + 
      dup56*XbglC*UCSbgl[1, 3] + dup55*Xbgl*UCSbglC[1, 3]), 
  se["HHHH"] -> 
   (Alfa1L*((SB*((3*dup79*MBgl2*(A0q12 + 2*B0q16*MBgl2))/8 - 
          (3*MBgl*(B0q18*(-(dup104*MBgl) + 2*dup79*MBgl2)*(MBgl - U2s14) + 
             B0q17*(dup104*MBgl + 2*dup79*MBgl2)*(MBgl + U2s14) - 
             B0q14*MBgl*(dup100*U2c12 + dup96*Conjugate[U2c12])))/16) + 
        (3*(A0q16 + A0q17)*MBgl2*(-4*dCB1*SB2 + 
           CB*(S2B*(2*dA1gl[zM] + dZHiggs1gl[1, 2]) + 
             2*SB2*dZHiggs1gl[2, 2])))/16)/(CB*CB2) + 
      (3*(2*MTgl2*(A0q13 + 2*B0q19*MTgl2)*(4*CB2*dSB1 + 
           CB*(SB2*(4*dA1gl[zM] + 2*dZHiggs1gl[1, 2]) - 
             S2B*dZHiggs1gl[2, 2])) - CB*((A0q14 + A0q15)*dup78*MTgl2 + 
           MTgl*(B0q20*(MTgl - U2s13)*(2*dup53*MTgl2 + MTgl*
                (8*SB2*(MTgl - U2s11)*dA1gl[zM] - (dup57*SB2*XtglC + 
                   dup51*YtglC - 2*S2B*Conjugate[dAt])*UCStgl[1, 3] - 
                 dup93*UCStglC[1, 3])) + B0q21*(MTgl + U2s13)*
              (2*dup53*MTgl2 + MTgl*(8*SB2*(MTgl + U2s11)*dA1gl[zM] + 
                 (dup57*SB2*XtglC + dup51*YtglC - 2*S2B*Conjugate[dAt])*
                  UCStgl[1, 3] + dup93*UCStglC[1, 3])) + B0q15*MTgl*
              (Conjugate[U2c14]*((dup51*YtglC - 2*S2B*Conjugate[dAt])*
                  UCStgl[3, 3] + SB2*(8*U2c13*dA1gl[zM] + dup57*XtglC*
                    UCStgl[3, 3]) + dup93*UCStgl[3, 4]) + U2c14*
                (dup93*UCStglC[3, 3] + (dup51*YtglC - 2*S2B*Conjugate[dAt])*
                  UCStglC[3, 4] + SB2*(8*Conjugate[U2c13]*dA1gl[zM] + 
                   dup57*XtglC*UCStglC[3, 4])))))))/(16*SB*SB2)))/
    (MW2*Pi*SW2), DebugLine[1, se["HHHH"]], U2c25 -> U2c2[UCSbgl, Ybgl], 
  U2c26 -> U2c2[UCSbgl, Xbgl], U2c27 -> U2c2[UCStgl, Ytgl], 
  U2c28 -> U2c2[UCStgl, Xtgl], dup66 -> S2B*Conjugate[dAt] + 
    SB2*XtglC*dZHiggs1gl[3, 4], dup92 -> -(CB2*dup59*Xbgl) + 4*dCB1*SB*Ybgl + 
    S2B*(-2*dTB1*MUEC - Ybgl*dZHiggs1gl[2, 2]), 
  dup109 -> (2*I)*CB2*(-2*MBgl2*dZHiggs1gl[1, 3] + 
      MBgl*(4*U2s28*dBn1gl[zM] + dup60*XbglC*UCSbgl[1, 3])) + 
    MBgl*((dup54*S2B - 4*dCB1*SB*YbglC)*UCSbgl[1, 3] + dup92*UCSbglC[1, 3]), 
  dup110 -> (-2*I)*CB2*(2*MBgl2*dZHiggs1gl[1, 3] + 
      MBgl*(4*U2s28*dBn1gl[zM] + dup60*XbglC*UCSbgl[1, 3])) - 
    MBgl*((dup54*S2B - 4*dCB1*SB*YbglC)*UCSbgl[1, 3] + dup92*UCSbglC[1, 3]), 
  se["A0A0"] -> 
   (-3*Alfa1L*(CB2*MTgl2*(2*(A0q14 + A0q15)*CB*(4*CB2*dSB1 + 
          CB*(-(S2B*dZHiggs1gl[2, 2]) + SB2*(4*dBn1gl[zM] - 2*dZHiggs1gl[3, 
                4]))) - 2*A0q13*(8*CB*CB2*dSB1 + 
          S2B*(-2*CB2*dZHiggs1gl[2, 2] + S2B*(2*dBn1gl[zM] - 
              dZHiggs1gl[3, 4])))) + A0q16*MBgl2*SB2*
       (8*dCB1*SB*SB2 - 2*S2B^2*dBn1gl[zM] + S2B*(-2*SB2*dZHiggs1gl[2, 2] + 
          S2B*dZHiggs1gl[3, 4])) + 
      SB*((2*A0q12 - A0q17)*MBgl2*(-8*dCB1*SB2^2 + 
          SB*(dup63*S2B + 2*S2B^2*dBn1gl[zM])) + MBgl*SB*SB2*
         ((2*I)*(B0q18*dup109 - B0q17*dup110)*U2s26 - 2*B0q14*MBgl*
           (Conjugate[U2c25]*(8*CB2*U2c26*dBn1gl[zM] - dup65*UCSbgl[3, 3] - 
              dup64*UCSbgl[3, 4]) + U2c25*(-(dup92*UCSbglC[3, 3]) + 
              (dup52*YbglC - 2*MUE*S2B*Conjugate[dTB1])*UCSbglC[3, 4] + 
              CB2*(8*Conjugate[U2c26]*dBn1gl[zM] + dup62*XbglC*UCSbglC[3, 
                  4]))))) - (2*I)*CB2^2*MTgl*((B0q21*dup116 + B0q20*dup117)*
         SB2*U2s25 + I*B0q15*MTgl*(Conjugate[U2c27]*(8*SB2*U2c28*dBn1gl[zM] + 
            dup89*UCStgl[3, 3] - dup94*UCStgl[3, 4]) + 
          U2c27*(8*SB2*Conjugate[U2c28]*dBn1gl[zM] + dup94*UCStglC[3, 3] + 
            (2*dup66 - dup51*YtglC + (2*I)*SB2*XtglC*dZHiggs1gl[1, 3])*
             UCStglC[3, 4])))))/(4*MW2*Pi*S2B^3*SW2), 
  DebugLine[1, se["A0A0"]], B0q22 -> B0q[0, MBgl2, MTgl2, Q2], 
  B0q23 -> B0q[0, MStgl2[1], MSbgl2[1], Q2], 
  B0q24 -> B0q[0, MStgl2[1], MSbgl2[2], Q2], 
  B0q25 -> B0q[0, MSbgl2[1], MStgl2[2], Q2], 
  B0q26 -> B0q[0, MStgl2[2], MSbgl2[2], Q2], 
  dup75 -> dup72*MBgl2 + MTgl2*(8*CB*CB2*dSB1 - dup71*S2B), 
  dup76 -> MTgl2*(8*CB*CB2*dSB1 - dup68*S2B) + 
    MBgl2*(-(dup67*S2B) + 8*dCB1*SB*SB2), 
  dup90 -> 16*dCB1*SB*SB2 - S2B*(dup70*S2B + 4*SB2*dZHiggs1gl[2, 2]), 
  tmp2 -> B0q23*((MBgl*USbglC[1, 2]*((-(dup69*S2B) + dup39*SB2)*UStgl[1, 1] + 
         2*dup77*MTgl*UStgl[1, 2]) + USbglC[1, 1]*(dup76*UStgl[1, 1] + 
         MTgl*(-(dup95*S2B) + 8*CB*CB2*dSB1*Ytgl)*UStgl[1, 2]) + 
       2*S2B^2*dBc1gl[zM]*((MTgl2*USbglC[1, 1] - MBgl*(MBgl*USbglC[1, 1] + 
             XbglC*USbglC[1, 2]))*UStgl[1, 1] + MTgl*Xtgl*USbglC[1, 1]*
          UStgl[1, 2]))*(((CB2*MTgl2 + MBgl2*SB2)*USbgl[1, 1] + 
         MBgl*SB2*Ybgl*USbgl[1, 2])*UStglC[1, 1] + 
       MTgl*(CB2*YtglC*USbgl[1, 1] + MBgl*USbgl[1, 2])*UStglC[1, 2]) + 
     (((CB2*MTgl2 + MBgl2*SB2)*USbglC[1, 1] + MBgl*SB2*YbglC*USbglC[1, 2])*
        UStgl[1, 1] + MTgl*(CB2*Ytgl*USbglC[1, 1] + MBgl*USbglC[1, 2])*
        UStgl[1, 2])*((dup75*USbgl[1, 1] + dup73*MBgl*USbgl[1, 2])*
        UStglC[1, 1] + MTgl*((8*CB*CB2*dSB1*YtglC - S2B^2*XtglC*
            dZHiggs1gl[6, 5])*USbgl[1, 1] + 8*dCB1*MBgl*SB*SB2*USbgl[1, 2])*
        UStglC[1, 2] - 2*(MBgl*USbgl[1, 2]*(MUEC*S2B*SB2*dZHiggs1gl[6, 5]*
            UStglC[1, 1] - 4*CB*CB2*dSB1*MTgl*UStglC[1, 2]) + 
         S2B*(2*CB2*MTgl*Conjugate[dAt]*USbgl[1, 1]*UStglC[1, 2] + 
           MBgl*USbgl[1, 2]*(2*dTB1*MUEC*SB2*UStglC[1, 1] + 
             MTgl*dZHiggs1gl[2, 2]*UStglC[1, 2]) + S2B*dBn1gl[zM]*
            (((MBgl2 - MTgl2)*USbgl[1, 1] + MBgl*Xbgl*USbgl[1, 2])*
              UStglC[1, 1] - MTgl*XtglC*USbgl[1, 1]*UStglC[1, 2]) + 
           dZHiggs1gl[2, 2]*(SB2*(MBgl2*USbgl[1, 1] + MBgl*Ybgl*USbgl[1, 2])*
              UStglC[1, 1] + CB2*MTgl*YtglC*USbgl[1, 1]*UStglC[1, 2]))))), 
  RuleAdd[tmp2, 
   B0q24*((MBgl*USbglC[2, 2]*((-(dup69*S2B) + dup39*SB2)*UStgl[1, 1] + 
         2*dup77*MTgl*UStgl[1, 2]) + USbglC[2, 1]*(dup76*UStgl[1, 1] + 
         MTgl*(-(dup95*S2B) + 8*CB*CB2*dSB1*Ytgl)*UStgl[1, 2]) + 
       2*S2B^2*dBc1gl[zM]*((MTgl2*USbglC[2, 1] - MBgl*(MBgl*USbglC[2, 1] + 
             XbglC*USbglC[2, 2]))*UStgl[1, 1] + MTgl*Xtgl*USbglC[2, 1]*
          UStgl[1, 2]))*(((CB2*MTgl2 + MBgl2*SB2)*USbgl[2, 1] + 
         MBgl*SB2*Ybgl*USbgl[2, 2])*UStglC[1, 1] + 
       MTgl*(CB2*YtglC*USbgl[2, 1] + MBgl*USbgl[2, 2])*UStglC[1, 2]) + 
     (((CB2*MTgl2 + MBgl2*SB2)*USbglC[2, 1] + MBgl*SB2*YbglC*USbglC[2, 2])*
        UStgl[1, 1] + MTgl*(CB2*Ytgl*USbglC[2, 1] + MBgl*USbglC[2, 2])*
        UStgl[1, 2])*((dup75*USbgl[2, 1] + dup73*MBgl*USbgl[2, 2])*
        UStglC[1, 1] + MTgl*((8*CB*CB2*dSB1*YtglC - S2B^2*XtglC*
            dZHiggs1gl[6, 5])*USbgl[2, 1] + 8*dCB1*MBgl*SB*SB2*USbgl[2, 2])*
        UStglC[1, 2] - 2*(MBgl*USbgl[2, 2]*(MUEC*S2B*SB2*dZHiggs1gl[6, 5]*
            UStglC[1, 1] - 4*CB*CB2*dSB1*MTgl*UStglC[1, 2]) + 
         S2B*(2*CB2*MTgl*Conjugate[dAt]*USbgl[2, 1]*UStglC[1, 2] + 
           MBgl*USbgl[2, 2]*(2*dTB1*MUEC*SB2*UStglC[1, 1] + 
             MTgl*dZHiggs1gl[2, 2]*UStglC[1, 2]) + S2B*dBn1gl[zM]*
            (((MBgl2 - MTgl2)*USbgl[2, 1] + MBgl*Xbgl*USbgl[2, 2])*
              UStglC[1, 1] - MTgl*XtglC*USbgl[2, 1]*UStglC[1, 2]) + 
           dZHiggs1gl[2, 2]*(SB2*(MBgl2*USbgl[2, 1] + MBgl*Ybgl*USbgl[2, 2])*
              UStglC[1, 1] + CB2*MTgl*YtglC*USbgl[2, 1]*UStglC[1, 2])))))], 
  RuleAdd[tmp2, 
   B0q25*((MBgl*USbglC[1, 2]*((-(dup69*S2B) + dup39*SB2)*UStgl[2, 1] + 
         2*dup77*MTgl*UStgl[2, 2]) + USbglC[1, 1]*(dup76*UStgl[2, 1] + 
         MTgl*(-(dup95*S2B) + 8*CB*CB2*dSB1*Ytgl)*UStgl[2, 2]) + 
       2*S2B^2*dBc1gl[zM]*((MTgl2*USbglC[1, 1] - MBgl*(MBgl*USbglC[1, 1] + 
             XbglC*USbglC[1, 2]))*UStgl[2, 1] + MTgl*Xtgl*USbglC[1, 1]*
          UStgl[2, 2]))*(((CB2*MTgl2 + MBgl2*SB2)*USbgl[1, 1] + 
         MBgl*SB2*Ybgl*USbgl[1, 2])*UStglC[2, 1] + 
       MTgl*(CB2*YtglC*USbgl[1, 1] + MBgl*USbgl[1, 2])*UStglC[2, 2]) + 
     (((CB2*MTgl2 + MBgl2*SB2)*USbglC[1, 1] + MBgl*SB2*YbglC*USbglC[1, 2])*
        UStgl[2, 1] + MTgl*(CB2*Ytgl*USbglC[1, 1] + MBgl*USbglC[1, 2])*
        UStgl[2, 2])*((dup75*USbgl[1, 1] + dup73*MBgl*USbgl[1, 2])*
        UStglC[2, 1] + MTgl*((8*CB*CB2*dSB1*YtglC - S2B^2*XtglC*
            dZHiggs1gl[6, 5])*USbgl[1, 1] + 8*dCB1*MBgl*SB*SB2*USbgl[1, 2])*
        UStglC[2, 2] - 2*(MBgl*USbgl[1, 2]*(MUEC*S2B*SB2*dZHiggs1gl[6, 5]*
            UStglC[2, 1] - 4*CB*CB2*dSB1*MTgl*UStglC[2, 2]) + 
         S2B*(2*CB2*MTgl*Conjugate[dAt]*USbgl[1, 1]*UStglC[2, 2] + 
           MBgl*USbgl[1, 2]*(2*dTB1*MUEC*SB2*UStglC[2, 1] + 
             MTgl*dZHiggs1gl[2, 2]*UStglC[2, 2]) + S2B*dBn1gl[zM]*
            (((MBgl2 - MTgl2)*USbgl[1, 1] + MBgl*Xbgl*USbgl[1, 2])*
              UStglC[2, 1] - MTgl*XtglC*USbgl[1, 1]*UStglC[2, 2]) + 
           dZHiggs1gl[2, 2]*(SB2*(MBgl2*USbgl[1, 1] + MBgl*Ybgl*USbgl[1, 2])*
              UStglC[2, 1] + CB2*MTgl*YtglC*USbgl[1, 1]*UStglC[2, 2])))))], 
  RuleAdd[tmp2, 
   B0q26*((MBgl*USbglC[2, 2]*((-(dup69*S2B) + dup39*SB2)*UStgl[2, 1] + 
         2*dup77*MTgl*UStgl[2, 2]) + USbglC[2, 1]*(dup76*UStgl[2, 1] + 
         MTgl*(-(dup95*S2B) + 8*CB*CB2*dSB1*Ytgl)*UStgl[2, 2]) + 
       2*S2B^2*dBc1gl[zM]*((MTgl2*USbglC[2, 1] - MBgl*(MBgl*USbglC[2, 1] + 
             XbglC*USbglC[2, 2]))*UStgl[2, 1] + MTgl*Xtgl*USbglC[2, 1]*
          UStgl[2, 2]))*(((CB2*MTgl2 + MBgl2*SB2)*USbgl[2, 1] + 
         MBgl*SB2*Ybgl*USbgl[2, 2])*UStglC[2, 1] + 
       MTgl*(CB2*YtglC*USbgl[2, 1] + MBgl*USbgl[2, 2])*UStglC[2, 2]) + 
     (((CB2*MTgl2 + MBgl2*SB2)*USbglC[2, 1] + MBgl*SB2*YbglC*USbglC[2, 2])*
        UStgl[2, 1] + MTgl*(CB2*Ytgl*USbglC[2, 1] + MBgl*USbglC[2, 2])*
        UStgl[2, 2])*((dup75*USbgl[2, 1] + dup73*MBgl*USbgl[2, 2])*
        UStglC[2, 1] + MTgl*((8*CB*CB2*dSB1*YtglC - S2B^2*XtglC*
            dZHiggs1gl[6, 5])*USbgl[2, 1] + 8*dCB1*MBgl*SB*SB2*USbgl[2, 2])*
        UStglC[2, 2] - 2*(MBgl*USbgl[2, 2]*(MUEC*S2B*SB2*dZHiggs1gl[6, 5]*
            UStglC[2, 1] - 4*CB*CB2*dSB1*MTgl*UStglC[2, 2]) + 
         S2B*(2*CB2*MTgl*Conjugate[dAt]*USbgl[2, 1]*UStglC[2, 2] + 
           MBgl*USbgl[2, 2]*(2*dTB1*MUEC*SB2*UStglC[2, 1] + 
             MTgl*dZHiggs1gl[2, 2]*UStglC[2, 2]) + S2B*dBn1gl[zM]*
            (((MBgl2 - MTgl2)*USbgl[2, 1] + MBgl*Xbgl*USbgl[2, 2])*
              UStglC[2, 1] - MTgl*XtglC*USbgl[2, 1]*UStglC[2, 2]) + 
           dZHiggs1gl[2, 2]*(SB2*(MBgl2*USbgl[2, 1] + MBgl*Ybgl*USbgl[2, 2])*
              UStglC[2, 1] + CB2*MTgl*YtglC*USbgl[2, 1]*UStglC[2, 2])))))], 
  se["HmHp"] -> 
   (3*Alfa1L*(-tmp2 + A0q12*(4*CB2^2*MTgl2*(8*CB*dSB1 - 2*dup80 + 
          4*SB2*(dBc1gl[zM] + dBn1gl[zM])) - 8*MBgl2*SB2^2*
         (CB2*dup70 - 4*dCB1*SB + S2B*dZHiggs1gl[2, 2])) + 
      2*B0q22*(CB2*MTgl2^2*(4*CB2*dup51 + dup70*S2B^2) + 
        MBgl2*MTgl2*S2B^2*(4*CB*dSB1 - CB2*dup70 + 4*dCB1*SB*(1 + TB2) - 
          S2B*(2 + TB2)*dZHiggs1gl[2, 2])) - 4*A0q17*CB2^2*MTgl2*
       (4*CB*dSB1 + dup70*SB2 - S2B*dZHiggs1gl[2, 2])*USbgl2[1, 1] + 
      4*A0q16*MBgl2*SB2^2*(CB2*dup70 - 4*dCB1*SB + S2B*dZHiggs1gl[2, 2])*
       USbgl2[1, 1] - 4*A0q16*CB2^2*MTgl2*(4*CB*dSB1 + dup70*SB2 - 
        S2B*dZHiggs1gl[2, 2])*USbgl2[1, 2] + 4*A0q17*MBgl2*SB2^2*
       (CB2*dup70 - 4*dCB1*SB + S2B*dZHiggs1gl[2, 2])*USbgl2[1, 2] - 
      A0q14*(4*CB2^2*MTgl2*(4*CB*dSB1 + dup70*SB2 - S2B*dZHiggs1gl[2, 2])*
         UStgl2[1, 1] + dup90*MBgl2*SB2*UStgl2[1, 2]) - 
      A0q15*(dup90*MBgl2*SB2*UStgl2[1, 1] + 4*CB2^2*MTgl2*
         (4*CB*dSB1 + dup70*SB2 - S2B*dZHiggs1gl[2, 2])*UStgl2[1, 2])))/
    (8*MW2*Pi*S2B^3*SW2), DebugLine[1, se["HmHp"]], 
  dup44 -> 8*CB*dCB1 - 4*CB2*dZHiggs1gl[1, 1] - 2*S2B*dZHiggs1gl[1, 2], 
  dup58 -> -2*S2B*Conjugate[dAt] + YtglC*(4*CB*dSB1 - S2B*dZHiggs1gl[2, 2]) + 
    SB2*XtglC*(2*dZHiggs1gl[1, 2] - (2*I)*dZHiggs1gl[2, 4]), 
  dup82 -> 4*(MBgl - U2s14)*dA1gl[zM] + dup48*YbglC*UCSbgl[1, 3] + 
    dup50*Ybgl*UCSbglC[1, 3], 
  dup99 -> 4*(S2B*(MBgl + U2s14)*dA1gl[zM] + CB*dCB1*XbglC*UCSbgl[1, 3]) + 
    2*dup40*Xbgl*UCSbglC[1, 3] - 2*CB2*(dup42*UCSbgl[1, 3] + 
      2*dTB1*MUEC*UCSbglC[1, 3]) - S2B*(dup48*YbglC*UCSbgl[1, 3] + 
      dup50*Ybgl*UCSbglC[1, 3]), 
  dup105 -> dup44*MBgl2 + MBgl*(dup82*S2B + (2*CB2*dup42 - 4*CB*dCB1*XbglC)*
       UCSbgl[1, 3] + (4*CB2*dTB1*MUEC - 2*dup40*Xbgl)*UCSbglC[1, 3]), 
  dup107 -> 2*dup53*MTgl2 + MTgl*(8*SB2*(MTgl - U2s11)*dA1gl[zM] - 
      dup58*UCStgl[1, 3] - dup93*UCStglC[1, 3]), 
  dup111 -> 2*dup53*MTgl2 + MTgl*(8*SB2*(MTgl + U2s11)*dA1gl[zM] + 
      (dup57*SB2*XtglC + dup51*YtglC - 2*S2B*Conjugate[dAt])*UCStgl[1, 3] + 
      dup93*UCStglC[1, 3]), 
  dup112 -> 8*CB*MTgl*SB2*(MTgl - U2s13)*dA1gl[zM] + 
    SB*(-(dup41*MTgl2) - dup102*S2B + MTgl*(2*dup86*UCStgl[1, 3] + 
        dup84*UCStglC[1, 3])), 
  dup114 -> 8*CB*MTgl*SB2*(MTgl + U2s13)*dA1gl[zM] - 
    SB*(dup41*MTgl2 + dup103*S2B + MTgl*(2*dup86*UCStgl[1, 3] + 
        dup84*UCStglC[1, 3])), se["h0HH"] -> 
   (-3*Alfa1L*
     ((-2*(MBgl*(-(B0q18*(CB*(dup98*MBgl + 2*dup79*MBgl2)*(MBgl - U2s12) + 
              dup105*SB*(MBgl - U2s14))) + 
           B0q17*(CB*(dup98*MBgl - 2*dup79*MBgl2)*(MBgl + U2s12) - 
             (dup99*MBgl + dup44*MBgl2)*SB*(MBgl + U2s14))) + 
         CB*MBgl2*((A0q16 + A0q17)*(-8*dCB1*SB + (4*CB2 - 4*SB2)*dA1gl[zM] + 
             2*dZHiggs1gl[1, 2] + S2B*(dZHiggs1gl[1, 1] + dZHiggs1gl[2, 
                2])) - 2*(A0q12 + 2*B0q16*MBgl2)*(-8*dCB1*SB + 
             2*CB2*(2*dA1gl[zM] + dZHiggs1gl[1, 2]) + 
             SB2*(-4*dA1gl[zM] + 2*dZHiggs1gl[1, 2]) + 
             S2B*(dZHiggs1gl[1, 1] + dZHiggs1gl[2, 2]))) + 
         B0q14*MBgl^2*(SB*(dup101*U2c12 - dup97*Conjugate[U2c12]) + 
           CB*(Conjugate[U2c11]*(16*CB2*U2c11*dA1gl[zM] + (2*dup46 - 
                 dup52*YbglC + (2*I)*CB2*XbglC*dZHiggs1gl[2, 4])*
                UCSbgl[3, 3] + (2*dup45 - dup52*Ybgl - (2*I)*CB2*Xbgl*
                  dZHiggs1gl[2, 4])*UCSbgl[3, 4]) + 
             U2c11*((2*dup45 - dup52*Ybgl - (2*I)*CB2*Xbgl*dZHiggs1gl[2, 4])*
                UCSbglC[3, 3] + (2*dup46 - dup52*YbglC + (2*I)*CB2*XbglC*
                  dZHiggs1gl[2, 4])*UCSbglC[3, 4])))))/(CB*CB2) + 
      (2*MTgl2*(2*(A0q13 + 2*B0q19*MTgl2)*SB2*
           (4*(dSB1*S2B + SB*SB2*dA1gl[zM]) + SB2*(-2*CB*dZHiggs1gl[1, 1] + 
              2*SB*dZHiggs1gl[1, 2]) + SB*(2*CB2*(-2*dA1gl[zM] + 
                dZHiggs1gl[1, 2]) - S2B*dZHiggs1gl[2, 2])) + 
          (A0q14 + A0q15)*(SB*(S2B^2 - 4*SB2^2)*dA1gl[zM] + 
            SB2*(-2*SB*dZHiggs1gl[1, 2] + S2B*(-4*dSB1 + 
                SB*(dZHiggs1gl[1, 1] + dZHiggs1gl[2, 2]))))) - 
        2*MTgl*SB*(-(CB*(B0q20*dup112*(MTgl - U2s13) + B0q21*dup114*
              (MTgl + U2s13) + B0q15*MTgl*(Conjugate[U2c14]*
                (16*CB*SB2*U2c14*dA1gl[zM] - SB*((-(dup43*SB2) + 4*dSB1*SB*
                      XtglC + dup50*S2B*YtglC)*UCStgl[3, 3] + dup85*
                    UCStgl[3, 4])) - SB*U2c14*(dup85*UCStglC[3, 3] + 
                 (-(dup43*SB2) + 4*dSB1*SB*XtglC + dup50*S2B*YtglC)*
                  UCStglC[3, 4])))) + SB2*(B0q20*dup107*(MTgl - U2s11) + 
            B0q21*dup111*(MTgl + U2s11) + B0q15*MTgl*(Conjugate[U2c13]*(
                (dup51*YtglC - 2*S2B*Conjugate[dAt])*UCStgl[3, 3] + 
                SB2*(16*U2c13*dA1gl[zM] + dup57*XtglC*UCStgl[3, 3]) + 
                dup93*UCStgl[3, 4]) + U2c13*(dup93*UCStglC[3, 3] + 
                (dup57*SB2*XtglC + dup51*YtglC - 2*S2B*Conjugate[dAt])*
                 UCStglC[3, 4])))))/(SB*SB2^2)))/(64*MW2*Pi*SW2), 
  DebugLine[1, se["h0HH"]], se["h0A0"] -> 
   (3*Alfa1L*(CB*MTgl2*S2B^2*(SB*(A0q14 + A0q15 - 4*B0q19*MTgl2*SB2) - 
        A0q13*(CB*S2B + 2*SB*SB2))*dZHiggs1gl[1, 3] - 
      2*SB*SB2^2*(2*CB*MBgl2*(2*A0q12 - A0q16 - A0q17 + 4*B0q16*CB2*MBgl2)*
         dZHiggs1gl[1, 3] - I*MBgl*(B0q17*(CB*dup110*(MBgl + U2s12) + 
            I*(dup99*MBgl + dup44*MBgl2)*SB*U2s26) - 
          B0q18*(I*dup105*SB*U2s26 + CB*(MBgl - U2s12)*
             ((4*I)*CB2*MBgl2*dZHiggs1gl[1, 3] - MBgl*(dup54*S2B*UCSbgl[1, 
                  3] - 4*dCB1*SB*YbglC*UCSbgl[1, 3] + (2*I)*CB2*
                 (4*U2s28*dBn1gl[zM] + dup60*XbglC*UCSbgl[1, 3]) + 
                dup92*UCSbglC[1, 3]))) - B0q14*MBgl*
           (SB*(dup101*U2c25 + dup97*Conjugate[U2c25]) - 8*CB*CB2*U2c11*
             Conjugate[U2c26]*dBn1gl[zM] + CB*Conjugate[U2c11]*
             (8*CB2*U2c26*dBn1gl[zM] - dup65*UCSbgl[3, 3] - dup64*UCSbgl[3, 
                4]) + CB*dup92*U2c11*UCSbglC[3, 3] - CB*U2c11*
             (CB2*dup62*XbglC + dup52*YbglC - 2*MUE*S2B*Conjugate[dTB1])*
             UCSbglC[3, 4]))) - (I/2)*CB2*MTgl*S2B*
       (B0q20*((-2*I)*CB*dup112*U2s25 + 2*SB2*(MTgl - U2s11)*
           (SB2*((8*I)*MTgl*U2s27*dBn1gl[zM] + (4*I)*MTgl2*dZHiggs1gl[1, 
                3]) + MTgl*(dup89*UCStgl[1, 3] - dup94*UCStglC[1, 3]))) + 
        B0q21*((2*I)*CB*dup114*U2s25 + 2*SB2*(MTgl + U2s11)*
           (SB2*((-8*I)*MTgl*U2s27*dBn1gl[zM] + (4*I)*MTgl2*dZHiggs1gl[1, 
                3]) + MTgl*(-(dup89*UCStgl[1, 3]) + dup94*UCStglC[1, 3]))) - 
        B0q15*MTgl*(CB*(-2*dup108*U2c27 + 2*dup106*Conjugate[U2c27]) + 
          SB2*(2*Conjugate[U2c13]*(8*SB2*U2c28*dBn1gl[zM] + dup89*UCStgl[3, 
                3] - dup94*UCStgl[3, 4]) - 2*U2c13*(8*SB2*Conjugate[U2c28]*
               dBn1gl[zM] + dup94*UCStglC[3, 3] + (2*dup66 - dup51*YtglC + 
                (2*I)*SB2*XtglC*dZHiggs1gl[1, 3])*UCStglC[3, 4]))))))/
    (32*CB2*MW2*Pi*S2B*SB2^2*SW2), DebugLine[1, se["h0A0"]], 
  se["HHA0"] -> 
   (-3*Alfa1L*((8*(MBgl2*TB2*(2*B0q16*MBgl2*dZHiggs1gl[1, 3] + 
           A0q12*(dZHiggs1gl[1, 3] - dZHiggs1gl[2, 4])) + 
         A0q13*MTgl2*(-dZHiggs1gl[1, 3] + dZHiggs1gl[2, 4])))/TB - 
      (SB*(4*(A0q16 + A0q17)*CB2*MBgl2*(dZHiggs1gl[1, 3] - 
           dZHiggs1gl[2, 4]) + (2*I)*MBgl*(B0q17*(dup110*(MBgl + U2s14) - 
             I*(dup98*MBgl - 2*dup79*MBgl2)*U2s26) + 
           B0q18*(dup109*(MBgl - U2s14) - I*(dup98*MBgl + 2*dup79*MBgl2)*
              U2s26) + B0q14*MBgl*(-(dup100*U2c25) + dup96*Conjugate[U2c25] + 
             8*CB2*U2c12*Conjugate[U2c26]*dBn1gl[zM] - Conjugate[U2c12]*
              (8*CB2*U2c26*dBn1gl[zM] - dup65*UCSbgl[3, 3] - dup64*
                UCSbgl[3, 4]) - dup92*U2c12*UCSbglC[3, 3] + 
             U2c12*(CB2*dup62*XbglC + dup52*YbglC - 2*MUE*S2B*Conjugate[
                 dTB1])*UCSbglC[3, 4]))))/(CB*CB2) + 
      (CB*(4*MTgl2*SB2*(-4*B0q19*MTgl2*dZHiggs1gl[1, 3] + 
           (A0q14 + A0q15)*(dZHiggs1gl[1, 3] - dZHiggs1gl[2, 4])) - 
         (2*I)*MTgl*((I*B0q20*dup107 - I*B0q21*dup111)*U2s25 + 
           B0q20*(MTgl - U2s13)*(SB2*((8*I)*MTgl*U2s27*dBn1gl[zM] + (4*I)*
                MTgl2*dZHiggs1gl[1, 3]) + MTgl*(dup89*UCStgl[1, 3] - dup94*
                UCStglC[1, 3])) + B0q21*(MTgl + U2s13)*
            (SB2*((-8*I)*MTgl*U2s27*dBn1gl[zM] + (4*I)*MTgl2*dZHiggs1gl[1, 
                 3]) + MTgl*(-(dup89*UCStgl[1, 3]) + dup94*UCStglC[1, 3])) + 
           B0q15*MTgl*(Conjugate[U2c27]*(8*SB2*U2c13*dA1gl[zM] + dup57*SB2*
                XtglC*UCStgl[3, 3] + dup51*YtglC*UCStgl[3, 3] - 2*S2B*
                Conjugate[dAt]*UCStgl[3, 3] + dup93*UCStgl[3, 4]) - 
             Conjugate[U2c14]*(8*SB2*U2c28*dBn1gl[zM] + dup89*UCStgl[3, 3] - 
               dup94*UCStgl[3, 4]) - U2c27*(8*SB2*Conjugate[U2c13]*
                dA1gl[zM] + dup93*UCStglC[3, 3] + dup57*SB2*XtglC*
                UCStglC[3, 4] + dup51*YtglC*UCStglC[3, 4] - 2*S2B*
                Conjugate[dAt]*UCStglC[3, 4]) + U2c14*(8*SB2*Conjugate[U2c28]*
                dBn1gl[zM] + dup94*UCStglC[3, 3] + (2*dup66 - dup51*YtglC + 
                 (2*I)*SB2*XtglC*dZHiggs1gl[1, 3])*UCStglC[3, 4])))))/
       (SB*SB2)))/(64*MW2*Pi*SW2), DebugLine[1, se["HHA0"]], 
  seshift["h0h0"] -> -dMHiggsZ2gl[1, 1] + se["h0h0"], 
  DebugLine[1, seshift["h0h0"]], seshift["HHHH"] -> 
   -dMHiggsZ2gl[2, 2] + se["HHHH"], DebugLine[1, seshift["HHHH"]], 
  seshift["A0A0"] -> -dMHiggsZ2gl[3, 3] + se["A0A0"], 
  DebugLine[1, seshift["A0A0"]], seshift["HmHp"] -> 
   -dMHiggsZ2gl[5, 5] + se["HmHp"], DebugLine[1, seshift["HmHp"]], 
  seshift["h0HH"] -> -dMHiggsZ2gl[1, 2] + se["h0HH"], 
  DebugLine[1, seshift["h0HH"]], seshift["h0A0"] -> 
   -dMHiggsZ2gl[1, 3] + se["h0A0"], DebugLine[1, seshift["h0A0"]], 
  seshift["HHA0"] -> -dMHiggsZ2gl[2, 3] + se["HHA0"], 
  DebugLine[1, seshift["HHA0"]], "\tend\n"}]
