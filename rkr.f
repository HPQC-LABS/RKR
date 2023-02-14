c=======================================================================
c|..optional...Register...optional...Register...optional...Register...|c
c|--------------------------------------------------------------------|c
c| You have choosen to download the attached source code for my       |c
c| Fortran program RKR1. I would appreciate it if you would please go |c
c| to the www address  http://scienide2.uwaterloo.ca/~rleroy/RKR16/   |c
c| fill in the registration form there if you wish to be accessible   |c
c| so that I can send you possible future updates and/or corrections  |c
c| for this code.  This address list will be held securely by me and  |c
c| used for no other purpose................. Robert J. Le Roy .......|c
c|..Register...optional....Register...optional...Register...optional..|c
c=======================================================================

c**********************************************************************
c** R.J. Le Roy's program 'RKR1' for calculating RKR turning points in
c  either simple first order or first-order Kaiser approximation.
c** 'G(v)' and 'B(v)' may be generated using conventional (Dunham) 
c  polynomials in (v+1/2), using Near-Dissociation expansions (NDE's), 
c  or using Tellinguisen's MXS 'mixed' Dunha,-at-low-v and NDE-at-high-v
c  functional representations.
c** If desired (see manual), this program will also extrapolate up the 
c  inner wall past a chosen point with an exponential fitted to the 3 
c  preceeding points, & adjust 'RT(outer)' accordingly.
c******************** Version of  8 January 2016  **********************
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 2003-2016  by  Robert J. Le Roy             +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the explicit written permission of the author.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER MXTP,MXLEV,MXDUN
      PARAMETER (MXLEV=500,MXTP=1001,MXDUN=25)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c** NOTE: Current dimensioning allows calculation of pairs of turning 
c  points for up to MXLEV 500 vibrational levels with a total of 
c  MXTP=2*MXLEV+1 turning points, and for Dunham and NDE 
c  polynomial expansions of defined by up to MXDUN= 25 parameters. 
c  These limits may may be changed by the user.
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CHARACTER*78 TITLE
      CHARACTER*6 PTYPE(6)
      CHARACTER*2 NAME1,NAME2,CCDC(0:1)
      INTEGER IAN1,IAN2,IMN1,IMN2,GEL1,GEL2,GNS1,GNS2,CHARGE,
     1 KORDR,MORSE,NGP,Kaiser,IBIS,NBIS,IDIV,NDIV,NSV,NTP,
     2 IVIB,NVIB,IVB,I634,I638,I,J,K,L,NN,I1,I2,I1P1,I1P2,I640,I644,
     3 IR1,IDR
c
      REAL*8 CU,SQCU,MASS1,MASS2,ZMASE,ABUND1,ABUND2,TOLER,ZMU,Req,Be,
     1 AE,De,WE,WEXE,BETA,V00,DV0,Sw,SwLR,VEXT,F,FF,FB,G,GG,GB,VST,VFN,
     2 RANGE,TSTF,TSTG,TSTFB,TSTGB,VUP,GUP,BUP,DGUP,RMIN,RMAX,RMIN2,
     3 RMIN3,VRAT,CEXT,DCEXT,ADCEXT,BDCEXT,EXP1,EXP2,EXP3,FUN,
     4 DFUN,AEXT,BEXT,REXT,DR1,DEI,DEIB,DRI,DRIB,D1,D2,TT,HEL,
     5 V1(9),V2(9),DV(9),V(MXLEV),RT(-4:MXTP),
     4 ET(-4:MXTP),VX(MXDUN),BV(MXDUN),GV(MXDUN),DGDV(MXDUN),WW(16)
c** Common block for Dunham & MXS function parameters
      INTEGER LMAXGv,LMAXBv,NDEGv,NDEBv
      REAL*8 VS,DVS,Y00,YL0(0:MXDUN),YL1(0:MXDUN)
      COMMON /DUNPRM/Y00,VS,DVS,YL0,YL1,LMAXGv,LMAXBv,NDEGv,NDEBv
c** Common block for NDE function parameters
      INTEGER NLR,ITYPE,ITYPB,IZP0,IZQ0,IZP1,IZQ1,NP0,NQ0,NP1,NQ1
      REAL*8 VD,DLIM,XCN0,XCN1,P0(MXDUN),Q0(MXDUN),P1(MXDUN),Q1(MXDUN)
      COMMON /NDEPRM/VD,DLIM,XCN0,XCN1,P0,Q0,P1,Q1,NLR,ITYPE,ITYPB,
     1  IZP0,IZQ0,IZP1,IZQ1,NP0,NQ0,NP1,NQ1
c** Common block for quadrature weights & points
      REAL*8 XG(32),WG(32),X2(16),W2(16)
      COMMON /GWGHT/XG,WG,X2,W2
c
      DATA CCDC/'Gv','Bv'/
      DATA PTYPE/' OUTER',' INNER','Expone',' Pade ',' Pade ','ntial '/
      DATA ZMASE /5.4857990946D-04/  !! 2010 physical constants d:mohr12
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Quadratures performed using NGP point Gaussian integrations.
      NGP= 16
      CALL WGHT(NGP)
c** Quadrature convergence criterion is  TOLER
      TOLER= 1.d-10
c** Use up to NBIS interval bisections in the NGP-point gaussian
c  integration for the "f" and "g" integrals (typically set  NBIS=3-5).
      NBIS= 5
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Begin by reading in the (integer) atomic numbers and mass numbers
c  defining the effective reduced mass of the system considered.
c** IAN1 & IAM2, and IMN1 & IMN2 are, respectively, the atomic numbers
c    and the mass numbers identifying the atoms forming the molecule.
c    Their masses are extracted from data subroutine MASSES and used
c    to generate the the reduced mass ZMU.
c** If  IMN1  or  IMN2  lie outside the range of mass numbers for normal
c  stable isotopes of that species, subroutine MASSES returns the
c  average atomic mass based on the natural isotope abundance.
c** If the read-in value of IAN1 and/or IAN2 is .LE.0, then instead of
c  using the MASS table, read an actual particle mass for it/them.
c** CHARGE (integer) is the charge on the molecule (=0 for neutral). If
c   (CHARGE.ne.0)  generate & use Watson's charge-adjusted reduced mass.
c** NDEGv(s) & NDEBv(s) specify whether Gv & Bv are represented by:
c      (a)  pure Dunham expansions    when    NDEXv(s) = 0
c      (b)  pure NDE expressions      when    NDEXv(s) = 1
c      (c)  MXS mixed NDE/Dunham expressions when    NDEXv(s) > 1
c!! SPECIAL CASE: if have NO v-dependent rotational data and wish to 
c   define inner potential wall by Morse function, set:  NDEBv =-1 !!
c!! NOTE: Program does NOT allow cases with  NDEBv > NDEGv
c----------------------------------------------------------------------
    2 READ(5,*,END=99) IAN1, IMN1, IAN2, IMN2, CHARGE, NDEGv, NDEBv
c----------------------------------------------------------------------
      Y00= 0.d0
      V00= -0.5d0
      MORSE= 0
      IF(NDEBv.GT.NDEGv) THEN
          WRITE(6,598) NDEGv, NDEBv
          STOP
          ENDIF
      I640= 0
      I644= 0
c** Subroutine MASSES returns the names of the atoms NAMEi,ground
c  electronic state degeneracy GELi, nuclear spin degeneracy GNSi,
c  mass MASSi, and isotopic abundance ABUNDi for a given atomic isotope.
      IF((IAN1.GT.0).AND.(IAN1.LE.109)) THEN
          CALL MASSES(IAN1,IMN1,NAME1,GEL1,GNS1,MASS1,ABUND1)
        ELSE
c** If particle-i is not a normal atomic isotope, read a 2-character
c   NAME (enclosed between ', as in 'mu') and its actual mass.
c----------------------------------------------------------------------
          READ(5,*) NAME1, MASS1
c----------------------------------------------------------------------
        ENDIF
      IF((IAN2.GT.0).AND.(IAN2.LE.109)) THEN
          CALL MASSES(IAN2,IMN2,NAME2,GEL2,GNS2,MASS2,ABUND2)
        ELSE
c----------------------------------------------------------------------
          READ(5,*) NAME2, MASS2
c----------------------------------------------------------------------
        ENDIF
      ZMU= (MASS1*MASS2)/(MASS1+MASS2 - DBLE(CHARGE)*ZMASE)
c** Numerical factor  16.857629206 (+/- 0.000,000,013) based on Compton
c  wavelength of proton & proton mass (u) from 2011 physical constants.
      CU= 16.857629206d0/ZMU
      SQCU= DSQRT(CU)
c=======================================================================
c TITLE is a title or output header of up to 78 characters, read on a
c   single line enclosed between single quotes: e.g.  'title of problem'
c=======================================================================
      READ(5,*) TITLE
c-----------------------------------------------------------------------
      WRITE(6,600) TITLE,NAME1,IMN1,NAME2,IMN2,CHARGE,ZMU,CU,MASS1,
     1                                                           MASS2
      WRITE(6,602) TOLER,NBIS,NGP
c=======================================================================
      IF((NDEGv.LE.0).OR.(NDEGv.GE.2)) THEN
c** For Dunham or MXS vibrational energy function, read in order LMAXGv
c    and values  YL0(i) (i=1,LMAXGv) of Dunham vibrational coefficients
c-----------------------------------------------------------------------
          READ(5,*) LMAXGv
          READ(5,*) (YL0(L),L= 1,LMAXGv)
c-----------------------------------------------------------------------
          WE= YL0(1)
          WEXE= -YL0(2)
c=======================================================================
c** If using Tellinghuisen-style Mixed Representations ... read
c  VS   the v value (floating point) where Dunham switches to NDE, 
c  DVS  the width parameter on the switching function, and
c  DLIM  the well depth, or energy at the asymptote assuming an energy
c      zero at  v= -1/2
c=======================================================================
          IF(NDEGv.GE.2) READ(5,*) VS, DVS, DLIM
c-----------------------------------------------------------------------
          IF(NDEGv.GE.2) WRITE(6,604) CCDC(0),LMAXGv,VS,DVS,DLIM
          WRITE(6,606) LMAXGv,CCDC(0),(YL0(I),I= 1,LMAXGv)
          ENDIF
      IF(NDEGv.GE.1) THEN
c
c** If use NDE or MXS expression for vibrational energies, read NDE 
c  control parameters and expansion constants here.
c* For an "outer" Pade expansion, ITYPE=1 ;  for an "inner" Pade,
c  ITYPE=2 ;  if ITYPE=3, use an exponential polynomial NDE .
c* Expansion variable is  (vD-v), and the leading non-zero contribution
c   to the NP-term numerator polynomial is the power  IZP0 of (vD-v), 
c   while the corresponding leading term in the NQ0-term denominator 
c   polynomial is  (vD-v)**IZQ0
c-----------------------------------------------------------------------
          READ(5,*) NLR, ITYPE, IZP0, IZQ0, NP0, NQ0, VD, XCN0
          IF(NP0.GT.0) READ(5,*) (P0(I),I= 1,NP0)
          IF(NQ0.GT.0) READ(5,*) (Q0(I),I= 1,NQ0)
c-----------------------------------------------------------------------
          IF(NDEGv.EQ.1) DLIM= 0.d0
          IZP0= IZP0- 1
          IZQ0= IZQ0- 1
          KORDR= 0
          FF= V00 - 1.d-5
          DO  I= 1,3
              CALL NDEDKM(FF,KORDR,GV(I),DGDV(I),NLR,XCN0,DLIM,VD,
     1                                  IZP0,IZQ0,ITYPE,NP0,NQ0,P0,Q0)
              FF= FF+ 1.d-5
              ENDDO
          IF(NDEGv.EQ.1) DLIM= -GV(2)
          WRITE(6,608) CCDC(0),NP0,NQ0,(PTYPE(I),I= ITYPE,6,3),0,NLR,
     1                                      XCN0,IZP0+1,IZQ0+1,VD,DLIM
          IF(NP0.GT.0) WRITE(6,610) (P0(I),I= 1,NP0)
          IF(NQ0.GT.0) WRITE(6,612) (Q0(I),I= 1,NQ0)
          WE= DGDV(2)
          WEXE= (DGDV(1)-DGDV(3))/4.D-5
          IF(NDEGv.GE.2) THEN
              SwLR= DEXP(-(VS+ 0.5d0)/DVS)
              Sw= 1.d0/(1.d0+ SwLR)
              SwLR= SwLR/(1.d0 + SwLR)
              WE= SwLR*WE + Sw*YL0(1)
              WEXE= SwLR*WEXE - Sw*YL0(2)
              ENDIF
          ENDIF
c
      IF(NDEBv.LT.0) THEN
c** If have NO v-dependent rotational data and wish to define inner 
c  potential wall by Morse function, read  Req  value to define position
c  of potential minimum.
c-----------------------------------------------------------------------
          READ(5,*) Req
c-----------------------------------------------------------------------
c** Use vibrational constants to determine Morse parameters at minimum.
          BE= CU/Req**2
          AE= (SQRT(WEXE/BE) - 1.d0)*6.d0*BE**2/WE
          De= 0.25d0*WE**2/WEXE
          BETA= DSQRT(DABS(WEXE)/CU)
          WRITE(6,614) Req,WE,WEXE,De,BETA
          MORSE= 1
          ENDIF
      IF((NDEBv.EQ.0).OR.(NDEBv.GE.2)) THEN
c** For Dunham or MXS  Bv  function, read in order LMAXBv of Dunham 
c  expansion in (v+1/2) and coefficients YL1(i) (i=0,LMAXBv)
c-----------------------------------------------------------------------
          READ(5,*) LMAXBv
          IF(LMAXBv.GE.0) READ(5,*) (YL1(L),L= 0,LMAXBv)
c-----------------------------------------------------------------------
c** Use vibrational constants to determine Morse parameters at minimum.
          BE= YL1(0)
          AE= -YL1(1)
          IF(NDEBv.GE.2) WRITE(6,604) CCDC(1),LMAXBv,VS
          WRITE(6,606) LMAXBv+1,CCDC(1),(YL1(I),I= 0,LMAXBv)
          ENDIF
      IF(NDEBv.GE.1) THEN
c** If use NDE or MXS expression for Bv values, read NDE control 
c  parameters and expansion constants which are defined in exactly the
c  same way as those for the vibrational NDE.
c-----------------------------------------------------------------------
          READ(5,*) ITYPB, IZP1, IZQ1, NP1, NQ1, XCN1
          IF(NP1.GT.0) READ(5,*) (P1(I),I= 1,NP1)
          IF(NQ1.GT.0) READ(5,*) (Q1(I),I= 1,NQ1)
c-----------------------------------------------------------------------
          WRITE(6,608) CCDC(1),NP1,NQ1,(PTYPE(I),I= ITYPB,6,3),1,NLR,
     1                                                  XCN1,IZP1,IZQ1
          IF(NP1.GT.0) WRITE(6,610) (P1(I),I= 1,NP1)
          IF(NQ1.GT.0) WRITE(6,612) (Q1(I),I= 1,NQ1)
          IZP1= IZP1- 1
          IZQ1= IZQ1- 1
          KORDR= 1
          FF= V00
          CALL NDEDKM(FF,KORDR,BE,AE,NLR,XCN1,DLIM,VD,
     1                                  IZP1,IZQ1,ITYPB,NP1,NQ1,P1,Q1)
          AE= -AE
          IF(NDEBv.GE.2) THEN
              IF(MAX(DABS(BE/YL1(0)),DABS(AE/YL1(1))).GT.0.01d0/SwLR)
     1                       WRITE(6,609) FF,SwLR,YL1(0),BE,-YL1(1),AE
              BE= SwLR*BE + Sw*YL1(0)
              AE= SwLR*AE - Sw*YL1(1)
              ENDIF
          ENDIF
c** Kaiser (integer) controls sophistication of RKR calculation. 
c  If(Kaiser.le.0) do simple first order with Y00=0  
c  If(Kaiser > 0) apply "Kaiser" correction by setting lower bound of
c    integration as  v(min)= -1/2 - Y00/Y10  where Dunham Y00 & Y10 
c    calculated from energy derivatives evaluated at v=-1/2.
c** NSV is the number of different v-increments to be used in defining
c  the set of v-values for which turning points are to be calculated.
c** If(VEXT.le.0) perform RKR calculation with no inner wall smoothing
c  but calculate & print exponent coefficients cEXT Cext of exponential
c  functions approximately fitted to the 3 preceeding inner turning pts.
c** If(VEXT.gt.0) extrapolate inner wall above  v=VEXT  with an
c  exponential exactly fitted to last 3 inner T.P. with  v.le.VEXT
c-----------------------------------------------------------------------
      READ(5,*)  Kaiser, NSV, VEXT
c-----------------------------------------------------------------------
      IF(MORSE.EQ.1) THEN
          VEXT= -1.d0
          Kaiser= 0
          ENDIF
      IF(Kaiser.GT.0) then
c** For (Kaiser > 0) calculate Y00 and vib quantum no. at minimum V00 and
c  in NDE case, adjust D value to actual well depth  De= G(v=-1/2)+Y00
          Y00= AE*WE/(BE*12.d0)
          Y00= (BE- WEXE)/4.d0 + Y00 + Y00**2/BE
          IF(NDEGv.GE.1) DLIM= DLIM + Y00
          DV0= -Y00/WE
          V00= V00+ DV0
c* Iterate to correct for linear Y00/WE approximation
          I= 1
          VX(1)= V00
          CALL GVBV(I,VX,GV,DGDV,BV)
          DV0= -GV(1)/WE
          V00= V00+ DV0
          DV0= DV0- Y00/WE
          WRITE(6,618) Y00,DV0,V00,WE,WEXE,BE,AE
          IF(NDEGv.GE.1) WRITE(6,620) DLIM
          ENDIF
      I=1
      VX(I)= V00
      CALL GVBV(I,VX,GV,DGDV,BV)
      IF(MORSE.NE.1) BE= BV(1)
      REQ= DSQRT(CU/BE)
      WRITE(6,622) V00,GV(1),DGDV(1),-WEXE,BE,REQ,AE
      IF(VEXT.GT.0.d0) WRITE(6,624) VEXT
      DO  I= 1, NSV
c** Read and generate v-s at which turning points desired
c** For each of the NSV cases, generate v values ranging from V1(j)
c   to V2(j) with (positive) increment DV(j) (for j= 1 to NSV).
c*  The resulting values should increase monotonically.
c-----------------------------------------------------------------------
      READ(5,*,end=99) V1(I), DV(I), V2(I)
c-----------------------------------------------------------------------
          ENDDO
      IF((NDEGv.GE.1).AND.(V2(NSV).GT.VD)) V2(NSV)= DINT(VD)
      I= 0
      DO  J= 1,NSV
          NN= (V2(J)- V1(J)+ 1.d-4)/DV(J)
          IF(J.LT.NSV) NN= (V1(J+1)-V1(J))/DV(J)- 1
          I= I+1
          V(I)= V1(J)
          DO  K= 1,NN
              I= I+1
              V(I)= V(I-1) + DV(J)
              ENDDO
          IF(DABS(V(I)).LT.1.D-12) V(I)= 0.d0
          ENDDO
      NVIB= I
      NTP= 2*NVIB+1
      WRITE(6,626) NVIB,(V(I),I= 1,NVIB)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++ Commence actual turning point calculation here +++++++++++++++++++++
      WRITE(6,630)
      I638= 0
      DEI= 0.d0
      DRI= 0.d0
c** Commence loop over the NVIB values of vib quantum no. v ********
      DO  IVB= 1,NVIB
          I634= 0
          VUP= V(IVB)
          VX(1)= VUP
          I= 1
          CALL GVBV(I,VX,GV,DGDV,BV)
          GUP= GV(1)
          BUP= BV(1)
          DGUP= DGDV(1)
          IF(DGUP.LE.0.d0) THEN
              WRITE(6,632) VUP,GUP,DGUP
              IVIB= IVB-1
              NTP= 2*IVIB+ 1
              I1= NVIB- IVIB+ 1
              GO TO 72
              ENDIF
          F= 0.d0
          G= 0.d0
          NDIV= 1
          TSTF= 9.d99
          TSTG= 9.d99
c** In attempting to achieve integral convergence, consider up to NBIS
c  bisections of the interval.
          DO 50 IBIS= 0,NBIS
              TSTFB= TSTF
              TSTGB= TSTG
              VFN= V00
              FF= 0.d0
              GG= 0.d0
              RANGE= (VUP-V00)/DBLE(NDIV)
c** Sum over contributions from the NDIV segments of the interval
              DO  IDIV= 1,NDIV
                  VST= VFN
                  VFN= VST+RANGE
                  IF(IDIV.NE.NDIV) THEN
c** Quadrature weights and points for range with no singularity
                      D1= 0.5d0*(VFN+ VST)
                      D2= 0.5d0*(VFN- VST)
                      DO  I= 1,NGP
                          WW(I)= WG(I)
                          VX(I)= D1+D2*XG(I)
                          ENDDO
                    ELSE
                      VFN= VUP
                      FF= FF/2.d0
                      GG= GG/2.d0
c** Weights and points for range with 1/sqrt(E-V) singularity
                      D2= RANGE
                      DO  I= 1,NGP
                          WW(I)= W2(I)
                          VX(I)= VST+X2(I)*D2
                          ENDDO
                    ENDIF
                  CALL GVBV(NGP,VX,GV,DGDV,BV)
c** Actual quadrature loop starts here
                  DO  I= 1,NGP
                      TT= WW(I)/DSQRT(GUP-GV(I))
                      FF= FF+TT
                      GG= GG+TT*BV(I)
                      ENDDO
                  ENDDO
c** End of quadrature - now test for convergence
              FB= F
              GB= G
              F= FF*RANGE*SQCU
              G= GG*RANGE/SQCU
              TSTF= DABS(1.d0- FB/F)
              TSTG= 0.d0
              IF(((NDEGv.NE.0).OR.(LMAXBv.GE.1)).AND.(VUP.le.VEXT)) 
     1                                          TSTG= DABS(1.d0- GB/G)
              IF(MORSE.EQ.1) TSTG= 0.d0
              IF(DMAX1(TSTF,TSTG).LT.TOLER) GO TO 54
              IF((TSTF.GT.TSTFB).OR.(TSTG.GT.TSTGB)) THEN
                  I634= I634+ 1
                  IF((I634.GE.1).AND.(NDIV.GT.4)) THEN
                      WRITE(6,634) NDIV, TSTF,TSTFB, TSTG,TSTGB
                      GOTO 54   
                      ENDIF
                  ENDIF
              NDIV= 2*NDIV
   50         CONTINUE
          NDIV= NDIV/2
          WRITE(6,636) NDIV,TSTF,TSTG,TOLER
   54     IF(MORSE.EQ.1) THEN
c** Use Morse inner turning point when no rotational data input
              RMIN= Req - DLOG(1.d0 + DSQRT(GUP/De))/BETA
              RMAX= RMIN + F + F
              BUP= BE
            ELSE
              HEL= DSQRT(F/G+F*F)
              RMIN= HEL-F
              RMAX= HEL+F
            ENDIF
          I1= NVIB+1-IVB
          I2= NVIB+1+IVB
          IF(IVB.LE.1) GO TO 60
          DEIB= DEI
          DRIB= DRI
          I1P1= I1+1
          RMIN2= RT(I1P1)
          DEI= GUP- ET(I1P1)
          DRI= RMIN- RMIN2
          IF(IVB.LE.2) GO TO 60
          IF((VEXT.GT.0.d0).AND.(VUP.GT.VEXT)) GO TO 56
          IF(RMIN.GE.RT(I1P1)) THEN
c** If inner wall turns over, print warning
              IF(I638.LE.0) THEN
                  I638= I1P1
                  WRITE(6,638) VUP
                  ENDIF
              GOTO 60
              ENDIF
          IF((MORSE.GT.0).OR.(I638.GT.0)) GO TO 60
c** Use differences to get rough estimate of exponent coefficient CEXT
c  for local exponential fit to inner wall.
          I1P2= I1+2
          RMIN3= RT(I1P2)
          VRAT= DEI/DEIB
          CEXT= 2.d0*(DRI-VRAT*DRIB)/(DRI**2 + VRAT*DRIB**2)
c** Iteratively converge on exact (to machine precision) value of CEXT
          IF((DABS(CEXT*RMIN3).GT.70.D0).AND.(VUP.GT.0.d0)) THEN
              IF(I640.LE.0) WRITE(6,640) 
              I640= 1
              GO TO 55
              ENDIF
          VRAT= (GUP-ET(I1P1))/(GUP-ET(I1P2))
          ADCEXT= 1.d99
          DO  I= 1,15
              BDCEXT= ADCEXT
              EXP1= 1.d0
              EXP2= DEXP(-CEXT*(RT(I1P1)-RMIN))
              EXP3= DEXP(-CEXT*(RT(I1P2)-RMIN))
              FUN= (EXP1-EXP2)/(EXP1-EXP3)-VRAT
              DFUN= ((RMIN*EXP1-RT(I1P1)*EXP2) - (RMIN*EXP1-
     1             RT(I1P2)*EXP3)*(EXP1-EXP2)/(EXP1-EXP3))/(EXP1-EXP3)
              DCEXT= FUN/DFUN
              ADCEXT= DABS(DCEXT)
              IF((ADCEXT.GE.BDCEXT).AND.(ADCEXT.LT.1.d-10)) GO TO 230
              IF(ADCEXT.LE.0.d0) GO TO 230
              CEXT= CEXT+DCEXT
              ENDDO
          WRITE(6,642) DCEXT
  230     CONTINUE
          IF(CEXT.LE.0.d0) THEN
              I644= I644+ 1
              IF(I644.EQ.1) WRITE(6,644) VUP
              ENDIF
          EXP1= DEXP(-CEXT*RMIN)
          EXP2= EXP2*EXP1
          BEXT= (GUP-ET(I1P1))/(EXP1-EXP2)
          AEXT= GUP-BEXT*EXP1
c** Parameters for possible inner wall extrapolation now determined.
c
   55     WRITE(6,646) VUP,GUP,DGUP,BUP,RMIN,RMAX,NDIV,TSTF,TSTG,CEXT
          GO TO 64
c** Apply option to smoothly extrapolate inner wall above  v=VEXT
c  using a simple exponential
   56     REXT= DLOG(BEXT/(GUP-AEXT))/CEXT
          DR1= REXT-RMIN
          RMIN= REXT
          RMAX= RMAX+DR1
          WRITE(6,646)VUP,GUP,DGUP,BUP,RMIN,RMAX,NDIV,TSTF,TSTG,CEXT,DR1
          GO TO 64
   60     WRITE(6,646) VUP,GUP,DGUP,BUP,RMIN,RMAX,NDIV,TSTF,TSTG
   64     RT(I1)= RMIN
          RT(I2)= RMAX
          ET(I1)= GUP
          ET(I2)= GUP
          ENDDO
      I1= 1
      I2= NTP
   72 IF(VEXT.GT.0.d0) WRITE(6,648) VEXT,AEXT,BEXT,CEXT
c** Write ordered turning points compactly on channel 7 
c** To facilitate use of resulting potential array, add 5 extrapolated
c   inner-wall points in the output file  
      IF((CEXT.GT.0.d0).AND.(CEXT.LT.20.d0)) THEN
          IDR= -INT(1.0d3*ET(I1)*(RT(I1+1)-RT(I1))/(ET(I1+1)-ET(I1)))
          IR1= INT(1.d4*RT(I1))
          DO  I= 1,5
              RT(I1-I)= 1.D-4*(IR1- I*IDR)
              ET(I1-I)= AEXT + BEXT*DEXP(-CEXT*RT(I1-I))
              ENDDO
          I1= I1 - 5
          I2= NTP + I1 + 4
          ENDIF
      RT(NVIB+1)= REQ
      ET(NVIB+1)= 0.d0
      IF(I1.LT.0) WRITE(7,700) TITLE,NTP+5,ZMU
      IF(I1.GT.0) WRITE(7,700) TITLE,NTP,ZMU
      WRITE(7,702) (RT(I),ET(I),I= I1,I2)
      WRITE(6,650)
      GO TO 2
   99 STOP
  598 FORMAT(/' *** Input ERROR ***   NDEBv=',i3,' .GT. NDEGv=',i3,
     1  '  is NOT allowed!')
  600 FORMAT(1x,A78/1x,35('**')/' RKR potential for ',A2,'(',I3,')-',
     1 A2,'(',I3,')   with   Charge=',I2/' Reduced mass   ZMU=',
     2 F15.11,'   and constant   C_u/ZMU =',F16.12/5x,'from atomic masse
     3s:',f15.10,'  & ',F15.10,'(u)')
  602 FORMAT(/' Seek relative quadrature convergence',1PD8.1,
     1 '.   Bisect interval up to',i3,' times.'/5x,'performing',i3,
     2 '-point Gaussian quadrature in each segment')
  604 FORMAT(/' Represent  ',A2,"'s  by Tellinghuisen-type MXS mixed rep
     1resentation:"/ 1x,16('==')/I5,"'th order Dunham for  v .le. VS   &
     2   NDE for  v > VS,   with   VS=",F8.4:/4x,'with switching functio
     3n   F_s = 1/[1 + exp{(v-VS)/DVS}]   with   DVS=',f7.4/5x,'and a
     4sympotote energy (dissociation limit)   DLIM=',F11.4,' [cm-1]')
  606 FORMAT(/' The',i3,' Dunham  ',A2,'  expansion coefficients are'/
     1  (6X,1P4D18.10:))
  608 FORMAT(/' NDE for  ',A2,'  is an  (NP=',I2,'/NQ=',I2,') ',2A6,
     1 ' expansion in  (vD-v)  with'/8x,'X',I1,'(n=',I1,')=',1PD14.7:
     2 '  and leading num. and denom. powers',I3,'  & ',I3:/8X,'vD=',
     3 0PF12.6:'    D-G(v=-1/2)=',F14.6)
  609 FORMAT(/' *** CAUTION *** DVS may be too large.  At  v=',f7.3,
     1 '  where  Sw(LR)=',1PD9.2/6x,'Be(Dun)=',d8.1,'   Be(NDE)=',
     2 d8.1,'   ae(Dun)=',d8.1,'   ae(NDE)=',d8.1)
  610 FORMAT(5X,'Numerator coefficients are:',2(1PD20.12:)/
     1  (12X,3D20.12:))
  612 FORMAT(5X,'Denominator coefficients : ',2(1PD20.12:)/
     1  (12X,3D20.12:))
  614 FORMAT(/' NO rotational constants input, so inner wall of potentia
     1l is Morse function.'/'   Input   Req=',f9.6,'(Angst)   plus   we=
     2',f9.3,'  &  wexe=',f9.5,' [cm-1]'/'   yields Morse with  De=',
     3 f10.3,' [cm-1]   and   beta=',f9.6,' [1/Angst.]')
  618 FORMAT(/' Calculate   Y00=',F13.9,4X,'v(cor)=',F14.10,4x,'v(min)='
     1 ,F14.10/5x,'using    we=',f9.4,'    wexe=',f9.6,'    Be=',f10.6,
     2 '    ae=',f11.8)
  620 FORMAT(5x,'and corrected effective   De=',F13.6,' (after adding Y0
     10)')
  622 FORMAT(/' At  v00=',F9.5,'   Gv=',F12.8,'   dG/dv=',F9.4,
     1 '   (1/2)d2G/dv2=',F10.6/21x,'Bv=',F12.8,'   { ==>   Req=',F12.9,
     2 '(A) }'/21x,'alpha_e =',F12.9)
  624 FORMAT(/' Above   v =',f8.3,'   extrapolate inner wall with expone
     1ntial'/9x,'fitted to last 3 points ( & shift  RMAX  accordingly)')
  626 FORMAT(/' Calculate turning points at the',i4,' v-values'/
     1  (1X,11F7.2:))
  630 FORMAT(/' Resulting Turning Points:'/
     1 5X,'v',8X,'E(v)',4X,'dE(v)/dv',7X,'B(v)',10X,'Rmin(v)',8X,
     2 'Rmax(v)',3X,'NDIV  tst(f)   tst(g)',4X,'C(exp)',7X,'d(RMIN)'/
     3 1X,61('**'))
  632 format(/' STOP calculation at   v=',f6.2,'  where  E(v)=',f8.2,
     1  '  &  dE(v)/dv =',G14.7)
  634 FORMAT(' *** STOP ITERATION:  At  NDIV=',i3,
     1 '   tst(f)/(previous)=',1Pd8.1,'/',d7.1,'   tst(g)/(previous)=',
     2 d8.1,'/',d7.1)
  636 FORMAT(' *** CAUTION:',i3,' interval incomplete convergence:  tst(
     1f) & tst(g)=',1P2d8.1,'  while  TOLER=',d8.1)
  638 FORMAT(' *** WARNING ***  inner wall becomes unstable at   v =',
     1  F6.2,'  where  RMIN  turns outward!')
  640 FORMAT(' *** CAUTION *** inner wall exponent parameter becomes ver
     1y large so skip converging it.')
  642 FORMAT(' *** CAUTION *** FAIL to converge C(exp) after 15 tries.',
     1 '  Last step =',1PD10.2)
  644 FORMAT(' *** CAUTION ***  Inner potential wall has negative curvat
     1ure and requires smoothing for   VEXT .ge.',f7.2)
  646 FORMAT(F7.3,F12.4,F11.4,F15.10,2F15.10,I4,1P2D9.1,0PF11.6,F15.10)
  648 FORMAT(1X,61('**')//' For   v .GE.',f6.2,'  inner wall extrapolate
     1d as:   V(R) =',F13.4,' +',D15.8,'*exp(-',f12.8,'*R)')
  650 FORMAT(1x,61('**')/////)
  700 FORMAT(/1x,A78/' NTP=',I5,'   RKR turning points for  mu=',f14.10)
  702 FORMAT((F20.14,F19.11))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE GVBV(NL,VX,GV,DG,BV)
c** At the NL values of the vibrational quantum number VX, generate
c  values of Bv, Gv and their first derivatives w.r.t. v from 
c   (i) pure Dunham (NDEXv= 0),  (ii) pure NDE (NDEXv= 1)   or 
c   (iii) Tellinghuisen 'mixed' MXS function NDEXv .ge. 2   where
c   require   NDEBv .le. NDEGv
c
      INTEGER MXDUN           !! max polynomial order for vib expansions
      PARAMETER (MXDUN=25)  
c*** define types for local variables
      INTEGER I,L,KORDR,NL
      REAL*8 VX(NL),GV(NL),DG(NL),BV(NL),KM,DKM,Sw,SwLR,VPH,YY,DY
c** Common block for Dunham & MXS function parameters
      INTEGER LMAXGv,LMAXBv,NDEGv,NDEBv
      REAL*8 VS,DVS,Y00,YL0(0:MXDUN),YL1(0:MXDUN)
      COMMON /DUNPRM/Y00,VS,DVS,YL0,YL1,LMAXGv,LMAXBv,NDEGv,NDEBv
c** Common block for NDE function parameters
      INTEGER NLR,ITYPE,ITYPB,IZP0,IZQ0,IZP1,IZQ1,NP0,NQ0,NP1,NQ1
      REAL*8 VD,DLIM,XCN0,XCN1,P0(MXDUN),Q0(MXDUN),P1(MXDUN),Q1(MXDUN)
      COMMON /NDEPRM/VD,DLIM,XCN0,XCN1,P0,Q0,P1,Q1,NLR,ITYPE,ITYPB,
     1  IZP0,IZQ0,IZP1,IZQ1,NP0,NQ0,NP1,NQ1
c-----------------------------------------------------------------------
      DO  I= 1, NL
c** Loop over all levels, calculating energies, Bv's and their first
c  derivatives w.r.t. v
          VPH= VX(I)+ 0.5d0
c** First - for Dunham vib energy (or Dunham part of MXS energy)
          IF(NDEGv.NE.1) THEN
              YY= YL0(LMAXGv)
              DY= LMAXGv*YY
              DO  L= LMAXGv-1, 1, -1
                  YY= YY*VPH + YL0(L)
                  dY= DY*VPH + L*YL0(L)
                  ENDDO
              YY= YY*VPH + Y00
              GV(I)= YY
              DG(I)= DY
              ENDIF
          IF(NDEGv.GE.1) THEN
c** For NDE or MXS vibrational energy .....
              KORDR= 0
              CALL NDEDKM(VX(I),KORDR,KM,DKM,NLR,XCN0,DLIM,VD,
     1                                  IZP0,IZQ0,ITYPE,NP0,NQ0,P0,Q0)
              IF(NDEGv.EQ.1) THEN
                  GV(I)= KM
                  DG(I)= DKM
                ELSE
                  SwLR= DEXP((VX(i)- VS)/DVS)
                  Sw= 1.d0/(1.d0+ SwLR)
                  SwLR= SwLR * Sw
                  GV(I)= Sw*YY + SwLR*KM
                  DG(I)= Sw*DY + SwLR*DKM - (YY - KM)*Sw*SwLR/DVS
                ENDIF
              ENDIF 
          IF(NDEBv.NE.1) THEN
c
c** First - for Dunham Bv value (or Dunham part of MXS function for Bv)
              DY= 1.d0
              YY= 0.d0
              IF(LMAXBv.GE.0) THEN
                  YY= YL1(LMAXBv)
                  DY= LMAXBv*YY
                  DO  L= LMAXBv-1, 0, -1
                      YY= YY*VPH + YL1(L)
                      DY= DY*VPH + L*YL1(L)
                      ENDDO
                  BV(I)= YY
                  ENDIF
              ENDIF
          IF(NDEBv.GE.1) THEN
c** For NDE or MXS function for Bv value ...
              KORDR= 1
              CALL NDEDKm(VX(I),KORDR,KM,DKM,NLR,XCN1,DLIM,VD,
     1                                  IZP1,IZQ1,ITYPB,NP1,NQ1,P1,Q1)
              IF(NDEBv.EQ.1) THEN
                  BV(I)= KM
                ELSE
                  BV(I)= Sw*YY + SwLR*KM
                ENDIF
              ENDIF
          ENDDO
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c**********************************************************************
      SUBROUTINE NDEDKM(VV,KORDR,KM,DKM,NLR,XCN,DLIM,VD,IZP,IZQ,ITYP,
     1 NP,NQ,P,Q)
c** Subroutine which uses Near-Dissociation Expansions to generate
c  predicted value KM and vibrational first derivative DKM of band
c  constant of rotational order KORDR (KORDR= 0 for Gv, =1 for Bv, 
c  =2 for -Dv, =3 for Hv, ... etc., for vibrational level (integer) IV
c-----------------------------------------------------------------------
      INTEGER KORDR,IZP,IZQ,ITYP,IPW,J,NLR,NP,NQ
      REAL*8 KM,DKM,XCN,DLIM,VD,P(NP),Q(NQ), PW,PAB,DV,DVP,SMN,SMD,DSMN,
     1  DSMD,VV
      IPW= 2*NLR/(NLR-2)
      PW= 2.d0*NLR/(DBLE(NLR)- 2.d0)
      IF(KORDR.GE.1) THEN
          IPW= IPW - 2*KORDR
          PW= PW - 2.d0*KORDR
          ENDIF
      PAB= 1.d0
      IF(ITYP.EQ.2) PAB= PW
      DV= VD- VV
c** Evaluate the NDE numerator & denominator polynomials
      SMN= 1.d0
      SMD= 1.d0
      DSMN= 0.d0
      DSMD= 0.d0
      IF(NP.GT.0) THEN
c ... numerator polynomial ...
          DVP= 1.d0
          IF(IZP.GT.0) DVP= DV**IZP
          DO  J= 1,NP
              DSMN= DSMN+(J+IZP)*P(J)*DVP
              DVP= DVP*DV
              SMN= SMN+P(J)*DVP
              ENDDO
          ENDIF
      IF(NQ.GT.0) THEN
c ... denominator polynomial ...
          DVP= 1.d0
          IF(IZQ.GT.0) DVP= DV**IZQ
          DO  J= 1,NQ
              DSMD= DSMD+(J+IZQ)*Q(J)*DVP
              DVP= DVP*DV
              SMD= SMD+Q(J)*DVP
              ENDDO
          ENDIF
      IF(ITYP.EQ.3) THEN
          KM= XCN*DV**PW *DEXP(SMN- 1.d0)
          DKM= -KM*(PW/DV + DSMN)
          IF(KORDR.EQ.0) THEN
              KM= DLIM - KM
              DKM= -DKM
              ENDIF 
        ELSE
          IF(ITYP.EQ.1) THEN
              KM= XCN*SMN/SMD
              DKM= DV
              ENDIF
          IF(ITYP.EQ.2) THEN
              KM= XCN
              DKM= DV*SMN/SMD
              ENDIF
          IF(NLR.EQ.5) KM= KM*DABS(DKM)**PW
          IF(NLR.NE.5) KM= KM*DKM**IPW
          DKM= - KM*(PW/DV + PAB*(DSMN/SMN- DSMD/SMD))
          IF(KORDR.EQ.0) THEN
              KM= DLIM - KM
              DKM= -DKM
              ENDIF
        ENDIF
      RETURN
      END
c=======================================================================

c***********************************************************************
      SUBROUTINE MASSES(IAN,IMN,NAME,GELGS,DGNS,MASS,ABUND)
c***********************************************************************
c** For isotope with (input) atomic number IAN and mass number IMN,
c  return (output):  (i) as the right-adjusted 2-character variable NAME
c  the alphabetic symbol for that element,  (ii) the ground state
c  electronic degeneracy GELGS, (iii) the nuclear spin degeneracy DGNS,
c  (iv) the atomic mass MASS [amu], and  (v) the natural isotopic
c  abundance ABUND [in percent].   GELGS values based on atomic states
c  in Moore's "Atomic Energy Level" tables, the isotope masses are taken
c  from the 2012 mass table [Wang, Audi, Wapstra, Kondev, MacCormick, Xu
c  & Pfeiffer, Chin.Phys.C 36, 1603-2014 (2012)] ,the proton, deuteron,
c  and triton masses are taken from the 2010 fundamental constants table 
c  [Mohr, Taylor, & Newell, Rev. Mod. Phys. 84, 1587-1591 (2012)] and other
c  quantities from Tables 6.2 and 6.3 of "Quantities, Units and Symbols in
c  Physical Chemistry", by Mills et al.(Blackwell,2'nd Edition, Oxford,1993).
c** If the input value of IMN does not equal one of the tabulated values
c  for atomic species IAN, return the abundance-averaged standard atomic
c  weight of that atom and set DGNS=-1 and ABUND=-1.
c** For Atomic number IAN=0 and isotope mass numbers IMN=1-3,  return the
c    masses of the proton, deuteron, and triton, p,d & t, respectively
c Masses and properties of selected Halo nuclei an unstable nuclei included
c                 COPYRIGHT 2005-2015  :  last  updated 10 January 2016
c** By R.J. Le Roy, with assistance from 
c                 G.T. Kraemer, J.Y. Seto and K.V. Slaughter.
c***********************************************************************
      REAL*8 zm(0:123,0:15),mass,ab(0:123,15),abund
      INTEGER i,ian,imn,gel(0:123),nmn(0:123),mn(0:123,15),
     1                                        gns(0:123,15),DGNS,gelgs
      CHARACTER*2 NAME,AT(0:123)
cc
      DATA  at(0),gel(0),nmn(0),(mn(0,i),i=1,3)/' p',1,3,1,2,3/
      DATA  (zm(0,i),i=0,3)/1.008d0,1.007276466812d0,2.013553212712d0,
     2                 3.0155007134d0/
      DATA  (gns(0,i),i=1,3)/2,3,2/
      DATA  (ab(0,i),i=1,3)/0.d0, 0.d0, 0.d0/
c
      DATA  at(1),gel(1),nmn(1),(mn(1,i),i=1,3)/' H',2,3,1,2,3/
      DATA  (zm(1,i),i=0,3)/1.00794d0, 1.00782503223d0, 2.01410177812d0,
     1                 3.0160492779d0/
      DATA  (gns(1,i),i=1,3)/2,3,2/
      DATA  (ab(1,i),i=1,3)/99.985d0,0.015d0,0.d0/
c
      DATA  at(2),gel(2),nmn(2),(mn(2,i),i=1,4)/'He',1,4,3,4,6,8/
      DATA  (zm(2,i),i=0,4)/4.002602d0, 3.0160293201d0, 4.00260325413d0,
     1                                        6.0188891d0, 8.033922d0/
      DATA  (gns(2,i),i=1,4)/2,1,1,1/
      DATA  (ab(2,i),i=1,4)/0.000137d0,99.999863d0, 2*0.d0/
c
      DATA  at(3),gel(3),nmn(3),(mn(3,i),i=1,6)/'Li',2,6,6,7,8,9,11,12/
      DATA  (zm(3,i),i=0,6)/6.941d0, 6.0151228874d0, 7.016003437d0,
     1     8.02248736d0,9.0267895d0,11.043798d0,12.05378d0/
      DATA  (gns(3,i),i=1,6)/3,4,5,4,4,1/
      DATA  (ab(3,i),i=1,6)/7.5d0, 92.5d0, 4*0.d0/
c
      DATA  at(4),gel(4),nmn(4),(mn(4,i),i=1,8)/'Be',1,8,7,9,10,11,12,
     1                                                       14,15,16/
      DATA  (zm(4,i),i=0,8)/9.012182d0, 7.01692983d0, 9.01218307d0,
     1 10.0135338d0, 11.021658d0, 12.026921d0, 14.04289d0, 15.05346d0,
     2 16.06192d0/
      DATA  (gns(4,i),i=1,8)/4,4,3,2,1,1,2,1/
      DATA  (ab(4,i),i=1,8)/0.d0, 100.d0, 6*0.d0/
c
      DATA at(5),gel(5),nmn(5),(mn(5,i),i=1,10)/' B',2,10,8,10,11,12,
     1                                              13,14,15,17,18,19/
      DATA (zm(5,i),i=0,10)/10.811d0, 8.0246072d0, 10.0129369d0, 
     1          11.0093054d0, 12.0143521d0, 13.0177802d0, 14.025404d0,
     2          15.031103d0, 17.04699d0, 18.05617d0,19.06373d0/
      DATA  (gns(5,i),i=1,10)/5,7,4,3,4,5,4,4,1,4/
      DATA  (ab(5,i),i=1,10)/0.d0, 19.9d0,80.1d0, 7*0.d0/
c
      DATA at(6),gel(6),nmn(6),(mn(6,i),i=1,14)/' C',1,14,9,10,11,12,13,
     1               14,15,16,17,18,19,20,21,22/
      DATA (zm(6,i),i=0,14)/12.011d0, 9.0310367d0, 10.0168532d0,
     1          11.0114336d0, 12.d0, 13.00335483507d0, 14.003241989d0, 
     1  15.0105993d0, 16.014701d0, 17.022586d0, 18.02676d0, 19.03481d0,
     2  20.04032d0, 21.04934d0, 22.05720d0/
      DATA  (gns(6,i),i=1,14)/4,1,4,1,2,1,2,1,4,1,2,1,2,1/
      DATA  (ab(6,i),i=1,14)/3*0.d0, 98.90d0,1.10d0, 9*0.d0/
c
      DATA at(7),gel(7),nmn(7),(mn(7,i),i=1,2)/' N',4,2,14,15/
      DATA (zm(7,i),i=0,2)/14.00674d0, 14.00307400443d0,15.0001088989d0/
      DATA (gns(7,i),i=1,2)/3,2/
      DATA (ab(7,i),i=1,2)/99.634d0,0.366d0/
c
      DATA at(8),gel(8),nmn(8),(mn(8,i),i=1,3)/' O',5,3,16,17,18/
      DATA (zm(8,i),i=0,3)/15.9994d0, 15.99491461957d0, 16.9991317565d0,
     1                      17.9991596129d0/
      DATA (gns(8,i),i=1,3)/1,6,1/
      DATA (ab(8,i),i=1,3)/99.762d0, 0.038d0, 0.200d0/
c
      DATA at(9),gel(9),nmn(9),(mn(9,i),i=1,1)/' F',4,1,19/
      DATA (zm(9,i),i=0,1)/18.9984032d0, 18.9984031627d0/
      DATA (gns(9,i),i=1,1)/2/
      DATA (ab(9,i),i=1,1)/100.d0/
c
      DATA at(10),gel(10),nmn(10),(mn(10,i),i=1,4)/'Ne',1,4,17,20,21,22/
      DATA (zm(10,i),i=0,4)/20.1797d0, 17.017672d0, 19.9924401762d0, 
     1                                   20.99384669d0,21.991385115d0/
      DATA (gns(10,i),i=1,4)/2,1,4,1/
      DATA (ab(10,i),i=1,4)/0.d0, 90.48d0, 0.27d0, 9.25d0/
c
      DATA at(11),gel(11),nmn(11),(mn(11,i),i=1,1)/'Na',2,1,23/
      DATA (zm(11,i),i=0,1)/22.989768d0, 22.9897692820d0/
      DATA (gns(11,i),i=1,1)/4/
      DATA (ab(11,i),i=1,1)/100.d0/
c
      DATA at(12),gel(12),nmn(12),(mn(12,i),i=1,3)/'Mg',1,3,24,25,26/
      DATA (zm(12,i),i=0,3)/24.3050d0, 23.985041698d0, 24.98583698d0,
     1                       25.98259297d0/
      DATA (gns(12,i),i=1,3)/1,6,1/
      DATA (ab(12,i),i=1,3)/78.99d0, 10.00d0, 11.01d0/
c
      DATA at(13),gel(13),nmn(13),(mn(13,i),i=1,1)/'Al',2,1,27/
      DATA (zm(13,i),i=0,1)/26.981539d0, 26.98153853d0/
      DATA (gns(13,i),i=1,1)/6/
      DATA (ab(13,i),i=1,1)/100.d0/
c
      DATA at(14),gel(14),nmn(14),(mn(14,i),i=1,3)/'Si',1,3,28,29,30/
      DATA (zm(14,i),i=0,3)/28.0855d0, 27.9769265346d0, 28.9764946649d0,
     1                       29.973770136d0/
      DATA (gns(14,i),i=1,3)/1,2,1/
      DATA (ab(14,i),i=1,3)/92.23d0, 4.67d0, 3.10d0/
 
      DATA at(15),gel(15),nmn(15),(mn(15,i),i=1,2)/' P',4,2,26,31/
      DATA (zm(15,i),i=0,2)/30.973762d0, 26.01178d0, 30.9737619984d0/
      DATA (gns(15,i),i=1,2)/15,2/
      DATA (ab(15,i),i=1,2)/0.d0, 100.d0/
c
      DATA at(16),gel(16),nmn(16),(mn(16,i),i=1,5)/' S',5,5,27,32,33,
     1                                                          34,36/
      DATA (zm(16,i),i=0,5)/32.066d0, 27.01883d0, 31.9720711744d0,
     1                   32.9714589098d0,33.96786700d0, 35.96708071d0/
      DATA (gns(16,i),i=1,5)/6,1,4,1,1/
      DATA (ab(16,i),i=1,5)/0.d0, 95.02d0, 0.75d0, 4.21d0, 0.02d0/
c
      DATA at(17),gel(17),nmn(17),(mn(17,i),i=1,2)/'Cl',4,2,35,37/
      DATA (zm(17,i),i=0,2)/35.4527d0, 34.96885268d0, 36.96590260d0/
      DATA (gns(17,i),i=1,2)/4,4/
      DATA (ab(17,i),i=1,2)/75.77d0, 24.23d0/
c
      DATA at(18),gel(18),nmn(18),(mn(18,i),i=1,3)/'Ar',1,3,36,38,40/
      DATA (zm(18,i),i=0,3)/39.948d0, 35.967545105d0, 37.96273211d0,
     1                       39.9623831237d0/
      DATA (gns(18,i),i=1,3)/1,1,1/
      DATA (ab(18,i),i=1,3)/0.337d0, 0.063d0, 99.600d0/
c
      DATA at(19),gel(19),nmn(19),(mn(19,i),i=1,3)/' K',2,3,39,40,41/
      DATA (zm(19,i),i=0,3)/39.0983d0, 38.963706486d0, 39.96399817d0,
     1                       40.961825258d0/
      DATA (gns(19,i),i=1,3)/4,9,4/
      DATA (ab(19,i),i=1,3)/93.2581d0, 0.0117d0, 6.7302d0/
 
      DATA at(20),gel(20),nmn(20),(mn(20,i),i=1,6)/'Ca',1,6,40,42,43,44,
     1                                              46,48/
      DATA (zm(20,i),i=0,6)/40.078d0, 39.962590864d0, 41.95861783d0,
     1         42.95876644d0, 43.9554816d0, 45.9536890d0, 47.95252277d0/
      DATA (gns(20,i),i=1,6)/1,1,8,1,1,1/
      DATA (ab(20,i),i=1,6)/96.941d0, 0.647d0, 0.135d0, 2.086d0,
     1                      0.004d0, 0.187d0/
c
      DATA at(21),gel(21),nmn(21),(mn(21,i),i=1,1)/'Sc',4,1,45/
      DATA (zm(21,i),i=0,1)/44.955910d0, 44.9559083d0/
      DATA (gns(21,i),i=1,1)/8/
      DATA (ab(21,i),i=1,1)/100.d0/
c
      DATA at(22),gel(22),nmn(22),(mn(22,i),i=1,5)/'Ti',5,5,46,47,48,49,
     1                                              50/
      DATA (zm(22,i),i=0,5)/47.88d0, 45.9526277d0, 46.9517588d0,
     1         47.9479420d0, 48.9478657d0, 49.9447869d0/
      DATA (gns(22,i),i=1,5)/1,6,1,8,1/
      DATA (ab(22,i),i=1,5)/8.0d0, 7.3d0, 73.8d0, 5.5d0, 5.4d0/
c
      DATA at(23),gel(23),nmn(23),(mn(23,i),i=1,2)/' V',4,2,50,51/
      DATA (zm(23,i),i=0,2)/50.9415d0, 49.9471560d0, 50.9439570d0/
      DATA (gns(23,i),i=1,2)/13,8/
      DATA (ab(23,i),i=1,2)/0.250d0, 99.750d0/
c
      DATA at(24),gel(24),nmn(24),(mn(24,i),i=1,4)/'Cr',7,4,50,52,53,54/
      DATA (zm(24,i),i=0,4)/51.9961d0, 49.9460418d0, 51.9405062d0,
     1                       52.9406481d0, 53.9388792d0/
      DATA (gns(24,i),i=1,4)/1,1,4,1/
      DATA (ab(24,i),i=1,4)/4.345d0, 83.789d0, 9.501d0, 2.365d0/
c
      DATA at(25),gel(25),nmn(25),(mn(25,i),i=1,1)/'Mn',6,1,55/
      DATA (zm(25,i),i=0,1)/54.93805d0, 54.938049d0/
      DATA (gns(25,i),i=1,1)/6/
      DATA (ab(25,i),i=1,1)/100.d0/
c
      DATA at(26),gel(26),nmn(26),(mn(26,i),i=1,4)/'Fe',9,4,54,56,57,58/
      DATA (zm(26,i),i=0,4)/55.847d0, 53.9396090d0, 55.9349363d0,
     1                       56.9353928d0, 57.9332744d0/
      DATA (gns(26,i),i=1,4)/1,1,2,1/
      DATA (ab(26,i),i=1,4)/5.8d0, 91.72d0, 2.2d0, 0.28d0/
c
      DATA at(27),gel(27),nmn(27),(mn(27,i),i=1,1)/'Co',10,1,59/
      DATA (zm(27,i),i=0,1)/58.93320d0, 58.9331943d0/
      DATA (gns(27,i),i=1,1)/8/
      DATA (ab(27,i),i=1,1)/100.d0/
c
      DATA at(28),gel(28),nmn(28),(mn(28,i),i=1,5)/'Ni',9,5,58,60,61,62,
     1                                              64/
      DATA (zm(28,i),i=0,5)/58.69d0, 57.9353424d0, 59.9307859d0,
     1         60.9310556d0, 61.9283454d0, 63.9279668d0/
      DATA (gns(28,i),i=1,5)/1,1,4,1,1/
      DATA (ab(28,i),i=1,5)/68.077d0,26.223d0,1.140d0,3.634d0,0.926d0/
c
      DATA at(29),gel(29),nmn(29),(mn(29,i),i=1,2)/'Cu',2,2,63,65/
      DATA (zm(29,i),i=0,2)/63.546d0, 62.9295977d0,64.9277897d0/
      DATA (gns(29,i),i=1,2)/4,4/
      DATA (ab(29,i),i=1,2)/69.17d0, 30.83d0/
c
      DATA at(30),gel(30),nmn(30),(mn(30,i),i=1,5)/'Zn',1,5,64,66,67,68,
     1                                              70/
      DATA (zm(30,i),i=0,5)/65.40d0, 63.9291420d0, 65.9260338d0,
     1         66.9271277d0, 67.9248446d0, 69.9253192d0/
      DATA (gns(30,i),i=1,5)/1,1,6,1,1/
      DATA (ab(30,i),i=1,5)/48.6d0, 27.9d0, 4.1d0, 18.8d0, 0.6d0/
c
      DATA at(31),gel(31),nmn(31),(mn(31,i),i=1,2)/'Ga',2,2,69,71/
      DATA (zm(31,i),i=0,2)/69.723d0, 68.9255735d0, 70.9247026d0/
      DATA (gns(31,i),i=1,2)/4,4/
      DATA (ab(31,i),i=1,2)/60.108d0, 39.892d0/
c
      DATA at(32),gel(32),nmn(32),(mn(32,i),i=1,5)/'Ge',1,5,70,72,73,74,
     1                                              76/
      DATA (zm(32,i),i=0,5)/72.61d0, 69.9242488d0, 71.92207583d0,
     1         72.92345896d0, 73.921177762d0, 75.921402726d0/
      DATA (gns(32,i),i=1,5)/1,1,10,1,1/
      DATA (ab(32,i),i=1,5)/21.23d0, 27.66d0, 7.73d0, 35.94d0, 7.44d0/
c
      DATA at(33),gel(33),nmn(33),(mn(33,i),i=1,1)/'As',4,1,75/
      DATA (zm(33,i),i=0,1)/74.92159d0, 74.9215946d0/
      DATA (gns(33,i),i=1,1)/4/
      DATA (ab(33,i),i=1,1)/100.d0/
c
      DATA at(34),gel(34),nmn(34),(mn(34,i),i=1,6)/'Se',5,6,74,76,77,78,
     1                                              80,82/
      DATA (zm(34,i),i=0,6)/78.96d0, 73.922475935d0, 75.919213704d0,
     1         76.91991415d0, 77.91730928d0, 79.9165218d0, 81.9166995d0/
      DATA (gns(34,i),i=1,6)/1,1,2,1,1,1/
      DATA (ab(34,i),i=1,6)/0.89d0, 9.36d0, 7.63d0, 23.78d0, 49.61d0,
     1                      8.73d0/
c
      DATA at(35),gel(35),nmn(35),(mn(35,i),i=1,2)/'Br',4,2,79,81/
      DATA (zm(35,i),i=0,2)/79.904d0, 78.9183376d0, 80.9162897d0/
      DATA (gns(35,i),i=1,2)/4,4/
      DATA (ab(35,i),i=1,2)/50.69d0, 49.31d0/
c
      DATA at(36),gel(36),nmn(36),(mn(36,i),i=1,6)/'Kr',1,6,78,80,82,83,
     1                                              84,86/
      DATA (zm(36,i),i=0,6)/83.80d0, 77.9203649d0, 79.9163781d0,
     1     81.9134827d0, 82.9141272d0, 83.911497728d0, 85.910610627d0/
      DATA (gns(36,i),i=1,6)/1,1,1,10,1,1/
      DATA (ab(36,i),i=1,6)/0.35d0, 2.25d0, 11.6d0, 11.5d0, 57.0d0,
     1                      17.3d0/
c
      DATA at(37),gel(37),nmn(37),(mn(37,i),i=1,2)/'Rb',2,2,85,87/
      DATA (zm(37,i),i=0,2)/85.4678d0, 84.911789738d0, 86.909180532d0/
      DATA (gns(37,i),i=1,2)/6,4/
      DATA (ab(37,i),i=1,2)/72.165d0, 27.835d0/
c
      DATA at(38),gel(38),nmn(38),(mn(38,i),i=1,4)/'Sr',1,4,84,86,87,88/
      DATA (zm(38,i),i=0,4)/87.62d0, 83.9134191d0, 85.9092606d0,
     1                      86.9088775d0, 87.9056125d0/
      DATA (gns(38,i),i=1,4)/1,1,10,1/
      DATA (ab(38,i),i=1,4)/0.56d0, 9.86d0, 7.00d0, 82.58d0/
c
      DATA at(39),gel(39),nmn(39),(mn(39,i),i=1,1)/' Y',4,1,89/
      DATA (zm(39,i),i=0,1)/88.90585d0, 88.9058403d0/
      DATA (gns(39,i),i=1,1)/2/
      DATA (ab(39,i),i=1,1)/100.d0/
c
      DATA at(40),gel(40),nmn(40),(mn(40,i),i=1,5)/'Zr',5,5,90,91,92,94,
     1                                              96/
      DATA (zm(40,i),i=0,5)/91.224d0, 89.9046977d0, 90.9056396d0,
     1                      91.9050347d0, 93.9063108d0, 95.9082714d0/
      DATA (gns(40,i),i=1,5)/1,6,1,1,1/
      DATA (ab(40,i),i=1,5)/51.45d0, 11.22d0, 17.15d0, 17.38d0, 2.80d0/
c
      DATA at(41),gel(41),nmn(41),(mn(41,i),i=1,1)/'Nb',2,1,93/
      DATA (zm(41,i),i=0,1)/92.90638d0, 92.9063730d0/
      DATA (gns(41,i),i=1,1)/10/
      DATA (ab(41,i),i=1,1)/100.d0/
c
      DATA at(42),gel(42),nmn(42),(mn(42,i),i=1,7)/'Mo',7,7,92,94,95,96,
     1                                              97,98,100/
      DATA (zm(42,i),i=0,7)/95.94d0, 91.9068080d0, 93.9050849d0,
     1        94.9058388d0, 95.9046761d0, 96.9060181d0, 97.9054048d0,
     2        99.9074718d0/
      DATA (gns(42,i),i=1,7)/1,1,6,1,6,1,1/
      DATA (ab(42,i),i=1,7)/14.84d0, 9.25d0, 15.92d0, 16.68d0, 9.55d0,
     1                      24.13d0, 9.63d0/
c
      DATA at(43),gel(43),nmn(43),(mn(43,i),i=1,1)/'Tc',6,1,98/
      DATA (zm(43,i),i=0,1)/97.907215d0, 97.907212d0/
      DATA (gns(43,i),i=1,1)/13/
      DATA (ab(43,i),i=1,1)/100.d0/
c
      DATA at(44),gel(44),nmn(44),(mn(44,i),i=1,7)/'Ru',11,7,96,98,99,
     1                                              100,101,102,104/
      DATA (zm(44,i),i=0,7)/101.07d0, 95.9075903d0, 97.905287d0,
     1     98.9059341d0, 99.9042143d0, 100.9055769d0, 101.9043441d0,
     2     103.9054275d0/
      DATA (gns(44,i),i=1,7)/1,1,6,1,6,1,1/
      DATA (ab(44,i),i=1,7)/5.52d0, 1.88d0, 12.7d0, 12.6d0, 17.0d0,
     1                      31.6d0, 18.7d0/
c
      DATA at(45),gel(45),nmn(45),(mn(45,i),i=1,1)/'Rh',10,1,103/
      DATA (zm(45,i),i=0,1)/102.90550d0, 102.9054980d0/
      DATA (gns(45,i),i=1,1)/2/
      DATA (ab(45,i),i=1,1)/100.d0/
c
      DATA at(46),gel(46),nmn(46),(mn(46,i),i=1,6)/'Pd',1,6,102,104,105,
     1                                              106,108,110/
      DATA (zm(46,i),i=0,6)/106.42d0, 101.9056022d0, 103.9040305d0,
     1       104.9050796d0, 105.9034804d0, 107.9038916d0, 109.9051722d0/
      DATA (gns(46,i),i=1,6)/1,1,6,1,1,1/
      DATA (ab(46,i),i=1,6)/1.02d0, 11.14d0, 22.33d0, 27.33d0, 26.46d0,
     1                      11.72d0/
c
      DATA at(47),gel(47),nmn(47),(mn(47,i),i=1,2)/'Ag',2,2,107,109/
      DATA (zm(47,i),i=0,2)/107.8682d0, 106.9050916d0, 108.9047553d0/
      DATA (gns(47,i),i=1,2)/2,2/
      DATA (ab(47,i),i=1,2)/51.839d0, 48.161d0/
c
      DATA at(48),gel(48),nmn(48),(mn(48,i),i=1,8)/'Cd',1,8,106,108,110,
     1                                             111,112,113,114,116/ 
      DATA (zm(48,i),i=0,8)/112.411d0, 105.9064599d0, 107.9041834d0, 
     1       109.9030066d0, 110.9041829d0, 111.9027629d0, 112.9044081d0,
     2       113.9033651d0, 115.90476315d0/
      DATA (gns(48,i),i=1,8)/1,1,1,2,1,2,1,1/
      DATA (ab(48,i),i=1,8)/1.25d0, 0.89d0, 12.49d0, 12.80d0, 24.13d0,
     1                      12.22d0, 28.73d0, 7.49d0/
c
      DATA at(49),gel(49),nmn(49),(mn(49,i),i=1,2)/'In',2,2,113,115/
      DATA (zm(49,i),i=0,2)/114.818d0, 112.9040618d0, 114.903878776d0/
      DATA  (gns(49,i),i=1,2)/10,10/
      DATA (ab(49,i),i=1,2)/4.3d0, 95.7d0/
c
      DATA at(50),gel(50),nmn(50),(mn(50,i),i=1,10)/'Sn',1,10,112,114,
     1                                 115,116,117,118,119,120,122,124/
      DATA (zm(50,i),i=0,10)/118.710d0, 111.9048239d0, 113.9027827d0,
     1    114.903344699d0, 115.90174280d0, 116.9029540d0, 117.9016066d0,
     2    118.9033112d0, 119.9022016d0, 121.9034438d0, 123.9052766d0/
      DATA (gns(50,i),i=1,10)/1,1,2,1,2,1,2,1,1,1/
      DATA (ab(50,i),i=1,10)/0.97d0, 0.65d0, 0.34d0, 14.53d0, 7.68d0,
     1                       24.23d0, 8.59d0, 32.59d0, 4.63d0, 5.79d0/
c
      DATA at(51),gel(51),nmn(51),(mn(51,i),i=1,2)/'Sb',4,2,121,123/
      DATA (zm(51,i),i=0,2)/121.757d0, 120.903812d0, 122.9042132d0/
      DATA (gns(51,i),i=1,2)/6,8/
      DATA (ab(51,i),i=1,2)/57.36d0, 42.64d0/
c
      DATA at(52),gel(52),nmn(52),(mn(52,i),i=1,8)/'Te',5,8,120,122,123,
     1                                             124,125,126,128,130/
      DATA (zm(52,i),i=0,8)/127.60d0, 119.904059d0, 121.9030435d0,
     1    122.9042698d0, 123.9028171d0, 124.9044299d0, 125.9033109d0,
     2    127.9044613d0, 129.906222749d0/
      DATA (gns(52,i),i=1,8)/1,1,2,1,2,1,1,1/
      DATA (ab(52,i),i=1,8)/0.096d0, 2.603d0, 0.908d0, 4.816d0,
     1                      7.139d0, 18.95d0, 31.69d0, 33.80d0/
c
      DATA at(53),gel(53),nmn(53),(mn(53,i),i=1,2)/' I',4,2,127,129/
      DATA (zm(53,i),i=0,2)/126.90447d0, 126.904472d0, 128.904984d0/
      DATA (gns(53,i),i=1,2)/6,8/
      DATA (ab(53,i),i=1,2)/100.d0,0.d0/
c
      DATA at(54),gel(54),nmn(54),(mn(54,i),i=1,9)/'Xe',1,9,124,126,128,
     1                                          129,130,131,132,134,136/
      DATA (zm(54,i),i=0,9)/131.29d0, 123.9058920d0, 125.904298d0,
     1    127.9035310d0, 128.904780861d0,129.903509350d0,130.90508406d0,
     2    131.904155086d0, 133.9053947d0, 135.907214484d0/
      DATA (gns(54,i),i=1,9)/1,1,1,2,1,4,1,1,1/
      DATA (ab(54,i),i=1,9)/0.10d0, 0.09d0, 1.91d0, 26.4d0, 4.1d0,
     1                      21.2d0, 26.9d0, 10.4d0, 8.9d0/
c
      DATA at(55),gel(55),nmn(55),(mn(55,i),i=1,1)/'Cs',2,1,133/
      DATA (zm(55,i),i=0,1)/132.90543d0, 132.905451961d0/
      DATA (gns(55,i),i=1,1)/8/
      DATA (ab(55,i),i=1,1)/100.d0/
c
      DATA at(56),gel(56),nmn(56),(mn(56,i),i=1,7)/'Ba',1,7,130,132,134,
     1                                             135,136,137,138/
      DATA (zm(56,i),i=0,7)/137.327d0, 129.9063207d0, 131.9050611d0,
     1    133.90450818d0, 134.90568838d0, 135.90457573d0, 136.9058271d0,
     2    137.9052470d0/
      DATA (gns(56,i),i=1,7)/1,1,1,4,1,4,1/
      DATA (ab(56,i),i=1,7)/0.106d0, 0.101d0, 2.417d0, 6.592d0, 
     1                      7.854d0, 11.23d0, 71.70d0/
c
      DATA at(57),gel(57),nmn(57),(mn(57,i),i=1,2)/'La',4,2,138,139/
      DATA (zm(57,i),i=0,2)/138.9055d0, 137.907115d0, 138.9063563d0/
      DATA (gns(57,i),i=1,2)/11,8/ 
      DATA (ab(57,i),i=1,2)/0.0902d0, 99.9098d0/
c
      DATA at(58),gel(58),nmn(58),(mn(58,i),i=1,4)/'Ce',9,4,136,138,140,
     1                                             142/
      DATA (zm(58,i),i=0,4)/140.115d0, 135.9071292d0, 137.905991d0,
     1    139.9054431d0, 141.9092504d0/
      DATA (gns(58,i),i=1,4)/1,1,1,1/
      DATA (ab(58,i),i=1,4)/0.19d0, 0.25d0, 88.48d0, 11.08d0/
c
      DATA at(59),gel(59),nmn(59),(mn(59,i),i=1,1)/'Pr',10,1,141/
      DATA (zm(59,i),i=0,1)/140.90765d0, 140.9076576d0/
      DATA (gns(59,i),i=1,1)/6/
      DATA (ab(59,i),i=1,1)/100.d0/
c
      DATA at(60),gel(60),nmn(60),(mn(60,i),i=1,7)/'Nd',9,7,142,143,144,
     1                                             145,146,148,150/
      DATA (zm(60,i),i=0,7)/144.24d0, 141.9077290d0, 142.9098200d0,
     1    143.9100930d0, 144.9125793d0, 145.9131226d0, 147.9168993d0,
     2    149.9209022d0/
      DATA (gns(60,i),i=1,7)/1,8,1,8,1,1,1/
      DATA (ab(60,i),i=1,7)/27.13d0, 12.18d0, 23.80d0, 8.30d0, 17.19d0,
     1                       5.76d0, 5.64d0/
c
      DATA at(61),gel(61),nmn(61),(mn(61,i),i=1,1)/'Pm',6,1,145/
      DATA (zm(61,i),i=0,1)/144.912743d0, 144.912756d0/
      DATA (gns(61,i),i=1,1)/6/
      DATA (ab(61,i),i=1,1)/100.d0/
c
      DATA at(62),gel(62),nmn(62),(mn(62,i),i=1,7)/'Sm',1,7,144,147,148,
     1                                             149,150,152,154/
      DATA (zm(62,i),i=0,7)/150.36d0, 143.9120065d0, 146.9149044d0,
     1    147.9148292d0, 148.9171921d0, 149.9172829d0, 151.9197397d0,
     2    153.9222169d0/
      DATA (gns(62,i),i=1,7)/1,8,1,8,1,1,1/
      DATA (ab(62,i),i=1,7)/3.1d0, 15.0d0, 11.3d0, 13.8d0, 7.4d0,
     1                      26.7d0, 22.7d0/
c
      DATA at(63),gel(63),nmn(63),(mn(63,i),i=1,2)/'Eu',8,2,151,153/
      DATA (zm(63,i),i=0,2)/151.965d0, 150.9198578d0, 152.9212380d0/
      DATA (gns(63,i),i=1,2)/6,6/
      DATA (ab(63,i),i=1,2)/47.8d0, 52.2d0/
c
      DATA at(64),gel(64),nmn(64),(mn(64,i),i=1,7)/'Gd',5,7,152,154,155,
     1                                              156,157,158,160/
      DATA (zm(64,i),i=0,7)/157.25d0, 151.9197995d0, 153.9208741d0,
     1    154.9226305d0, 155.9221312d0, 156.9239686d0, 157.9241123d0,
     2    159.9270624d0/
      DATA (gns(64,i),i=1,7)/1,1,4,1,4,1,1/
      DATA (ab(64,i),i=1,7)/0.20d0, 2.18d0, 14.80d0, 20.47d0, 15.65d0,
     1                      24.84d0, 21.86d0/
c
      DATA at(65),gel(65),nmn(65),(mn(65,i),i=1,1)/'Tb',16,1,159/
      DATA (zm(65,i),i=0,1)/158.92534d0, 158.9253547d0/
      DATA (gns(65,i),i=1,1)/4/
      DATA (ab(65,i),i=1,1)/100.d0/
c
      DATA at(66),gel(66),nmn(66),(mn(66,i),i=1,7)/'Dy',17,7,156,158,
     1                                           160,161,162,163,164/
      DATA (zm(66,i),i=0,7)/162.50d0, 155.9242847d0, 157.924416d0,
     1    159.9252046d0, 160.9269405d0, 161.9268056d0, 162.9287383d0,
     2    163.9291819d0/
      DATA (gns(66,i),i=1,7)/1,1,1,6,1,6,1/
      DATA (ab(66,i),i=1,7)/0.06d0, 0.10d0, 2.34d0, 18.9d0, 25.5d0,
     1                      24.9d0, 28.2d0/
c
      DATA at(67),gel(67),nmn(67),(mn(67,i),i=1,1)/'Ho',16,1,165/
      DATA (zm(67,i),i=0,1)/164.93032d0, 164.9303288d0/
      DATA (gns(67,i),i=1,1)/8/
      DATA (ab(67,i),i=1,1)/100.d0/
     
      DATA at(68),gel(68),nmn(68),(mn(68,i),i=1,6)/'Er',13,6,162,164,
     1                                            166,167,168,170/
      DATA (zm(68,i),i=0,6)/167.26d0, 161.9287884d0, 163.9292088d0,
     1    165.9302995d0, 166.9320546d0, 167.9323767d0, 169.9354702d0/
      DATA (gns(68,i),i=1,6)/1,1,1,8,1,1/
      DATA (ab(68,i),i=1,6)/0.14d0, 1.61d0, 33.6d0, 22.95d0, 26.8d0,
     1                      14.9d0/
c
      DATA at(69),gel(69),nmn(69),(mn(69,i),i=1,1)/'Tm',8,1,169/  
      DATA (zm(69,i),i=0,1)/168.93421d0, 168.9342179d0/
      DATA (gns(69,i),i=1,1)/2/
      DATA (ab(69,i),i=1,1)/100.d0/
c
      DATA at(70),gel(70),nmn(70),(mn(70,i),i=1,7)/'Yb',1,7,168,170,171,
     1                                            172,173,174,176/
      DATA (zm(70,i),i=0,7)/173.04d0, 167.9338896d0, 169.9347664d0,
     1    170.9363302d0, 171.9363859d0, 172.9382151d0, 173.9388664d0,
     2    175.9425764d0/
      DATA (gns(70,i),i=1,7)/1,1,2,1,6,1,1/
      DATA (ab(70,i),i=1,7)/0.13d0, 3.05d0, 14.3d0, 21.9d0, 16.12d0,
     1                      31.8d0, 12.7d0/
c
      DATA at(71),gel(71),nmn(71),(mn(71,i),i=1,2)/'Lu',4,2,175,176/
      DATA (zm(71,i),i=0,2)/174.967d0, 174.9407752d0, 175.9426897d0/
      DATA (gns(71,i),i=1,2)/6,15/
      DATA (ab(71,i),i=1,2)/97.41d0, 2.59d0/
c
      DATA at(72),gel(72),nmn(72),(mn(72,i),i=1,6)/'Hf',5,6,174,176,177,
     1                                             178,179,180/
      DATA (zm(72,i),i=0,6)/178.49d0, 173.9400461d0, 175.9414076d0,
     1    176.9432277d0, 177.9437058d0, 178.9458232d0, 179.9465570d0/
      DATA (gns(72,i),i=1,6)/1,1,8,1,10,1/
      DATA (ab(72,i),i=1,6)/0.162d0, 5.206d0, 18.606d0, 27.297d0,
     1                      13.629d0, 35.100d0/
c
      DATA at(73),gel(73),nmn(73),(mn(73,i),i=1,2)/'Ta',4,2,180,181/
      DATA (zm(73,i),i=0,2)/180.9479d0, 179.9474648d0, 180.9479958d0/
      DATA (gns(73,i),i=1,2)/17,8/
      DATA (ab(73,i),i=1,2)/0.012d0, 99.988d0/
c
      DATA at(74),gel(74),nmn(74),(mn(74,i),i=1,5)/' W',1,5,180,182,183,
     1                                             184,186/
      DATA (zm(74,i),i=0,5)/183.84d0, 179.9467108d0, 181.9482039d0,
     1    182.9502227d0, 183.9509309d0, 185.9543628d0/
      DATA (gns(74,i),i=1,5)/1,1,2,1,1/
      DATA (ab(74,i),i=1,5)/0.13d0, 26.3d0, 14.3d0, 30.67d0, 28.6d0/
c
      DATA at(75),gel(75),nmn(75),(mn(75,i),i=1,2)/'Re',6,2,185,187/
      DATA (zm(75,i),i=0,2)/186.207d0, 184.9529545d0, 186.9557501d0/
      DATA (gns(75,i),i=1,2)/6,6/
      DATA (ab(75,i),i=1,2)/37.40d0, 62.60d0/
c
      DATA at(76),gel(76),nmn(76),(mn(76,i),i=1,7)/'Os',9,7,184,186,187,
     1                                             188,189,190,192/
      DATA (zm(76,i),i=0,7)/190.23d0, 183.9524885d0, 185.9538350d0,
     1    186.9557474d0, 187.9558352d0, 188.9581442d0, 189.9584437d0,
     2    191.9614770d0/
      DATA (gns(76,i),i=1,7)/1,1,2,1,4,1,1/
      DATA (ab(76,i),i=1,7)/0.02d0, 1.58d0, 1.6d0, 13.3d0, 16.1d0,
     1                      26.4d0, 41.0d0/
c
      DATA at(77),gel(77),nmn(77),(mn(77,i),i=1,2)/'Ir',10,2,191,193/
      DATA (zm(77,i),i=0,2)/192.22d0, 190.9605893d0, 192.9629216d0/
      DATA (gns(77,i),i=1,2)/4,4/
      DATA (ab(77,i),i=1,2)/37.3d0, 62.7d0/
c
c
      DATA at(78),gel(78),nmn(78),(mn(78,i),i=1,6)/'Pt',7,6,190,192,194,
     1                                            195,196,198/
      DATA (zm(78,i),i=0,6)/195.08d0, 189.959930d0, 191.961039d0,
     1    193.9626809d0, 194.9647917d0, 195.9649521d0, 197.9678949d0/
      DATA (gns(78,i),i=1,6)/1,1,1,2,1,1/
      DATA (ab(78,i),i=1,6)/0.01d0,0.79d0,32.9d0,33.8d0,25.3d0,7.2d0/
c
      DATA at(79),gel(79),nmn(79),(mn(79,i),i=1,1)/'Au',2,1,197/
      DATA (zm(79,i),i=0,1)/196.96654d0, 196.9665688d0/
      DATA (gns(79,i),i=1,1)/4/
      DATA (ab(79,i),i=1,1)/100.d0/
c
      DATA at(80),gel(80),nmn(80),(mn(80,i),i=1,7)/'Hg',1,7,196,198,199,
     1                                            200,201,202,204/
      DATA (zm(80,i),i=0,7)/200.59d0, 195.965833d0, 197.9667686d0,
     1    198.9682806d0, 199.9683266d0, 200.9703028d0, 201.9706434d0,
     2    203.9734940d0/
      DATA (gns(80,i),i=1,7)/1,1,2,1,4,1,1/
      DATA (ab(80,i),i=1,7)/0.15d0, 9.97d0, 16.87d0, 23.10d0, 13.18d0,
     1                      29.86d0, 6.87d0/
c
      DATA at(81),gel(81),nmn(81),(mn(81,i),i=1,2)/'Tl',2,2,203,205/
      DATA (zm(81,i),i=0,2)/204.3833d0, 202.9723446d0, 204.9744278d0/
      DATA (gns(81,i),i=1,2)/2,2/
      DATA (ab(81,i),i=1,2)/29.524d0, 70.476d0/
c
      DATA at(82),gel(82),nmn(82),(mn(82,i),i=1,4)/'Pb',1,4,204,206,207,
     1                                             208/
      DATA (zm(82,i),i=0,4)/207.2d0, 203.9730440d0, 205.9744657d0,
     1    206.9758973d0, 207.9766525d0/
      DATA (gns(82,i),i=1,4)/1,1,2,1/
      DATA (ab(82,i),i=1,4)/1.4d0, 24.1d0, 22.1d0, 52.4d0/
c
      DATA at(83),gel(83),nmn(83),(mn(83,i),i=1,1)/'Bi',4,1,209/
      DATA (zm(83,i),i=0,1)/208.98037d0, 208.9803991d0/
      DATA (gns(83,i),i=1,1)/10/
      DATA (ab(83,i),i=1,1)/100.d0/
c
      DATA at(84),gel(84),nmn(84),(mn(84,i),i=1,1)/'Po',5,1,209/
      DATA (zm(84,i),i=0,1)/208.982404d0, 208.9824308d0/
      DATA (gns(84,i),i=1,1)/2/
      DATA (ab(84,i),i=1,1)/100.d0/
c
      DATA at(85),gel(85),nmn(85),(mn(85,i),i=1,1)/'At',-1,1,210/
      DATA (zm(85,i),i=0,1)/209.987126d0, 209.987148d0/
      DATA (gns(85,i),i=1,1)/11/
      DATA (ab(85,i),i=1,1)/100.d0/
c
      DATA at(86),gel(86),nmn(86),(mn(86,i),i=1,1)/'Rn',1,1,222/
      DATA (zm(86,i),i=0,1)/222.017571d0, 222.0175782d0/
      DATA (gns(86,i),i=1,1)/1/
      DATA (ab(86,i),i=1,1)/100.d0/
c
      DATA at(87),gel(87),nmn(87),(mn(87,i),i=1,1)/'Fr',-1,1,223/
      DATA (zm(87,i),i=0,1)/223.019733d0, 223.0197360d0/
      DATA (gns(87,i),i=1,1)/4/
      DATA (ab(87,i),i=1,1)/100.d0/
c
      DATA at(88),gel(88),nmn(88),(mn(88,i),i=1,1)/'Ra',1,1,226/
      DATA (zm(88,i),i=0,1)/226.025403d0, 226.0254103d0/
      DATA (gns(88,i),i=1,1)/1/
      DATA (ab(88,i),i=1,1)/100.d0/
c
      DATA at(89),gel(89),nmn(89),(mn(89,i),i=1,1)/'Ac',4,1,227/
      DATA (zm(89,i),i=0,1)/227.027750d0, 227.0277523d0/
      DATA (gns(89,i),i=1,1)/4/
      DATA (ab(89,i),i=1,1)/100.d0/
c
      DATA at(90),gel(90),nmn(90),(mn(90,i),i=1,1)/'Th',-1,1,232/
      DATA (zm(90,i),i=0,1)/232.038d0, 232.0380558d0/
      DATA (gns(90,i),i=1,1)/1/
      DATA (ab(90,i),i=1,1)/100.d0/
c
      DATA at(91),gel(91),nmn(91),(mn(91,i),i=1,1)/'Pa',-1,1,231/
      DATA (zm(91,i),i=0,1)/231.03588d0, 231.0358842d0/
      DATA (gns(91,i),i=1,1)/4/
      DATA (ab(91,i),i=1,1)/100.d0/
c
      DATA at(92),gel(92),nmn(92),(mn(92,i),i=1,4)/' U',-1,4,233,234,
     1                                             235,238/
      DATA (zm(92,i),i=0,4)/238.0289d0, 233.0396355d0, 234.0409523d0,
     1    235.0439301d0, 238.0507884d0/
      DATA (gns(92,i),i=1,4)/6,1,8,1/
      DATA (ab(92,i),i=1,4)/0.d0, 0.0055d0, 0.7200d0, 99.2745d0/
c
      DATA at(93),gel(93),nmn(93),(mn(93,i),i=1,1)/'Np',-1,1,237/
      DATA (zm(93,i),i=0,1)/237.0481678d0, 237.0481736d0/
      DATA (gns(93,i),i=1,1)/6/
      DATA (ab(93,i),i=1,1)/100.d0/
c
      DATA at(94),gel(94),nmn(94),(mn(94,i),i=1,1)/'Pu',-1,1,244/
      DATA (zm(94,i),i=0,1)/244.064199d0, 244.064205d0/
      DATA (gns(94,i),i=1,1)/1/
      DATA (ab(94,i),i=1,1)/100.d0/
c
      DATA at(95),gel(95),nmn(95),(mn(95,i),i=1,1)/'Am',-1,1,243/
      DATA (zm(95,i),i=0,1)/243.061375d0, 243.0613815d0/
      DATA (gns(95,i),i=1,1)/6/
      DATA (ab(95,i),i=1,1)/100.d0/
c
      DATA at(96),gel(96),nmn(96),(mn(96,i),i=1,1)/'Cm',-1,1,247/
      DATA (zm(96,i),i=0,1)/247.070347d0, 247.070354d0/
      DATA (gns(96,i),i=1,1)/10/
      DATA (ab(96,i),i=1,1)/100.d0/
c
      DATA at(97),gel(97),nmn(97),(mn(97,i),i=1,1)/'Bk',-1,1,247/
      DATA (zm(97,i),i=0,1)/247.070300d0, 247.070307d0/
      DATA (gns(97,i),i=1,1)/4/
      DATA (ab(97,i),i=1,1)/100.d0/
c
      DATA at(98),gel(98),nmn(98),(mn(98,i),i=1,1)/'Cf',-1,1,251/
      DATA (zm(98,i),i=0,1)/251.079580d0, 251.079589d0/
      DATA (gns(98,i),i=1,1)/2/
      DATA (ab(98,i),i=1,1)/100.d0/
c
      DATA at(99),gel(99),nmn(99),(mn(99,i),i=1,1)/'Es',-1,1,252/
      DATA (zm(99,i),i=0,1)/252.082944d0, 252.082980d0/
      DATA (gns(99,i),i=1,1)/11/
      DATA (ab(99,i),i=1,1)/100.d0/
c
      DATA at(100),gel(100),nmn(100),(mn(100,i),i=1,1)/'Fm',-1,1,257/
      DATA (zm(100,i),i=0,1)/257.095099d0, 257.095106d0/
      DATA (gns(100,i),i=1,1)/10/
      DATA (ab(100,i),i=1,1)/100.d0/
c
      DATA at(101),gel(101),nmn(101),(mn(101,i),i=1,1)/'Md',-1,1,258/
      DATA (zm(101,i),i=0,1)/258.09857d0, 258.098431d0/
      DATA (gns(101,i),i=1,1)/17/
      DATA (ab(101,i),i=1,1)/100.d0/
c
      DATA at(102),gel(102),nmn(102),(mn(102,i),i=1,1)/'No',-1,1,259/
      DATA (zm(102,i),i=0,1)/259.100931d0, 259.101030d0/
      DATA (gns(102,i),i=1,1)/10/
      DATA (ab(102,i),i=1,1)/100.d0/
c
      DATA at(103),gel(103),nmn(103),(mn(103,i),i=1,1)/'Lr',-1,1,260/
      DATA (zm(103,i),i=0,1)/260.105320d0, 260.105510d0/
      DATA (gns(103,i),i=1,1)/-1/
      DATA (ab(103,i),i=1,1)/100.d0/
c
      DATA at(104),gel(104),nmn(104),(mn(104,i),i=1,1)/'Rf',-1,1,261/
      DATA (zm(104,i),i=0,1)/261.10869d0, 261.108770d0/
      DATA (gns(104,i),i=1,1)/-1/
      DATA (ab(104,i),i=1,1)/100.d0/
c
      DATA at(105),gel(105),nmn(105),(mn(105,i),i=1,1)/'Db',-1,1,262/
      DATA (zm(105,i),i=0,1)/262.11376d0, 262.114070d0/
      DATA (gns(105,i),i=1,1)/-1/
      DATA (ab(105,i),i=1,1)/100.d0/
c
      DATA at(106),gel(106),nmn(106),(mn(106,i),i=1,1)/'Sg',-1,1,263/
      DATA (zm(106,i),i=0,1)/263.11822d0, 263.118290d0/
      DATA (gns(106,i),i=1,1)/-1/
      DATA (ab(106,i),i=1,1)/100.d0/
c
      DATA at(107),gel(107),nmn(107),(mn(107,i),i=1,1)/'Bh',-1,1,262/
      DATA (zm(107,i),i=0,1)/262.12293d0, 262.122970d0/
      DATA (gns(107,i),i=1,1)/-1/
      DATA (ab(107,i),i=1,1)/100.d0/
c
      DATA at(108),gel(108),nmn(108),(mn(108,i),i=1,1)/'Hs',-1,1,265/
      DATA (zm(108,i),i=0,1)/265.13016d0, 265.129793d0/
      DATA (gns(108,i),i=1,1)/-1/
      DATA (ab(108,i),i=1,1)/100.d0/
c
      DATA at(109),gel(109),nmn(109),(mn(109,i),i=1,1)/'Mt',-1,1,266/
      DATA (zm(109,i),i=0,1)/266.13764d0, 266.137370d0/
      DATA (gns(109,i),i=1,1)/-1/
      DATA (ab(109,i),i=1,1)/100.d0/
c
      IF((IAN.LT.0).OR.(IAN.GT.109)) THEN
          MASS= 0.d0
          NAME= 'XX'
          IMN= 0
          WRITE(6,601) IAN
          RETURN
        ELSE
          NAME= AT(IAN)
        ENDIF
      IF((IAN.EQ.1).AND.(IMN.GT.1)) THEN
c** Special case: insert common name for deuterium or tritium
          IF(IMN.EQ.2) NAME=' D'
          IF(IMN.EQ.3) NAME=' T'
          ENDIF
      IF((IAN.EQ.0).AND.(IMN.GT.1)) THEN
          IF(IMN.EQ.2) NAME=' d'
          IF(IMN.EQ.3) NAME=' t'
          ENDIF
      GELGS= GEL(IAN)
      MASS= -1.d0
      DGNS= -1
      ABUND = -1.d0
      DO  I= 1,NMN(IAN)
          if(i.gt.15)  write(6,606) ian,imn,nmn(ian)
          IF(IMN.EQ.MN(IAN,I)) THEN
              MASS= ZM(IAN,I)
              DGNS= gns(IAN,I)
              ABUND = AB(IAN,I)
              ENDIF
          ENDDO
      IF(MASS.LT.0.d0) THEN
          MASS= ZM(IAN,0)
          IF(IMN.NE.0) WRITE(6,602) AT(IAN),IMN
          IMN= 0
          ENDIF
      RETURN
  601 FORMAT(' *** MASSES Data base does not include Atomic Number=',i4)
  602 FORMAT(' *** MASSES Data base does not include ',A2,'(',i3,
     1 '), so use average atomic mass.')
  606  format(/' *** ERROR *** called MASSES for atom with  AN=',I4,
     1  '  MN=',I4,'n(MN)=',I4)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE WGHT(NGP)
c** Subroutine to generate points (XG & X2) and weights (WG & W2) for
c  both regular (XG & WG) and singular  1/sqrt(1-X) (X2 & W2) 
c  Gaussian integration for  NGP=8 and 16 (for both) and for NGP=7 & 32
c  (for regular Gaussian only)
c-----------------------------------------------------------------------
      INTEGER I,J,NGP,NGPH
      REAL*8 x7(4),w7(4),aa(4),bb(4),A(8),B(8),XX(16),WX(16)
c** Common block for quadrature weights & points
      REAL*8 XG(32),WG(32),X2(16),W2(16)
      COMMON /GWGHT/XG,WG,X2,W2
c
      data x7/0.949107912342759d0, 0.741531185599394d0,
     1        0.405845151377397d0, 0.d0/,
     2     w7/0.129484966168870d0, 0.279705391489277d0,
     3        0.381830550505119d0, 0.417959183673469d0/
      data aa/0.960289856497536d0, 0.796666477413627d0,
     1        0.525532409916329d0, 0.183434642495650d0/,
     2     bb/0.101228536290376d0, 0.222381034453374d0,
     3        0.313706645877887d0, 0.362683783378362d0/
      DATA A/0.989400934991649932596154D0,0.944575023073232576077988D0,
     1       0.865631202387831743880468D0,0.755404408355003033895101D0,
     2       0.617876244402643748446672D0,0.458016777657227386342420D0,
     3       0.281603550779258913230461D0,0.095012509837637440185319D0/,
     4     B/0.02715245941175409485178D0,0.06225352393864789286284D0,
     5       0.09515851168249278480993D0,0.12462897125553387205248D0,
     6       0.14959598881657673208150D0,0.16915651939500253818931D0,
     7       0.18260341504492358886676D0,0.18945061045506849628540D0/
      DATA XX/0.997263861849481563544981D0,0.985611511545268335400175D0,
     1        0.964762255587506430773812D0,0.934906075937739689170919D0,
     2        0.896321155766052123965307D0,0.849367613732569970133693D0,
     3        0.794483795967942406963097D0,0.732182118740289680387427D0,
     4        0.663044266930215200975115D0,0.587715757240762329040746D0,
     5        0.506899908932229390023747D0,0.421351276130635345364120D0,
     6        0.331868602282127649779917D0,0.239287362252137074544603D0,
     7        0.144471961582796493485186D0,0.048307665687738316234813D0/
     8    ,WX/0.00701861000947009660041D0,0.01627439473090567060517D0,
     9        0.02539206530926205945575D0,0.03427386291302143310269D0,
     1        0.04283589802222668065688D0,0.05099805926237617619616D0,
     2        0.05868409347853554714528D0,0.06582222277636184683765D0,
     3        0.07234579410884850622540D0,0.07819389578707030647174D0,
     4        0.08331192422694675522220D0,0.08765209300440381114277D0,
     5        0.09117387869576388471287D0,0.09384439908080456563918D0,
     6        0.09563872007927485941908D0,0.09654008851472780056676D0/
c** For simple Gaussian case
      NGPH= (NGP+1)/2
      J= NGP+1
      if(ngp.eq.7) then
          DO  I= 1,4
              J= J-1
              XG(I)= -x7(I)
              WG(I)= w7(I)
              XG(J)= x7(I)
              WG(J)= w7(I)
              ENDDO
c%%   WRITE(6,601) NGP,(XG(I),WG(I),I= 1,NGP)
c*** For  1/SQRT(E-V)  case
c%%   WRITE(6,603) (XX(I),WX(I),I= 1,NGP)
          write(6,606) ngp
          go to 99
          endif
      if(ngp.eq.8) then
          DO  I= 1,4
              J= J-1
              XG(I)= -aa(I)
              WG(I)= bb(I)
              XG(J)= aa(I)
              WG(J)= bb(I)
              ENDDO
c%%   WRITE(6,601) NGP,(XG(I),WG(I),I= 1,NGP)
c*** For  1/SQRT(E-V)  case
c%%   WRITE(6,603) (XX(I),WX(I),I= 1,NGP)
          DO    I= 1,8
              X2(I)= 1.d0-A(I)**2
c** Include SQRT(1-X) factor in weight, not integrand.
              W2(I)= 2.d0*B(I)*A(I)
              ENDDO
c%%       WRITE(6,604) NGP,(X2(I),W2(I),I=1,NGP)
          go to 99
          ENDIF
      if(ngp.eq.16) then
          DO  I=1,8
              J=J-1
              XG(I)=-A(I)
              WG(I)=B(I)
              XG(J)=A(I)
              WG(J)=B(I)
              ENDDO
c%%   WRITE(6,601) NGP,(XG(I),WG(I),I=1,NGP)
c*** For  1/SQRT(E-V)  case
c%%   WRITE(6,603) (XX(I),WX(I),I=1,NGP)
          DO  I=1,16
              X2(I)=1.d0-XX(I)**2
c** Include SQRT(1-X) factor in weight, not integrand.
              W2(I)=2.d0*WX(I)*XX(I)
              ENDDO
c%%       WRITE(6,604) NGP,(X2(I),W2(I),I=1,NGP)
          go to 99
          endif
      if(ngp.eq.32) then
          DO  I=1,NGPH
              J=J-1
              XG(I)=-XX(I)
              WG(I)=WX(I)
              XG(J)=XX(I)
              WG(J)=WX(I)
              ENDDO
c%%   WRITE(6,601) NGP,(XG(I),WG(I),I=1,NGP)
          write(6,606) ngp
          go to 99
          endif
      write(6,608) ngp
   99 RETURN
c 601 FORMAT(//2X,'Points and weights for simple Gaussian',i3,
c    1 '-point quadrature'/(1x,2(F25.20,F23.20)))
c 603 FORMAT(/' Simple Gaussian points and weights used for generating G
c    1auss-Mehler points and weights:'/(5X,2(F25.20,F23.20)))
c 604     FORMAT(/'  Points and weights for singular integrand',i3,
c    1     '-point quadratures'/(1x,2(F25.20,F23.20)))
  606 format(/'  NOTE: cannot generate singular integrand Gaussian point
     1s and weights for   NGP =',i3)
  608 format(/' NOTE:  cannot generate Gaussian points & weights for',
     1   '   NGP =',i3)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
