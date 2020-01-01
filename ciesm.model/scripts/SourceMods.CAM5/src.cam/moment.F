      SUBROUTINE MOMENT(DUDT,DVDT,dpdx,dpdy,
     1                  DP,ED,EU,MC,MD,MU,MB,
     1                  P,QD,QU,QHAT,SD,SU,SHAT,T,U,V,Z,ZF,
     2                  DSUBCLD,JB,JD,JT,
     3                  MSG,DT,GRAV,CPRES,RGAS,ILEV,IL1G,IL2G,ILG)

C     * JULY 17/92. - GUANG JUN ZHANG, M.LAZARE.
C
C     * PERFORMS MOMENTUM MIXING DUE TO CUMULUS PARAMETRIZATION.
C
!songxl 20060901---------------------
       use shr_kind_mod, only: r8=>shr_kind_r8   
       implicit none
       integer, parameter :: NBF = 10 
!songxl      PARAMETER(NBF=10)
!songxl      PARAMETER(JLG=128, JLEV=18)
!songxl 20060901---------------------
     
      real(r8) DUDT(ILG,ILEV),  DVDT(ILG,ILEV), 
     1     dpdx(ILG,ILEV),  dpdy(ILG,ILEV),
     1     ALPHA(ILG,ILEV), DP(ILG,ILEV),
     2     ED(ILG,ILEV),    EU(ILG,ILEV),    MB(ILG),
     3     MC(ILG,ILEV),    MD(ILG,ILEV),    MU(ILG,ILEV),
     4     P(ILG,ILEV),     QD(ILG,ILEV),    QU(ILG,ILEV),
     5     QHAT(ILG,ILEV),  SD(ILG,ILEV),    SU(ILG,ILEV),
     6     SHAT(ILG,ILEV),  T(ILG,ILEV),     U(ILG,ILEV), 
     7     V(ILG,ILEV),     Z(ILG,ILEV),     ZF(ILG,ILEV+1)    !songxl 
!songxl     7     V(ILG,ILEV),     Z(ILG,ILEV),     ZF(ILG,ILEV) 
C    
      real(r8) DSUBCLD(ILG)
!songxl 20060901-----------------
      real(r8) DT, GRAV, CPRES, RGAS
      integer MSG, ILEV, IL1G, IL2G, ILG
!songxl 20060901-----------------
C
      INTEGER   JB(ILG),    JD(ILG),         JT(ILG) 
C
C     * INTERNAL WORK FIELDS.
C
      real(r8) AC(ILG,ILEV),    AD(ILG,ILEV),    AU(ILG,ILEV),
     1     ACFL(ILG,ILEV),  ADFL(ILG,ILEV),  AUFL(ILG,ILEV),
     2     B1(ILG,ILEV),    B1FL(ILG,ILEV),  BD(ILG,ILEV),
     3     BU(ILG,ILEV),    D0(ILG,ILEV),    D0HAT(ILG,ILEV), 
     4     DELPX(ILG,ILEV), DELPY(ILG,ILEV), DZ(ILG,ILEV),
     5     DZF(ILG,ILEV),   EC(ILG,ILEV),    E(ILG,ILEV), 
     6     EHAT(ILG,ILEV),  RHO(ILG,ILEV),   RHOHAT(ILG,ILEV),
     7     UHAT(ILG,ILEV),  VHAT(ILG,ILEV),  UC(ILG,ILEV), 
     8     VC(ILG,ILEV),    W0(ILG,ILEV),    W1(ILG,ILEV),
     9     W0FL(ILG,ILEV),  W1FL(ILG,ILEV) 
C
      real(r8) COSA(ILG),    SINA(ILG),    CCX(ILG),     CCY(ILG),
     1     UMN(ILG),     VMN(ILG),     DEP(ILG)  
C
      real(r8) FX(ILG,ILEV,NBF), DX(ILG,ILEV,NBF),
     1     AA(ILG,ILEV,NBF), BB(ILG,ILEV,NBF),  CC(ILG,ILEV,NBF), 
     2     FY(ILG,ILEV,NBF), DY(ILG,ILEV,NBF) 
C  
!songxl      COMMON/CONST/ALFA1(NBF),ALFA2(NBF),BSJ0(NBF),BSJ1(NBF), 
!songxl     1             BSJ0CHI,FACTOR 
!songxl      COMMON/TAU/TAU(NBF) 
!songxl 20060901-----------------
      real(r8) ALFA1(NBF),ALFA2(NBF),BSJ0(NBF),BSJ1(NBF),BSJ0CHI,FACTOR,TAU(NBF)
!songxl 20060901-----------------
      DATA ALFA1/4.33675e-05_r8, 7.88429e-05_r8, 1.13065e-04_r8, 1.45836e-04_r8,
     *           1.76823e-04_r8, 2.05660e-04_r8, 2.32101e-04_r8, 2.55829e-04_r8,
     *           2.76672e-04_r8, 2.94347e-04_r8/
      DATA ALFA2/-8.64731e-05_r8,-1.56092e-04_r8,-2.21315e-04_r8,-2.80995e-04_r8,
     *           -3.33799e-04_r8,-3.78455e-04_r8,-4.14020e-04_r8,-4.39600e-04_r8,
     *           -4.54652e-04_r8,-4.58761e-04_r8/
      DATA BSJ0/-4.02759e-01_r8, 3.00116e-01_r8, -2.49705e-01_r8, 2.18359e-01_r8,
     *          -1.96465e-01_r8, 1.80062e-01_r8, -1.67183e-01_r8, 1.56722e-01_r8,
     *          -1.48011e-01_r8, 1.40605e-01_r8/
      DATA BSJ1/1.14193e-01_r8, 2.05841e-01_r8, 2.91209e-01_r8, 3.68619e-01_r8,
     *          4.36201e-01_r8, 4.92233e-01_r8, 5.35486e-01_r8, 5.64886e-01_r8,
     *          5.79879e-01_r8, 5.80195e-01_r8/
      DATA BSJ0CHI/-0.402759_r8/
      DATA FACTOR/0.179503_r8/
      DATA TAU/ 3.83170e+00_r8, 7.01560e+00_r8, 1.01735e+01_r8, 1.33237e+01_r8,
     *          1.64705e+01_r8, 1.96120e+01_r8, 2.27560e+01_r8, 2.58980e+01_r8,
     *          2.90480e+01_r8, 3.21926e+01_r8/
!songxl 20060901------------------------
      integer J, IL, N ,jj
      real(r8) RC, RMAX, CHI, SHAPE, RNU, WINDMAG
!songxl 20060901-----------------------

C----------------------------------------------------------------------
C***********************************************************************
CCCCC INITIALIZE RELEVANT INTERFACIAL AND MIDLAYER VARIABLES     CCCCC
C***********************************************************************
      RC=3000._r8         ! cloud radius (m)
      RMAX=50000._r8      ! distance where perturb. vanishes (m)  
      CHI=3.8317_r8/RC 
      SHAPE=0.8_r8 
C
      DO 5 J=MSG+1,ILEV 
      DO 5 IL=IL1G,IL2G 
        BU(IL,J)=0._r8 
        ALPHA(IL,J)=0.5_r8 
        BD(IL,J)=0._r8 
        DELPX(IL,J)=0._r8 
        DELPY(IL,J)=0._r8 
        EC(IL,J)=EU(IL,J)+ED(IL,J)     ! unit: 1/s
        IF(T(IL,J).GT.0._r8)        THEN
            RHO(IL,J)=100._r8*P(IL,J)/(RGAS*T(IL,J)) 
        ENDIF 
    5 CONTINUE
C
      DO 10 N=1,NBF 
      DO 10 J=MSG+1,ILEV
      DO 10 IL=IL1G,IL2G
        FX(IL,J,N)=0._r8 
        FY(IL,J,N)=0._r8 
   10 CONTINUE
C
      DO 15 J=MSG+2,ILEV
      DO 15 IL=IL1G,IL2G
        BU(IL,J)=(   SU(IL,J)-SHAT(IL,J)+0.608_r8*( QU(IL,J) 
     1          *(SU(IL,J)-GRAV/CPRES*ZF(IL,J))-QHAT(IL,J)
     2          *(SHAT(IL,J)-GRAV/CPRES*ZF(IL,J)) )  )
     3          /( (SHAT(IL,J)-GRAV/CPRES*ZF(IL,J)) 
     4          *(1._r8+.608_r8*QHAT(IL,J)) )*GRAV
        BD(IL,J)=(   SD(IL,J)-SHAT(IL,J)+0.608_r8*( QD(IL,J) 
     1          *(SD(IL,J)-GRAV/CPRES*ZF(IL,J))-QHAT(IL,J)
     2          *(SHAT(IL,J)-GRAV/CPRES*ZF(IL,J)) )  )
     3          /( (SHAT(IL,J)-GRAV/CPRES*ZF(IL,J)) 
     4          *(1._r8+.608_r8*QHAT(IL,J)) )*GRAV
   15 CONTINUE
C
      DO 20 IL=IL1G,IL2G
        UMN(IL)=0._r8
        VMN(IL)=0._r8
C       DEP(IL)=0.
        COSA(IL)=0._r8
        SINA(IL)=0._r8  
        CCX(IL)=0._r8
        CCY(IL)=0._r8
   20 CONTINUE
C                                                                
      DO 25 J=1,ILEV
      DO 25 IL=IL1G,IL2G
        IF(J.GE.JT(IL) .AND. J.LT.JB(IL))     THEN
          UMN(IL)=UMN(IL)+U(IL,J)*P(IL,J)/T(IL,J)*
     1            (ZF(IL,J)-ZF(IL,J+1))                
          VMN(IL)=VMN(IL)+V(IL,J)*P(IL,J)/T(IL,J)*
     1            (ZF(IL,J)-ZF(IL,J+1))                
C         DEP(IL)=DEP(IL)+P(IL,J)/T(IL,J)*(ZF(IL,J)-ZF(IL,J+1))   
        ENDIF                      
   25 CONTINUE  
C
      DO 30 IL=IL1G,IL2G
        WINDMAG=SQRT(UMN(IL)**2+VMN(IL)**2)
C       IF(DEP(IL).NE.0. .AND. WINDMAG.NE.0.)   THEN
        IF(WINDMAG.NE.0._r8)                       THEN
          COSA(IL)=UMN(IL)/WINDMAG
          SINA(IL)=VMN(IL)/WINDMAG
C         CCX(IL)=0.*UMN(IL)/DEP(IL)
C         CCY(IL)=0.*VMN(IL)/DEP(IL)
        ENDIF
   30 CONTINUE 
C
      DO 35 J=MSG+1,ILEV
      DO 35 IL=IL1G,IL2G
        IF(J.GT.MSG+1)                          THEN
         RHOHAT(IL,J)=ALPHA(IL,J)*RHO(IL,J-1)+(1._r8-ALPHA(IL,J))*RHO(IL,J) 
         DZF(IL,J)=Z(IL,J-1)-Z(IL,J) 
         UHAT(IL,J)=ALPHA(IL,J)*U(IL,J-1)+(1._r8-ALPHA(IL,J))*U(IL,J) 
         VHAT(IL,J)=ALPHA(IL,J)*V(IL,J-1)+(1._r8-ALPHA(IL,J))*V(IL,J) 
        ELSE
         RHOHAT(IL,J)=RHO(IL,J)
         DZF(IL,J)=ZF(IL,J)-Z(IL,J)
         UHAT(IL,J)=U(IL,J)
         VHAT(IL,J)=V(IL,J)
        ENDIF 
        UC(IL,J)=UHAT(IL,J) 
        VC(IL,J)=VHAT(IL,J) 
   35 CONTINUE
C
      DO 40 J=MSG+1,ILEV
      DO 40 IL=IL1G,IL2G
        IF(J.LT.ILEV)                       THEN
           DZ(IL,J)=ZF(IL,J)-ZF(IL,J+1) 
        ELSE
c          DZ(IL,ILEV)=ZF(IL,ILEV)
           DZ(IL,ILEV)=DP(IL,ILEV)*100._r8/(RHO(IL,ILEV)*GRAV) !m
        ENDIF 
   40 CONTINUE
C     ******************************************************************
C     AU, AD,AC ARE ACTUAL FRACTIONAL CLOUD AREA TIMES GRAV 
C     ******************************************************************
      DO 75 J=MSG+1,ILEV
      DO 75 IL=IL1G,IL2G
        IF(J.GE.JT(IL))                     THEN 
           AU(IL,J)=MU(IL,JB(IL))*100._r8/GRAV/RHOHAT(IL,JB(IL)) 
        ELSE
           AU(IL,J)=0._r8
        ENDIF 
        IF(J.GE.JD(IL))                     THEN
           AD(IL,J)=-MD(IL,JD(IL))*100._r8/GRAV/RHOHAT(IL,JD(IL)) 
        ELSE
           AD(IL,J)=0._r8
        ENDIF 
        AC(IL,J)=AU(IL,J)+AD(IL,J)
        IF(AC(IL,J).GT.0._r8)                  THEN 
           W0(IL,J)=MC(IL,J)*100._r8/GRAV/(RHOHAT(IL,J)*AC(IL,J))
        ELSE
           W0(IL,J)=0._r8
        ENDIF 
   75 CONTINUE
C
      DO 80 J=MSG+1,ILEV
      DO 80 IL=IL1G,IL2G
        IF(J.LT.ILEV)                       THEN
           ACFL(IL,J)=ALPHA(IL,J)*AC(IL,J)+(1._r8-ALPHA(IL,J))*AC(IL,J+1)
           AUFL(IL,J)=ALPHA(IL,J)*AU(IL,J)+(1._r8-ALPHA(IL,J))*AU(IL,J+1)
           ADFL(IL,J)=ALPHA(IL,J)*AD(IL,J)+(1._r8-ALPHA(IL,J))*AD(IL,J+1)
           W0FL(IL,J)=ALPHA(IL,J)*W0(IL,J)+(1._r8-ALPHA(IL,J))*W0(IL,J+1)
        ELSE
           ACFL(IL,J)=ALPHA(IL,J)*AC(IL,J)
           AUFL(IL,J)=ALPHA(IL,J)*AU(IL,J)
           ADFL(IL,J)=ALPHA(IL,J)*AD(IL,J)
           W0FL(IL,J)=ALPHA(IL,J)*W0(IL,J)
        ENDIF 
   80 CONTINUE
C
      DO 250 J=MSG+1,ILEV 
      DO 250 IL=IL1G,IL2G 
        IF(J.LT.JB(IL) .AND. J.GE.JT(IL) .AND. ACFL(IL,J).GT.0._r8) THEN 
           D0(IL,J)=( MC(IL,J+1)-MC(IL,J) )*100._r8/GRAV 
     1              /(RHO(IL,J)*DZ(IL,J)*ACFL(IL,J))
        ELSE
           IF(J.EQ.JB(IL) .AND. ACFL(IL,J).GT.0._r8 )               THEN 
              D0(IL,J)= -MC(IL,J)*100._r8/GRAV
     1              /(RHO(IL,J)*DZ(IL,J)*ACFL(IL,J)) 
           ELSE 
              D0(IL,J)=0._r8 
           ENDIF
        ENDIF 
        IF(J.LT.JB(IL) .AND. J.GE.JT(IL) .AND. ACFL(IL,J).GT.0._r8 
     1     .AND. ADFL(IL,J).GT.0._r8 .AND. AUFL(IL,J).GT.0._r8)        THEN 
           E(IL,J)=1._r8/CHI**2*ADFL(IL,J)/(ACFL(IL,J)*RHO(IL,J))
     1             *( 1._r8/(AUFL(IL,J)*DZ(IL,J))
     2             *(MU(IL,J)-MU(IL,J+1))*100._r8/GRAV -
     2                1._r8/(ADFL(IL,J)*DZ(IL,J))
     3             *(MD(IL,J)-MD(IL,J+1))*100._r8/GRAV ) 
     4       *SHAPE/FACTOR
        ELSE
           E(IL,J)=0._r8 
        ENDIF 
  250 CONTINUE
C
      DO 275 J=MSG+1,ILEV 
      DO 275 IL=IL1G,IL2G 
        IF(J.GT.MSG+1)                                           THEN
          D0HAT(IL,J)=ALPHA(IL,J)*D0(IL,J-1)+(1._r8-ALPHA(IL,J))*D0(IL,J)
          EHAT(IL,J)=ALPHA(IL,J)*E(IL,J-1)+(1._r8-ALPHA(IL,J))*E(IL,J) 
        ELSE
          D0HAT(IL,J)=(1._r8-ALPHA(IL,J))*D0(IL,J) 
          EHAT(IL,J)=(1._r8-ALPHA(IL,J))*E(IL,J) 
        ENDIF 
  275 CONTINUE
C     **********************************************
C     CALCULATE FIRST HARMONICS IN THERMODYNAMICS 
C     **********************************************
      DO 325 J=MSG+1,ILEV 
      DO 325 IL=IL1G,IL2G 
        IF( J.GE.JD(IL) .AND. J.LE.JB(IL) .AND. AU(IL,J).GT.0._r8
     1      .AND. AD(IL,J).GT.0._r8 .AND. AC(IL,J).GT.0._r8 )          THEN
          W1(IL,J)=SHAPE/FACTOR*(MU(IL,J)/AU(IL,J)-MD(IL,J)/AD(IL,J)) 
     1             *100._r8/GRAV*AD(IL,J)/(AC(IL,J)*RHOHAT(IL,J))
          B1(IL,J)=SHAPE/FACTOR*(BU(IL,J)-BD(IL,J))*AD(IL,J)/AC(IL,J) 
        ELSE
          W1(IL,J)=0._r8 
          B1(IL,J)=0._r8 
        ENDIF 
  325 CONTINUE
C
      DO 350 J=MSG+1,ILEV 
      DO 350 IL=IL1G,IL2G 
        IF(J.LT.ILEV)                                      THEN 
          W1FL(IL,J)=ALPHA(IL,J)*W1(IL,J)+(1._r8-ALPHA(IL,J))*W1(IL,J+1) 
          B1FL(IL,J)=ALPHA(IL,J)*B1(IL,J)+(1._r8-ALPHA(IL,J))*B1(IL,J+1) 
        ELSE
          W1FL(IL,J)=ALPHA(IL,J)*W1(IL,J) 
          B1FL(IL,J)=ALPHA(IL,J)*B1(IL,J) 
        ENDIF 
  350 CONTINUE
C
      DO 500 N=1,NBF
      DO 500 J=MSG+1,ILEV 
      DO 500 IL=IL1G,IL2G 
        IF(J.LT.JB(IL) .AND. J.GE.JT(IL) .AND. MC(IL,J).GT.0._r8 
     1     .AND. MC(IL,J+1).GT.0._r8)                          THEN
           FX(IL,J,N)=2._r8*RHO(IL,J)/BSJ0(N)**2
     1         *( D0(IL,J)*E(IL,J)*CHI**2*COSA(IL)*ALFA1(N) 
     2         +2._r8*W0FL(IL,J)*RC/RMAX**2*BSJ1(N)
     3         *( (EHAT(IL,J)-EHAT(IL,J+1))/DZ(IL,J)*CHI*COSA(IL) 
     4         *BSJ0CHI+(UHAT(IL,J)-UHAT(IL,J+1))/DZ(IL,J) )
     5         -(D0HAT(IL,J)-D0HAT(IL,J+1))/DZ(IL,J)*W1FL(IL,J)*COSA(IL)
     6         *(ALFA2(N)-ALFA1(N))    -2._r8*( (W0(IL,J)-W0(IL,J+1))
     7         *(W1(IL,J)-W1(IL,J+1))/DZ(IL,J)**2 
     8         -W0FL(IL,J)*W1FL(IL,J) 
     9         *( LOG(RHO(IL,J-1))/(DZ(IL,J)*DZF(IL,J))
     A         -(1._r8/DZF(IL,J)+1._r8/DZF(IL,J+1))/DZ(IL,J)
     B         *LOG(RHO(IL,J)) 
     C         +LOG(RHO(IL,J+1))/(DZ(IL,J)*DZF(IL,J+1)) ) )
     D         *COSA(IL)*ALFA1(N) 
     E         +(RHOHAT(IL,J)*B1(IL,J)-RHOHAT(IL,J+1)*B1(IL,J+1)) 
     F         /(DZ(IL,J)*RHO(IL,J))*COSA(IL)*ALFA1(N) )
           FY(IL,J,N)=2._r8*RHO(IL,J)/BSJ0(N)**2
     1         *( D0(IL,J)*E(IL,J)*CHI**2*SINA(IL)*ALFA1(N) 
     2         +2._r8*W0FL(IL,J)*RC/RMAX**2*BSJ1(N)
     3         *( (EHAT(IL,J)-EHAT(IL,J+1))/DZ(IL,J)*CHI*SINA(IL) 
     4         *BSJ0CHI+(VHAT(IL,J)-VHAT(IL,J+1))/DZ(IL,J) )
     5         -(D0HAT(IL,J)-D0HAT(IL,J+1))/DZ(IL,J)*W1FL(IL,J)*SINA(IL)
     6         *(ALFA2(N)-ALFA1(N))    -2._r8*( (W0(IL,J)-W0(IL,J+1))
     7         *(W1(IL,J)-W1(IL,J+1))/DZ(IL,J)**2 
     8         -W0FL(IL,J)*W1FL(IL,J) 
     9         *( LOG(RHO(IL,J-1))/(DZ(IL,J)*DZF(IL,J))
     A         -(1._r8/DZF(IL,J)+1._r8/DZF(IL,J+1))/DZ(IL,J)
     B         *LOG(RHO(IL,J)) 
     C         +LOG(RHO(IL,J+1))/(DZ(IL,J)*DZF(IL,J+1)) ) )
     D         *SINA(IL)*ALFA1(N) 
     E         +(RHOHAT(IL,J)*B1(IL,J)-RHOHAT(IL,J+1)*B1(IL,J+1)) 
     F         /(DZ(IL,J)*RHO(IL,J))*SINA(IL)*ALFA1(N) )
        ENDIF 
  500 CONTINUE
C  
      DO 525 N=1,NBF
      DO 525 IL=IL1G,IL2G 
        AA(IL,MSG+1,N)=0._r8 
        BB(IL,MSG+1,N)=1._r8/ DZF(IL,MSG+2)
        CC(IL,MSG+1,N)=-1._r8/ DZF(IL,MSG+2) 
  
        AA(IL,ILEV,N)=1._r8/ DZF(IL,ILEV)
        BB(IL,ILEV,N)=-1._r8/ DZF(IL,ILEV) 
        CC(IL,ILEV,N)=0._r8
  
        DX(IL,MSG+1,N)=B1(IL,MSG+2)*ALFA1(N)*2._r8
     1               *COSA(IL)/BSJ0(N)**2 
        DX(IL,ILEV,N)=B1(IL,ILEV)*ALFA1(N)*2._r8
     1              *COSA(IL)/BSJ0(N)**2
  
        DY(IL,MSG+1,N)=B1(IL,MSG+2)*ALFA1(N)*2._r8
     1               *SINA(IL)/BSJ0(N)**2 
        DY(IL,ILEV,N)=B1(IL,ILEV)*ALFA1(N)*2._r8
     1              *SINA(IL)/BSJ0(N)**2
  525 CONTINUE
C
      DO 550 N=1,NBF
      DO 550 J=MSG+2,ILEV-1 
      DO 550 IL=IL1G,IL2G 
        AA(IL,J,N)=1._r8/( DZ(IL,J)*DZF(IL,J) )
        BB(IL,J,N)=-( 1._r8/DZ(IL,J)*(1._r8/DZF(IL,J)+1._r8/DZF(IL,J+1)) 
     1           +(TAU(N)/RMAX)**2 )
        CC(IL,J,N)=1._r8/( DZ(IL,J)*DZF(IL,J+1) )
  
        DX(IL,J,N)=FX(IL,J,N) 
        DY(IL,J,N)=FY(IL,J,N) 
  550 CONTINUE
C
      DO 575 N=1,NBF
      DO 575 IL=IL1G,IL2G 
        CC(IL,MSG+1,N)=CC(IL,MSG+1,N)/BB(IL,MSG+1,N)
        DX(IL,MSG+1,N)=DX(IL,MSG+1,N)/BB(IL,MSG+1,N)
        DY(IL,MSG+1,N)=DY(IL,MSG+1,N)/BB(IL,MSG+1,N)
  575 CONTINUE
C
      DO 600 N=1,NBF
      DO 600 J=MSG+2,ILEV 
      DO 600 IL=IL1G,IL2G 
        CC(IL,J,N)=CC(IL,J,N)/(BB(IL,J,N)-AA(IL,J,N)*CC(IL,J-1,N))
        DX(IL,J,N)=(DX(IL,J,N)-AA(IL,J,N)*DX(IL,J-1,N)) 
     1           /(BB(IL,J,N)-AA(IL,J,N)*CC(IL,J-1,N))
        DY(IL,J,N)=(DY(IL,J,N)-AA(IL,J,N)*DY(IL,J-1,N)) 
     1           /(BB(IL,J,N)-AA(IL,J,N)*CC(IL,J-1,N))
  600 CONTINUE
C
      DO 650 N=1,NBF
      DO 650 J=ILEV-1,MSG+1,-1
      DO 650 IL=IL1G,IL2G 
        DX(IL,J,N)=DX(IL,J,N)-CC(IL,J,N)*DX(IL,J+1,N) 
        DY(IL,J,N)=DY(IL,J,N)-CC(IL,J,N)*DY(IL,J+1,N) 
  650 CONTINUE
C
      DO 700 N=1,NBF
      DO 700 J=MSG+1,ILEV 
      DO 700 IL=IL1G,IL2G 
        DELPX(IL,J)=DELPX(IL,J)+DX(IL,J,N)*BSJ1(N)   !kg/m/s**2
        DELPY(IL,J)=DELPY(IL,J)+DY(IL,J,N)*BSJ1(N)   !kg/m/s**2
  700 CONTINUE
C
      DO 850 J=MSG+1,ILEV 
      DO 850 IL=IL1G,IL2G 
        DELPX(IL,J)=DELPX(IL,J)/RC             !kg/m**2/s**2
        DELPY(IL,J)=DELPY(IL,J)/RC             !kg/m**2/s**2
c to get cloud-scale pressure gradient by multiplying cloud fraction  
        DELPX(IL,J)=ACFL(IL,J)*DELPX(IL,J)/RHO(IL,J)  !m/s/s 
        DELPY(IL,J)=ACFL(IL,J)*DELPY(IL,J)/RHO(IL,J)  !m/s/s
        dpdx(IL,J)=DELPX(IL,J)
        dpdy(IL,J)=DELPY(IL,J)
  850 CONTINUE
C     ************************************
C     CALCULATE THE CLOUD MEAN WIND 
C     ************************************
      DO 875 J=ILEV-1,MSG+1,-1
      DO 875 IL=IL1G,IL2G 
        IF(MC(IL,J).GT.0._r8 .and. MC(IL,J+1) .GT. 0._r8 
     1         .AND. J.GT.JT(IL) .AND.J.LT.JB(IL))      THEN 
       UC(IL,J)=UC(IL,J+1) + RHO(IL,J)*DZ(IL,J)
     1                /((MC(IL,J)+MC(IL,J+1))*0.5_r8*100._r8/GRAV)
     1                *( EC(IL,J)*(U(IL,J)-UC(IL,J+1))-DELPX(IL,J) )
       VC(IL,J)=VC(IL,J+1) + RHO(IL,J)*DZ(IL,J)
     1                /((MC(IL,J)+MC(IL,J+1))*0.5_r8*100._r8/GRAV)
     1                *( EC(IL,J)*(V(IL,J)-VC(IL,J+1))-DELPY(IL,J) )
        ENDIF 
  875 CONTINUE 
C 
c     RNU=1.
      RNU=0._r8
      DO 950 J=MSG+1,ILEV
      DO 950 IL=IL1G,IL2G
        IF( J.GE.JT(IL) .AND. J.LE.JB(IL) )                 THEN
           UHAT(IL,J)=UHAT(IL,J)+RNU*ALPHA(IL,J)*DT*DUDT(IL,J-1)
           VHAT(IL,J)=VHAT(IL,J)+RNU*ALPHA(IL,J)*DT*DVDT(IL,J-1)
           IF(J.EQ.JT(IL))                                     THEN
             DUDT(IL,J)=1._r8/(1._r8+RNU*ALPHA(IL,J+1)*MC(IL,J+1)*DT/DP(IL,J))*
     1                  MC(IL,J+1)*(UC(IL,J+1)-UHAT(IL,J+1))/DP(IL,J) 
             DVDT(IL,J)=1._r8/(1._r8+RNU*ALPHA(IL,J+1)*MC(IL,J+1)*DT/DP(IL,J))*
     1                  MC(IL,J+1)*(VC(IL,J+1)-VHAT(IL,J+1))/DP(IL,J) 
           ELSE IF(J.LT.JB(IL))                                THEN
             DUDT(IL,J)=1._r8/(1._r8+RNU*ALPHA(IL,J+1)*MC(IL,J+1)*DT/DP(IL,J))
     1                 *(MC(IL,J+1)*(UC(IL,J+1)-UHAT(IL,J+1)) 
     2                 -MC(IL,J)*(UC(IL,J)-UHAT(IL,J)))/DP(IL,J)
             DVDT(IL,J)=1._r8/(1._r8+RNU*ALPHA(IL,J+1)*MC(IL,J+1)*DT/DP(IL,J))
     1                 *(MC(IL,J+1)*(VC(IL,J+1)-VHAT(IL,J+1)) 
     2                 -MC(IL,J)*(VC(IL,J)-VHAT(IL,J)))/DP(IL,J)
           ELSE 
             DUDT(IL,J)=-1._r8/DSUBCLD(IL)*MC(IL,J)*(UC(IL,J)-UHAT(IL,J)) 
             DVDT(IL,J)=-1._r8/DSUBCLD(IL)*MC(IL,J)*(VC(IL,J)-VHAT(IL,J)) 
           ENDIF
        ENDIF 
c dudt and dvdt (m/s/s)

      if (abs(dudt(il,j)).gt.5.0e-2_r8.or.abs(dvdt(il,j)).gt.5.0e-2_r8) then
       print*,'moment',il,j,dt,jt(il),jb(il)
     $       ,dudt(il,j),dvdt(il,j),dp(il,j)
     $       ,MC(IL,J+1),UC(IL,J+1),VC(IL,J+1),uhat(il,j+1),vhat(il,j+1)
     $     ,DSUBCLD(IL),MC(IL,J),UC(IL,J),VC(IL,J),UHAT(IL,J),VHAT(IL,J)
       print*,'mb,msg,ilev=',mb(il),msg,ilev
       print*,'uc=',(uc(il,jj),jj=msg+1,ilev)
       print*,'vc=',(vc(il,jj),jj=msg+1,ilev)
       print*,'mc=',(mc(il,jj),jj=msg+1,ilev)
       print*,'mu=',(mu(il,jj),jj=msg+1,ilev)
       print*,'md=',(md(il,jj),jj=msg+1,ilev)
       print*,'u=',(u(il,jj),jj=msg+1,ilev)
       print*,'v=',(v(il,jj),jj=msg+1,ilev)
       print*,'ec=',(ec(il,jj),jj=msg+1,ilev)
       print*,'delpx=',(delpx(il,jj),jj=msg+1,ilev)
       print*,'delpy=',(delpy(il,jj),jj=msg+1,ilev)
       print*,'RHO=',(RHO(il,jj),jj=msg+1,ilev)
       print*,'dz=',(dz(il,jj),jj=msg+1,ilev)
      endif

c     if (abs(dudt(il,j)).gt.1.0e-2.or.abs(dvdt(il,j)).gt.1.0e-2) then
c      do i9=IL1G,IL2G
c      do j9=MSG+1,ILEV
c      print* ,'bad',i9,j9,dt,jt(i9),jb(i9)
c    $       ,dudt(i9,j9),dvdt(i9,j9),dp(i9,j9),dz(i9,j9)
c    $       ,mc(i9,j9+1),uc(i9,j9+1),vc(i9,j9+1),ec(i9,j9+1)
c    $       ,du(i9,j9+1)
c    $       ,uhat(i9,j9+1),vhat(i9,j9+1),delpx(i9,j9),delpx(i9,j9+1)
c    $       ,delpy(i9,j9),delpy(i9,j9+1)
c      enddo
c      enddo
c     endif
  950 CONTINUE
C  
      RETURN
      END 
