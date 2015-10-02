! **************************************************
! 												
!   This is a module for some basic subroutines 
!    
!    ETI+ POXY + DPOT + DINP
!    SPFUNC6 + SPFUNC8 + TRIPOL
!    WAVECK+DISPER+ BSP+DBSP 
!    INVERSE+RLUDCMP + RLUBKSB +
!	
!	                                        
!       March. 29, 2005           by   Bin Teng		
!
! ===================================================
!   
!
        MODULE MFUNC_mod

       CONTAINS  
C
C
C  Incident wave 
C
        FUNCTION ETI(X,Y)
        USE MVAR_MOD
        IMPLICIT  NONE
!
	  REAL*8,INTENT(IN):: X,Y
        Real*8  ETI
!
         ETI=Amp*DCos(WK*(X*DCOS(BETA)+Y*DSIN(BETA)) -W1*TimeRK)
	   ETI=RAMPF*ETI

        RETURN
        END FUNCTION ETI


C
C  Incident potential 
C
        FUNCTION POXY(X,Y,Z)
        USE MVAR_MOD
        
        IMPLICIT  NONE
!
	  REAL*8,INTENT(IN):: X,Y,Z
        Real*8  POXY
	  REAL*8  WKX,DUM
C
        IF(H.GT.0.0d0) THEN
          DUM=AMP*G/W1*DCOSH(WK*(Z+H))/DCOSH(WK*H)
        ELSE
          DUM=AMP*G/W1*DEXP(WK*Z)
        ENDIF
         
         WKX=WK*( X*DCOS(BETA)+Y*DSIN(BETA))
         POXY=DUM*DSIN(WKX-W1*TimeRK)
            

C             POXY = DCOS(X)*DSIN(Y)*DEXP(DSQRT(2.0D0)*Z)*DCOS(-W1*TimeRK)
C

! H   : water depth
! Amp : amplitude of incident waves
! BETA: incident angle
! W1  : angular frequency of incident waves
! V   : wave number in deep water
! WK  : wave number
! TPER: wave period
! Tstep: time step
! Time : simulation time at each time step
! TimeRK: simulation time at each RK step
! RampF: ramp function for incident potential
! RampV: ramp function for damping
!	   POXY=RAMPF*POXY

        RETURN
        END FUNCTION POXY
!
! ----------------------------------------------

C Time Derivatives of incident wave potential
C   
C
        FUNCTION DPOT(X,Y,Z)
        USE MVAR_MOD
        IMPLICIT  NONE
!
	  REAL*8,INTENT(IN):: X,Y,Z
        Real*8  DPOT
	  REAL*8  WKX,DUM
C C
        IF(H.GT.0.0d0) THEN
          DUM=-AMP*G*DCOSH(WK*(Z+H))/DCOSH(WK*H)
        ELSE
          DUM=-AMP*G*DEXP(WK*Z)
        ENDIF
!
         WKX=WK*( X*DCOS(BETA)+Y*DSIN(BETA))
         DPOT=DUM*DCOS(WKX-W1*TimeRK)
C
!	   DPOT=RAMPF*DPOT
C         DPOT = DCOS(X)*DSIN(Y)*DEXP(DSQRT(2.0D0)*Z)*DCOS(-W1*TimeRK)*-W1

        RETURN
        END FUNCTION DPOT


!

C -------------------------------------------
C Spacial Derivatives of incident wave potential
C IORDER=0: for a current
C IORDER=1: for the first order potential
C IORDER=2: for the second order potential 
C
        SUBROUTINE  DINP(X,Y,Z,DPOX,DPOY,DPOZ)
        USE MVAR_MOD
	  IMPLICIT    NONE
	  
	  REAL*8,INTENT(IN)::   X,Y,Z
        REAL*8,INTENT(OUT)::  DPOX,DPOY,DPOZ

	  REAL*8 DUM,WKX
C
         DUM=AMP*G/W1
         WKX=WK*(X*DCOS(BETA)+Y*DSIN(BETA))
		  
		IF (H.GT.0.0d0) THEN           
         DPOX= DUM*WK*DCOS(BETA)*
	1		   DCOSH(WK*(Z+H))/DCOSH(WK*H)*DCOS(WKX-W1*TimeRK)
         DPOY= DUM*WK*DSIN(BETA)*
	1		   DCOSH(WK*(Z+H))/DCOSH(WK*H)*DCOS(WKX-W1*TimeRK)
         DPOZ= DUM*WK*
	1		   DSINH(WK*(Z+H))/DCOSH(WK*H)*DSIN(WKX-W1*TimeRK)

        ELSE
        
         DPOX= DUM*WK*DCOS(BETA)*DEXP(WK*Z)*DCOS(WKX-W1*TimeRK)
         DPOY= DUM*WK*DSIN(BETA)*DEXP(WK*Z)*DCOS(WKX-W1*TimeRK)
         DPOZ= DUM*WK*DEXP(WK*Z)*DSIN(WKX-W1*TimeRK)
        
        ENDIF
!
!	   DPOX=RAMPF*DPOX
!	   DPOY=RAMPF*DPOY
!	   DPOZ=RAMPF*DPOZ

C             DPOX = -DSIN(X)*DSIN(Y)
C             DPOX = DPOX*DEXP(DSQRT(2.0D0)*Z)*DCOS(-W1*TimeRK)
C             DPOY = DCOS(X)*DCOS(Y)*DEXP(DSQRT(2.0D0)*Z)*DCOS(-W1*TimeRK)

C             DPOZ = DSQRT(2.0D0)
C             DPOZ = DPOZ * DCOS(X)*DSIN(Y)*DEXP(DSQRT(2.0D0)*Z)
C             DPOZ = DPOZ*DCOS(-W1*TimeRK)
C       DPOX = 0
C       DPOY = 0
C       DPOZ = 0
!
        RETURN
        END SUBROUTINE  DINP 
        
        
!
! ======================================================
! Spacial potential and derivatives of fiction potential
!  Kind=0: for constant potential phi=1.0 
!  Kind=1: for a uniform current in the x-direction phi=x 
!  Kind=2: for a uniform current in the y-direction phi=y 
!  Kind=3: for a uniform current in the z-direction phi=z 
!
!
        SUBROUTINE  DINP0(KIND,X,Y,Z,PHI,DPOX,DPOY,DPOZ)
	  IMPLICIT    NONE
	  
	  INTEGER KIND
	  REAL*8,INTENT(IN)::   X,Y,Z
        REAL*8,INTENT(OUT)::  PHI,DPOX,DPOY,DPOZ

		  
		IF (KIND==0) THEN
		  PHI=1.0d0           
          DPOX=0.0d0
          DPOY=0.0d0
          DPOZ=0.0d0
        ELSE IF (KIND==1) THEN
		  PHI=X       
          DPOX=1.0d0
          DPOY=0.0d0
          DPOZ=0.0d0
        ELSE IF (KIND==2) THEN
		  PHI=Y           
          DPOX=0.0d0
          DPOY=1.0d0
          DPOZ=0.0d0
        ELSE IF (KIND==3) THEN
		  PHI=Z           
          DPOX=0.0d0
          DPOY=0.0d0
          DPOZ=1.0d0          
        ENDIF
!
        RETURN
        END SUBROUTINE  DINP0         
        
C
C *********************************************************
C *                                                       *
C * CALCULATE THE SHAPE FUNCTIONS AND THEIR DERIVATIVES   *
C * FOR 6 NODE TRIANGULAR ELEMENTS                        *
C *                                                       *
C *********************************************************     
C
        SUBROUTINE SPFUNC6(SI,ETA,SF,DSF)
      	IMPLICIT  NONE

	  REAL*8,INTENT(IN) :: SI,ETA
	  REAL*8,INTENT(OUT):: SF(8),DSF(2,8)
C
        SF(1)=(1.0D0-SI-ETA)*(1.0D0-2.0D0*SI-2.0D0*ETA)
        SF(2)=SI *(2.0D0*SI -1.0D0)
        SF(3)=ETA*(2.0D0*ETA-1.0D0)
        SF(4)=4.0D0*SI*(1.0D0-SI-ETA)
        SF(5)=4.0D0*SI*ETA
        SF(6)=4.0D0*ETA*(1.0D0-SI-ETA)
c
        DSF(1,1)=4.0D0*SI+4.0*ETA-3.0D0
        DSF(1,2)=4.0D0*SI-1.0D0
        DSF(1,3)=0.0D0
        DSF(1,4)=4.0D0-8.0D0*SI-4.0D0*ETA
        DSF(1,5)=4.0D0*ETA
        DSF(1,6)=-4.0*ETA
C
        DSF(2,1)=4.0D0*SI+4.0D0*ETA-3.0D0
        DSF(2,2)=0.0D0
        DSF(2,3)=4.0D0*ETA-1.0D0
        DSF(2,4)=-4.0*SI
        DSF(2,5)=4.0*SI
        DSF(2,6)=4.0-8.0D0*ETA-4.0*SI
 
C
	RETURN
	END SUBROUTINE SPFUNC6
C
C
C
C
C *********************************************************
C *                                                       *
C * CALCULATE THE SHAPE FUNCTIONS AND THEIR DERIVATIVES   *
C * FOR 8 NODE QUADRILATERIAL ELEMENTS                    *
C *                                                       *
C *********************************************************     
C
        SUBROUTINE SPFUNC8(SI,ETA,SF,DSF)
        IMPLICIT  NONE

	  REAL*8,INTENT(IN) :: SI,ETA
	  REAL*8,INTENT(OUT):: SF(8),DSF(2,8)
C
        SF(1)=-0.25D0*(1.0D0-SI)*(1.0D0-ETA)*(1.0D0+SI+ETA)
        SF(2)= 0.5D0 *(1.0D0-SI*SI)*(1.0-ETA)
        SF(3)= 0.25D0*(1.0D0+SI)*(1.0D0-ETA)*(SI-ETA-1.0D0)
        SF(4)= 0.5D0 *(1.0D0-ETA*ETA)*(1.0D0+SI)
        SF(5)= 0.25D0*(1.0D0+SI)*(1.0D0+ETA)*(SI+ETA-1.0D0)
        SF(6)= 0.5D0 *(1.0-SI*SI)*(1.0D0+ETA)
        SF(7)=-0.25D0*(1.0D0-SI)*(1.0D0+ETA)*(SI-ETA+1.0D0)
        SF(8)= 0.5D0 *(1.0D0-ETA*ETA)*(1.0-SI)            
      
        DSF(1,1)= 0.25D0*(2.0D0*SI+ETA)*(1.0D0-ETA)
        DSF(1,2)=-SI*(1.0D0-ETA)     	
        DSF(1,3)= 0.25D0*(2.0D0*SI-ETA)*(1.0D0-ETA)
        DSF(1,4)= 0.5D0*(1.0D0-ETA*ETA)
        DSF(1,5)= 0.25D0*(2.0*SI+ETA)*(1.0D0+ETA)
        DSF(1,6)=-SI*(1.0D0+ETA)     
        DSF(1,7)= 0.25D0*(2.0D0*SI-ETA)*(1.0D0+ETA)
        DSF(1,8)=-0.5D0*(1.0-ETA*ETA)
C
        DSF(2,1)= 0.25D0*(SI+2.0D0*ETA)*(1.0D0-SI)
        DSF(2,2)=-0.5D0*(1.0-SI*SI)
        DSF(2,3)= 0.25D0*(1.0D0+SI)*(2.0D0*ETA-SI)
        DSF(2,4)=-(1.0D0+SI)*ETA
        DSF(2,5)= 0.25*(1.0D0+SI)*(SI+2.0D0*ETA)
        DSF(2,6)= 0.5D0*(1.0D0-SI*SI)
        DSF(2,7)=-0.25D0*(1.0D0-SI)*(SI-2.0D0*ETA)
        DSF(2,8)=-(1.0D0-SI)*ETA   
C
	RETURN
	END SUBROUTINE SPFUNC8


C
C *********************************************************
C *                                                       *
C * CALCULATE THE SHAPE FUNCTIONS AND THEIR DERIVATIVES   *
C * FOR 6 NODE TRIANGULAR ELEMENTS                        *
C *                                                       *
C *********************************************************     
C
        SUBROUTINE SPFUNC6_1(SI,ETA,SF,DSF,DDSF)
      	IMPLICIT  NONE

	  REAL*8,INTENT(IN) :: SI,ETA
	  REAL*8,INTENT(OUT):: SF(8),DSF(2,8),DDSF(3,8)
C
        SF(1)=(1.0D0-SI-ETA)*(1.0D0-2.0D0*SI-2.0D0*ETA)
        SF(2)=SI *(2.0D0*SI -1.0D0)
        SF(3)=ETA*(2.0D0*ETA-1.0D0)
        SF(4)=4.0D0*SI*(1.0D0-SI-ETA)
        SF(5)=4.0D0*SI*ETA
        SF(6)=4.0D0*ETA*(1.0D0-SI-ETA)
! ------------------------------------------
!
        DSF(1,1)=4.0D0*SI+4.0*ETA-3.0D0
        DSF(1,2)=4.0D0*SI-1.0D0
        DSF(1,3)=0.0D0
        DSF(1,4)=4.0D0-8.0D0*SI-4.0D0*ETA
        DSF(1,5)=4.0D0*ETA
        DSF(1,6)=-4.0*ETA
C
        DSF(2,1)=4.0D0*SI+4.0D0*ETA-3.0D0
        DSF(2,2)=0.0D0
        DSF(2,3)=4.0D0*ETA-1.0D0
        DSF(2,4)=-4.0*SI
        DSF(2,5)=4.0*SI
        DSF(2,6)=4.0-8.0D0*ETA-4.0*SI

!  ---------------------------
        DDSF(1,1)= 4.D0
        DDSF(1,2)= 4.0D0     	
        DDSF(1,3)= 0.0D0
        DDSF(1,4)=-8.0D0
        DDSF(1,5)= 0.0D0
        DDSF(1,6)= 0.0D0     

        DDSF(2,1)= 4.0D0
        DDSF(2,2)= 0.0D0
        DDSF(2,3)= 4.0D0
        DDSF(2,4)= 0.0D0
        DDSF(2,5)= 0.0D0
        DDSF(2,6)=-8.0D0

        DDSF(3,1)= 4.0D0
        DDSF(3,2)= 0.0D0   	
        DDSF(3,3)= 0.0D0
        DDSF(3,4)=-4.0d0
        DDSF(3,5)= 4.0D0
        DDSF(3,6)=-4.0d0     

	  RETURN
	 END SUBROUTINE SPFUNC6_1
C



        SUBROUTINE SPFUNC8_1(SI,ETA,SF,DSF,DDSF)
        IMPLICIT  NONE

        REAL*8,INTENT(IN) :: SI,ETA
        REAL*8,INTENT(OUT):: SF(8),DSF(2,8),DDSF(3,8)
!
        SF(1)=-0.25D0*(1.0D0-SI)*(1.0D0-ETA)*(1.0D0+SI+ETA)
        SF(2)= 0.5D0 *(1.0D0-SI*SI)*(1.0-ETA)
        SF(3)= 0.25D0*(1.0D0+SI)*(1.0D0-ETA)*(SI-ETA-1.0D0)
        SF(4)= 0.5D0 *(1.0D0-ETA*ETA)*(1.0D0+SI)
        SF(5)= 0.25D0*(1.0D0+SI)*(1.0D0+ETA)*(SI+ETA-1.0D0)
        SF(6)= 0.5D0 *(1.0-SI*SI)*(1.0D0+ETA)
        SF(7)=-0.25D0*(1.0D0-SI)*(1.0D0+ETA)*(SI-ETA+1.0D0)
        SF(8)= 0.5D0 *(1.0D0-ETA*ETA)*(1.0-SI)            
! ------------------------------------------------------------
      
        DSF(1,1)= 0.25D0*(2.0D0*SI+ETA)*(1.0D0-ETA)
        DSF(1,2)=-SI*(1.0D0-ETA)     	
        DSF(1,3)= 0.25D0*(2.0D0*SI-ETA)*(1.0D0-ETA)
        DSF(1,4)= 0.5D0*(1.0D0-ETA*ETA)
        DSF(1,5)= 0.25D0*(2.0*SI+ETA)*(1.0D0+ETA)
        DSF(1,6)=-SI*(1.0D0+ETA)     
        DSF(1,7)= 0.25D0*(2.0D0*SI-ETA)*(1.0D0+ETA)
        DSF(1,8)=-0.5D0*(1.0-ETA*ETA)
!
        DSF(2,1)= 0.25D0*(SI+2.0D0*ETA)*(1.0D0-SI)
        DSF(2,2)=-0.5D0*(1.0-SI*SI)
        DSF(2,3)= 0.25D0*(1.0D0+SI)*(2.0D0*ETA-SI)
        DSF(2,4)=-(1.0D0+SI)*ETA
        DSF(2,5)= 0.25*(1.0D0+SI)*(SI+2.0D0*ETA)
        DSF(2,6)= 0.5D0*(1.0D0-SI*SI)
        DSF(2,7)=-0.25D0*(1.0D0-SI)*(SI-2.0D0*ETA)
        DSF(2,8)=-(1.0D0-SI)*ETA   
! --------------------------------------------------------------

        DDSF(1,1)=0.5D0*(1.0D0-ETA)
        DDSF(1,2)=-(1.0D0-ETA)     	
        DDSF(1,3)= 0.5D0*(1.0D0-ETA)
        DDSF(1,4)= 0.0D0
        DDSF(1,5)= 0.5D0*(1.0D0+ETA)
        DDSF(1,6)=-(1.0D0+ETA)     
        DDSF(1,7)= 0.5D0*(1.0D0+ETA)
        DDSF(1,8)= 0.0D0

        DDSF(2,1)= 0.50D0*(1.0D0-SI)
        DDSF(2,2)= 0.0D0
        DDSF(2,3)= 0.50D0*(1.0D0+SI)
        DDSF(2,4)=-(1.0D0+SI)
        DDSF(2,5)= 0.50D0*(1.0D0+SI)
        DDSF(2,6)= 0.0D0
        DDSF(2,7)= 0.50D0*(1.0D0-SI)
        DDSF(2,8)=-(1.0D0-SI) 

        DDSF(3,1)= 0.25D0*(-2.0D0*SI-ETA*2.0D0+1.0D0)
        DDSF(3,2)= SI     	
        DDSF(3,3)= 0.25D0*(-2.0D0*SI+ETA*2.0D0-1.0D0)
        DDSF(3,4)= -ETA
        DDSF(3,5)= 0.25D0*(2.0D0*SI+ETA*2.0D0+1.0D0)
        DDSF(3,6)=-SI     
        DDSF(3,7)= 0.25D0*(2.0D0*SI-ETA*2.0D0-1.0D0)
        DDSF(3,8)= ETA

        RETURN
	  END SUBROUTINE SPFUNC8_1
!C
!!C
! ===================================================
!C
!!C Tri-pole transformation for removing singularity
!C

!	SUBROUTINE TRIPOL(NODNUM,NCNE,XYZT,DXYZT)
!	USE MVAR_MOD
!
!	USE PVAR_MOD
!	USE TRVAR_MOD
!	IMPLICIT NONE
!!C
!	 INTEGER NCNE,I,J,LK,LJ,LI,IXRHO,IYRHO,NOTRI,NODNUM
!	REAL*8 XIQ(4),WIQ(4)
!	REAL*8 SF(8),DSF(2,8),XJ(3,3)
!	REAL*8 SITRI(3,3),ETATRI(3,3)
!	REAL*8 XYZT(3,8),DXYZT(3,8)
!	REAL*8 X,Y,Z,DX,DY,DZ,SI,ETA,XRHO,YRHO,AG,
!	1	DET,DET1,DET2,DET3,DUM,DUM1,DUM2,DUM3 
!!C
!!C
!!C
!!!C** matrix XIQ store the Gauss-Legendre sampling points (4) for the 
!!!C   quadrilateral element 
!!C
!      DATA XIQ/-0.86113 63115 94053D0,-0.33998 10435 84856D0,
!     1          0.33998 10435 84856D0, 0.86113 63115 94053D0/
!!C
!!!C** matrix WIQ store the Gauss-Legendre weighting fa!Ctors for the 
!!!C   quadrilateral element 
!!C
!      DATA WIQ/0.34785 48451 37454D0,0.65214 51548 62546D0,
!     1         0.65214 51548 62546D0,0.34785 48451 37454D0/
!!C
!!C
!!!C** identity the local nodal position of INODE IN IELEM
!!C
!      	IF(NCNE.EQ.8) THEN
!          IF(NODNUM.EQ.1.OR.NODNUM.EQ.3.OR.
!     1     NODNUM.EQ.5.OR.NODNUM.EQ.7) THEN
!           NOTRI=2
!          ELSE
!           NOTRI=3
!          ENDIF
!      	ELSE
!          IF(NODNUM.LE.3)THEN
!            NOTRI=1
!          ELSE
!            NOTRI=2
!          ENDIF
!      	ENDIF
!!C
!!!C** identity the local coordinate(SI,ETA) of the trianges
!!C
!      IF(NCNE.EQ.8) THEN
!!C
!!!C** quadrilateral element
!!C
!      IF(NODNUM.EQ.1) THEN
!!C
!      SITRI (1,1)=-1.0D0
!      SITRI (1,2)= 1.0D0
!      SITRI (1,3)=-1.0D0
!!C
!      SITRI (2,1)=-1.0D0
!      SITRI (2,2)= 1.0D0
!      SITRI (2,3)= 1.0D0
!!C
!      SITRI (3,1)= 0.0D0
!      SITRI (3,2)= 0.0D0
!      SITRI (3,3)= 0.0D0
!!C
!      ETATRI(1,1)=-1.0D0
!      ETATRI(1,2)= 1.0D0
!      ETATRI(1,3)= 1.0D0
!!C
!      ETATRI(2,1)=-1.0D0
!      ETATRI(2,2)=-1.0D0
!      ETATRI(2,3)= 1.0D0
!!C
!      ETATRI(3,1)= 0.0D0
!      ETATRI(3,2)= 0.0D0
!      ETATRI(3,3)= 0.0D0
!!C
!      ELSE IF(NODNUM.EQ.3) THEN
!!C
!      SITRI (1,1)= 1.0D0
!      SITRI (1,2)=-1.0D0
!      SITRI (1,3)=-1.0D0
!!C
!      SITRI (2,1)= 1.0D0
!      SITRI (2,2)= 1.0D0
!      SITRI (2,3)=-1.0D0
!!C
!      SITRI (3,1)= 0.0D0
!      SITRI (3,2)= 0.0D0
!      SITRI (3,3)= 0.0D0
!!C
!      ETATRI(1,1)=-1.0D0
!      ETATRI(1,2)= 1.0D0
!      ETATRI(1,3)=-1.0D0
!!C
!      ETATRI(2,1)=-1.0D0
!      ETATRI(2,2)= 1.0D0
!      ETATRI(2,3)= 1.0D0
!!C
!      ETATRI(3,1)= 0.0D0
!      ETATRI(3,2)= 0.0D0
!      ETATRI(3,3)= 0.0D0
!!C
!      ELSE IF(NODNUM.EQ.5) THEN
!!C
!      SITRI (1,1)= 1.0D0
!      SITRI (1,2)=-1.0D0
!      SITRI (1,3)=-1.0D0
!!C
!      SITRI (2,1)= 1.0D0
!      SITRI (2,2)=-1.0D0
!      SITRI (2,3)= 1.0D0
!!C
!      SITRI (3,1)= 0.0D0
!      SITRI (3,2)= 0.0D0
!      SITRI (3,3)= 0.0D0
!!C
!      ETATRI(1,1)= 1.0D0
!      ETATRI(1,2)= 1.0D0
!      ETATRI(1,3)=-1.0D0
!!C
!      ETATRI(2,1)= 1.0D0
!      ETATRI(2,2)=-1.0D0
!      ETATRI(2,3)=-1.0D0
!!C
!      ETATRI(3,1)= 0.0D0
!      ETATRI(3,2)= 0.0D0
!      ETATRI(3,3)= 0.0D0
!!C
!      ELSE IF(NODNUM.EQ.7) THEN
!!C
!      SITRI (1,1)=-1.0D0
!      SITRI (1,2)=-1.0D0
!      SITRI (1,3)= 1.0D0
!!C
!      SITRI (2,1)=-1.0D0
!      SITRI (2,2)= 1.0D0
!      SITRI (2,3)= 1.0D0
!!C
!      SITRI (3,1)= 0.0D0
!      SITRI (3,2)= 0.0D0
!      SITRI (3,3)= 0.0D0
!!C
!      ETATRI(1,1)= 1.0D0
!      ETATRI(1,2)=-1.0D0
!      ETATRI(1,3)=-1.0D0
!!C
!      ETATRI(2,1)= 1.0D0
!      ETATRI(2,2)=-1.0D0
!      ETATRI(2,3)= 1.0D0
!!C
!      ETATRI(3,1)= 0.0D0
!      ETATRI(3,2)= 0.0D0
!      ETATRI(3,3)= 0.0D0
!!C
!      ELSE IF(NODNUM.EQ.2) THEN
!!C
!      SITRI (1,1)= 0.0D0
!      SITRI (1,2)=-1.0D0
!      SITRI (1,3)=-1.0D0
!!C
!      SITRI (2,1)= 0.0D0
!      SITRI (2,2)= 1.0D0
!      SITRI (2,3)=-1.0D0
!!C
!      SITRI (3,1)= 0.0D0
!      SITRI (3,2)= 1.0D0
!      SITRI (3,3)= 1.0D0
!!C
!      ETATRI(1,1)=-1.0D0
!      ETATRI(1,2)= 1.0D0
!      ETATRI(1,3)=-1.0D0
!!C
!      ETATRI(2,1)=-1.0D0
!      ETATRI(2,2)= 1.0D0
!      ETATRI(2,3)= 1.0D0
!!C
!      ETATRI(3,1)=-1.0D0
!      ETATRI(3,2)=-1.0D0
!      ETATRI(3,3)= 1.0D0
!!C
!      ELSE IF(NODNUM.EQ.4) THEN
!!C
!      SITRI (1,1)= 1.0D0
!      SITRI (1,2)= 1.0D0
!      SITRI (1,3)=-1.0D0
!!C
!      SITRI (2,1)= 1.0D0
!      SITRI (2,2)=-1.0D0
!      SITRI (2,3)=-1.0D0
!!C
!      SITRI (3,1)= 1.0D0
!      SITRI (3,2)=-1.0D0
!      SITRI (3,3)= 1.0D0
!!C
!      ETATRI(1,1)= 0.0D0
!      ETATRI(1,2)= 1.0D0
!      ETATRI(1,3)= 1.0D0
!!C
!      ETATRI(2,1)= 0.0D0
!      ETATRI(2,2)= 1.0D0
!      ETATRI(2,3)=-1.0D0
!!C
!      ETATRI(3,1)= 0.0D0
!      ETATRI(3,2)=-1.0D0
!      ETATRI(3,3)=-1.0D0
!!C
!      ELSE IF(NODNUM.EQ.6) THEN
!!C
!      SITRI (1,1)= 0.0D0
!      SITRI (1,2)=-1.0D0
!      SITRI (1,3)=-1.0D0
!!C
!      SITRI (2,1)= 0.0D0
!      SITRI (2,2)=-1.0D0
!      SITRI (2,3)= 1.0D0
!!C
!      SITRI (3,1)= 0.0D0
!      SITRI (3,2)= 1.0D0
!      SITRI (3,3)= 1.0D0
!!C
!      ETATRI(1,1)= 1.0D0
!      ETATRI(1,2)= 1.0D0
!      ETATRI(1,3)=-1.0D0
!!C
!      ETATRI(2,1)= 1.0D0
!      ETATRI(2,2)=-1.0D0
!      ETATRI(2,3)=-1.0D0
!!C
!      ETATRI(3,1)= 1.0D0
!      ETATRI(3,2)=-1.0D0
!      ETATRI(3,3)= 1.0D0
!!C
!      ELSE IF(NODNUM.EQ.8) THEN
!!C
!      SITRI (1,1)=-1.0D0
!      SITRI (1,2)= 1.0D0
!      SITRI (1,3)=-1.0D0
!!C
!      SITRI (2,1)=-1.0D0
!      SITRI (2,2)= 1.0D0
!      SITRI (2,3)= 1.0D0
!!C
!      SITRI (3,1)=-1.0D0
!      SITRI (3,2)=-1.0D0
!      SITRI (3,3)= 1.0D0
!!C
!      ETATRI(1,1)= 0.0D0
!      ETATRI(1,2)= 1.0D0
!      ETATRI(1,3)= 1.0D0
!!C
!      ETATRI(2,1)= 0.0D0
!      ETATRI(2,2)=-1.0D0
!      ETATRI(2,3)= 1.0D0
!!C
!      ETATRI(3,1)= 0.0D0
!      ETATRI(3,2)=-1.0D0
!      ETATRI(3,3)=-1.0D0
!!C
!      ENDIF
!!C
!      ELSE
!!C
!!!C** triangular element ***************
!!!C         
!      IF(NODNUM.EQ.1) THEN
!!C
!      SITRI (1,1)= 0.0D0
!      SITRI (1,2)= 1.0D0
!      SITRI (1,3)= 0.0D0
!!C
!      SITRI (2,1)= 0.0D0
!      SITRI (2,2)= 0.0D0
!      SITRI (2,3)= 0.0D0
!!C
!      SITRI (3,1)= 0.0D0
!      SITRI (3,2)= 0.0D0
!      SITRI (3,3)= 0.0D0
!!C
!      ETATRI(1,1)= 0.0D0
!      ETATRI(1,2)= 0.0D0
!      ETATRI(1,3)= 1.0D0
!!C
!      ETATRI(2,1)= 0.0D0
!      ETATRI(2,2)= 0.0D0
!      ETATRI(2,3)= 0.0D0
!!C
!      ETATRI(3,1)= 0.0D0
!      ETATRI(3,2)= 0.0D0
!      ETATRI(3,3)= 0.0D0
!!C
!      ELSE IF(NODNUM.EQ.2) THEN
!!C
!      SITRI (1,1)= 1.0D0
!      SITRI (1,2)= 0.0D0
!      SITRI (1,3)= 0.0D0
!!C
!      SITRI (2,1)= 0.0D0
!      SITRI (2,2)= 0.0D0
!      SITRI (2,3)= 0.0D0
!C
!      SITRI (3,1)= 0.0D0
!      SITRI (3,2)= 0.0D0
!      SITRI (3,3)= 0.0D0
!!C
!      ETATRI(1,1)= 0.0D0
!      ETATRI(1,2)= 1.0D0
!      ETATRI(1,3)= 0.0D0
!!C
!      ETATRI(2,1)= 0.0D0
!      ETATRI(2,2)= 0.0D0
!      ETATRI(2,3)= 0.0D0
!!C
!      ETATRI(3,1)= 0.0D0
!      ETATRI(3,2)= 0.0D0
!      ETATRI(3,3)= 0.0D0
!!C
!      ELSE IF(NODNUM.EQ.3) THEN
!!C
!      SITRI (1,1)= 0.0D0
!      SITRI (1,2)= 0.0D0
!      SITRI (1,3)= 1.0D0
!!C
!      SITRI (2,1)= 0.0D0
!      SITRI (2,2)= 0.0D0
!      SITRI (2,3)= 0.0D0
!!C
!      SITRI (3,1)= 0.0D0
!      SITRI (3,2)= 0.0D0
!      SITRI (3,3)= 0.0D0
!!C
!      ETATRI(1,1)= 1.0D0
!      ETATRI(1,2)= 0.0D0
!      ETATRI(1,3)= 0.0D0
!!C
!      ETATRI(2,1)= 0.0D0
!      ETATRI(2,2)= 0.0D0
!      ETATRI(2,3)= 0.0D0
!!C
!      ETATRI(3,1)= 0.0D0
!      ETATRI(3,2)= 0.0D0
!      ETATRI(3,3)= 0.0D0
!!C
!      ELSE IF(NODNUM.EQ.4) THEN
!!C
!      SITRI (1,1)= 0.5D0
!      SITRI (1,2)= 0.0D0
!      SITRI (1,3)= 0.0D0
!!C
!      SITRI (2,1)= 0.5D0
!      SITRI (2,2)= 1.0D0
!      SITRI (2,3)= 0.0D0
!!C
!      SITRI (3,1)= 0.0D0
!      SITRI (3,2)= 0.0D0
!      SITRI (3,3)= 0.0D0
!!C
!      ETATRI(1,1)= 0.0D0
!      ETATRI(1,2)= 1.0D0
!      ETATRI(1,3)= 0.0D0
!!C
!      ETATRI(2,1)= 0.0D0
!      ETATRI(2,2)= 0.0D0
!      ETATRI(2,3)= 1.0D0
!!C
!      ETATRI(3,1)= 0.0D0
!      ETATRI(3,2)= 0.0D0
!      ETATRI(3,3)= 0.0D0
!!C
!      ELSE IF(NODNUM.EQ.5) THEN
!!C
!      SITRI (1,1)= 0.5D0
!      SITRI (1,2)= 0.0D0
!      SITRI (1,3)= 1.0D0
!!C
!      SITRI (2,1)= 0.5D0
!      SITRI (2,2)= 0.0D0
!      SITRI (2,3)= 0.0D0
!!C
!      SITRI (3,1)= 0.0D0
!      SITRI (3,2)= 0.0D0
!      SITRI (3,3)= 0.0D0
!!C
!      ETATRI(1,1)= 0.5D0
!      ETATRI(1,2)= 0.0D0
!      ETATRI(1,3)= 0.0D0
!!C
!      ETATRI(2,1)= 0.5D0
!      ETATRI(2,2)= 1.0D0
!      ETATRI(2,3)= 0.0D0
!!C
!      ETATRI(3,1)= 0.0D0
!      ETATRI(3,2)= 0.0D0
!      ETATRI(3,3)= 0.0D0
!!C
!      ELSE IF(NODNUM.EQ.6) THEN
!!C
!      SITRI (1,1)= 0.0D0
!      SITRI (1,2)= 0.0D0
!      SITRI (1,3)= 1.0D0
!!C
!      SITRI (2,1)= 0.0D0
!      SITRI (2,2)= 1.0D0
!      SITRI (2,3)= 0.0D0
!!C
!      SITRI (3,1)= 0.0D0
!      SITRI (3,2)= 0.0D0
!      SITRI (3,3)= 0.0D0
!!C
!      ETATRI(1,1)= 0.5D0
!      ETATRI(1,2)= 0.0D0
!      ETATRI(1,3)= 0.0D0
!!C
!      ETATRI(2,1)= 0.5D0
!      ETATRI(2,2)= 0.0D0
!      ETATRI(2,3)= 1.0D0
!!C
!      ETATRI(3,1)= 0.0D0
!      ETATRI(3,2)= 0.0D0
!      ETATRI(3,3)= 0.0D0
!!C
!      ENDIF
!!C
!      ENDIF
!!C
!      NOSAMP=0
!      DO 100 I=1,NOTRI
!      IF(NCNE.EQ.8) THEN
!      AG=2.0D0
!         IF(NOTRI.EQ.3) THEN
!            IF(I.EQ.1.OR.I.EQ.3) AG=1.0D0
!         ENDIF
!      ELSE
!         AG=0.5D0
!            IF(NOTRI.EQ.2) AG=0.25D0
!      ENDIF
!      DO 110 IXRHO=1,4
!      DO 120 IYRHO=1,4
!      NOSAMP=NOSAMP+1
!!C
!!!C**  calculate the shape function at the sampling points(XRHO,YRHO)
!!C
!      XRHO = XIQ(IXRHO)
!      YRHO = XIQ(IYRHO)
!!C
!      DUM1=0.5D0*(1.0D0-XRHO)
!      DUM2=0.25D0*(1.0D0+XRHO)*(1.0D0-YRHO)
!      DUM3=0.25D0*(1.0D0+XRHO)*(1.0D0+YRHO)
!!C
!      SI =DUM1*SITRI (I,1)+DUM2*SITRI (I,2)+DUM3*SITRI (I,3)
!      ETA=DUM1*ETATRI(I,1)+DUM2*ETATRI(I,2)+DUM3*ETATRI(I,3)
!!C
!!C
!      IF(NCNE.EQ.8) THEN
! 	CALL SPFUNC8(SI,ETA,SF,DSF)
!      ELSE
! 	CALL SPFUNC6(SI,ETA,SF,DSF)
!      ENDIF
!!C
!!!C** evaluate the Jacobian matrix at (SI,ETA)
!!C
!      	DO 130 LI=1,2
!      	DO 130 LJ=1,3
!      	DUM=0.0D0
!      	DO 140 LK=1,NCNE
!      	DUM=DUM+DSF(LI,LK)*XYZT(LJ,LK)
!140   	CONTINUE
!130   	XJ(LI,LJ)=DUM
!!C
!!!C** compute the determinant of the Jacobian maxtix at (XRHO,YRHO)
!!C
!      	DET1=XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2) 
!      	DET2=XJ(1,1)*XJ(2,3)-XJ(1,3)*XJ(2,1) 
!      	DET3=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1) 
!      	DET=DSQRT(DET1*DET1+DET2*DET2+DET3*DET3)
!!C
!!C
!      	X=0.0D0
!      	Y=0.0D0
!      	Z=0.0D0
!!C
!      	DX=0.0D0
!      	DY=0.0D0
!      	DZ=0.0D0
!!C
!      	DO 150 LK=1, NCNE
!      	X=X+SF(LK)*XYZT(1,LK)
!      	Y=Y+SF(LK)*XYZT(2,LK)
!      	Z=Z+SF(LK)*XYZT(3,LK)
!      	DX=DX+SF(LK)*DXYZT(1,LK)
!      	DY=DY+SF(LK)*DXYZT(2,LK)
!      	DZ=DZ+SF(LK)*DXYZT(3,LK)
!150   	CONTINUE
!!C
!      	XYNOD(1,NOSAMP)=X
!      	XYNOD(2,NOSAMP)=Y
!      	XYNOD(3,NOSAMP)=Z
!      	DXYNOD(1,NOSAMP)=DX
!      	DXYNOD(2,NOSAMP)=DY
!      	DXYNOD(3,NOSAMP)=DZ 
!        DXYNOD(4,NOSAMP)=(Y-YC)*DZ-(Z-ZC)*DY
!        DXYNOD(5,NOSAMP)=(Z-ZC)*DX-(X-XC)*DZ
!        DXYNOD(6,NOSAMP)=(X-XC)*DY-(Y-YC)*DX
!!C
!      	SAMNOD(NOSAMP,0)=0.25D0*AG*(1.0+XRHO)*WIQ(IXRHO)*WIQ(IYRHO)*DET
!!C
!      	DO 160 J=1,NCNE
!      	SAMNOD(NOSAMP,J)=SF(J)*SAMNOD(NOSAMP,0)
!160   	CONTINUE
!!C
!!120   	CONTINUE
!110   	CONTINUE
!100   	CONTINUE
!!C
!!C
!      	RETURN
!      	END SUBROUTINE TRIPOL
!

!
!C
C
C      ************************************************
C      *                                              *
C      *    CALCULATING THE Bj,4(X).                  *
C      *                                              *
C      ************************************************
C
       SUBROUTINE BSP1(X,BJ4,MSTA,MEND)
C 
C -------------------------------------------------
C   
       USE BVAR_mod
	 IMPLICIT NONE
C
	 INTEGER,INTENT(OUT):: MSTA,MEND
       REAL*8,INTENT(OUT):: BJ4(NAA)
	 REAL*8,INTENT(IN):: X
	 INTEGER I,I1,I2,I3,I4,LN,KK
	 REAL*8 D,E
C
	 BJ4=0.0D0
	 B  =0.0D0
C
	 DO I=1+KBSPL, NBU+KBSPL
	   I1=I
	   I2=I+1 
	   IF(X.GE.UTJ(I1).AND.X.LT.UTJ(I2)) THEN
	     B(I,1)=1 
	     LN=I
	   ELSE	IF(I.EQ.(NBU+KBSPL).AND. ABS(X-UTJ(I2))
     1                	       .LE.1E-6) THEN
	     B(I,1)=1 
	     LN=I
	   ELSE
	     B(I,1)=0
	   END IF
	 END DO
C
	 DO KK=2,KBSPL+1
	   DO I=1,NBU+KBSPL 	  	  
	     I1=I
		 I2=I+KK-1   
	     I3=I+KK
	     I4=I+1
	     IF(ABS(B(I1,KK-1)).LE.1E-9) THEN
	       D=0
	     ELSE
	       D=(X-UTJ(I1))*B(I1,KK-1)/(UTJ(I2)
     1        		   -UTJ(I1))
	     END IF
           IF(ABS(B(I4,KK-1)).LE.1E-9) THEN
	       E=0
	     ELSE
	       E=(UTJ(I3)-X)*B(I4,KK-1)/(UTJ(I3)
     1        		   -UTJ(I4))
	     END IF
	     B(I,KK)=D+E
	   END DO
       END DO
C

C	  
       DO I=LN-KBSPL, LN
         BJ4(I)=B(I,KBSPL+1)
	 END DO
C
	 MSTA=LN-KBSPL
	 MEND=LN
C
       END SUBROUTINE BSP1


C   
C      ************************************************
C      *                                              *
C      *    CALCULATING THE DBj,4(X).                 *
C      *                                              *
C      ************************************************
C 
       SUBROUTINE DBSP1(X,BJ4,DBJ4,MSTA,MEND)
C 
  
       USE BVAR_mod
       IMPLICIT NONE
C
	 INTEGER,INTENT(OUT):: MSTA,MEND
       REAL*8,INTENT(OUT):: BJ4(NAA),DBJ4(NAA)
	 REAL*8,INTENT(IN):: X
	 INTEGER I,I1,I2,I3,I4,LN,KK,J
	 REAL*8 D,E
C
	 BJ4=0.0D0
	 B  =0.0D0
C
C
	 DO I=1+KBSPL, NBU+KBSPL
	   I1=I
	   I2=I+1 
	   IF(X.GE.UTJ(I1).AND.X.LT.UTJ(I2)) THEN
	     B(I,1)=1 
	     LN=I
	   ELSE	IF(I.EQ.(NBU+KBSPL).AND. DABS(X-UTJ(I2))
     1                	       .LE.1E-6) THEN
	     B(I,1)=1 
	     LN=I
	   ELSE
	     B(I,1)=0
	   END IF
	 END DO
C
	 DO KK=2,KBSPL+1
	   DO I=1,NBU+KBSPL 	  	  
	     I1=I
		 I2=I+KK-1   
	     I3=I+KK
	     I4=I+1
	     IF( DABS(B(I1,KK-1)).LE.1E-9) THEN
	       D=0
	     ELSE
	       D=(X-UTJ(I1))*B(I1,KK-1)/(UTJ(I2)
     1        		   -UTJ(I1))
	     END IF
           IF( DABS(B(I4,KK-1)).LE.1E-9) THEN
	       E=0
	     ELSE
	       E=(UTJ(I3)-X)*B(I4,KK-1)/(UTJ(I3)
     1        		   -UTJ(I4))
	     END IF
	     B(I,KK)=D+E
	   END DO
       END DO
C
C	  
       DO I=LN-KBSPL,LN
         BJ4(I)=B(I,KBSPL+1)
	 END DO
C
	 MSTA=LN-KBSPL
	 MEND=LN
C
	 DO 400 J=MSTA, MEND  	  
	     I1=J
	     I2=J+ KBSPL   
	     I3=J+ KBSPL+1
	     I4=J+ 1
C
	   IF( DABS(B(J, KBSPL)) .LE. 1E-8) THEN
	     D=0.0D0
	   ELSE
	     D=B(J,KBSPL)/(UTJ(I2)-UTJ(I1))
	   END IF
C           
	   IF( DABS(B(J+1, KBSPL)).LE.1E-8) THEN
	     E=0.0D0
	   ELSE
	     E=B(J+1,KBSPL)/(UTJ(I3)-UTJ(I4)) 
	   END IF
C
         DBJ4(J)=KBSPL*( D - E )
400	 CONTINUE
C
       END SUBROUTINE DBSP1
!
!
!    
! ---------------------------------------- 
! A   :
! N   :
! NP  :
! IP  :
! NSYS:
!
! B   :
! LI  :
! NMOD:
!
! INDX:
! ----------------------------------------
!

! 
! ********************************************* 
! * SOLUTION OF LINEAR EQUATIONS  [A]{X}= {B} * 
! *           LU DECOMPOSITION                * 
! * FROM 'NUMERICAL RECIPES'     pp. 35-37    * 
! ********************************************* 
!
!
! ---------------------------------------- 
! NSYS: 矩阵A[:,:]的个数，用于开数组 
! IP  : 本次计算的矩阵序号
! A   :
! N   :
! NP  :
! LI  :
! NMOD:
! INDX:
! B   :
! ---------------------------------------- 
! 

        SUBROUTINE RLUDCMP(IP,A,N,NP,NSYS,INDX,D)           
        IMPLICIT REAL*8(A-H,O-Z)  
        PARAMETER (NMAX=5000, TINY=1.0E-20) 
        DIMENSION INDX(NP,NSYS) 
        REAL*8 SUM,DUM,AAMAX 
        REAL*8 A(NP,NP,NSYS),VV(NMAX) 
C 
        D=1. 
        DO 12 I=1, N 
        AAMAX=(0., 0.) 
          DO 11 J=1, N 
          IF ( DABS(A(I,J,IP)).GT. DABS(AAMAX) )  
     1                        AAMAX=A(I,J,IP) 
11      CONTINUE 
!
        IF (DABS(AAMAX) .EQ. 0.0)  THEN	  
	    Print  *, ' SINGULAR MATRIX   inside RLUDCMP' 
	    Print  *, ' IP=',IP,' I=',I
		PAUSE
        ENDIF
!
	  VV(I)=1./AAMAX 
12      CONTINUE 
        DO 19 J=1, N 
          DO 14 I=1, J-1 
            SUM=A(I,J,IP) 
            DO 13 K=1, I-1 
              SUM=SUM-A(I,K,IP)*A(K,J,IP) 
13            CONTINUE 
            A(I,J,IP)=SUM 
14          CONTINUE 
          AAMAX=(0., 0.) 
          DO 16 I=J, N                          
            SUM=A(I,J,IP) 
            DO 15 K=1, J-1 
              SUM=SUM-A(I,K,IP)*A(K,J,IP) 
15            CONTINUE 
            A(I,J,IP)=SUM 
            DUM=VV(I)*SUM 
            IF (DABS(DUM) .GE. DABS(AAMAX)) THEN 
            IMAX=I 
            AAMAX=DUM 
            END IF 
16          CONTINUE 
         IF(J .NE. IMAX) THEN 
         DO 17 K=1, N 
           DUM=A(IMAX,K,IP) 
           A(IMAX,K,IP)=A(J,K,IP) 
           A(J,K,IP)=DUM 
17         CONTINUE 
          D=-D 
          VV(IMAX)=VV(J) 
          END IF 
          INDX(J,IP)=IMAX 
          IF(A(J,J,IP).EQ.0.) A(J,J,IP)=TINY 
          IF(J.NE.N) THEN 
          DUM=1./A(J,J,IP) 
          DO 18 I=J+1,N 
            A(I,J,IP)=A(I,J,IP)*DUM 
18          CONTINUE 
          END IF 
19        CONTINUE 
        RETURN 
	 END SUBROUTINE RLUDCMP





!    
! ---------------------------------------- 
! A   :
! N   :
! NP  :
! IP  :
! NSYS:
!
! B   :
! LI  :
! NMOD:
!
! INDX:
! ----------------------------------------
!
        SUBROUTINE RLUBKSB(IP,A,N,NP,LI,NSYS,NMOD,INDX,B)  
        IMPLICIT REAL*8(A-H,O-Z)  
        DIMENSION INDX(NP,NSYS) 
        REAL*8 SUM,A(NP,NP,NSYS),B(NP,NMOD,NSYS) 
C                                     
        II=0 
        DO 12 I=1,N 
         LL=INDX(I,IP) 
        SUM=B(LL,LI,IP) 
        B(LL,LI,IP)=B(I,LI,IP) 
        IF(II.NE.0) THEN 
        DO 11 J=II, I-1 
        SUM=SUM-A(I,J,IP)*B(J,LI,IP) 
11      CONTINUE 
        ELSE IF (DABS(SUM).NE. 0.) THEN 
        II=I 
        END IF 
        B(I,LI,IP)=SUM 
12      CONTINUE 
        DO 14 I=N, 1, -1 
        SUM=B(I,LI,IP) 
        DO 13 J=I+1, N 
        SUM=SUM-A(I,J,IP)*B(J,LI,IP) 
13      CONTINUE 
        B(I,LI,IP)=SUM/A(I,I,IP) 
14      CONTINUE 
        RETURN 
	  END SUBROUTINE RLUBKSB

!
!
!
C 
C ********************************************* 
C * Inverse of a matric  [A] by               * 
C *           LU decomposition                * 
C * From 'NUMERICAL RECIPES'     pp. 35-37    * 
C ********************************************* 
C 
       SUBROUTINE INVERSE(a,y,n,np)
       INTEGER n,np,indx(n)
       REAL*8 a(np,np),b(np),y(np,np),d

	 CALL ludcmp(a,n,np,indx,d)
	  y=0.0d0
       do 10 i=1, n
	  y(i,i)=1.0d0
10	 continue
c
	 do 100 j=1, n
	  do 20 i=1, n
	   b(i)=y(j,i) !i,j
20	  continue
c
        call  lubksb(a,n,np,indx,b)
c
	  do 30 i=1, n
	   y(j,i)=b(i) !j,i
30	  continue

100	 continue


      END	 SUBROUTINE INVERSE


!
C 
C ********************************************* 
C * SOLUTION OF LINEAR EQUATIONS  [A]{X}= {B} * 
C *           LU DECOMPOSITION                * 
C * FROM 'NUMERICAL RECIPES'     pp. 35-37    * 
C ********************************************* 
C 
      SUBROUTINE ludcmp(a,n,np,indx,d)

      INTEGER n,np,indx(n),NMAX
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (DAbs(a(i,j)).gt.aamax) aamax=DAbs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*DAbs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax) then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n) then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END	 SUBROUTINE ludcmp




      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL*8 sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0) then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12     continue
!       
	 do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14     continue
       return
       END	 SUBROUTINE lubksb


C 
C
C ***********************************************************
C         From Chapter  2.6 of Numerical Recipes
C       
C
C ---------------------------------------------------------
C   Given a matrix A, with logical dimensions M by N and physical 
C dimensions MP by NP, this routine computes its singular value 
C decomposition, A=U.W.V^T. The matrix U replaces A on output. 
C The diagonal matrix of singular values W is output as a vector W. 
C The matrix V (not the transpose V^T) is output as V. M must be 
C greater or equal to N; if it is amaller, then A should be filled 
C up to square with zero rows. 
C **************************************************************
C
!      SUBROUTINE svdcmp(a,m,n,mp,np,w,v)  
!      IMPLICIT  REAL*4(A-H, O-Z)
!      REAL*8    a(mp,np),v(np,np),w(np)
!      PARAMETER (NMAX=5000)
!C
!      DIMENSION rv1(NMAX)
!C
!      g=0.0
!      scale=0.0
!      anorm=0.0
!      do 25 i=1,n
!        l=i+1
!        rv1(i)=scale*g
!        g=0.0
!        s=0.0
!        scale=0.0
!        if(i.le.m)then
!          do 11 k=i,m
!            scale=scale+Abs(a(k,i))
!11        continue
!          if(scale.ne.0.0)then
!            do 12 k=i,m
!              a(k,i)=a(k,i)/scale
!              s=s+a(k,i)*a(k,i)
!12          continue
!            f=a(i,i)
!            g=-sign(Sqrt(s),f)
!            h=f*g-s
!            a(i,i)=f-g
!            do 15 j=l,n
!              s=0.0
!              do 13 k=i,m
!                s=s+a(k,i)*a(k,j)
!13            continue
!              f=s/h
!              do 14 k=i,m
!                a(k,j)=a(k,j)+f*a(k,i)
!14            continue
!15          continue
!            do 16 k=i,m
!              a(k,i)=scale*a(k,i)
!16          continue
!          endif
!        endif
!        w(i)=scale *g
!        g=0.0
!        s=0.0
!        scale=0.0
!        if((i.le.m).and.(i.ne.n))then
!          do 17 k=l,n
!            scale=scale+Abs(a(i,k))
!17        continue
!          if(scale.ne.0.0)then
!            do 18 k=l,n
!              a(i,k)=a(i,k)/scale
!              s=s+a(i,k)*a(i,k)
!18          continue
!            f=a(i,l)
!            g=-sign(Sqrt(s),f)
!            h=f*g-s
!            a(i,l)=f-g
!            do 19 k=l,n
!              rv1(k)=a(i,k)/h
!19          continue
!            do 23 j=l,m
!              s=0.0
!              do 21 k=l,n
!                s=s+a(j,k)*a(i,k)
!21            continue
!              do 22 k=l,n
!                a(j,k)=a(j,k)+s*rv1(k)
!22            continue
!23          continue
!            do 24 k=l,n
!              a(i,k)=scale*a(i,k)
!24          continue
!          endif
!        endif
!        anorm=max(anorm,(Abs(w(i))+Abs(rv1(i))))
!25    continue
!      do 32 i=n,1,-1
!        if(i.lt.n)then
!          if(g.ne.0.0)then
!            do 26 j=l,n
!              v(j,i)=(a(i,j)/a(i,l))/g
!26          continue
!            do 29 j=l,n
!              s=0.0
!              do 27 k=l,n
!                s=s+a(i,k)*v(k,j)
!27            continue
!              do 28 k=l,n
!                v(k,j)=v(k,j)+s*v(k,i)
!28            continue
!29          continue
!          endif
!          do 31 j=l,n
!            v(i,j)=0.0
!            v(j,i)=0.0
!31        continue
!        endif
!        v(i,i)=1.0
!        g=rv1(i)
!        l=i
!32    continue
!      do 39 i=min(m,n),1,-1
!        l=i+1
!        g=w(i)
!        do 33 j=l,n
!          a(i,j)=0.0
!33      continue
!        if(g.ne.0.0)then
!          g=1.0/g
!          do 36 j=l,n
!            s=0.0
!            do 34 k=l,m
!              s=s+a(k,i)*a(k,j)
!34          continue
!            f=(s/a(i,i))*g
!            do 35 k=i,m
!              a(k,j)=a(k,j)+f*a(k,i)
!35          continue
!36        continue
!          do 37 j=i,m
!            a(j,i)=a(j,i)*g
!37        continue
!        else
!          do 38 j= i,m
!            a(j,i)=0.0
!38        continue
!        endif
!        a(i,i)=a(i,i)+1.0
!39    continue
!      do 49 k=n,1,-1
!        do 48 its=1,30
!          do 41 l=k,1,-1
!            nm=l-1
!            if((Abs(rv1(l))+anorm).eq.anorm)  goto 2
!            if((Abs(w(nm))+anorm).eq.anorm)  goto 1
!41        continue
!1         c=0.0
!          s=1.0
!          do 43 i=l,k
!            f=s*rv1(i)
!            rv1(i)=c*rv1(i)
!            if((Abs(f)+anorm).eq.anorm) goto 2
!            g=w(i)
!            h=Sqrt(f*f+g*g)
!            w(i)=h
!            h=1.0/h
!            c= (g*h)
!            s=-(f*h)
!            do 42 j=1,m
!              y=a(j,nm)
!              z=a(j,i)
!              a(j,nm)=(y*c)+(z*s)
!              a(j,i)=-(y*s)+(z*c)
!42          continue
!43        continue
!2         z=w(k)
!          if(l.eq.k)then
!            if(z.lt.0.0)then
!              w(k)=-z
!              do 44 j=1,n
!                v(j,k)=-v(j,k)
!44            continue
!            endif
!            goto 3
!          endif
!          if(its.eq.30) pause 'no convergence in svdcmp'
!          x=w(l)
!          nm=k-1
!          y=w(nm)
!          g=rv1(nm)
!          h=rv1(k)
!          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
!          g=Sqrt(f*f+1.0)
!          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
!          c=1.0
!          s=1.0
!          do 47 j=l,nm
!            i=j+1
!            g=rv1(i)
!            y=w(i)
!            h=s*g
!            g=c*g
!            z=Sqrt(f*f+h*h)
!            rv1(j)=z
!            c=f/z
!            s=h/z
!            f= (x*c)+(g*s)
!            g=-(x*s)+(g*c)
!            h=y*s
!            y=y*c
!            do 45 jj=1,n
!              x=v(jj,j)
!              z=v(jj,i)
!              v(jj,j)= (x*c)+(z*s)
!              v(jj,i)=-(x*s)+(z*c)
!45          continue
!            z=Sqrt(f*f+h*h)
!            w(j)=z
!            if(z.ne.0.0)then
!              z=1.0/z
!              c=f*z
!              s=h*z
!            endif
!            f= (c*g)+(s*y)
!            x=-(s*g)+(c*y)
!            do 46 jj=1,m
!              y=a(jj,j)
!              z=a(jj,i)
!              a(jj,j)= (y*c)+(z*s)
!              a(jj,i)=-(y*s)+(z*c)
!46          continue
!47        continue
!          rv1(l)=0.0
!          rv1(k)=f
!          w(k)=x
!48      continue
!3       continue
!49    continue
!      return
!      END	SUBROUTINE svdcmp
C
C -----------------------------------------------------------------
C  Solve A.X=B for a vector X, where A is specified by the arrays 
C  U, W, V as returned by SVDCMP. M and N are the logical dimensions 
C  of A, and will be equal for square matrices. MP and NP are the 
C  physical dimensions of A. B is the input right-hand side. X is 
C  the output solution vector. No input quantities are destroyed, 
C  so the routine may be called sequentially with different B's. 
C  M must be greater or equal to N; see SVDCMP
C
      SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
      IMPLICIT  REAL*8(A-H, O-Z)

      REAL*8    b(mp),u(mp,np),v(np,np),w(np),x(np)
      PARAMETER (NMAX=5000)
      DIMENSION tmp(NMAX) 
C      
      do 12 j=1,n
        s=0.
        if(w(j).ne.0.)then
          do 11 i=1,m
            s=s+u(i,j)*b(i)
11        continue
          s=s/w(j)
        endif
        tmp(j)=s
12    continue
      do 14 j=1,n
        s=0.
        do 13 jj=1,n
          s=s+v(j,jj)*tmp(jj)
13      continue
        x(j)=s
14    continue
      return
      END SUBROUTINE svbksb

       END MODULE MFUNC_mod



