C  BODMASS 
C
C ******************************************************
C *                                                    *
C * Forming Mass and Hydro-restoring Matrices          *
C *                                                    *
C ******************************************************
C
        SUBROUTINE BODMASS
	  USE PVAR_MOD
	  USE MFUNC_MOD
        IMPLICIT   NONE  

        CHARACTER  CHAR
	  INTEGER I,J,IFL_MASS
	  REAL*8    XMS,XMC,YMC,ZMC
C               
         RMAS(:,:)=0.0   

C 
C  XC, YC, ZC   : Rotation centre
C  XMC, YMC, ZMC: Mass centre
C  RMAS : mass matrix         CRS  : hydro-restoring matrix
C  STKM : stiffness matrix    CDAP : damping matrix     
C          
        READ(7,20)   CHAR
        READ(7,*)    IFL_MASS 
C          
        READ(7,20)   CHAR
        READ(7,*)    XC,YC,ZC
!          
        READ(7,20)   CHAR
        READ(7,*)    XMC,YMC,ZMC
!
        READ(7,20)   CHAR
        READ(7,*)    ARE,XF,YF,XK2,YK2,XCF
C
        READ(7,20)   CHAR
        READ(7,*)    VOLM, XB,   YB,   ZB
!
! ----------------------------------------
!                     
!  Mass matrix
!
	  IF(IFL_MASS .EQ. 0) THEN
         READ(7,20)    CHAR 
         READ(7,*)     XMS
         DO I=1, 3
          READ(7,*)  (XIA(I,J),  J=1,3)
        END DO
!
	  RMAS=0.0d0
!
        RMAS(1,1)= XMS
        RMAS(1,5)= XMS*(ZMC-ZC)
        RMAS(1,6)=-XMS*(YMC-YC)
!
        RMAS(2,2)= XMS
        RMAS(2,4)=-XMS*(ZMC-ZC)
        RMAS(2,6)= XMS*(XMC-XC)
!
        RMAS(3,3)= XMS
        RMAS(3,4)= XMS*(YMC-YC)
        RMAS(3,5)=-XMS*(XMC-XC)
!
        RMAS(4,2)=-XMS*(ZMC-ZC)
        RMAS(4,3)= XMS*(YMC-YC)
        RMAS(4,4)= XIA(2,2)+XIA(3,3) +
	1	         XMS*(YMC-YC)*(YMC-YC)+XMS*(ZMC-ZC)*(ZMC-ZC)   
        RMAS(4,5)=-XIA(1,2)-XMS*(XMC-XC)*(YMC-YC)
        RMAS(4,6)=-XIA(1,3)-XMS*(XMC-XC)*(ZMC-ZC)
!
        RMAS(5,1)= XMS*(ZMC-ZC)
        RMAS(5,3)=-XMS*(XMC-XC)
        RMAS(5,4)=-XIA(2,1)-XMS*(YMC-YC)*(XMC-XC)          
        RMAS(5,5)= XIA(1,1)+XIA(3,3) +	
	1	         XMS*(XMC-XC)*(XMC-XC)+XMS*(ZMC-ZC)*(ZMC-ZC)   
        RMAS(5,6)=-XIA(2,3)-XMS*(YMC-YC)*(ZMC-ZC)
!
        RMAS(6,1)=-XMS*(YMC-YC)
        RMAS(6,2)= XMS*(XMC-XC)
        RMAS(6,4)=-XIA(3,1)-XMS*(ZMC-ZC)*(XMC-XC)
        RMAS(6,5)=-XIA(3,2)-XMS*(ZMC-ZC)*(YMC-YC)
        RMAS(6,6)= XIA(1,1)+XIA(2,2)+
	1	         XMS*(XMC-XC)*(XMC-XC)+XMS*(YMC-YC)*(YMC-YC)   
!
	 ELSE IF(IFL_MASS .EQ. 1) THEN
        READ(7,20)    CHAR
        DO I=1, 6
          READ(7,*)  (RMAS(I,J),   J=1,6)
        END DO
	 ENDIF
C 
        READ(7,20)    CHAR
        DO I=1, 6
          READ(7,*)  (CRS(I,J),   J=1,6)
        END DO
C
        READ(7,20)    CHAR
        DO I=1, 6
          READ(7,*)  (STKM(I,J),   J=1,6)
        END DO
C
        READ(7,20)    CHAR
        DO I=1, 6
          READ(7,*)  (CDMP(I,J),  J=1,6)
        END DO

C   
! -------------------------------------
!
	  
	 CALL INVERSE(RMAS,TRMAS,6,6)
c
        WRITE(6,143)
        WRITE(9,143) 
        DO   I=1,   6
          WRITE(6,102) I, (TRMAS(I,J),  J=1, 6) 
          WRITE(9,102) I, (TRMAS(I,J),  J=1, 6)
        END  DO     

C
20      FORMAT(A80)          
102     FORMAT('M',I1,'J',3x,6(2x,E10.3))  
122     FORMAT(10x,' CENTRE OF BODY MASS AND ROTATION',/,
     1          12X,'XC          YC          ZC ',/,
     2          5X,3F12.4 )
132     FORMAT(/,10x,'      CENTRE OF BUOYANCY',/,
     1          12X,'XB          YB          ZB ',/,
     2          5X,3F12.4 )
142     FORMAT(/,10x,'        MATRIX of MASS')
143     FORMAT(/,10x,'  INVERSE  MATRIX of MASS')
                                 
C


        RETURN                                      

        END 


