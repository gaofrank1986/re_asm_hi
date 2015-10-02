      
!  NORM_ELE0+TSING0
!  NORM_INT0+SING_INT0
! ======================================================
!
!   Integration on an element without source in itself
!   nor its symmetrical ones
!
! ======================================================
        include './add_on/ElemIntgl1.f90'
!
        SUBROUTINE NORM_ELE0(IELEM,XP,YP,ZP,AVAL,BVAL)
        
	    USE MVAR_MOD
        use tripole_mod
        IMPLICIT   NONE 
	  
	    INTEGER  IS,ND,J,NP,IELEM,NCNE
	    REAL*8  XP,YP,ZP 
        REAL*8  AVAL(4,8),BVAL(4,8)
!   
!            
	      NCNE=NCN(IELEM)
	                 
!	 Print *,' In TINBOD   PI=',PI             
!
          AVAL=0.0D0
          BVAL=0.0D0
!
        DO 100   IS=1,   NSYS  
          CALL NORM_INT0(IS,IELEM,NCNE,XP,YP,ZP,AVAL,BVAL)
100     CONTINUE
!
        RETURN
        END           

!
! ======================================================
!   Integration on an element with source in itself
!   or its mirror ones about any symmetrical axis
!
! ======================================================
!
       SUBROUTINE SING_ELE0(INODE,IELEM,NUMQUA,XP,YP,ZP,AVAL,BVAL)
	   USE MVAR_MOD
	   USE MFUNC_mod
!
        IMPLICIT   NONE  
!
	  INTEGER I,J,IS,IELEM,INODE,NODNUM,ND,NP,NUMQUA
        REAL*8  XP,YP,ZP,XYZT(3,8),DXYZT(3,8)
        REAL*8 AVAL(4,8),BVAL(4,8)
!
          AVAL= 0.0d0 
          BVAL= 0.0d0 
!

        DO 5     I=1,  NCN(IELEM)
        XYZT(1, I)  =  XYZE(1, I, IELEM)  
        XYZT(2, I)  =  XYZE(2, I, IELEM)  
        XYZT(3, I)  =  XYZE(3, I, IELEM)
	    
        DXYZT(1, I) = DXYZE(1, I, IELEM)  
        DXYZT(2, I) = DXYZE(2, I, IELEM)  
        DXYZT(3, I) = DXYZE(3, I, IELEM)

        IF(INODE.EQ.NCON(IELEM,I)) NODNUM=I
5       CONTINUE
!
        CALL TRIPOL(NODNUM,NCN(IELEM),XYZT,DXYZT)
!	  PRINT *,' AFTER TRIPOL'
!
!
        !write(*,*) 'NUMQUA=',NUMQUA,'nsys=',nsys
        !write(*,*) 'NODNUM=',nodnum
        !pause
        IF(NUMQUA.EQ.0)       THEN
         DO 100 IS=1,  NSYS
          IF(IS.EQ.1) THEN 
            CALL SING_INT0(IS,IELEM,XP,YP,ZP,AVAL,BVAL)
          ELSE IF(IS.NE.1 ) THEN   
            CALL NORM_INT0(IS,IELEM,NCN(IELEM),XP,YP,ZP,AVAL,BVAL)
          END IF
100      CONTINUE
!
        ELSE IF(NUMQUA.EQ.2) THEN
         DO 200 IS=1,NSYS     
          IF(IS.EQ.1.OR.IS.EQ.2) THEN
            CALL SING_INT0(IS,IELEM,XP,YP,ZP,AVAL,BVAL)
          ELSE IF(IS.EQ.3.OR.IS.EQ.4) THEN  
            CALL NORM_INT0(IS,IELEM,NCN(IELEM),XP,YP,ZP,AVAL,BVAL)
          END IF
!
200      CONTINUE  
!
        ELSE IF(NUMQUA.EQ.4) THEN
         DO 300  IS=1,  NSYS
          IF(IS.EQ.1.OR.IS.EQ.4) THEN   
            CALL SING_INT0(IS,IELEM,XP,YP,ZP,AVAL,BVAL)
          ELSE IF(IS.EQ.2.OR.IS.EQ.3) THEN  
            CALL NORM_INT0(IS,IELEM,NCN(IELEM),XP,YP,ZP,AVAL,BVAL)
          END IF
300      CONTINUE
!
        ELSE IF(NUMQUA.EQ.5) THEN
         DO 400 IS=1, NSYS  
            CALL SING_INT0(IS,IELEM,XP,YP,ZP,AVAL,BVAL)
400      CONTINUE
        ENDIF
!
        RETURN
        END
                           
!C ======================================================
!C
!C Integration on an element without source point
!C 
!C ======================================================
!C                      
      SUBROUTINE NORM_INT0(IS,IELEM,NCNE,XP,YP,ZP,AVAL,BVAL)
	  USE MVAR_MOD
	    USE MFUNC_mod
      IMPLICIT   NONE  
 
	  INTEGER IS,IELEM,N,NSAMB,NCNE,J,IP
      REAL*8  XP,YP,ZP,EX(4,4),EY(4,4)
	  REAL*8  X,Y,Z,X0,Y0,Z0   

      REAL*8  NX,NY,NZ,DGN
      REAL*8  AVAL(4,8),BVAL(4,8),GXF(4)
!   
	  DATA EX/  1.0d0,  1.0d0, -1.0d0, -1.0d0,    &
                1.0d0,  1.0d0, -1.0d0, -1.0d0,    &
               -1.0d0, -1.0d0,  1.0d0,  1.0d0,    &
               -1.0d0, -1.0d0,  1.0d0,  1.0d0/
!                                                  
	  DATA EY/  1.0d0, -1.0d0, -1.0d0,  1.0d0,    &
               -1.0d0,  1.0d0,  1.0d0, -1.0d0,    &
               -1.0d0,  1.0d0,  1.0d0, -1.0d0,    &
                1.0d0, -1.0d0, -1.0d0,  1.0d0/
!
!
!	 PRINT *,' IN  NORM_INT0'
!

	   X0=EX(IS,1)*XP
	   Y0=EY(IS,1)*YP
	   Z0= ZP

        NSAMB=16
        IF(NCNE.EQ.6)   NSAMB=4

        DO 100    N=1,   NSAMB     

	   X =SAMBXY(IELEM,N,1)
	   Y =SAMBXY(IELEM,N,2)
	   Z =SAMBXY(IELEM,N,3)

	   CALL	TGRN (H,X,X0,Y,Y0,Z,Z0,GXF) 
! 
         NX=DSAMB(IELEM,N,1)
         NY=DSAMB(IELEM,N,2)
         NZ=DSAMB(IELEM,N,3)
                                           
         DGN=GXF(2)*Nx+GXF(3)*Ny+GXF(4)*Nz           
!
        DO  J=1,   NCNE
         BVAL(IS,J)=BVAL(IS,J)+GXF(1)*SAMB(IELEM,N,J)
         AVAL(IS,J)=AVAL(IS,J)+DGN*SAMB(IELEM,N,J)
        ENDDO
!
!
100     CONTINUE
!
        RETURN
        END
!
! ===========================================================
! Integration on an element with source point
!
       SUBROUTINE SING_INT0(IS,IELEM,XP,YP,ZP,AVAL,BVAL)
	   USE MVAR_MOD
	    USE MFUNC_mod
	   USE TRVAR_MOD
       IMPLICIT   NONE  
!
	   INTEGER IS,IELEM,N,J,IP       
	   REAL*8  XP,YP,ZP,EX(4,4),EY(4,4)
	   REAL*8  X,Y,Z,X0,Y0,Z0      
	   REAL*8  NX,NY,NZ,DGN
       REAL*8  AVAL(4,8),BVAL(4,8),GXF(4)
!
	  DATA EX/  1.0d0,  1.0d0, -1.0d0, -1.0d0,    &
                1.0d0,  1.0d0, -1.0d0, -1.0d0,    &
               -1.0d0, -1.0d0,  1.0d0,  1.0d0,    &
               -1.0d0, -1.0d0,  1.0d0,  1.0d0/
!                                                  
	  DATA EY/  1.0d0, -1.0d0, -1.0d0,  1.0d0,    &
               -1.0d0,  1.0d0,  1.0d0, -1.0d0,    &
               -1.0d0,  1.0d0,  1.0d0, -1.0d0,    &
                1.0d0, -1.0d0, -1.0d0,  1.0d0/
!
!   
201	  FORMAT(' In SING_INT0    XP, YP, ZP=',3F12.4)
!
 
	  X0=EX(IS,1)*XP
	  Y0=EY(IS,1)*YP
	  Z0= ZP

      DO 130 N=1, NOSAMP
	
	  X =XYNOD(1,N)
	  Y =XYNOD(2,N)
	  Z =XYNOD(3,N)

	  CALL	TGRN (H,X,X0,Y,Y0,Z,Z0,GXF) 
!
        DGN=GXF(2)*DXYNOD(1,N)+GXF(3)*DXYNOD(2,N)+   &
              GXF(4)*DXYNOD(3,N)             

!	  write(6,*) ' DGN=',DGN

        DO J=1, NCN(IELEM)
         AVAL(IS,J)=AVAL(IS,J)+DGN*SAMNOD(N,J)
         BVAL(IS,J)=BVAL(IS,J)+GXF(1)*SAMNOD(N,J)
       ENDDO
!

130     CONTINUE
!
        RETURN
        END


