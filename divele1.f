!
!*******************************************************************
!动态加密高斯点

      RECURSIVE SUBROUTINE Gauss1(NCNE,NSAM,EL_MIN,AL,
	1			               XYZC,DXYZC,SIC,ETC,XYZF) 
	IMPLICIT NONE
!
	INTEGER,INTENT(IN)::  NCNE
	INTEGER  NSAM,IE

	REAL*8,INTENT(IN):: XYZC(3,8),DXYZC(3,8),SIC(8),ETC(8)
	REAL*8 XYZC1(3,8),DXYZC1(3,8),XYZC2(3,8,4),DXYZC2(3,8,4)
	REAL*8 SIC1(8),ETC1(8),SIC2(8,4),ETC2(8,4)
	
	REAL*8 XYZF(3),EL,EL_MIN,DIST,AL
!
! XYZF: Surce point
! AL: a factor
!
	
C 计算单元的特征长度EL和空间点（XP,YP,ZP）到单元各节点的最小距离DIST 

	CALL ELDIS(EL,DIST,NCNE,XYZF,XYZC)
	
	IF (DIST.GT.AL*EL .or. EL. LT. EL_MIN) THEN !最小距离大于特征长度
		CALL BODDIV(ncne,NSAM,XYZC,DXYZC,SIC,ETC)
	ELSE
		!细分单元
		CALL ELEDIV(ncne,XYZC,XYZC2,DXYZC,DXYZC2,
	1		 	    SIC,ETC,SIC2,ETC2)
			                
		DO IE=1, 4
		XYZC1(:,:) = XYZC2(:,:,IE)
		DXYZC1(:,:)=DXYZC2(:,:,IE)

		SIC1(:)=SIC2(:,IE)
		ETC1(:)=ETC2(:,IE)

		CALL gauss1(NCNE,NSAM,EL_MIN,AL,
	1		    	XYZC1,DXYZC1,SIC1,ETC1,XYZF)
		END DO
	END IF
	
	END



!
!
!
C **************************************************************
C *					                						 *
C *    Compute the characteristic length of an element         *
C *											                 *
C **************************************************************
C
	 SUBROUTINE EL_LENGTH(EL,NN,XYZC)
			
	 INTEGER NN,N
	 REAL*8  XYZC(3,NN)
	 REAL*8  EL,ELN(NN)
	 
C
	 IF (NN .EQ. 8) THEN
c	  
	  ELN(1)=DSQRT((XYZC(1,1)-XYZC(1,2))**2+
	1	           (XYZC(2,1)-XYZC(2,2))**2+
     2		       (XYZC(3,1)-XYZC(3,2))**2 )
	  ELN(2)=DSQRT((XYZC(1,2)-XYZC(1,3))**2+
	1	           (XYZC(2,2)-XYZC(2,3))**2+
     2		       (XYZC(3,2)-XYZC(3,3))**2 )
	  ELN(3)=DSQRT((XYZC(1,3)-XYZC(1,4))**2+
	1	           (XYZC(2,3)-XYZC(2,4))**2+
     2	 	       (XYZC(3,3)-XYZC(3,4))**2 )
	  ELN(4)=DSQRT((XYZC(1,4)-XYZC(1,5))**2+
	1	           (XYZC(2,4)-XYZC(2,5))**2+
     2		       (XYZC(3,4)-XYZC(3,5))**2 )
	  ELN(5)=DSQRT((XYZC(1,5)-XYZC(1,6))**2+
	1	           (XYZC(2,5)-XYZC(2,6))**2+
     2		       (XYZC(3,5)-XYZC(3,6))**2 )
	  ELN(6)=DSQRT((XYZC(1,6)-XYZC(1,7))**2+
	1	           (XYZC(2,6)-XYZC(2,7))**2+
     2	 	       (XYZC(3,6)-XYZC(3,7))**2 )
	  ELN(7)=DSQRT((XYZC(1,7)-XYZC(1,8))**2+
	1	           (XYZC(2,7)-XYZC(2,8))**2+
     2		       (XYZC(3,7)-XYZC(3,8))**2 )
	  ELN(8)=DSQRT((XYZC(1,8)-XYZC(1,1))**2+
	1	           (XYZC(2,8)-XYZC(2,1))**2+
     2		       (XYZC(3,8)-XYZC(3,1))**2 )

	ELSE IF(NN .EQ. 6) THEN
        ELN(1)=DSQRT((XYZC(1,1)-XYZC(1,4))**2+
	1	           (XYZC(2,1)-XYZC(2,4))**2+
     2	           (XYZC(3,1)-XYZC(3,4))**2 )
	  ELN(2)=DSQRT((XYZC(1,4)-XYZC(1,2))**2+
	1	           (XYZC(2,4)-XYZC(2,2))**2+
     2		       (XYZC(3,4)-XYZC(3,2))**2 )
	  ELN(3)=DSQRT((XYZC(1,2)-XYZC(1,5))**2+
	1	           (XYZC(2,2)-XYZC(2,5))**2+
     2	 	       (XYZC(3,2)-XYZC(3,5))**2 )
	  ELN(4)=DSQRT((XYZC(1,5)-XYZC(1,3))**2+
	1	           (XYZC(2,5)-XYZC(2,3))**2+
     2		       (XYZC(3,5)-XYZC(3,3))**2 )
	  ELN(5)=DSQRT((XYZC(1,3)-XYZC(1,6))**2+
	1	           (XYZC(2,3)-XYZC(2,6))**2+
     2		       (XYZC(3,3)-XYZC(3,6))**2 )
	  ELN(6)=DSQRT((XYZC(1,6)-XYZC(1,1))**2+
	1	           (XYZC(2,6)-XYZC(2,1))**2+
     2	 	       (XYZC(3,6)-XYZC(3,1))**2 )
	 END IF
C
	 el=0.0
	 DO N=1, NN
          el=el+ELN(n)
	 ENDDO
!
	 el=el/nn*2.0
	 	
		
	END


C *********************************************
C *											*
C *   Compute the characteristic length of 	*
C * an element, and the distances from the 	*
C * point to the (corner) nodes 		        *
C *											*
C *********************************************
C
	 SUBROUTINE ELDIS(EL,DIST,NN,XYZF,XYZC)
	IMPLICIT NONE
			
	 INTEGER NN,N
	 REAL*8  XYZF(3),XYZC(3,NN)
	 REAL*8  EL,ELN(NN),DIST,DISTN(NN)
	 
C
	 IF (NN .EQ. 8) THEN
c	  
	  ELN(1)=DSQRT((XYZC(1,1)-XYZC(1,2))**2+
	1	           (XYZC(2,1)-XYZC(2,2))**2+
     2		       (XYZC(3,1)-XYZC(3,2))**2 )
	  ELN(2)=DSQRT((XYZC(1,2)-XYZC(1,3))**2+
	1	           (XYZC(2,2)-XYZC(2,3))**2+
     2		       (XYZC(3,2)-XYZC(3,3))**2 )
	  ELN(3)=DSQRT((XYZC(1,3)-XYZC(1,4))**2+
	1	           (XYZC(2,3)-XYZC(2,4))**2+
     2	 	       (XYZC(3,3)-XYZC(3,4))**2 )
	  ELN(4)=DSQRT((XYZC(1,4)-XYZC(1,5))**2+
	1	           (XYZC(2,4)-XYZC(2,5))**2+
     2		       (XYZC(3,4)-XYZC(3,5))**2 )
	  ELN(5)=DSQRT((XYZC(1,5)-XYZC(1,6))**2+
	1	           (XYZC(2,5)-XYZC(2,6))**2+
     2		       (XYZC(3,5)-XYZC(3,6))**2 )
	  ELN(6)=DSQRT((XYZC(1,6)-XYZC(1,7))**2+
	1	           (XYZC(2,6)-XYZC(2,7))**2+
     2	 	       (XYZC(3,6)-XYZC(3,7))**2 )
	  ELN(7)=DSQRT((XYZC(1,7)-XYZC(1,8))**2+
	1	           (XYZC(2,7)-XYZC(2,8))**2+
     2		       (XYZC(3,7)-XYZC(3,8))**2 )
	  ELN(8)=DSQRT((XYZC(1,8)-XYZC(1,1))**2+
	1	           (XYZC(2,8)-XYZC(2,1))**2+
     2		       (XYZC(3,8)-XYZC(3,1))**2 )

	ELSE IF(NN .EQ. 6) THEN
        ELN(1)=DSQRT((XYZC(1,1)-XYZC(1,4))**2+
	1	           (XYZC(2,1)-XYZC(2,4))**2+
     2	           (XYZC(3,1)-XYZC(3,4))**2 )
	  ELN(2)=DSQRT((XYZC(1,4)-XYZC(1,2))**2+
	1	           (XYZC(2,4)-XYZC(2,2))**2+
     2		       (XYZC(3,4)-XYZC(3,2))**2 )
	  ELN(3)=DSQRT((XYZC(1,2)-XYZC(1,5))**2+
	1	           (XYZC(2,2)-XYZC(2,5))**2+
     2	 	       (XYZC(3,2)-XYZC(3,5))**2 )
	  ELN(4)=DSQRT((XYZC(1,5)-XYZC(1,3))**2+
	1	           (XYZC(2,5)-XYZC(2,3))**2+
     2		       (XYZC(3,5)-XYZC(3,3))**2 )
	  ELN(5)=DSQRT((XYZC(1,3)-XYZC(1,6))**2+
	1	           (XYZC(2,3)-XYZC(2,6))**2+
     2		       (XYZC(3,3)-XYZC(3,6))**2 )
	  ELN(6)=DSQRT((XYZC(1,6)-XYZC(1,1))**2+
	1	           (XYZC(2,6)-XYZC(2,1))**2+
     2	 	       (XYZC(3,6)-XYZC(3,1))**2 )
	 END IF
C
	 el=0.0
	 DO 120 N=1, NN
	    DISTN(N)=DSQRT((XYZF(1)-XYZC(1,N))**2+
	1		           (XYZF(2)-XYZC(2,N))**2+
     2                   (XYZF(3)-XYZC(3,N))**2 ) 
       el=el+ELN(n)
120	 CONTINUE
	 el=el/nn*2.0
	 	
	DIST=DISTN(1)		
	DO N=2, NN
	IF(DISTN(N).LE.DIST) DIST=DISTN(N)
	END DO
		
	END


C ****************************************
C Dividing the element into smaller ones
C 细分单元
C
C ****************************************
C
      SUBROUTINE ELEDIV(NCNE,XYZC,XYZC2,DXYZC,DXYZC2,
	1           SIC,ETC,SIC2,ETC2)
C
	USE MFUNC_mod   
      IMPLICIT NONE
	INTEGER NCNE,I,N,IE,LJ,LK,IP
      REAL*8  SI,ETA,XIN(4,8),ETN(4,8),TRXIN(4,6),TRETN(4,6)
      REAL*8  SF(8),DSF(2,8)
	REAL*8  XYZC(3,8),XYZC2(3,8,4),DXYZC(3,8),DXYZC2(3,8,4)
	REAL*8  SIC(8),ETC(8),SIC2(8,4),ETC2(8,4)

C  NCNE: number of nodes in the element
C
C  Coordinates of nodes of divided elements
       ! 四边形单元
       DATA (XIN(1,I),I=1,8)/-1.0D0,-0.5D0, 0.0D0, 0.0D0,
	1		                0.0D0,-0.5D0,-1.0D0,-1.0D0  /
       DATA (ETN(1,I),I=1,8)/-1.0D0,-1.0D0,-1.0D0,-0.5D0,
	1		                0.0D0, 0.0D0, 0.0D0,-0.5D0  /
c
	 DATA (XIN(2,I),I=1,8)/ 0.0D0, 0.5D0, 1.0D0, 1.0D0,
	1		                1.0D0, 0.5D0, 0.0D0, 0.0D0  /
	 DATA (ETN(2,I),I=1,8)/-1.0D0,-1.0D0,-1.0D0,-0.5D0,
 	1		                0.0D0, 0.0D0, 0.0D0,-0.5D0  /
c
	 DATA (XIN(3,I),I=1,8)/ 0.0D0, 0.5D0, 1.0D0, 1.0D0,
	1		                1.0D0, 0.5D0, 0.0D0, 0.0D0  /
	 DATA (ETN(3,I),I=1,8)/ 0.0D0, 0.0D0, 0.0D0, 0.5D0,
	1	                    1.0D0, 1.0D0, 1.0D0, 0.5D0  /
c
	 DATA (ETN(4,I),I=1,8)/ 0.0D0, 0.0D0, 0.0D0, 0.5D0,
 	1	                    1.0D0, 1.0D0, 1.0D0, 0.5D0  /
C
	 DATA (XIN(4,I),I=1,8)/-1.0D0,-0.5D0, 0.0D0, 0.0D0,
	1		                0.0D0,-0.5D0,-1.0D0,-1.0D0  /
c
       ! 三角形单元
       DATA (TRXIN(1,I),I=1,6)/0.0D0, 0.5D0, 0.0D0,
	1		                 0.25D0,0.25D0,0.0D0  /
       DATA (TRETN(1,I),I=1,6)/0.5D0, 0.5D0, 1.0D0,
	1		                 0.5D0, 0.75D0,0.75D0  /
c
	 DATA (TRXIN(2,I),I=1,6)/0.0D0, 0.50D0, 0.0D0,
	1		                 0.25D0,0.25D0, 0.0D0  /
	 DATA (TRETN(2,I),I=1,6)/0.0D0, 0.00D0, 0.5D0,
 	1		                 0.0D0, 0.25D0, 0.25D0  /
c
	 DATA (TRXIN(3,I),I=1,6)/0.5D0,  0.5D0,  0.0D0,
	1		                 0.5D0,  0.25D0, 0.25D0 /
	 DATA (TRETN(3,I),I=1,6)/0.0D0,  0.5D0,  0.5D0,
	1	                     0.25D0, 0.5D0,  0.25D0  /
c
	 DATA (TRXIN(4,I),I=1,6)/0.5D0,  1.0D0,  0.5D0,
	1		                 0.75D0, 0.75D0, 0.5D0  /
	 DATA (TRETN(4,I),I=1,6)/0.0D0,  0.0D0,  0.5D0,
 	1	                     0.0D0,  0.25D0, 0.25D0  /
c
C =====================================================
C
       SIC2=0.0d0
       ETC2=0.0d0


	 IF(NCNE .EQ. 8) THEN

	 DO 200 IE=1, 4
       DO 200 N=1, NCNE
C
C **  Calculate the shape function at the nodes of new elements
C
       SI =XIN(IE,N)
       ETA=ETN(IE,N)
	
C	 PRINT *,' SI=',SI,' ETA=',ETA
C                                      
	CALL SPFUNC8(SI,ETA,SF,DSF)
C
C -----------------------------
C
      DO 40 LJ=1,3
        XYZC2(LJ,N,IE)=0.0D0
       DO 20 LK=1,  NCNE
         XYZC2(LJ,N,IE)=XYZC2(LJ,N,IE)+SF(LK)*XYZC(LJ, LK)
20    CONTINUE
40    CONTINUE

C	 WRITE(6,41) XYZC2(1,N,IE),XYZC2(2,N,IE),XYZC2(3,N,IE)
C41	 FORMAT(3F14.5)

C
C ----------------------
C
      DO 80 LJ=1,3
        DXYZC2(LJ,N,IE)=0.0D0
       DO 60 LK=1,  NCNE
         DXYZC2(LJ,N,IE)=DXYZC2(LJ,N,IE)+SF(LK)*DXYZC(LJ, LK)
60    CONTINUE
80    CONTINUE
C
C -------------------
C
       DO 90 LK=1,  NCNE
         SIC2(N,IE)=SIC2(N,IE)+SF(LK)*SIC(LK)
         ETC2(N,IE)=ETC2(N,IE)+SF(LK)*ETC(LK)
90    CONTINUE

200	 CONTINUE
C
C =============================================================
C
	 ELSE IF(NCNE .EQ. 6) THEN

	 DO 400 IE=1, 4
C
       DO 400 N=1, NCNE
C
C **  calculate the shape function at the sampling points
C
       SI =TRXIN(IE,N)
       ETA=TRETN(IE,N)
C                                      
	 CALL SPFUNC6(SI,ETA,SF,DSF)
C
C -----------------------------
C
      DO 240 LJ=1,3
        XYZC2(LJ,N,IE)=0.0D0
       DO 220 LK=1,  NCNE
         XYZC2(LJ,N,IE)=XYZC2(LJ,N,IE)+SF(LK)*XYZC(LJ, LK)
220    CONTINUE
240    CONTINUE
C
C ----------------------
C
      DO 280 LJ=1,3
        DXYZC2(LJ,N,IE)=0.0D0
       DO 260 LK=1,  NCNE
         DXYZC2(LJ,N,IE)=DXYZC2(LJ,N,IE)+SF(LK)*DXYZC(LJ, LK)
260    CONTINUE
280    CONTINUE
C
C -------------------
C
       DO 290 LK=1,  NCNE
         SIC2(N,IE)=SIC2(N,IE)+SF(LK)*SIC(LK)
         ETC2(N,IE)=ETC2(N,IE)+SF(LK)*ETC(LK)
290    CONTINUE

C
400	 CONTINUE
C
	 END IF
C
	 END

C
C ****************************************
! Compute Gauss 
! 计算Gauss样点，Jaccobi行列式等
C
C ****************************************
C
      SUBROUTINE BODDIV(NCNE,NSAM,XYZC,DXYZC,SIC,ETC)	  
	USE MFUNC_mod         
      IMPLICIT NONE
C
	 INTEGER  NCNE,NSAM
	 REAL*8  XYZC(3,8),DXYZC(3,8),SIC(8),ETC(8)

	 INTEGER ISI,IETA,LI,LJ,LK,IP,J
       REAL*8  XITSI(7),XITETA(7),WIT(7),XIQ(4),WIQ(4)
	 REAL*8  SI,ETA,SIT,ETAT,DUMM,DET1,DET2,DET3,DET
C
       REAL*8 SF(8),DSF(2,8),SF1(8),DSF1(2,8),XJ(3,3)
       REAL*8 SAMXYZ(3,4000),DSAXYZ(3,4000),SAMBF(4000,0:8)
C
       COMMON/DIVIDE/SAMXYZ,DSAXYZ,SAMBF
C
C ** XITSI, XITETA store the Gauss-Legendre sampling points(7) for the 
C    triangular element
C
C
C ** matrix WIT store the Gauss-Legendre weighting factors for the 
C    triangular element
C
      DATA XITSI/ 0.101286507323456D0, 0.797426985353087D0,
     1            0.101286507323456D0, 0.470142064105115D0,
     1            0.470142064105115D0, 0.059715871789770D0,
     1            0.333333333333333D0/

      DATA XITETA/0.101286507323456D0, 0.101286507323456D0,
     1            0.797426985353087D0, 0.059715871789770D0,
     1            0.470142064105115D0, 0.470142064105115D0,
     1            0.333333333333333D0/

      DATA WIT/   0.062969590272414D0, 0.062969590272414D0,
     1            0.062969590272414D0, 0.066197076394253D0,
     1            0.066197076394253D0, 0.066197076394253D0,
     1            0.112500000000000D0/

C
C ** matrix XIQ store the Gauss-Legendre sampling points for the 
C    quadrilateral elements
C
       DATA XIQ/ 0.861136311594052D+00, 0.339981043584856D+00,
     +          -0.339981043584856D+00,-0.861136311594052D+00/
C
C ** matrix WIQ store the Gauss-Legendre weighting factors for the 
C     quadrilateral elements
C
       DATA WIQ/ 0.347854845137454D+00, 0.652145154862546D+00,
     +           0.652145154862546D+00, 0.347854845137454D+00/
C

	 IF(NCNE .EQ. 8) THEN
C
       DO 200 ISI =1, 4
       DO 200 IETA=1, 4
C
	 NSAM=NSAM+1
C
C **  calculate the shape function at the sampling points
C
      SI =XIQ(ISI)
      ETA=XIQ(IETA)
C                                      
	CALL SPFUNC8(SI,ETA,SF,DSF)
C
C ** evaluate the Jacobian matrix at (SI,ETA)
C
      DO 130 LI=1,2
      DO 130 LJ=1,3
       DUMM=0.0D0
       DO 140 LK=1,  NCNE
        DUMM=DUMM+DSF(LI,LK)*XYZC(LJ, LK)
140    CONTINUE
130    XJ(LI,LJ)=DUMM
C
C ** compute the determinant of the Jacobian maxtix at (SI,ETA)
C
      DET1=XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2) 
      DET2=XJ(1,1)*XJ(2,3)-XJ(1,3)*XJ(2,1) 
      DET3=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1) 
      DET=DSQRT(DET1*DET1+DET2*DET2+DET3*DET3)
C
C ** transform the local coordinates of the sampling points to 
C    global coordinates
C
      DO 160 LJ=1,  3
      SAMXYZ(LJ,NSAM)=0.0D0
      DO 155 LK=1, NCNE
       SAMXYZ(LJ,NSAM)=SAMXYZ(LJ,NSAM)+SF(LK)*XYZC(LJ,LK)
155   CONTINUE
160   CONTINUE
C
      DO 180 LJ=1,3
      DSAXYZ(LJ,NSAM)=0.0D0
      DO 175 LK=1, NCNE
       DSAXYZ(LJ,NSAM)=DSAXYZ(LJ,NSAM)+SF(LK)*DXYZC(LJ,LK)
175   CONTINUE
180   CONTINUE
C              
C ** integration weighting
C
	 SAMBF(NSAM,0)=WIQ(ISI)*WIQ(IETA)*DET

	 SIT=0.0d0
	 ETAT=0.0d0
       DO  LK=1, NCNE
        SIT=SIT+SF(LK)*SIC(LK)
        ETAT=ETAT+SF(LK)*ETC(LK)
	 ENDDO

	 CALL SPFUNC8(SIT,ETAT,SF1,DSF1)

       DO  LK=1, NCNE
	 SAMBF(NSAM,LK)=SAMBF(NSAM,0)*SF1(LK)
	 ENDDO

!	 print *,' WIQ(ISI) =',WIQ(ISI)
!	 print *,' WIQ(IETA)=',WIQ(IETA)
!	 Print *,' NSAM=',NSAM,' DET=',DET
	  
C
200	 CONTINUE
C
C --------------------------------
C
	 ELSE IF(NCNE .EQ. 6) THEN

       DO 400 J =1, 7

	 NSAM=NSAM+1
C
C **  calculate the shape function at the sampling points
C
      SI =XITSI(J)
      ETA=XITETA(J)
C
      CALL SPFUNC6(SI,ETA,SF,DSF)
C                                      
C ** evaluate the Jacobian matrix at (SI,ETA)
C
      DO 330 LI=1,2
      DO 330 LJ=1,3
       DUMM=0.0D0
       DO 320 LK=1,  NCNE
        DUMM=DUMM+DSF(LI,LK)*XYZC(LJ, LK)
320    CONTINUE
330    XJ(LI,LJ)=DUMM
C
C ** compute the determinant of the Jacobian maxtix at (SI,ETA)
C
      DET1=XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2) 
      DET2=XJ(1,1)*XJ(2,3)-XJ(1,3)*XJ(2,1) 
      DET3=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1) 
      DET=DSQRT(DET1*DET1+DET2*DET2+DET3*DET3)
C
C ** transform the local coordinates of the sampling points to 
C    global coordinates
C
      DO 370 LJ=1,3
      SAMXYZ(LJ,NSAM)=0.0D0
      DO 360 LK=1, NCNE
       SAMXYZ(LJ,NSAM)=SAMXYZ(LJ,NSAM)+SF(LK)*XYZC(LJ,LK)
360   CONTINUE
370   CONTINUE
C
      DO 380 LJ=1,3
      DSAXYZ(LJ,NSAM)=0.0D0
      DO 375 LK=1, NCNE
       DSAXYZ(LJ,NSAM)=DSAXYZ(LJ,NSAM)+SF(LK)*DXYZC(LJ,LK)
375   CONTINUE
380   CONTINUE
C              

390	CONTINUE
C    
C ** calculate the free surface boundary condition*WIT*DET
C
C
	  SAMBF(NSAM,0)=WIT(J)*DET
C
C
	 SIT=0.0d0
	 ETAT=0.0d0
       DO  LK=1, NCNE
        SIT=SIT+SF(LK)*SIC(LK)
        ETAT=ETAT+SF(LK)*ETC(LK)
	 ENDDO

	 CALL SPFUNC6(SIT,ETAT,SF1,DSF1)

       DO  LK=1, NCNE
	 SAMBF(NSAM,LK)=SAMBF(NSAM,0)*SF1(LK)
	 ENDDO


400	 CONTINUE
C
	 END IF
C
c	 Pause

	 END


