C
C  Programs: CONVSB 
C
C  *****************************************************************
C  *                                                               *
C  *  The programme is used to converge the format of DATBDMS      *
C  *                                                               *
C  *****************************************************************
C       
	 SUBROUTINE CONVSB
C
       USE MVAR_mod
	 IMPLICIT NONE
C
       INTEGER IEB,IEL,IND,I,K,L,M,N,KK,NN,NND0
	 REAL(8) X,Y,Z,DX,DY,DZ,DR2
C
C
C ------------------------------------------------------------
C  NELEM : total number of elements 
C  NELEMF: total number of elements 
C  NELEMB: total number of elements 
C
C  NNODE : total number of nodes 
C  NNF   : total number of nodes on the free surface 
C  NNB   : total number of nodes on the body surface
C ------------------------------------------------------------
C    
C ** check nodes' situation 
C		  

	 WRITE(6,*) '  Inside CONVSB'
	 write(11,*) '  Inside CONVSB'
	 write(11,*) '  NELEMF=',NELEMF, '  NELEMB=',NELEMB
	 
	 
	 DO  IEL=1, NELEMF+NELEMB
	 write(11,1000) IEL,NCN(IEL)
	 write(11,1100) (XYZE(1,K,IEL),K=1,NCN(IEL))
	 write(11,1100) (XYZE(2,K,IEL),K=1,NCN(IEL))
	 write(11,1100) (XYZE(3,K,IEL),K=1,NCN(IEL))	 
	 ENDDO
	 
1100    FORMAT(8E15.6)
	 
	 
	  M=0
      IF (NELEMF .EQ. 0) THEN
        NNF=0
      
      
      ELSE
         DO 200 IEL=1,  NELEMF
           DO 180 I=1, NCN(IEL)
C
	     X=XYZE(1,I,IEL)
	     Y=XYZE(2,I,IEL)
	     Z=XYZE(3,I,IEL)
C
             DO 150 L=1, M
               IF( DABS( X-XYZTP(1,L) ).LT.1.0E-04  .AND.
     1             DABS( Y-XYZTP(2,L) ).LT.1.0E-04  .AND.
     2             DABS( Z-XYZTP(3,L) ).LT.1.0E-04 ) THEN
                 NCON(IEL, I)=L
                 GOTO 180
               END IF
 150         CONTINUE
C
C ** if the point (XPM, YPM) does not concide with any other
C    points had been checked,  give a new code
C
             M=M+1
             NCON(IEL, I)=M
             XYZTP(1,M)=X
             XYZTP(2,M)=Y
             XYZTP(3,M)=Z
C
180       CONTINUE
200     CONTINUE
C
       NNF=M     
	 write(11,  *) ' NNF=',NNF
       END IF    
!
! -----------------------------------------
!
         DO 300 IEB=1,  NELEMB
	      IEL=IEB+NELEMF
           DO 280 I=1, NCN(IEL)
C

	     X=XYZE(1,I,IEL)
	     Y=XYZE(2,I,IEL)
	     Z=XYZE(3,I,IEL)
C
             DO 250 L=1, M
               IF( DABS( X-XYZTP(1,L) ).LT.1.0E-04  .AND.
     1             DABS( Y-XYZTP(2,L) ).LT.1.0E-04  .AND.
     2             DABS( Z-XYZTP(3,L) ).LT.1.0E-04 ) THEN
                 NCON(IEL, I)=L
                 GOTO 280
               END IF
 250         CONTINUE
C
C ** if the point (XPM, YPM) does not concide with any other
C    points had been checked,  give a new code
C
             M=M+1
             NCON(IEL, I)=M
             XYZTP(1,M)=X
             XYZTP(2,M)=Y
             XYZTP(3,M)=Z
C
280       CONTINUE
300     CONTINUE
C
       NNODE=M
	 NNB=NNODE-NNF
	 write(11,*)  ' NNB=',NNB,' NNF=',NNF,' NNODE=',NNODE
	 write(11,*) 
!
! ===========================================================
!
!
       NNODED=NNODE
!
!  对每个节点先找出一个方向向量
!
       DO 350 IND=1,  NNODE
         DO 320 IEL=1,  NELEM
           DO 320 I=1,  NCN(IEL)
C
	      IF(NCON(IEL, I) .EQ. IND) THEN
		       DX=DXYZE(1,I,IEL)
	         DY=DXYZE(2,I,IEL)
			   DZ=DXYZE(3,I,IEL)
	         NCOND(IEL, I)=IND
	         DXYZTP(1,IND)=DX
	         DXYZTP(2,IND)=DY
	         DXYZTP(3,IND)=DZ
		     NNORMN(IND)=IND
   		     GOTO 350
		 ENDIF
320	  CONTINUE
350	 CONTINUE
!
!  对每个节点查找是否还有其他方向向量
!
	 DO 500 IND=1,  NNODE

	  NND0=NNODED

         DO 480 IEL=1,  NELEM
           DO 450 I=1,  NCN(IEL)
C
	      IF(NCON(IEL, I)  .EQ. IND) THEN	   
		    DX=DXYZE(1,I,IEL)
	        DY=DXYZE(2,I,IEL)
			DZ=DXYZE(3,I,IEL)
	      
		    DR2=(( DX-DXYZTP(1,IND))**2+
	1		     ( DY-DXYZTP(2,IND))**2+
     2		     ( DZ-DXYZTP(3,IND))**2)

	       IF( DR2 .LT. 1.0E-6)  THEN 
	         NCOND(IEL, I)=IND
		     GOTO 450
	       ENDIF
	      
		   DO KK=NND0+1,  NNODED    ! 检查是否与已定义的法向量相同
	         DR2=(( DX-DXYZTP(1,KK))**2+
	1	          ( DY-DXYZTP(2,KK))**2+
     2	          ( DZ-DXYZTP(3,KK))**2)
	       
		    IF( DR2 .LT. 1.0E-6) THEN 
	         NCOND(IEL, I)=KK
			 GOTO 450
	        ENDIF
	          
		   ENDDO
!
	       
		      NNODED=NNODED+1
	        NCOND(IEL, I)=NNODED
	        DXYZTP(1,NNODED)=DX
	        DXYZTP(2,NNODED)=DY
	        DXYZTP(3,NNODED)=DZ

	        XYZTP(1,NNODED)=XYZTP(1,IND)
	        XYZTP(2,NNODED)=XYZTP(2,IND)
	        XYZTP(3,NNODED)=XYZTP(3,IND)
		      NNORMN(NNODED)=IND

	      ENDIF
C
!
450       CONTINUE
!	    

480      CONTINUE
500     CONTINUE
C
!
! 
! ======================================================
!                   
C
        write(11,*) ' NELEM =',NELEM, '  NNODE=', NNODE
        write(11,*) ' NELEMF', NELEMF,'  NNF=',NNF
        write(11,*) ' NELEMB=',NELEMB,'  NNB=',NNB
!
!
! ======================================================
!	  
        write(11,*) Isys
        write(11,*)  NELEM, NNODE,NNODED,1
        write(11,*) ' 1  0  0.0  0.0  0.0'

	 DO N=1, NNODE
	  write(11,1010) N, 1, XYZTP(1,N),XYZTP(2,N),XYZTP(3,N)
	 ENDDO
C
	 DO N=1, NNODED
	  write(11,1010) N, 1, DXYZTP(1,N),DXYZTP(2,N),DXYZTP(3,N)
	 ENDDO
!
	 DO IEL=1, NELEM
	  write(11,1001) IEL, NCN(IEL)
        write(11,1005) (NCON(IEL, I), I=1, NCN(IEL))
	 ENDDO
!
	 DO IEL=1, NELEM
	  write(11,1001) IEL, NCN(IEL)
        WRITE(11,1005) (NCOND(IEL, I), I=1, NCN(IEL))
	 ENDDO
!	  
! ======================================================	  
!
1000   FORMAT(1X,8I6)
1001   FORMAT(1X,2I6,3F14.6)
1010   FORMAT(1X,I6,I4,3F14.6)
1005   FORMAT(8(1X,I6))
1020   FORMAT(1X,I6,I4,6F14.6)
C  
       END
C
