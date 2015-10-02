!       MESHFS4 + MESHBD 
!
C *******************************************************************
C *                                                                 *
C *  Read in the data file of the mesh on the free surface          *
C *                                                                 *
C *                                                                 *
C *******************************************************************
C 
        SUBROUTINE MESHFS4

	  USE MVAR_MOD
        IMPLICIT   NONE  

	  INTEGER IE,J,M
C
	  XYZE(3,1:8, 1:NELEMF)=0.0d0
!
        DO 100 IE=1, NELEMF
	  IETYPE(IE)=2
	  READ(3, *)    M, NCN(IE)
        READ(3, *) (XYZE(1,J,IE), J=1, NCN(IE))
        READ(3, *) (XYZE(2,J,IE), J=1, NCN(IE))

	  dampe(:,:)=0.0d0
C
100    CONTINUE
!
!
!  数据文件中物面法方向指入流体为正
!
	  DXYZE(1, 1:8, 1:NELEMF)= 0.0d0
	  DXYZE(2, 1:8, 1:NELEMF)= 0.0d0
	  DXYZE(3, 1:8, 1:NELEMF)= 1.0d0 
!	  DXYZE(3, 1:8, 1:NELEMF)=-1.0d0 
C

      RETURN
      END



C ***************************************************************
C *                                                             *
C *  Generate nodal and element data on the body surface of     *
C *  an arbitrary body                                          *
C *                                                             *
C ***************************************************************
C 
        SUBROUTINE MESHBD(NCOR)

	  USE MVAR_MOD
        IMPLICIT   NONE  

!	  INTEGER NELEM,NELEMB,NELEMF,NNODE,NNODED,NNB,NNBD,NNF,NNTCH


	  INTEGER NCOR,IPL,IPOLAR(50)
	  INTEGER I,IE,M,INODE,NCNN,K,KK
	  REAL*8  A1,A2,R1,R2,Z1
        REAL*8  XOFSET(50),YOFSET(50),ZOFSET(50)
!
! --------------------------------------------
!
        DO 5 I=1, NCOR
        READ(2,*)  M, IPOLAR(I), XOFSET(I),YOFSET(I),ZOFSET(I)
5       CONTINUE
!
        DO 10 INODE=1, NNB
        READ(2,*) M, IPL, R1, A1, Z1
        IF(IPOLAR(IPL).EQ.0) THEN
         XYZB(1,INODE)=R1+XOFSET(IPL)
         XYZB(2,INODE)=A1+YOFSET(IPL)
        ELSE IF(IPOLAR(IPL).EQ.1) THEN
         XYZB(1,INODE)=R1*DCOS(A1*PI/180.0D0)+XOFSET(IPL)
         XYZB(2,INODE)=R1*DSIN(A1*PI/180.0D0)+YOFSET(IPL)
        ENDIF
	   XYZB(3,INODE)=Z1+ZOFSET(IPL)
10      CONTINUE
!       
        DO 11 INODE=1, NNBD
        READ(2,*) M,IPL,R2,A2,DXYZB(3,INODE)

        IF(IPOLAR(IPL).EQ.0) THEN
          DXYZB(1,INODE)= R2
          DXYZB(2,INODE)= A2
        ELSE IF(IPOLAR(IPL).NE.0) THEN
          DXYZB(1,INODE)=R2*DCOS( A2*PI/180.0D0 )
          DXYZB(2,INODE)=R2*DSIN( A2*PI/180.0D0 )
        ENDIF
11      CONTINUE
C
C      
        DO 20 I=1, NELEMB
	  IE=I+NELEMF
	  IETYPE(IE)=1
        READ(2,*) M, NCN(IE)    
        READ(2,*) (NCONB(I,K), K=1, NCN(IE))
20    CONTINUE
!
	 DO 30 I=1,  NELEMB
	  IE=I+NELEMF
        READ(2,*)  M, NCNN
        READ(2,*) (NCONDB(I,K), K=1, NCN(IE))
30     CONTINUE
!
! ===============================================================
!    转换网格数据格式，将水线上节点从物面上消除，归入到水面上
!
       DO 100 I=1, NELEMB
	  IE=I+NELEMF
	  DO 80 K=1, NCN(IE)
         XYZE(1,K,IE)=XYZB(1, NCONB(I,K))
         XYZE(2,K,IE)=XYZB(2, NCONB(I,K))
         XYZE(3,K,IE)=XYZB(3, NCONB(I,K))
80      CONTINUE
!
100    CONTINUE
!
       DO 200 I=1, NELEMB
	  IE=I+NELEMF
	  DO 180 K=1, NCN(IE)
         DXYZE(1,K,IE)=DXYZB(1,NCONDB(I,K))
         DXYZE(2,K,IE)=DXYZB(2,NCONDB(I,K))
         DXYZE(3,K,IE)=DXYZB(3,NCONDB(I,K))
180     CONTINUE
!
200    CONTINUE
!
       RETURN
       END


