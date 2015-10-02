!
! *****************************************************************
!
!    直接方法计算固角系数
!
!  对应的格林函数包含关于海底的像
! 可以考虑分别关于x轴和轴对称，以及同时关于x轴和y轴对称的情形
! *****************************************************************
!
       SUBROUTINE SOLIDANGLE(INODE,NNODE,NELEM,NCN,NCON,NODQUA,
     1                        H,XYZ,DXYZE,SANGLE)

        USE MFUNC_mod   
!	   
        IMPLICIT   NONE 
!	
       INTEGER,INTENT(IN):: INODE,NNODE,NELEM
       INTEGER,INTENT(IN):: NCN(NELEM),NCON(NELEM,8),NODQUA(NNODE)
       
       REAL*8,INTENT(IN)::  H,XYZ(3,NNODE),DXYZE(3,8,NELEM) 
       REAL*8,INTENT(OUT):: SANGLE
! 
	  INTEGER  IELEM,I,J,LK,LJ,LI,IEL,MINSIDE
	  INTEGER  MCELE,MCNT(100),MEND(100),MNOD(100)
	  INTEGER  NCONT(100),NSIDE(100,2)
	  INTEGER  NODS_EVEN(4),NODS_ODD(8),NODT_EVEN(3),NODT_ODD(6)

	  REAL*8   PI
	  REAL*8   DR,SI,ETA,DUM
	  REAL*8   SF(8),DSF(2,8),SIQ(12),ETAQ(12),DSIQ(12,2)

	  REAL*8   DXYZN(100,3),TXYZN(100,2,3),TXYZ(2,3),XJ(2,3)
!
       DATA PI/3.14159265359/

       DATA  SIQ/-1.0D0, 0.0D0, 0.0D0, 1.0D0, 1.0D0, 1.0D0,  
     1		      1.0D0, 0.0D0, 0.0D0,-1.0D0,-1.0D0,-1.0D0/ 
!
       DATA ETAQ/-1.0D0,-1.0D0,-1.0D0,-1.0D0, 0.0D0, 0.0D0,  
	1		      1.0D0, 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0/ 
!
	  DATA  DSIQ/ 1.0d0,-1.0d0, 1.0d0,-1.0d0,-1.0d0,-1.0d0, 
     1 	           -1.0d0, 1.0d0,-1.0d0, 1.0d0, 1.0d0, 1.0d0,      
     2              1.0d0, 1.0d0, 1.0d0, 1.0d0,-1.0d0, 1.0d0, 
     3 	           -1.0d0,-1.0d0,-1.0d0,-1.0d0, 1.0d0,-1.0d0/ 
!

       DATA NODS_EVEN/1,4,7,10/
       DATA NODS_ODD/2,3,5,6,8,9,11,12/
!
       DATA NODT_EVEN/1,2,3/
       DATA NODT_ODD/4,5,6,7,8,9/



      IF(INODE ==1) THEN
      Print *
      print *,'  Computing solid angle by direct method'
      print *,'The Program is developed by Prof. Bin Teng at DUT, China'
      Print *     
      ENDIF
       
!        WRITE(10,*)
!        WRITE(10,*) 'Inside SOLIDANGLE     INODE=',INODE  
!        WRITE(10,*) 'X, Y, Z=', XYZ(1,INODE),XYZ(2,INODE),XYZ(3,INODE)
!        WRITE(10,*) 'NODQUA(INODE)=',NODQUA(INODE)
!
! ============================================================
!
! MELE : Number of elements including the node      包含该节点的单元数量
! MCELE: Number of sub-elements including the node  包含该节点的子单元数量
! MCNT(100): codes the elements including the node  包含该节点的单元编号                       
! MEND(100): code the node in the element           该节点在单元中的旧序号（6、8节点）              
! MNOD(100): code the node in the element           该节点在单元中的新序号（9、12节点单元）             
! 

	 MCELE=0
	 DO 40 IELEM=1, NELEM
        DO 20  I=1, NCN(IELEM)
         IF(INODE .EQ. NCON(IELEM,I)) THEN
!
	     IF(MCELE .GT. 25) THEN
	      Print *,' There are ', MCELE,' elements which have the node'
	      Print *,' To inlarge the size of some arrays'
	      pause
	     ENDIF
!
	   IF(NCN(IELEM) .EQ. 8)  THEN

	     IF(I==1 .or. I==3 .or. I==5 .or. I==7)  THEN
	       MCELE=MCELE+1
	       MCNT(MCELE)=IELEM
	       MEND(MCELE)=I
	       MNOD(MCELE)=NODS_EVEN((I+1)/2)

	     ELSE IF(I==2 .or. I==4 .or. I ==6 .or. I==8)   THEN
	       MCELE=MCELE+1
	       MCNT(MCELE)=IELEM
	       MEND(MCELE)=I
	       MNOD(MCELE)=NODS_ODD(I-1)
!
	       MCELE=MCELE+1
	       MCNT(MCELE)=IELEM
	       MEND(MCELE)=I
	       MNOD(MCELE)=NODS_ODD(I)   	 
		   ENDIF
!
	   ELSE IF(NCN(IELEM) .EQ. 6)  THEN    

	     IF(I==1 .or. I==2 .or. I==3)  THEN
	       MCELE=MCELE+1
	       MCNT(MCELE)=IELEM
	       MEND(MCELE)=I
	       MNOD(MCELE)=NODT_EVEN(I)

	     ELSE IF(I==4 .or. I==5 .or. I ==6)   THEN
	       MCELE=MCELE+1
	       MCNT(MCELE)=IELEM
	       MEND(MCELE)=I
	       MNOD(MCELE)=NODT_ODD(2*(I-3)-1)
!
	       MCELE=MCELE+1
	       MCNT(MCELE)=IELEM
	       MEND(MCELE)=I
	       MNOD(MCELE)=NODT_ODD(2*(I-3))   	 
		   ENDIF
!	    
         ENDIF   ! NCN(IELEM) .EQ. 8 or 6
!
	   ENDIF
20      CONTINUE
40     CONTINUE
!

!  ================================================================
!
! DXYZN: 各单元在锲形体端点处的法向导数
! TXYZN: 各单元在锲形体端点处指出端点的切向量
!
		DXYZN=0.0d0
	  TXYZN=0.0d0
!
	  DO IEL=1, MCELE

	   IELEM=MCNT(IEL)
!
!        WRITE(16,*) 
!        WRITE(16,*) ' IEMEM=',IELEM,' IEL_NEW=',IEL,
!     1              ' MEND=',MEND(IEL),' MNOD=',MNOD(IEL)

		DXYZN(IEL,1)=DXYZE(1,MEND(IEL),IELEM)
		DXYZN(IEL,2)=DXYZE(2,MEND(IEL),IELEM)
		DXYZN(IEL,3)=DXYZE(3,MEND(IEL),IELEM)
!
         IF(NCN(IELEM).EQ.8) THEN
!
! ** Quadrilateral element **
!
         SI =SIQ(MNOD(IEL))
         ETA=ETAQ(MNOD(IEL))
!        
         CALL SPFUNC8(SI,ETA,SF,DSF)
  
        DO 140 LI=1,2
       DO 130 LJ=1,3
        DUM=0.0D0
        DO 120 LK=1, NCN(IELEM)
          DUM=DUM+DSF(LI,LK)*XYZ(LJ,NCON(IELEM,LK))
120     CONTINUE
	 XJ(LI,LJ)=DUM
130    CONTINUE
!
	 DR=SQRT(XJ(LI,1)*XJ(LI,1)+XJ(LI,2)*XJ(LI,2)+XJ(LI,3)*XJ(LI,3))

	  TXYZN(IEL,LI,1)=DSIQ(MNOD(IEL),LI)*XJ(LI,1)/DR
	  TXYZN(IEL,LI,2)=DSIQ(MNOD(IEL),LI)*XJ(LI,2)/DR
	  TXYZN(IEL,LI,3)=DSIQ(MNOD(IEL),LI)*XJ(LI,3)/DR
  
140    CONTINUE
!
       ELSE IF(NCN(IELEM).EQ.6) THEN

 
        CALL TANGEL6(IELEM,MNOD(IEL),NNODE,NELEM,NCON,XYZ,TXYZ)


       
         TXYZN(IEL,:,:)=TXYZ(:,:)
       ENDIF
      
      
 !        WRITE(16,112) (XYZ(1,NCON(IELEM,LK)), LK=1, NCN(IELEM))
 !        WRITE(16,114) (XYZ(2,NCON(IELEM,LK)), LK=1, NCN(IELEM))
 !        WRITE(16,116) (XYZ(3,NCON(IELEM,LK)), LK=1, NCN(IELEM))
!
112    FORMAT(' X=',8E13.5)
114    FORMAT(' Y=',8E13.5)
116    FORMAT(' Z=',8E13.5)

!           WRITE(16,*) '  TX=', TXYZN(IEL,1,1),
!     1      '  TY=', TXYZN(IEL,1,2),'  TZ=', TXYZN(IEL,1,3)
!           WRITE(16,*) '  TX=', TXYZN(IEL,2,1),
!     1      '  TY=', TXYZN(IEL,2,2),'  TZ=', TXYZN(IEL,2,3)
   
       ENDDO
       
!	   write(6,*) ' 50  MCELE=',MCELE
!         WRITE(6,*) '  Z=',XYZ(3,INODE),  '  H=',H!
!	   write(16,*) ' 50  MCELE=',MCELE
!         WRITE(16,*) '  Z=',XYZ(3,INODE),  '  H=',H
       

!    对于与水面处或海底相交的单元，补充关于z轴对称的单元
!
	 IF(ABS(XYZ(3,INODE)+H) .LT. 1.0D-6)  THEN
	   CALL XYSYM(DXYZN,TXYZN,MCELE)
	 ENDIF
!
!	   write(16,*) ' 60  MCELE=',MCELE
!	   write(6,*) ' 60  MCELE=',MCELE

!	   	 
	 IF(NODQUA(INODE) .EQ. 2) THEN
	   CALL ZXSYM(DXYZN,TXYZN,MCELE)
	 ELSE IF(NODQUA(INODE) .EQ. 4) THEN
	   CALL YZSYM(DXYZN,TXYZN,MCELE)
	 ELSE IF(NODQUA(INODE) .EQ. 5) THEN
	   CALL YZZXSYM(DXYZN,TXYZN,MCELE)
	 ENDIF

!	        write(16,*) ' 80  MCELE=',MCELE

	
	   CALL DOMAINFAS(DXYZN,TXYZN,MCELE,NCONT,NSIDE,MINSIDE)
   
!	           write(16,*) '  90  MCELE=',MCELE
	   
	 IF(MINSIDE .GT.0) THEN
	   CALL ANGLE0(DXYZN,TXYZN,MCELE,NCONT,NSIDE,SANGLE)

!	   WRITE(16,*)   '     Solid_ANGLE=',SANGLE
!	   WRITE(9,*)   '     Solid_ANGLE=',SANGLE

	   SANGLE=1.0d0-SANGLE/(4.0d0*PI)

!        write(10,*) '  DXYZN=',DXYZN
!        write(10,*) '  TXYZN=',TXYZN
!        write(10,*) '  MCELE=',MCELE
!        WRITE(10,*) '  NCONT=',NCONT
!        WRITE(10,*) '  NSIDE=',NSIDE
!	  WRITE(10,*)   '     Solid_ANGLE=',SANGLE
        
!        stop

       ELSE
         SANGLE=0.5
       ENDIF

!	  Print *,'     Solid_ANGLE=',SANGLE

!	  pause
!
       RETURN
       END           



!
! *****************************************************************
!
!   对于三角形单元，计算各节点处边的切向量
!
! *****************************************************************
!
        SUBROUTINE TANGEL6(IELEM,JNEW,NNODE,NELEM,NCON,XYZ,TXYZ)
!	  USE MVAR_MOD
!	   
        IMPLICIT   NONE 
!	  
       INTEGER,INTENT(IN):: IELEM,JNEW
       INTEGER,INTENT(IN):: NNODE,NELEM
       INTEGER,INTENT(IN):: NCON(NELEM,8)
      
       REAL*8,INTENT(IN)::  XYZ(3,NNODE) 
	 REAL*8,INTENT(OUT):: TXYZ(2,3)

	  REAL*8   DX,DY,DZ,DR
!            
!
! ** Trianglar element **
!
        IF (JNEW .EQ. 1)  THEN
          DX=XYZ(1,NCON(IELEM,4))-XYZ(1,NCON(IELEM,1))
          DY=XYZ(2,NCON(IELEM,4))-XYZ(2,NCON(IELEM,1))
          DZ=XYZ(3,NCON(IELEM,4))-XYZ(3,NCON(IELEM,1))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(1,1)=DX/DR
	    TXYZ(1,2)=DY/DR
	    TXYZ(1,3)=DZ/DR
!  
          DX=XYZ(1,NCON(IELEM,6))-XYZ(1,NCON(IELEM,1))
          DY=XYZ(2,NCON(IELEM,6))-XYZ(2,NCON(IELEM,1))
          DZ=XYZ(3,NCON(IELEM,6))-XYZ(3,NCON(IELEM,1))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(2,1)=DX/DR
	    TXYZ(2,2)=DY/DR
	    TXYZ(2,3)=DZ/DR        

        
        ELSE IF (JNEW .EQ. 2)  THEN
        
          DX=XYZ(1,NCON(IELEM,5))-XYZ(1,NCON(IELEM,2))
          DY=XYZ(2,NCON(IELEM,5))-XYZ(2,NCON(IELEM,2))
          DZ=XYZ(3,NCON(IELEM,5))-XYZ(3,NCON(IELEM,2))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(1,1)=DX/DR
	    TXYZ(1,2)=DY/DR
	    TXYZ(1,3)=DZ/DR
!  
          DX=XYZ(1,NCON(IELEM,4))-XYZ(1,NCON(IELEM,2))
          DY=XYZ(2,NCON(IELEM,4))-XYZ(2,NCON(IELEM,2))
          DZ=XYZ(3,NCON(IELEM,4))-XYZ(3,NCON(IELEM,2))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(2,1)=DX/DR
	    TXYZ(2,2)=DY/DR
	    TXYZ(2,3)=DZ/DR        

        ELSE IF (JNEW .EQ. 3)  THEN
          DX=XYZ(1,NCON(IELEM,6))-XYZ(1,NCON(IELEM,3))
          DY=XYZ(2,NCON(IELEM,6))-XYZ(2,NCON(IELEM,3))
          DZ=XYZ(3,NCON(IELEM,6))-XYZ(3,NCON(IELEM,3))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(1,1)=DX/DR
	    TXYZ(1,2)=DY/DR
	    TXYZ(1,3)=DZ/DR
!  
          DX=XYZ(1,NCON(IELEM,5))-XYZ(1,NCON(IELEM,3))
          DY=XYZ(2,NCON(IELEM,5))-XYZ(2,NCON(IELEM,3))
          DZ=XYZ(3,NCON(IELEM,5))-XYZ(3,NCON(IELEM,3))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(2,1)=DX/DR
	    TXYZ(2,2)=DY/DR
	    TXYZ(2,3)=DZ/DR        

        ELSE IF (JNEW .EQ. 4)  THEN
          DX=XYZ(1,NCON(IELEM,3))-XYZ(1,NCON(IELEM,4))
          DY=XYZ(2,NCON(IELEM,3))-XYZ(2,NCON(IELEM,4))
          DZ=XYZ(3,NCON(IELEM,3))-XYZ(3,NCON(IELEM,4))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(1,1)=DX/DR
	    TXYZ(1,2)=DY/DR
	    TXYZ(1,3)=DZ/DR
!  
          DX=XYZ(1,NCON(IELEM,1))-XYZ(1,NCON(IELEM,4))
          DY=XYZ(2,NCON(IELEM,1))-XYZ(2,NCON(IELEM,4))
          DZ=XYZ(3,NCON(IELEM,1))-XYZ(3,NCON(IELEM,4))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(2,1)=DX/DR
	    TXYZ(2,2)=DY/DR
	    TXYZ(2,3)=DZ/DR        
        ELSE IF (JNEW .EQ. 5)  THEN
          DX=XYZ(1,NCON(IELEM,2))-XYZ(1,NCON(IELEM,4))
          DY=XYZ(2,NCON(IELEM,2))-XYZ(2,NCON(IELEM,4))
          DZ=XYZ(3,NCON(IELEM,2))-XYZ(3,NCON(IELEM,4))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(1,1)=DX/DR
	    TXYZ(1,2)=DY/DR
	    TXYZ(1,3)=DZ/DR
!  
          DX=XYZ(1,NCON(IELEM,3))-XYZ(1,NCON(IELEM,4))
          DY=XYZ(2,NCON(IELEM,3))-XYZ(2,NCON(IELEM,4))
          DZ=XYZ(3,NCON(IELEM,3))-XYZ(3,NCON(IELEM,4))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(2,1)=DX/DR
	    TXYZ(2,2)=DY/DR
	    TXYZ(2,3)=DZ/DR        
        ELSE IF (JNEW .EQ. 6)  THEN
          DX=XYZ(1,NCON(IELEM,1))-XYZ(1,NCON(IELEM,5))
          DY=XYZ(2,NCON(IELEM,1))-XYZ(2,NCON(IELEM,5))
          DZ=XYZ(3,NCON(IELEM,1))-XYZ(3,NCON(IELEM,5))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(1,1)=DX/DR
	    TXYZ(1,2)=DY/DR
	    TXYZ(1,3)=DZ/DR
!  
          DX=XYZ(1,NCON(IELEM,2))-XYZ(1,NCON(IELEM,5))
          DY=XYZ(2,NCON(IELEM,2))-XYZ(2,NCON(IELEM,5))
          DZ=XYZ(3,NCON(IELEM,2))-XYZ(3,NCON(IELEM,5))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(2,1)=DX/DR
	    TXYZ(2,2)=DY/DR
	    TXYZ(2,3)=DZ/DR        
        ELSE IF (JNEW .EQ. 7)  THEN
          DX=XYZ(1,NCON(IELEM,3))-XYZ(1,NCON(IELEM,5))
          DY=XYZ(2,NCON(IELEM,3))-XYZ(2,NCON(IELEM,5))
          DZ=XYZ(3,NCON(IELEM,3))-XYZ(3,NCON(IELEM,5))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(1,1)=DX/DR
	    TXYZ(1,2)=DY/DR
	    TXYZ(1,3)=DZ/DR
!  
          DX=XYZ(1,NCON(IELEM,1))-XYZ(1,NCON(IELEM,5))
          DY=XYZ(2,NCON(IELEM,1))-XYZ(2,NCON(IELEM,5))
          DZ=XYZ(3,NCON(IELEM,1))-XYZ(3,NCON(IELEM,5))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(2,1)=DX/DR
	    TXYZ(2,2)=DY/DR
	    TXYZ(2,3)=DZ/DR        
        ELSE IF (JNEW .EQ. 8)  THEN
          DX=XYZ(1,NCON(IELEM,2))-XYZ(1,NCON(IELEM,6))
          DY=XYZ(2,NCON(IELEM,2))-XYZ(2,NCON(IELEM,6))
          DZ=XYZ(3,NCON(IELEM,2))-XYZ(3,NCON(IELEM,6))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(1,1)=DX/DR
	    TXYZ(1,2)=DY/DR
	    TXYZ(1,3)=DZ/DR
!  
          DX=XYZ(1,NCON(IELEM,3))-XYZ(1,NCON(IELEM,6))
          DY=XYZ(2,NCON(IELEM,3))-XYZ(2,NCON(IELEM,6))
          DZ=XYZ(3,NCON(IELEM,3))-XYZ(3,NCON(IELEM,6))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(2,1)=DX/DR
	    TXYZ(2,2)=DY/DR
	    TXYZ(2,3)=DZ/DR        
        ELSE IF (JNEW .EQ. 9)  THEN
          DX=XYZ(1,NCON(IELEM,1))-XYZ(1,NCON(IELEM,6))
          DY=XYZ(2,NCON(IELEM,1))-XYZ(2,NCON(IELEM,6))
          DZ=XYZ(3,NCON(IELEM,1))-XYZ(3,NCON(IELEM,6))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(1,1)=DX/DR
	    TXYZ(1,2)=DY/DR
	    TXYZ(1,3)=DZ/DR
!  
          DX=XYZ(1,NCON(IELEM,2))-XYZ(1,NCON(IELEM,6))
          DY=XYZ(2,NCON(IELEM,2))-XYZ(2,NCON(IELEM,6))
          DZ=XYZ(3,NCON(IELEM,2))-XYZ(3,NCON(IELEM,6))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(2,1)=DX/DR
	    TXYZ(2,2)=DY/DR
	    TXYZ(2,3)=DZ/DR              
!       
        ENDIF
        
!
!	  print *, TXYZN(IEL,LI,1),TXYZN(IEL,LI,2),TXYZN(IEL,LI,3)

       RETURN
       END           

!
!  *************************************************
!
!    对于域内单元，按某一旋转方向排列各单元
!
!  *************************************************
!  
        SUBROUTINE  DOMAINFAS(DXYZN,TXYZN,MCELE,NCONT,NSIDE,MIN)
        IMPLICIT   NONE 
	  
	  INTEGER  I,J,K,N,NA,NAXIS
	  INTEGER  MCELE,NCONT(100),NSIDE(100,2)

	  INTEGER  LCT(100),LSD(100),MIN

	  REAL*8   DXYZN(100,3),TXYZN(100,2,3)
	  REAL*8   VEC(3),DN(100),DMIN,DOT
!
! MCELE:       包含该点的子单元数量
! NCONT(N)  ： 按某一旋转方向排列的单元次序
! NSIDE(N,2)： 按某一旋转方向排列的单元中边的次序
!例：  NSIDE(5,1)=2,NSIDE(5,2)=1  第5个单元中第1个边为2号边，第2个边为1号边
!	
!
!       WRITE(10, *) ' Inside DOMAINFAS   MCEL=',MCELE
       
! 从第一个单元开始
!
	  I=1  
	  NCONT(I)=1

	   VEC(1)=TXYZN(I,1,2)*TXYZN(I,2,3)-TXYZN(I,1,3)*TXYZN(I,2,2)
	   VEC(2)=TXYZN(I,1,3)*TXYZN(I,2,1)-TXYZN(I,1,1)*TXYZN(I,2,3)
	   VEC(3)=TXYZN(I,1,1)*TXYZN(I,2,2)-TXYZN(I,1,2)*TXYZN(I,2,1)

	   DOT=VEC(1)*DXYZN(1,1)+VEC(2)*DXYZN(1,2)+VEC(3)*DXYZN(1,3)

	  IF(DOT .LT. 0.0D0) THEN
	    NSIDE(I,1)=1
	    NSIDE(I,2)=2
	  ELSE
	    NSIDE(I,1)=2
	    NSIDE(I,2)=1
	  ENDIF
!

	  DO 100 I=2, MCELE
!	  PRINT *,' I=',I
!
!   查找哪个单元中的边与刚查过单元的第二条边的方向最为接近
!	   
           NA=0
	   DO  40 J=2, MCELE
!	     PRINT *,'     J=',J
!
!  检查该单元是否已挑出
!	  
	  DO N=1, I-1
	    IF(NCONT(N).EQ.J)  GOTO  40
	  END DO

!	     PRINT *,'    J=',J,' NCONT(I-1)=',NCONT(I-1)
!	     Print *,' NSIDE(I-1,2)=',NSIDE(I-1,2)
	    DO 30 K=1, 2

	     NA=NA+1
!
	     DN(NA)=SQRT((TXYZN(J,K,1)-TXYZN(NCONT(I-1),NSIDE(I-1,2),1))**2+ 
     1         (TXYZN(J,K,2)-TXYZN(NCONT(I-1),NSIDE(I-1,2),2))**2+  
     2	         (TXYZN(J,K,3)-TXYZN(NCONT(I-1),NSIDE(I-1,2),3))**2 )	  
!
	     LCT(NA)=J
	     LSD(NA)=K
30	 CONTINUE

40	 CONTINUE


!	  WRITE(10, *) ' NA=',NA
	  MIN=0
	  DMIN=1.0d0
	  DO 80 N=1, NA
!	  WRITE(10, *) ' DN(N)=',DN(N)
	     
	     IF(DN(N) .LT. DMIN) THEN
		    DMIN=DN(N)
	      MIN=N
	     ENDIF
80	 CONTINUE

        IF (MIN .LE. 0) THEN
        
!          WRITE(6, *) '  There is a problem with the mesh.'
!          WRITE(6, *) '  Please consult Professor Bin Teng for it.'
!         STOP        
 
        ELSE IF (MIN .GT. 0) THEN   ! ???
!         WRITE(10,*)  ' MIN=',MIN
!         WRITE(10,*)  ' LCT(MIN)=',LCT(MIN)
         
	   NCONT(I)=LCT(MIN)

	   NSIDE(I,1)=LSD(MIN) 
	   IF(LSD(MIN) .EQ. 1) THEN
	    NSIDE(I,2)=2
	   ELSE IF(LSD(MIN) .EQ. 2) THEN
	    NSIDE(I,2)=1
	   ENDIF
!
        ENDIF
!
100	 CONTINUE
!
!	  Print *
!	  PAUSE
!
	  END



!
!  *************************************************
!
!     域内单元构成的固角
!
!  *************************************************
!  

        SUBROUTINE  ANGLE0(DXYZN,TXYZN,MCELE,NCONT,NSIDE,ANGLE)
	  IMPLICIT NONE
!
	  INTEGER  I,N,M
	  INTEGER  MCELE,MCNT(100),MNOD(100)
	  INTEGER  NCONT(100),NSIDE(100,2)
	  REAL*8   ANGLE,DXYZN(100,3),TXYZN(100,2,3)
	  REAL*8   ANG(100),VEC(3),DOT1,DOT2,SIGN,PI

        ANGLE=0.0
	 PI=4.0d0*DATAN(1.0d0)

	 DO I=1, MCELE-1
	   N=NCONT(I)
	   M=NCONT(I+1)
	   DOT1=DXYZN(N,1)*DXYZN(M,1)+DXYZN(N,2)*DXYZN(M,2)+ 
	1	    DXYZN(N,3)*DXYZN(M,3)
	   
	   VEC(1)=DXYZN(N,2)*DXYZN(M,3)-DXYZN(N,3)*DXYZN(M,2)
	   VEC(2)=DXYZN(N,3)*DXYZN(M,1)-DXYZN(N,1)*DXYZN(M,3)
	   VEC(3)=DXYZN(N,1)*DXYZN(M,2)-DXYZN(N,2)*DXYZN(M,1)
	   DOT2=VEC(1)*TXYZN(N,NSIDE(I,2),1)+VEC(2)*TXYZN(N,NSIDE(I,2),2)+ 
	1	    VEC(3)*TXYZN(N,NSIDE(I,2),3)

    
	   IF(ABS(DOT2) .LT. 1.0D-10) THEN
	    SIGN=1.0d0
	   ELSE
	    SIGN=DOT2/ABS(DOT2)
	   ENDIF
!
	 IF(DOT1  .GT. 1.0d0)   DOT1= 1.0d0
	 IF(DOT1  .LT.-1.0d0)   DOT1=-1.0d0

	   ANG(I)=SIGN*ACOS(DOT1)

	 ENDDO
	  
	   I=MCELE
	   N=NCONT(I)
	   M=NCONT(1)
	   DOT1=DXYZN(N,1)*DXYZN(M,1)+DXYZN(N,2)*DXYZN(M,2)+  
	1	    DXYZN(N,3)*DXYZN(M,3)
	   VEC(1)=DXYZN(N,2)*DXYZN(M,3)-DXYZN(N,3)*DXYZN(M,2)
	   VEC(2)=DXYZN(N,3)*DXYZN(M,1)-DXYZN(N,1)*DXYZN(M,3)
	   VEC(3)=DXYZN(N,1)*DXYZN(M,2)-DXYZN(N,2)*DXYZN(M,1)   
	   DOT2=VEC(1)*TXYZN(N,NSIDE(I,2),1)+VEC(2)*TXYZN(N,NSIDE(I,2),2)+  
	1	    VEC(3)*TXYZN(N,NSIDE(I,2),3)


	   IF(ABS(DOT2) .LT. 1.0D-10) THEN
	    SIGN=1.0d0
	   ELSE
	    SIGN=DOT2/ABS(DOT2)
	   ENDIF
!
	 IF(DOT1  .GT. 1.0d0)   DOT1= 1.0d0
	 IF(DOT1  .LT.-1.0d0)   DOT1=-1.0d0
!
	   ANG(I)=SIGN*ACOS(DOT1)
!
	   ANGLE=2.0d0*PI

	  DO I=1, MCELE
	  ANGLE=ANGLE-ANG(I)
	  ENDDO
!

	 END


!
!  *************************************************
!
!    补充关于x-y平面对称的单元
!
!  *************************************************
!  
        SUBROUTINE  XYSYM(DXYZN,TXYZN,MCELE)
        IMPLICIT   NONE 
	  
	  INTEGER  I,J
	  INTEGER  MCELE

	  REAL*8   DXYZN(100,3),TXYZN(100,2,3)
!
	  DO I=1, MCELE
	  J=I+MCELE
	   DXYZN(J,1)= DXYZN(I,1)
	   DXYZN(J,2)= DXYZN(I,2)
	   DXYZN(J,3)=-DXYZN(I,3)

	   TXYZN(J,1,1)= TXYZN(I,1,1)
	   TXYZN(J,1,2)= TXYZN(I,1,2)
	   TXYZN(J,1,3)=-TXYZN(I,1,3)

	   TXYZN(J,2,1)= TXYZN(I,2,1)
	   TXYZN(J,2,2)= TXYZN(I,2,2)
	   TXYZN(J,2,3)=-TXYZN(I,2,3)
	  ENDDO
	  MCELE=MCELE+MCELE

	  END

!
!  *************************************************
!
!    补充关于y-z平面对称的单元
!
!  *************************************************
!  
        SUBROUTINE  YZSYM(DXYZN,TXYZN,MCELE)
        IMPLICIT   NONE 
	  
	  INTEGER  I,J
	  INTEGER  MCELE

	  REAL*8   DXYZN(100,3),TXYZN(100,2,3)
!
	  DO I=1, MCELE
	  J=I+MCELE
	   DXYZN(J,1)=-DXYZN(I,1)
	   DXYZN(J,2)= DXYZN(I,2)
	   DXYZN(J,3)= DXYZN(I,3)

	   TXYZN(J,1,1)=-TXYZN(I,1,1)
	   TXYZN(J,1,2)= TXYZN(I,1,2)
	   TXYZN(J,1,3)= TXYZN(I,1,3)

	   TXYZN(J,2,1)=-TXYZN(I,2,1)
	   TXYZN(J,2,2)= TXYZN(I,2,2)
	   TXYZN(J,2,3)= TXYZN(I,2,3)
	  ENDDO
!
	  MCELE=MCELE+MCELE
!
	  END
!
!  *************************************************
!
!   补充关于z-x平面对称的单元
!
!  *************************************************
!  
        SUBROUTINE  ZXSYM(DXYZN,TXYZN,MCELE)
        IMPLICIT   NONE 
	  
	  INTEGER  I,J
	  INTEGER  MCELE

	  REAL*8   DXYZN(100,3),TXYZN(100,2,3)
!
	  DO I=1, MCELE
	  J=I+MCELE
	   DXYZN(J,1)= DXYZN(I,1)
	   DXYZN(J,2)=-DXYZN(I,2)
	   DXYZN(J,3)= DXYZN(I,3)

	   TXYZN(J,1,1)= TXYZN(I,1,1)
	   TXYZN(J,1,2)=-TXYZN(I,1,2)
	   TXYZN(J,1,3)= TXYZN(I,1,3)

	   TXYZN(J,2,1)= TXYZN(I,2,1)
	   TXYZN(J,2,2)=-TXYZN(I,2,2)
	   TXYZN(J,2,3)= TXYZN(I,2,3)
	  ENDDO
!
	  MCELE=MCELE+MCELE

	  END

!  *************************************************
!
!   补充关于 y-z, z-x 两平面对称的单元
!
!  *************************************************
!  
        SUBROUTINE  YZZXSYM(DXYZN,TXYZN,MCELE)
        IMPLICIT   NONE 
	  
	  INTEGER  I,J
	  INTEGER  MCELE

	  REAL*8   DXYZN(100,3),TXYZN(100,2,3)
!
	  DO I=1, MCELE
	  J=I+MCELE
	   DXYZN(J,1)= DXYZN(I,1)
	   DXYZN(J,2)=-DXYZN(I,2)
	   DXYZN(J,3)= DXYZN(I,3)

	   TXYZN(J,1,1)= TXYZN(I,1,1)
	   TXYZN(J,1,2)=-TXYZN(I,1,2)
	   TXYZN(J,1,3)= TXYZN(I,1,3)

	   TXYZN(J,2,1)= TXYZN(I,2,1)
	   TXYZN(J,2,2)=-TXYZN(I,2,2)
	   TXYZN(J,2,3)= TXYZN(I,2,3)
	  ENDDO
!
	  DO I=1, MCELE
	  J=I+2*MCELE
	   DXYZN(J,1)=-DXYZN(I,1)
	   DXYZN(J,2)=-DXYZN(I,2)
	   DXYZN(J,3)= DXYZN(I,3)

	   TXYZN(J,1,1)=-TXYZN(I,1,1)
	   TXYZN(J,1,2)=-TXYZN(I,1,2)
	   TXYZN(J,1,3)= TXYZN(I,1,3)

	   TXYZN(J,2,1)=-TXYZN(I,2,1)
	   TXYZN(J,2,2)=-TXYZN(I,2,2)
	   TXYZN(J,2,3)= TXYZN(I,2,3)
	  ENDDO
!
	  DO I=1, MCELE
	  J=I+3*MCELE
	   DXYZN(J,1)=-DXYZN(I,1)
	   DXYZN(J,2)= DXYZN(I,2)
	   DXYZN(J,3)= DXYZN(I,3)

	   TXYZN(J,1,1)=-TXYZN(I,1,1)
	   TXYZN(J,1,2)= TXYZN(I,1,2)
	   TXYZN(J,1,3)= TXYZN(I,1,3)

	   TXYZN(J,2,1)=-TXYZN(I,2,1)
	   TXYZN(J,2,2)= TXYZN(I,2,2)
	   TXYZN(J,2,3)= TXYZN(I,2,3)
	  ENDDO

	  MCELE=4*MCELE

	  END
!