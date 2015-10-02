
C  
C **************************************************************** 
C *                                                              * 
C *    Test program of a novel BEM method basded on hybrid       * 
C *   integral equations for internal problem.                   * 
C *							         *
C *  							         *
C * 	 December 1, 2014                                        *
C *                                                              * 
C **************************************************************** 
C 
	PROGRAM Hyper_Integral_Equation

	 USE MVAR_MOD
	  USE PVAR_MOD
	  use INTVar_mod
	  USE WMAT_MOD
      
          !use hi_SIEPPEM
          use hi_intg

	  IMPLICIT  NONE  
	  INTEGER I,IFWKO,NDRFA,IS,IND,M,IP,IPOL,K,NDMAX
	  INTEGER NTnum,IFLAG_T

	  real(8) ::  WL,AMFJ,R,EX(4),Alpha,WAVES  
	  real(8) ::  FAMPR(6),FAMPI(6),FORCER(6),FORCEI(6)
	  real(8) ::  PL_AMP(6),FORAMP
	  real(8) ::  FCD_AMR,  FCD_AMI
C !
C ! IFLAG_T,FCD_AMR, FCD_AMI: not used in this program
C !
C         DATA EX /  1.,  1., -1., -1./ 
C !
C ! FAMPR, FAMPI:   real and imaginary parts of body motion 
C !                 from frequency domain calculation
C ! FORCER, FORCEI: real and imaginary parts of wave force 
C !                 from frequency domain calculation
C ! PL_AMP:         amplitude of plotting variable
C
C ----------------------------------------
C Input data files

        OPEN(1, FILE='INPUT/DATIN.txt',      STATUS='OLD') 
        OPEN(2, FILE='INPUT/DATBDMS.txt',    STATUS='OLD') 
        OPEN(3, FILE='INPUT/DATWPMS.txt',    STATUS='OLD') 

C         OPEN(8, FILE='INPUT/SOLIDANGLE.txt',    STATUS='OLD') 
C         open(67, FILE='INPUT/M.txt',    STATUS='OLD')
C         read(67,*) wl
C         stop


C ! -----------------------------------
C !  Output data files
C !
        OPEN(6,  FILE='OUTPUT/OUTScreen.txt',    STATUS='UNKNOWN')
        OPEN(9,  FILE='OUTPUT/OUTPUT1.txt',    STATUS='UNKNOWN')
        OPEN(10, FILE='OUTPUT/OUTPUT.txt' ,    STATUS='UNKNOWN')
        OPEN(11, FILE='OUTPUT/mesh_track.txt' ,    STATUS='UNKNOWN')
        OPEN(12, FILE='OUTPUT/output_amatrix.txt' ,    STATUS='UNKNOWN')


        OPEN(101,FILE='OUTPUT/OUTAMT.txt' ,    STATUS='UNKNOWN')	  
        OPEN(102,FILE='OUTPUT/OUTBMT.txt' ,    STATUS='UNKNOWN')
        OPEN(103,FILE='OUTPUT/OUTCMT.txt' ,    STATUS='UNKNOWN')
        OPEN(109,FILE='OUTPUT/OUTPUT9.txt' ,    STATUS='UNKNOWN')
        OPEN(108,FILE='OUTPUT/OUTPUT7.txt' ,    STATUS='UNKNOWN')
        OPEN(110,FILE='OUTPUT/OUTPUT10.txt' ,    STATUS='UNKNOWN')    

C            OPEN(108,FILE='OUTPUT/sing1_debug_info.txt',STATUS='UNKNOWN')


C
C ----------------------------------------------- 
C
C
	MFREE = 1
	NBETA=0
	 
        WRITE(10, *) ' Test on huper-singular integration'

C
C ---------------------------------------- 
C
C  IFWKO=0 : input wave number, otherwise input wave frequency
C  H<0: infinite water depth
C  Nwave: Number of waves to simulate
C
        READ(1,*)      IFWKO

        IF (IFWKO .EQ. 0)  THEN
          READ(1,*)      H, AMP, WK, BETA
          IF (H .LE. 0.0D0) THEN
            W1=DSQRT(G*WK)
          ELSE
            W1=DSQRT(G*WK*DTANH(WK*H))
          END IF
        ELSE
          READ(1,*)      H, AMP, W1, BETA
          IF(H .LE. 0.0D0) THEN
            WK=W1**2/G
          ELSE
            CALL WAVECK(W1,H,WK)				 ! Compute wave number
          END IF
        END IF
C !
C !
C !  Waves: number of waves to simulate (may be 10 or 0.5)
C !  NTnum: number of time steps in one wave
C !  NPLOUT=1-5, Output the wave profile with step intervals 
C !  NDRFA=0, plot force; =1, plot response
C !
	  PI4=4.0*PI
        IORDER=1
   
        WRITE(10,*)  '  After 10'
!  
!  PI4，IORDER 的定义在mvar中 perturbation order of the problem 
! 
        TPER=2.*PI/W1 
        BETA=BETA*PI/180.0D0 
C C  beta化成弧度
        V=W1*W1/G 
	  WL=2.0D0*PI/WK
C 

        WRITE(10,*) 
        WRITE(9,*)
	  WRITE(10,*) '                   ================='
        WRITE(9,*) '                   ================='
C
        WRITE(10,1111)  H,AMP,WK,V,WL,W1,TPER,BETA*180./PI 
        WRITE(9,1111)  H,AMP,WK,V,WL,W1,TPER,BETA*180./PI 

C 1111    FORMAT(//,'  WATER DEPTH=',F9.3,'    WAVE AMPLITUDE=', F6.2,/,
C      1    '  WAVE NUMBER=',F9.5,'  K0=',F9.5,'  WAVE LENGTH=',F9.4,/, 
C      3    '  ANGULAR FREQU.=',F9.5,'   WAVE PERIOD=',F7.3,/,      
C      2    '  WAVE DIRECTION:',F7.3,'  Degree',/)
C !  ------------------------------------
C !bodmass在mass.f中

        XC=0
        YC=0
        ZC=0

C !    mvar中ISYS: number of symmetric planes
C !    NELEMB: number of elements on body surface
C !    NNB: number of nodes on the body surface according to coordinate
C !    NNBD: number of nodes on the body surface according to directives 
        READ(2,*) ISYS 
        READ(2,*)   NELEMB, NNB, NNBD, IPOL

   	 ALLOCATE (NCONB(NELEMB,8),NCONDB(NELEMB,8))




       ALLOCATE (XYZB(3,NNB),DXYZB(3,NNBD))

!	

        IF(ISYS.EQ.0) NSYS=1
        IF(ISYS.EQ.1) NSYS=2
        IF(ISYS.EQ.2) NSYS=4
C
C !    mvar中NELEMF:number of elements on the free surface
      

        READ(3,*)   NELEMF
       
!
	  WRITE(11,*) ' ISYS=',ISYS,' NSYS=',NSYS
	  WRITE(11,*) ' NELEMB=',NELEMB,' NELEMF=',NELEMF
C !    mvar中NELEM: number of total elements

	  NELEM=NELEMB+NELEMF

	  WRITE(11,*) ' NELEM=',NELEM,'  IOPL=',IPOL

        ALLOCATE(SAMB(NELEM,16,0:8),SAMBXY(NELEM,16,3),
	1		   DSAMB(NELEM,16,6),NCN(NELEM),NCON(NELEM,8),
     1		   NCOND(NELEM,8),IETYPE(NELEM),NNORMN(8*NELEM) )

        ALLOCATE( XYZE(3,8,NELEM),DXYZE(3,8,NELEM),DAMPE(8,NELEM))
        ALLOCATE( TXYZE(3,8,NELEM))
	  ALLOCATE( XYZTP(3,8*NELEM),DXYZTP(3,8*NELEM),DAMPTP(8*NELEM))

C ! mvar中的定义 NCN: number of nodes in the element
C ! IETYPE: type of the element; =1, on body surface; =2, on free surface
C ! SAMBXY: Coordinates of Gaussin points
C ! DSAMB:  Normal direvatives at Gaussian points
C ! XYZE  : Initial Coordinates of nodes of body mesh
C ! TXYZE : Coordinates of nodes of body mesh at the simulation time 
C ! --------------------------------------------
C ! MESHFS4和MESHBD在meshda4.f中

        call MESHFS4   		        ! Read in data on free surface mesh
        
        WRITE(11,*),'  After MESHFS4' 

        CALL MESHBD(IPOL) 		    ! Read in data on body mesh
        WRITE(11,*),'  After MESHBD' 

	   CLOSE(2)

     
        
C         OPEN(50, FILE='OUTPUT/DATBDMS.txt',    STATUS='UNKNOWN')

        WRITE(11,*),'  Before CONVSB' 
        WRITE(11,*) ' NNODE=',NNODE
	  
	      CALL CONVSB
        WRITE(11,*),'  After CONVSB' 
        
C ! mvr中 NNODE:  total number of nodes according to the coordinate
C ! NNODED: total number of nodes according to the normals
	  WRITE(11,*) ' NNODE=',NNODE

	  NDMAX=MAX(NNODE,NNODED)
       ALLOCATE(NODELE(NNODE,64),NODNOE(NNODE),NODELJ(NNODE,64),
     1         NODQUA(NNODE),NNORMC(NNODED),
     2		 XYZ(3,NDMAX),DXYZ(6,NDMAX),DAMPF(NNF))
        ALLOCATE(ANGLE(NNODE),FrA3(NNODE),
     1	        FrC31(NNODE),FrC32(NNODE),FrC33(NNODE))
!
! ---------------------------------------------------
!
	 WRITE(10,*) '  IND       X          Y          Z        DAMPing'
	 DO IND=1, NNF
	  DAMPF(IND)=DAMPTP(IND)
	  XYZ(1,IND)=XYZTP(1,IND)
	  XYZ(2,IND)=XYZTP(2,IND)
	  XYZ(3,IND)=XYZTP(3,IND)
	  !WRITE(10,111) IND, XYZ(1,IND),XYZ(2,IND),XYZ(3,IND),DAMPF(IND)
	 END DO
!	  	 
	 DO IND=NNF+1, NNODE
	  XYZ(1,IND)=XYZTP(1,IND)
	  XYZ(2,IND)=XYZTP(2,IND)
	  XYZ(3,IND)=XYZTP(3,IND)
	  !WRITE(10,111) IND, XYZ(1,IND),XYZ(2,IND),XYZ(3,IND)
	 END DO
!
	 DO IND=1, NNODED
	  DXYZ(1,IND)=DXYZTP(1,IND)
	  DXYZ(2,IND)=DXYZTP(2,IND)
	  DXYZ(3,IND)=DXYZTP(3,IND)
	  NNORMC(IND)=NNORMN(IND)
	 ! WRITE(10,111) IND, DXYZ(1,IND),DXYZ(2,IND),DXYZ(3,IND)
	 END DO
!
111	 FORMAT(I4,5(2x,F12.5))
!
	 DEALLOCATE(XYZTP,DXYZTP,NNORMN)
!
	 DAMPF(:)=W1*DAMPF(:)
!
!  --------------------------------------------

       ALLOCATE(AMATA(NNODE,NNODE,NSYS),
	1		  BMATA(NNODE,NSYS), INDX(NNODE,NSYS))
!       
       ALLOCATE(UNKN(NNODE,NSYS),  BKN(NNODED,NSYS),
	1		  UNKN_O(NNODE,NSYS),BKN_O(NNODED,NSYS),
	1	      DPDT(NNODE,NSYS))
	 ALLOCATE(DH(4,NNF,NSYS),DP(4,NNF,NSYS),Dposi(4,6))
!
!      
! Identify the Gaussian sampling points and evaluate the    
! Corresponding values for the shape functions and Jacobian matrices       
!
        CALL BODYFD                  
        WRITE(10,*)  'AFTER BODYFD' 
!
!  --------------------------------------------
!   
        print *,"read model"
        call read_model_from_WAVDUT()
C ! Assembling matrix and computing diffraction and radiation potentials
C !
         CALL TASSB0   
C 	  WRITE(10,*) '    AFTER TASSB0' 
! 
! =================
! 
C         do i = 1,NELEMF
C           write (200,999) i,ncon(i,1:8)
C           x = 0.5*(XYZ(1,NCON(2))+XYZ(1,NCON(6)))
C           Y = 0.5*(XYZ(2,NCON(2))+XYZ(2,NCON(6)))
C         end DO



!                              
1010    FORMAT(F7.3,1x,F7.3,1x,6E14.5) 
! 
1111    FORMAT(//,'  WATER DEPTH=',F9.3,'    WAVE AMPLITUDE=', F6.2,/,
     1    '  WAVE NUMBER=',F9.5,'  K0=',F9.5,'  WAVE LENGTH=',F9.4,/, 
     3    '  ANGULAR FREQU.=',F9.5,'   WAVE PERIOD=',F7.3,/,      
     2    '  WAVE DIRECTION:',F7.3,'  Degree',/)
! 
1115	  FORMAT(' I_time=',I5,'    at the time:',F10.5,' s')
1200	FORMAT(2x,I3,3F12.5,1x,2F13.5)
        print *,"=================== main program ends ==============="
        STOP 
        END      

 

C 
C
C
C  *************************************************
C  *    The subroutine computes wave number        *
C  *  by wave frequency and water depth.           *
C  *    For infinite water depth, give negative H. *
C  *************************************************

        SUBROUTINE WAVECK(SIGMA,H,WK)
        IMPLICIT  NONE

	  real(8) SIGMA,H,WK,B,G,A,Y,C
C
C  H: WATER DEPTH    SIGMA: WAVE FREQUENCY
C  C: WAVE CELERITY  Wk: wave number
C
        DATA G/9.807D0/
C
        IF( SIGMA .LE. 0.0D0 )  THEN
	    Print *,' IN WAVECK:  W=',SIGMA
		STOP
	
        ELSE 

          IF( H .GT. 0.0D0) THEN
           B = G * H
           Y = SIGMA * SIGMA *H/G
  12       A=1.0D0/(1.0D0+Y*(0.66667D0+Y*(0.35550D0+
     1       Y*(0.16084D0+Y*(0.063201D0+Y*(0.02174D0+
     2       Y*(0.00654D0+Y*(0.00171D0+
     2       Y*(0.00039D0+Y*0.00011D0)))))))))
  22       C=SQRT(B/(Y+A))
           WK=SIGMA/C
           ELSE IF( H .LE. 0.0D0) THEN
           WK=SIGMA**2/G
        END IF
        RETURN
C

	 END IF


        RETURN
        END

