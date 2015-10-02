C **************************************************
C 												
C   This is a module for declaring variables  
C                                                 
C       Nov. 3, 2000           by   Bin Teng		
C **************************************************
C
C  MODULES: MVAR_MOD + PVAR_MOD + TRVAR_MOD
C
C  MVAR_MOD : for constants, wave parameter and mesh
C  BVAR_MOD : for B-spline function
C  PVAR_MOD : for potentials, forces, and body mass, etc
C  TRVAR_MOD: for finer mesh from tri-pole transformation
C
C ===================================================
!C
!C ===================================================
!C  Varibles used in the pFFT method 
!C ===================================================
!C
      MODULE SEBSM_MOD
	REAL*8,ALLOCATABLE:: A_SEBSM(:,:),B_SEBSM(:)
      INTEGER N_SEBSM
      END MODULE SEBSM_MOD

      MODULE PFFT_mod
!C
	INTEGER Ncube,Nps,Npr,NNps
	INTEGER NCNcube_TOTAL,NCNcubeE_TOTAL
	INTEGER NcubeX,NcubeY,NcubeZ
	INTEGER NprX,NprY,NprZ
	INTEGER Npt

	REAL*8 Xmax,Xmin,Ymax,Ymin,Zmax,Zmin,Dlength,blap,Dsp

	INTEGER,ALLOCATABLE:: NCNcube(:,:,:),NCONcube(:,:,:,:)
	INTEGER,ALLOCATABLE:: NcubeNODX(:),NcubeNODY(:),NcubeNODZ(:) 
	 
	INTEGER,ALLOCATABLE:: NCNcubeE(:,:,:),NCONcubeE(:,:,:,:)
	INTEGER,ALLOCATABLE:: NcubeELEX(:),NcubeELEY(:),NcubeELEZ(:) 

	INTEGER,ALLOCATABLE:: IprX(:,:),IprY(:,:),IprZ(:,:)
	
	REAL*8,ALLOCATABLE:: Xcube_pFFT(:,:,:)
	REAL*8,ALLOCATABLE:: Ycube_pFFT(:,:,:)
	REAL*8,ALLOCATABLE:: Zcube_pFFT(:,:,:)
	REAL*8,ALLOCATABLE:: Xpr(:,:,:),Ypr(:,:,:),Zpr(:,:,:)
	REAL*8,ALLOCATABLE:: Xtest(:),Ytest(:),Ztest(:),RGTINV(:,:)
	

	 COMPLEX*16,ALLOCATABLE:: UUU(:,:)
	 COMPLEX*16,ALLOCATABLE:: HTPLZ(:,:,:),HHKL(:,:,:)


	 COMPLEX*16,ALLOCATABLE::  HHKLV(:,:,:,:)	
	 COMPLEX*16,ALLOCATABLE::  QVT(:,:,:,:),QVH(:,:,:,:)
	 COMPLEX*16,ALLOCATABLE::  QVV(:,:,:)

	 COMPLEX*16,ALLOCATABLE::  WPROL(:,:,:,:),WPROR(:,:,:,:)
	 COMPLEX*16,ALLOCATABLE::  HTPLZL(:,:,:,:),HTPLZR(:,:,:,:)


!	REAL*8,ALLOCATABLE:: WPQL(:,:,:),WPQR(:,:,:,:)

      END MODULE PFFT_mod
!
! ===================================================
!  Varibles used in the SUBROUTINEs PRECORRECTED_R and PRECORRECTED_L
! ===================================================
!
       MODULE  PRECOR_MOD

		 COMPLEX*16,ALLOCATABLE:: WPRO(:,:,:,:)
		 COMPLEX*16,ALLOCATABLE:: TPLZ(:,:,:,:)

       END MODULE PRECOR_MOD

       MODULE INTVar_mod

       INTEGER MFREE,NBETA   
	 END MODULE INTVar_mod  


!
!========================================================
C   Constants and Variables for body mesh
C
        MODULE MVar_mod
C
        INTEGER NTIME,ITIME,ISYS,NSYS,IORDER
	  INTEGER NELEM,NELEMB,NELEMF,NNODE,NNODED,NNB,NNBD,NNF,NNTCH

	  INTEGER NPLOUT
!
!  NPLOUT=1, Output the wave profile at every time step; =2, every 2 time steps;
!        =3, every 4 time step;  =4, every 8 time steps; =5, every 16 time steps.
!
! NTIME:  total time steps for simulation
! ITIME:  present time step 
! ISYS:   number of symmetric planes
! IORDER: perturbation order of the problem 
! NELEM:  number of total elements 
! NELEMB: number of elements on body surface
! NELEMF: number of elements on the free surface
! NNODE:  total number of nodes according to the coordinate
! NNODED: total number of nodes according to the normals
! NNB:    number of nodes on the body surface according to coordinate
! NNBD:   number of nodes on the body surface according to directives
! NNF:    number of nodes on the free surface
! NNTCH:  number of nodes and the centers of element (For plotting by Techplot)
! 
        REAL*8 G,RHO,PI,PI4     
        REAL*8 H,AMP,BETA,W1,V 
        REAL*8 WK,TPER
	  REAL*8 Tstep,TIME,TimeRK,RAMPF,RAMPV
        real*8 line_sum

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
!    
   	 INTEGER, ALLOCATABLE:: NCN(:),IETYPE(:),NCON(:,:),NCOND(:,:),
	1                        NNORMC(:)
       INTEGER, ALLOCATABLE:: NODELE(:,:),NODNOE(:),
	1			            NODELJ(:,:),NODQUA(:)
   	 INTEGER, ALLOCATABLE:: NCONB(:,:),NCONDB(:,:)

! NCN: number of nodes in the element
! IETYPE: type of the element; =1, on body surface; =2, on free surface
! NCON: 
! NCOND:
! NNORMC:
! NODELE:
! NODNOE:
! NODELJ:
! NODQUA:
!
       REAL*8,ALLOCATABLE:: XYZE(:,:,:),DXYZE(:,:,:),XYZ(:,:),DXYZ(:,:),
	1                      DAMPE(:,:),DAMPF(:)
       REAL*8,ALLOCATABLE:: XYZB(:,:),DXYZB(:,:)

       REAL*8,ALLOCATABLE:: TXYZE(:,:,:)
!
! XYZE  : Initial Coordinates of nodes of body mesh
! DXYZE:
! XYZ:
! DXYZ:
! DAMPE:
! DAMPF:
! TXYZE : Coordinates of nodes of body mesh at the simulation time 
!
	 REAL*8,ALLOCATABLE:: BKN(:,:),BKN_O(:,:),UNKN(:,:),UNKN_O(:,:)
	 REAL*8,ALLOCATABLE:: ET(:,:),ET_O(:,:),DPDT(:,:)
!	 REAL*4,ALLOCATABLE:: BKN(:,:),UNKN(:,:),BKN_O(:,:),UNKN_O(:,:,:),
!	1					  ET(:,:),ET_O(:,:),DPDT(:,:)
!
! BKN:  1:NNF, Potential; NNF+1:NNODE, normal derivative
! UNKN: 1:NNF, normal derivative; NNF+1:NNODE, potential 
! BKN_O: store one-step data
! UNKN_O: store two-step data
! ET: Wave profile
! ET_O: store one-step data
! DPDT: time difference of potential

!
       REAL*8,ALLOCATABLE:: SAMB(:,:,:),SAMBXY(:,:,:),DSAMB(:,:,:)
	 REAL*8,ALLOCATABLE:: ANGLE(:),FrA3(:),FrC31(:),FrC32(:),FrC33(:)
!
! SAMB: 
! SAMBXY: Coordinates of Gaussin points
! DSAMB:  Normal direvatives at Gaussian points
! ANGLE:  Solid angle
!
	 REAL*8,ALLOCATABLE:: DH(:,:,:),DP(:,:,:),Dposi(:,:)
!
! DH:
! DP:
! Dposi:
!
!  For linear equations
!
	 INTEGER, ALLOCATABLE:: INDX(:,:)
       REAL*8,ALLOCATABLE::   AMATA(:,:,:),BMATA(:,:)
!       REAL*8,ALLOCATABLE::   LEFT(:,:,:),RIGHT(:,:,:),RIGHT1(:,:)

       DATA G,PI,RHO/9.807,3.14159265359,
	1	               1023.0/  

! ------------------------------------
! Temporary arrayes
! 
   	 INTEGER, ALLOCATABLE:: NNORMN(:)
       REAL*8,  ALLOCATABLE:: XYZTP(:,:),DXYZTP(:,:),DAMPTP(:)
!
	   DATA AL/1.2/
	   DATA MAX_LEVEL/ 8/
       END MODULE MVar_mod


C
C Variable for potential, force and body motion
C =============================================
C
        MODULE PVar_mod
C
	  INTEGER  NNT
	  PARAMETER (NNT=2000)    

        REAL*8 XC,YC,ZC,XTC,YTC,ZTC
        REAL*8 ARE,XF,YF,XK2,YK2,XCF
        REAL*8 VOLM,XB,YB,ZB

        REAL*8 AMAS(6,6),BDMP(6,6)  
        REAL*8 RMAS(6,6),CDMP(6,6),CRS(6,6),STKM(6,6),XIA(3,3)
C
        REAL*8 FORCEW(6),FORCE0(6),FORSCD(6),AMPJ(6)
	  REAL*8 FORCE(6),FORCE_O(6)
C
        REAL*8 TRMAS(6,6),VISC(6,6)

        REAL*8  DISP(6),DSDT(6),DISP_O(6),
	2		  DSDT_O(6),DSDDTL(6),TRXYZ(3,3)   

        REAL*8  RESPR(6),RESPI(6)  


        END MODULE PVar_mod

C
C =========================================
C  Varibles used in the Tri-pole transform 
C							  
       MODULE TRVar_mod

	 INTEGER NOSAMP
	 REAL*8 XYNOD(3,50),DXYNOD(6,50),SAMNOD(50,0:8)

       END MODULE TRVar_mod

C
C  **************************************************
C  *												  *
C  *  This is a module for declaring variables      *
C  *      in B-spline expansion                     *
C  *                                                *
C  *      Aug. 10, 2002          by   Bin Teng      *
C  **************************************************				
C
       MODULE  BVAR_mod
C
       INTEGER KBSPL,NB,NT,NA,NAA,NBU,NA1
       INTEGER MUSTA,MUEND	 
C
C --------------------------------------------------------------------
C KBSPL : DEGREE (order-1) of the B-splines
C
C NB    : Number of given points 
C NT    : Maximum number of given points
C
C NBU   : Number of intervals in the u-direction 
C NA1	  : Number of control factors    (NBU+KBSPL)
C NAA	  : Number of control factors ??   (NBU+2*KBSPL+1)
C NA    : Maximum number of total control factors  
C -------------------------------------------------------------------
C
	 PARAMETER (NT=500, NA=100, KBSPL=3)
	 REAL*8 BXYZ(2,NT),BXYZ1(2,NT),UXYZ(NT)
C
C XYZ     : Catersian coordinates of given points
C UXYZ    : U coordinates of given points on the body surface
C
       REAL*8 BJU(NA),B(NA,4),DBJU(NA)
	 REAL*8 XYZCON(NA,2),UTJ(NA),UXYZNE(NA)
C
C UTJ    : U coordinates of nodes and co-located points
C XYZCON : factors for body geometry expansion
C UXYZNE : U coordinates of controlling points 
C
       END MODULE BVAR_mod
C


!
C  **************************************************
C  *												  *
C  *  variables for fenders and cables              *
C  *                                                *
C  *      April. 20, 2005        by   Bin Teng      *
C  **************************************************				
C
       MODULE FCVAR_mod
C
	  INTEGER NOCABLE,NCKIND,NOFender,NFender(50),JCPRO(50)
	  REAL*4  WTL
	  REAL*4  DLS(50),CLNT(50),CXYZ0(3,50),CXYZ1(3,50)
	  REAL*4  DLNT(50),DXYZ0(3,50),DXYZ1(3,50)
C 
C NCABLE: Number of cables;        
C NCKIND: Number of cable kinds;
C NOFender: Number of Fenders;
C NFender(KD): Number of data for tension curve for each Fender

!  WTL: water level
C    CLNT: length of each cable;     DLNT: length of each Fender
C  CXYZ0, CXYZ1: the coordinates of the two ends of each cable
C  DXYZ0, DXYZ1: the coordinates of the two ends of each Fender

        INTEGER NCABLEU,NFenderU,JDRMXY(50),JDRM(50),KindDrm
	  REAL*4 CABX1(50),CABX2(50),CABXYCON(200,2,50),CABUTJ(200,50)
	  REAL*4 DRMX1(50),DRMX2(50),DRMXYCON(200,2,50),DRMUTJ(200,50)
	  REAL*4 CEPSMAX(50), FCRMAX(50)
	  REAL*4 FDRMAX(6,50),EPSMAX(50)
C
C  JDRM    : The code of Fender character
C JDRMXY(N): the code for Fender direction, 
C          =1, force in x-direction; =2, force in y-direction
C KindDrm  : Number of Fender's kinds

       END MODULE FCVAR_mod

!
! ===================================================
!  Varibles used in the SUBROUTINE CGR
! ===================================================
!
       MODULE  GCR_MOD

		 COMPLEX*16,ALLOCATABLE:: R(:,:,:),P(:,:,:)	
	     COMPLEX*16,ALLOCATABLE:: AR(:,:,:)
	     COMPLEX*16,ALLOCATABLE:: AP(:,:,:),B(:,:)
	     COMPLEX*16,ALLOCATABLE:: X(:,:,:),XS(:,:)
	     COMPLEX*16,ALLOCATABLE:: AAP(:),AIT(:)
	     COMPLEX*16,ALLOCATABLE:: RAP(:),ERRR(:)

       END MODULE GCR_MOD

	 MODULE WMAT_MOD
       real*8,ALLOCATABLE:: WMAT_L(:,:,:),WMAT_R(:,:,:)
	 INTEGER COUNT_W_L,COUNT_W_R
	 END MODULE WMAT_MOD

      MODULE  GMRES_MOD
		   REAL*8,ALLOCATABLE:: work(:)	
      END MODULE GMRES_MOD