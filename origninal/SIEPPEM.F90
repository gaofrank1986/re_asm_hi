     
    
      
      PROGRAM SIEPPEM
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/RIM_COEF/PI,DLT(3,3),CNU
      DIMENSION XP(3)
      ALLOCATABLE CD(:,:),LNDB(:,:),NSEL(:),XIS(:,:),VINT(:)
	  DATA NGR,NGL/-10,-10/  !If |NGR| or |NGL| is bigger than 10, must change subroutine GAUSSV
      OPEN(5,FILE='SIEPPEM.DAT',STATUS='OLD')
      OPEN(7,FILE='SIEPPEM.OUT',STATUS='UNKNOWN')
      READ(5,*)NDIM,NTP,NBE,NODE,BETA,NF        ! Card set 1
      ALLOCATE(CD(NDIM,NTP),LNDB(NODE,NBE),NSEL(NBE),XIS(2,NBE),        &
     &         VINT(NF))
 !    Assign values to COMMON BLOCK variables
      DLT=RESHAPE((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))
      PI=4.D0*DATAN(1.D0); CNU=0.3; TOLGP=1.D-14
      IF(NDIM.EQ.2) NGL=1
 !    Input nodal coordinates and element connectivity
      DO 10 IP=1,NTP
  10  READ(5,*)NP,(CD(I,NP),I=1,NDIM)                 ! Card set 2
      DO 20 IE=1,NBE
  20  READ(5,*)NE,(LNDB(ID,NE),ID=1,NODE),NSEL(NE)    ! Card set 3
      READ(5,*)(XP(I),I=1,NDIM)                       ! Card set 4  
	DO 30 IE=1,NBE
30    IF(NSEL(IE).LT.0) READ(5,*)(XIS(I,IE),I=1,NDIM-1)  ! Card set 5
      !--------------------------------------------------------------------
!      CALL NASTRAN_MESH(32,'BOUNDARY MESH AND INTERNAL POINTS',NDIM,    &
!     &                  NDIM-1,NTP,NBE,CD, LNDB,NODE,1,1)
      !--------------------------------------------------------------------
 !    Evaluate boundary integrals 
      CALL RIM_ELEMS(NDIM,NODE,BETA,NBE,CD,LNDB,NSEL,XP,XIS,NGR,NGL,    &
     &               TOLGP,NF,VINT)
               write(*,*)'              Results :'
      WRITE(*,'(1P,8(3E14.6/,6X))')(VINT(I),I=1,NF)
      WRITE(7,'(1P,8(3E14.6/,6X))')(VINT(I),I=1,NF)
      DEALLOCATE (CD,LNDB,NSEL,XIS,VINT)
      STOP
      END

      SUBROUTINE RIM_ELEMS(NDIM,NODE,BETA,NBE,CD,LNDB,NSEL,XP,XIS,      &
     &                     NGR,NGL,TOLGP,NF,VINT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CD(NDIM,*),LNDB(NODE,NBE),NSEL(NBE),VINT(NF),XP(*),     &
     &          XIS(2,NBE),GPR(IABS(NGR)),GWR(IABS(NGR)),GPL(IABS(NGL)),&
     &          GWL(IABS(NGL)),CDL(18),NODEF(12),CP0(3),CK(3,NODE),     &
     &          V1E(NF)
      EXTERNAL INT_ELEM
      DATA CDL/-1.,-1.,1.,-1.,1.,1.,-1.,1.,0.,-1.,1.,0.,0.,1.,          &
     &         -1.,0.,0.,0./
      DATA NODEF/1,2,5, 2,3,6, 3,4,7, 4,1,8/,NPW/4/
      IF(NDIM.EQ.2) NPOWG=(NODE/3)*2               ! 0,  2
      IF(NDIM.EQ.3) NPOWG=NODE/2+(NODE/9)*2        ! 2,  4,  6 
      NBDM=NDIM-1
      NSUB=2*NBDM
      NDSID=2+NODE/8
      CP0=0.
      VINT=0.
      DO 20 IE=1,NBE
      DO 10 ID=1,NODE
10    CK(1:NDIM,ID)=CD(1:NDIM,LNDB(ID,IE))
      V1E=0.
      IF(NSEL(IE).EQ.0) THEN    ! EVALUATE INTEGRAL OVER REGULAR ELEMENT
       CALL ADAPTINT_ELEM(NDIM,NBDM,NODE,BETA,CP0,XP,CK,CDL,NF,V1E,     &
     &                    NGR,NGL,GPR,GWR,GPL,GWL,TOLGP,BETA,INT_ELEM)
      ELSE     ! EVALUATE INTEGRAL OVER SINGULAR ELEMENT
       CALL SINGULAR_ELEM(NDIM,NBDM,NODE,BETA,NPOWG,XP,CP0,CK,CDL,      &
     &                    TOLGP,IABS(NGL),GPL,GWL,NODEF,NSUB,NDSID,     &
     &                    NSEL(IE),XIS(1,IE),NPW,NF,V1E)
!       print *,NSEL
!       print *,XIS
      ENDIF
      VINT=VINT+V1E
20    CONTINUE
      END
    
      SUBROUTINE ADAPTINT_ELEM(NDIM,NDIMB,NODE,BETA,CP0,CP,CK,CDL,NF,   &
     &                         V1E,NGR,NGL,GPR,GWR,GPL,GWL,TOLGP,RLEVL, &
     &                         INT_COEF)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CP(3),CK(3,NODE),DGS(3),XI(3),XIC(3),AL(3),RI(3),       &
     &          NSUB(3),CSUB(3,8),CDL(*),NGSP(3),SHAP(NODE),V1E(NF),    &
     &          GPR(IABS(NGR)),GWR(IABS(NGR)),GPL(IABS(NGL)),CP0(3),    &
     &          GWL(IABS(NGL))
      DIMENSION COSB(NDIM),DN(2,NODE),GD(3,3)
      EXTERNAL INT_COEF
  !       Determine XIC, DGS and necessary number of sub-elements (cells)
      CALL CHOSEGP(NDIM,NDIMB,NODE,CP,CK,CDL,NGR,NGL,TOLGP,RLEVL,       &
     &             DGS,XIC,AL,DISL,1,NSUB,NGSP,GPR,GWR,GPL,GWL)
  !                      Non sub-division cases
      IF(NSUB(1)+NSUB(2)+NSUB(3).EQ.3) THEN
       DO IGL=1,NGSP(2)
       IF(NDIM.EQ.3)THEN; XI(2)=GPL(IGL); WL=GWL(IGL)
       ELSE; XI(2)=0.; WL=1.; ENDIF
       DO IGR=1,NGSP(1)
        XI(1)=GPR(IGR); WR=GWR(IGR)
        CALL INT_COEF(NDIM,NDIMB,NODE,BETA,CP,CK,XI,WR*WL,1.D0,CDL,NF,  &
     &                V1E)
       ENDDO
       ENDDO
       RETURN
      ENDIF 
  !                 Sub-division cases
      NCOR=2**NDIMB            !  Number of element corner nodes
      DO ID=1,NCOR; CSUB(2,ID)=-1.+(CDL(2*ID)+1.)/NSUB(2); ENDDO
      DO 70 ISUB3=1,NSUB(2)    !  No. of sub-cells in third direction
      DO ID=1,NCOR; CSUB(1,ID)=-1.+(CDL(2*ID-1)+1.)/NSUB(1); ENDDO
      DO 60 ISUB2=1,NSUB(1)    !  No. of sub-cells in second direction
  !         Find local proximate point XI within the sub-element
      CALL MINDISTD(NDIMB,NDIMB,NCOR,XIC,CSUB,1.D0,XI,DISL,COSB,DN,     &
     &              GD,CDL)
  !         Find intrinsic coordinates RI for  XI in the original element
      CALL SHAPEF(NDIMB,NCOR,XI,CSUB,CP0,RI,R2,CDL,SHAP)
  !     Find distance SQRT(R2)from source point to sub-element
      CALL SHAPEF(NDIM,NODE,RI,CK,CP,XI,R2,CDL,SHAP)
  !             Calculate Gauss orders for the sub-element
      CALL CHOSEGP(NDIM,NDIMB,NODE,CP,CK,CDL,NGR,NGL,TOLGP,RLEVL,       &
     &             DGS,XI,AL,DSQRT(R2),0,NSUB,NGSP,GPR,GWR,GPL,GWL)
  !                Integrating over the sub-element
      DO 40 IGL=1,NGSP(2)
      IF(NDIM.EQ.3)THEN; XIC(2)=GPL(IGL); WL=GWL(IGL)
      ELSE; XIC(2)=0.; WL=1.; ENDIF
      DO 40 IGR=1,NGSP(1)
      XIC(1)=GPR(IGR); WR=GWR(IGR)
  !      Find global intrinsic coordinates XI for the Gauss point IG
      CALL SHAPEF(NDIM,NCOR,XIC,CSUB,CP0,XI,R2,CDL,SHAP)
  !        Find local jacobian FJCBL for the sub-element (cell)
      CALL DSHAPE(NDIM,NDIMB,NCOR,XIC,CSUB,COSB,FJCBL,CDL,GD)
  !         Evaluate boundary integrals
      CALL INT_COEF(NDIM,NDIMB,NODE,BETA,CP,CK,XI,WL*WR,FJCBL,CDL,      &
     &              NF,V1E)     
  40  CONTINUE
  !          Compute intrinsic coordinates for next sub-element
  60  CSUB(1,1:NCOR)=CSUB(1,1:NCOR)+DGS(1)  ! Second intrinsic direction
  70  CSUB(2,1:NCOR)=CSUB(2,1:NCOR)+DGS(2)  ! Third intrinsic direction
      END
      
      SUBROUTINE INT_ELEM(NDIM,NBDM,NODE,BETA,XP,CK,XI,WXY,FJCBL,CDL,   &
     &                    NF,V1E)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XP(NDIM),XQ(NDIM),XI(NBDM),CK(3,*),RI(3),CDL(*),        &
     &          SHAP(NODE), COSN(NDIM),V1E(NF),GCD(NDIM,2),DRDX(NDIM),  &
     &          GD(3,3),FQ(NF)
      CALL SHAPEF(NDIM,NODE,XI,CK,XP,RI,RQ2,CDL,SHAP)
      CALL DSHAPE(NDIM,NBDM,NODE,XI,CK,COSN,FJCB,CDL,GD)
      RQ=DSQRT(RQ2); XQ=XP+RI; DRDX=RI/RQ
      DRDN=DOT_PRODUCT(COSN,DRDX)
      COMT=FJCB*WXY*FJCBL/RQ**BETA
      CALL F_BAR(NDIM,NBDM,DRDX,COSN,RQ,DRDN,XI,SHAP,XP,XQ,NF,FQ) 
      V1E=V1E+COMT*FQ
      ENDSUBROUTINE
    
      SUBROUTINE CHOSEGP(NDIM,NDIMB,NODE,CP,CF,CDL,NGR,NGL,TOLGP,RLEVL, &
     &                   DGS,XI,AL,DISL,NFLAG,NSUB,NGSP,GPR,GWR,GPL,GWL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CP(3),AL(3),XI(3),CF(3,NODE),DGS(3),NSUB(3),NGSP(3),    &
     &          GPR(*),GWR(*),GPL(*),GWL(*)
      DIMENSION DN(2,NODE),CDL(*),GD(3,3),COSB(NDIM),NG(2)
      NG(1)=NGR; NG(2)=NGL
      PB=DSQRT(RLEVL*2./3.+0.4); IF(NFLAG.EQ.0) GOTO 30
      NSUB(1:NDIMB)=1; NGSP(1:NDIMB)=IABS(NG(1:NDIMB))
      !NGSP(:)=IABS(NGL); IF(NDIM.EQ.2) NGSP(1)=IABS(NGR)
      DO I=NDIMB+1,3; DGS(I)=2.; NSUB(I)=1; NGSP(I)=1; ENDDO
      IF(NG(1).GT.0.AND.NG(2).GT.0) GOTO 90  ! No sub-division
  !       Calculate side lengths Li of an element using (4.34)
      CALL SETGAS(8,1,GPR,GWR,GPL,GWL,NGSS)  ! Set up Gauss points and weights
      AVL=0.; DO 20 ISID=1,NDIMB
      AL(ISID)=0.; DO I=1,NDIMB; IF(I.NE.ISID) XI(I)=0.0; ENDDO
      DO 10 IG=1,NGSS; XI(ISID)=GPR(IG)
      CALL DSHAPE(NDIM,NDIMB,NODE,XI,CF,COSB,FJCB,CDL,GD) 
  !          Find line JACOBIAN for calculation of element length
      FJCB=0.; DO I=1,NDIM; FJCB=FJCB+GD(I,ISID)*GD(I,ISID); ENDDO 
  10  AL(ISID)=AL(ISID)+DSQRT(FJCB)*GWR(IG)
  20  AVL=AVL+AL(ISID)/NDIMB
  !               Calculate minimum distance R
      CALL MINDISTD(NDIM,NDIMB,NODE,CP,CF,AVL,XI,DISL,COSB,DN,GD,CDL)
  !               Calculate Gauss integration orders
  30  WFA=PB*DLOG(DABS(TOLGP)/2.)
      DO 80 I=1,NDIMB
      IF(NG(I).GT.0) GOTO 80  
      IF(TOLGP.LT.0.) GOTO 40
      ALM=AL(I); IF(ALM.GT.3.9*DISL) ALM=3.9*DISL
      NGSP(I)=0.5*WFA/DLOG(ALM/DISL/4.)              ! Equation (4.28) of Gao's book (2002)
      GOTO 50
  40  NGSP(I)=-WFA/10.*((8./3.*AL(I)/DISL)**(3.D0/4.)+1.)    ! Eq.(4.31) of Gao's book (2002)
  50  IF(NGSP(I).LT.2) NGSP(I)=2
      IF(NFLAG.EQ.1.AND.NG(I).LT.0) GOTO 60  ! NAUTO  -->>  NG(I)
      IF(NGSP(I).GT.IABS(NG(I))) NGSP(I)=IABS(NG(I))  ! MGAUS  -->>  IABS(NG(I))
      GOTO 80
  60  IF(NGSP(I).LE.IABS(NG(I))) GOTO 70; NGSP(I)=IABS(NG(I)) ! MGAUS  -->>  IABS(NG(I))
  !       Calculating maximum length of sub-elements using (4.32)
      IF(TOLGP.LT.0.) ALI=3./8.*DISL*(-10.*NGSP(I)/WFA-1.)**(4.D0/3.)
  !       Calculating maximum length of sub-elements using (4.32) 
      IF(TOLGP.GT.0.) ALI=4.*DISL*(TOLGP/2.)**(0.5*PB/NGSP(I)) 
      NSUB(I)=AL(I)/ALI+0.95D0             ! No. of sub-elements 
      AL(I)=AL(I)/NSUB(I)                  ! Side length of sub-elements 
  70  DGS(I)=2.D0/NSUB(I)                  ! Step of intrinsic coordinates 
  80  CONTINUE
  !      Set up Gauss integration coordinates and weights
  90  CALL SETGAS(NGSP(1),NGSP(2),GPR,GWR,GPL,GWL,NGSS)
      ENDSUBROUTINE CHOSEGP
      
      SUBROUTINE SETGAS(NGSX,NGSY,POSGX,WEITX,POSGY,WEITY,NGSS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION POSGX(NGSX),POSGY(NGSY),WEITX(NGSX),WEITY(NGSY)
      CALL GAUSSV(NGSX,POSGX,WEITX)
      IF(NGSY.NE.1) CALL GAUSSV(NGSY,POSGY,WEITY)
      NGSS=NGSX*NGSY
      ENDSUBROUTINE   
   
      SUBROUTINE MINDISTD(NDIME,NDIMB,NODE,CP,CF,AVL,XIC,DISL,COSB,DN,  &
     &                    GD,CORDL) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GDT(3,3),RI(3),CP(*),CF(3,NODE),XI(3),XIC(*),REV(3,3),  &
     &          COSB(3),GD(3,3),DN(3,NODE),SHAP(NODE),CORDL(*)
!     MODIFIED ON 11/19/2002
      TOLD=1.D-11; TOLB=TOLD*10.
      NITER=50
      XI(1:NDIMB)=0.D0
      DISL=1.D12
      DO 60 ITER=1,NITER
      CALL DSHAPE(NDIME,NDIMB,NODE,XI,CF,COSB,FJCB,CORDL,GD)
      IF(DABS(FJCB).LT.TOLD) THEN    ! FOR DEGENERATED ELEMENTS
       XI(1:NDIMB)=0.382D0*XIC(1:NDIMB)+0.618D0*XI(1:NDIMB)
      CALL DSHAPE(NDIME,NDIMB,NODE,XI,CF,COSB,FJCB,CORDL,GD)
      ENDIF
!     RI(i)=X(i)-CP(i)
      CALL SHAPEF(NDIME,NODE,XI,CF,CP,RI,RQ2,CORDL,SHAP)
      RM=DSQRT(RQ2)
      IF(DABS(RM-DISL)/AVL.LT.TOLB) RETURN
!     Update the minimum value of the distance R. Once R increases, return
      IF(RM.LT.DISL) THEN
       DISL=RM
       XIC(1:NDIMB)=XI(1:NDIMB)
      ELSE
       XI(1:NDIMB)=0.382D0*XIC(1:NDIMB)+0.618D0*XI(1:NDIMB)
       GOTO 60
      ENDIF
      IF(NDIMB.EQ.NDIME) THEN
       CALL IVSNR123D(NDIME,GD,FJCB,GDT)
       GOTO 40
      ENDIF
      IF(NDIMB.EQ.1) THEN
       REV(1,1)=1.D0/FJCB/FJCB
       GOTO 20
      ENDIF
!     Multiplication of transpose of [dX(i)/dXI(K)] AND [dX(i)/dXI(K)]
      DO 10 I=1,NDIMB
      DO 10 J=1,NDIMB
      GDT(I,J)=0.D0
      DO K=1,NDIME
       GDT(I,J)=GDT(I,J)+GD(K,I)*GD(K,J)
      ENDDO
  10  CONTINUE
      FJCB=GDT(1,1)*GDT(2,2)-GDT(1,2)*GDT(2,1)
!     Inverse of matrix
      CALL IVSNR123D(NDIMB,GDT,FJCB,REV)
  20  DO 30 I=1,NDIMB
      DO 30 J=1,NDIME
      GDT(I,J)=0.D0
      DO K=1,NDIMB
       GDT(I,J)=GDT(I,J)+REV(I,K)*GD(J,K)
      ENDDO
  30  CONTINUE
  40  DO 50 I=1,NDIMB
!     Update intrinsic coordinates (XI=XI+dXI)
      TERM=0.D0
      DO K=1,NDIME
       TERM=TERM+GDT(I,K)*RI(K)
      ENDDO
  50  XI(I)=XI(I)-TERM
!     Scale intrinsic coordinates
      DO I=1,NDIMB
       IF(DABS(XI(I)).GT.1.D0) XI(I)=XI(I)/DABS(XI(I))
      ENDDO
  60  CONTINUE
      ENDSUBROUTINE MINDISTD
      
      SUBROUTINE IVSNR123D(NDIM,GD,FJCB,REV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GD(3,3),REV(3,3)
!     GD:     The original matrix to be inversed
!     REV:    Store the inverse of the matrix GD
!     NDIM:   The order of the matrix GD (=1, 2 or 3)
!     FJCB:   The determinant of the matrix GD
      GOTO(10,20,30),NDIM
  10  REV(1,1)=1.D0/FJCB
      RETURN
  20  REV(1,1)=GD(2,2)/FJCB; REV(2,2)=GD(1,1)/FJCB
      REV(1,2)=-GD(1,2)/FJCB; REV(2,1)=-GD(2,1)/FJCB
      RETURN
  30  REV(1,1)=(GD(2,2)*GD(3,3)-GD(3,2)*GD(2,3))/FJCB
      REV(1,2)=(GD(1,3)*GD(3,2)-GD(1,2)*GD(3,3))/FJCB
      REV(1,3)=(GD(1,2)*GD(2,3)-GD(2,2)*GD(1,3))/FJCB
      REV(2,1)=(GD(3,1)*GD(2,3)-GD(2,1)*GD(3,3))/FJCB
      REV(2,2)=(GD(1,1)*GD(3,3)-GD(1,3)*GD(3,1))/FJCB
      REV(2,3)=(GD(2,1)*GD(1,3)-GD(1,1)*GD(2,3))/FJCB
      REV(3,1)=(GD(2,1)*GD(3,2)-GD(3,1)*GD(2,2))/FJCB
      REV(3,2)=(GD(3,1)*GD(1,2)-GD(1,1)*GD(3,2))/FJCB
      REV(3,3)=(GD(1,1)*GD(2,2)-GD(1,2)*GD(2,1))/FJCB
      ENDSUBROUTINE IVSNR123D
    
    
      SUBROUTINE SINGULAR_ELEM(NDIM,NBDM,NODE,BETA,NPOWG,XP,CP0,CK,CDL, &
     &                         TOLGP,NGL,GPL,GWL,NODEF,NSUB,NDSID,NSP,  &
     &                         XIS,NPW,NF,V1E)
  !   THIS ROUTINE EVALUATES 2D AND 2D SINGULAR ELEMENT INTEGRALS
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XP(NDIM),XIS(NBDM),CK(3,*),CDL(*),NODEF(*),GPL(*),      &
     &          GWL(*),CP0(*),XIP(NBDM),XIQ(NBDM),CSUB(2,2),SHAP(NODE), &
     &          COEFG(0:NPOWG),COEFH(0:NPW),XI0(NBDM),KSB(4),V1E(NF),   &
     &          RINT(NF)
      TOL=1.D-5
      KSB=(/-2,1,2,-1/)


      IF(NSP.LT.0) THEN
       XIP=XIS
       CALL SHAPEF(NDIM,NODE,XIP,CK,CP0,XP,RQ2,CDL,SHAP)
       WRITE(*,*)' XP =',XP
      ENDIF
      IF(NSP.GT.0) XP(1:NDIM)=CK(1:NDIM,NSP) 

      print *,NSP
     
      PRINT *,XIP



      IF(NDIM.EQ.3) GOTO 30
  !   Evaluate 2D SINGULAR INTEGRALS
      IF(NSP.GT.0) XIP(1)=(-1)**NSP+NSP/3
      DO 20 ISUB=1,2
      XIQ(1)=(-1.D0)**ISUB
      IF(DABS(XIQ(1)-XIP(1)).LT.TOL) GOTO 20
  !   SET UP COEFFICIENTS Gn ADN Hn
      CALL COEFS_GH(NDIM,NBDM,NODE,BETA,NPOWG,CDL,CK,XP,XIP,XIQ,        &
     &              NPW,COEFG,COEFH,NF)
      CALL INT_RHO(NDIM,NBDM,NODE,BETA,NF,NPOWG,NPW,CDL,CK,CP0,         &
     &               XP,XIP,XIQ,SHAP,COEFG,COEFH,RINT)
      V1E=V1E+RINT   
  20  CONTINUE
      RETURN
  !   Evaluate 3D SINGULAR INTEGRALS    
30    CALL GAUSSV(NGL,GPL,GWL)
      IF(NSP.GT.0) XIP=CDL(2*NSP-1:2*NSP)
      WFA=DSQRT(BETA*2.D0/3.D0+0.4D0)*DLOG(DABS(TOLGP)/2.D0)
      FK=3.D0/8.D0*(-10.D0*DBLE(NGL)/WFA-1.D0)**(4.D0/3.D0)
      DO 80 ISID=1,NSUB

      KS=KSB(ISID)
      IF(DABS(XIP(IABS(KS))-DBLE(KS)/DABS(DBLE(KS))).LT.TOL) GOTO 80  
      
  !   Set up a sub-element from side ISID and the singular point NSP
      DO ID=1,2; JD=NODEF(3*(ISID-1)+ID)
       CSUB(1:2,ID)=CDL(2*JD-1:2*JD)
      ENDDO
      IXI=1; IF(ISID/2*2.EQ.ISID) IXI=2; ISWP=3-IXI
      SV=DSIGN(1.D0,CSUB(IXI,2)-CSUB(IXI,1))
      VLc2=(CSUB(ISWP,2)-XIP(ISWP))**2
      XI0=CSUB(:,1)
      XIQ(ISWP)=XI0(ISWP)
      DO 60 ISUB=1,500
      RS=DSQRT(DOT_PRODUCT(XI0-XIP,XI0-XIP))
      IF(SV*XI0(IXI)+1.D-8.GE.SV*CSUB(IXI,2)) GOTO 80
      IF(SV*XI0(IXI).GE.SV*XIP(IXI)) THEN
       XIEL=FK*RS
      ELSE
       SUM=VLc2*(1.-FK*FK)+(XIP(IXI)-XI0(IXI))**2
       IF(SUM.LT.0.D0) THEN
        XIEL=SV*(XIP(IXI)-XI0(IXI))
       ELSE
        XIEL=FK*(FK*SV*(XIP(IXI)-XI0(IXI))-DSQRT(SUM))/(FK*FK-1.)
       ENDIF
      ENDIF
      IF(SV*(XI0(IXI)+SV*XIEL)+1.D-8.GT.SV*CSUB(IXI,2))                 &
     &   XIEL=SV*(CSUB(IXI,2)-XI0(IXI))
      FJCBL=0.5D0*SV*XIEL
      DO 40 IGL=1,NGL
      XIQ(IXI)=XI0(IXI)+FJCBL*(1.D0+GPL(IGL))
      RHOQ=DSQRT(DOT_PRODUCT(XIQ-XIP,XIQ-XIP))
      DRDNP=DABS(XIQ(ISWP)-XIP(ISWP))/RHOQ
      CALL COEFS_GH(NDIM,NBDM,NODE,BETA,NPOWG,CDL,CK,XP,XIP,XIQ,        &
     &              NPW,COEFG,COEFH,NF)
      CALL INT_RHO(NDIM,NBDM,NODE,BETA,NF,NPOWG,NPW,CDL,CK,CP0,         &
     &             XP,XIP,XIQ,SHAP,COEFG,COEFH,RINT)
      COMT=DABS(FJCBL)*GWL(IGL)*DRDNP/RHOQ
  40  V1E=V1E+COMT*RINT
  60  XI0(IXI)=XI0(IXI)+SV*XIEL
  80  CONTINUE
      END

      SUBROUTINE INT_RHO(NDIM,NBDM,NODE,BETA,NF,NPOWG,NPW,CDL,CK,       &
     &                   CP0,XP,XIP,XIQ,SHAP,COEFG,COEFH,RINT)      
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XP(NDIM),XIP(NBDM),XIQ(NBDM),CK(3,*),CDL(*),CP0(*),     &
     &          DRDX(NDIM),SLOP(NBDM),COSN(NDIM),SHAP(*),RI(NDIM),      &
     &          GCD(NDIM,NBDM),COEFG(0:NPOWG),COEFH(0:NPW),             &
     &          COEFB(0:11,NF),XI(NBDM),XQ(NDIM),RINT(NF),FQ(NF)
      SLOP=XIQ-XIP
      RHOQ=DSQRT(DOT_PRODUCT(SLOP,SLOP))
      SLOP=SLOP/RHOQ
      NPOWF=3+2.1214*RHOQ    ! NPOWF IS FROM 3 TO 9
      IF(NPOWF.LT.BETA-NBDM) NPOWF=BETA-NBDM
      !   DETERMINE COEFFICIENTS Bn
      XI=XIP+RHOQ*SLOP
      CALL COEF_B(NDIM,NBDM,NODE,BETA,NPOWG,NPOWF,NF,CDL,CK,XP,XIP,     &
     &            XI,SLOP,RHOQ,RI,COSN,GCD,COEFG,COEFB)
      !   EVALUATE RESULTS USING Eqs.(3-6-44) AND (3-6-63)
      BETA1=BETA-NBDM+1
      RINT=0.D0  
      DO 40 K=0,BETA-NDIM
      PW=BETA-K-NBDM
!       print *,"beta = ",beta
!       print *,"k = ",k
!       print *,"pw = ",pw
  40  RINT=RINT-(1.D0/RHOQ**PW-COEFH(INT(PW)))/PW*COEFB(K,:)
      NBETA=BETA-NBDM  
      IF(NBETA.GE.0) RINT=RINT+COEFB(NBETA,:)*DLOG(RHOQ/COEFH(0))

      DO K=BETA-NBDM+1,NPOWF
       PW=K-BETA+NBDM
       RINT=RINT+RHOQ**PW/PW*COEFB(K,:)
      ENDDO
      END

      SUBROUTINE COEFS_GH(NDIM,NBDM,NODE,BETA,NPOWG,CDL,CK,XP,XIP,XIQ,  &
     &                    NPW,COEFG,COEFH,NF)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XP(NDIM),XIP(NBDM),XIQ(NBDM),CK(3,*),COSN(NDIM),        &
     &          GCD(NDIM,NBDM),XI(NDIM),RI(NDIM),CDL(*),SHAP(NODE),     &
     &          COEFG(0:NPOWG),COEFC(0:NPW),COEFH(0:NPW)
      ! DETERMINE COEFFICIENTS Gn
      COEFG=0.D0
      CALL COEF_G(NDIM,NBDM,NODE,NPOWG,CDL,CK,XP,XIP,XIQ, DRDN,COSN,    &
     &            GCD,RI,SHAP,COEFG)
      ! DETERMINE COEFFICIENTS Cn USING Eq. (3-6-37)
      COEFC(0)=DSQRT(COEFG(0))
      DO N=1,NPW
       COEFC(N)=0.D0
       DO I=1,N-1
        COEFC(N)=COEFC(N)-COEFC(I)*COEFC(N-I)     
       ENDDO
       IF(N.LE.NPOWG)COEFC(N)=COEFG(N)+COEFC(N)   ! Eq. (3-6-37)
       COEFC(N)=COEFC(N)/(2.D0*COEFC(0))
      ENDDO
      ! DETERMINE COEFFICIENTS Hi USING Eq.(3-6-47a)
      COEFC(1:NPW)=COEFC(1:NPW)/COEFC(0)    ! Eq.(3-6-47b)
      COEFH(0)=1.D0/COEFC(0)
      COEFH(1)=COEFC(1)
      COEFH(2)=2.*COEFC(2)-COEFC(1)*COEFC(1)
      COEFH(3)=3.*COEFC(3)-3.*COEFC(1)*COEFC(2)+COEFC(1)**3
      COEFH(4)=4.*COEFC(4)+4.*COEFC(1)*COEFC(1)*COEFC(2)                &
     &        -4.*COEFC(1)*COEFC(3)-2.*COEFC(2)*COEFC(2)                &
     &        -COEFC(1)**4
      END

      SUBROUTINE COEF_G(NDIM,NBDM,NODE,NPOWG,CDL,CK,XP,XIP,XIQ,         &
     &                  DRDN,COSN,GCD,RI,SHAP,COEFG)
      ! THIS ROUTINE DETERMINES COEFFICIENTS Gm USING Eqs.(3-6-58)~(3-6-60)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XP(NDIM),XIP(NBDM),XIQ(NBDM),CK(3,*),COSN(NDIM),        &
     &          GCD(NDIM,NBDM),XI(NDIM),SLOP(NBDM),RI(NDIM),CDL(*),     &
     &          SHAP(NODE),RMAT(NPOWG,NPOWG),COEFG(0:NPOWG)
      IF(NODE.EQ.2) THEN         ! Eq.(24)
       COEFG(0)=0.25*((CK(1,1)-CK(1,2))**2+(CK(2,1)-CK(2,2))**2)  
      ELSEIF(NODE.EQ.3) THEN             ! FORM Eq.(3-6-33)
       RI=CK(:,1)+CK(:,2)-2.*CK(:,3)     ! Xbar
       XI=CK(:,1)-CK(:,2)                ! X1i-X2i 
       XBAR2=RI(1)*RI(1)+RI(2)*RI(2)     ! Xbar^2
       XXBAR=XI(1)*RI(1)+XI(2)*RI(2)
       COEFG(0)=0.25*(XI(1)*XI(1)+XI(2)*XI(2))-XXBAR*XIP(1)             &
     &         +XBAR2*XIP(1)*XIP(1)
       COEFG(1)=(-0.5*XXBAR+XBAR2*XIP(1))*XIQ(1)
       COEFG(2)=0.25*XBAR2
      ELSE      ! QUADRILATERAL ELEMENTS
       SLOP=XIQ-XIP
       RHOQ=DSQRT(DOT_PRODUCT(SLOP,SLOP))
       SLOP=SLOP/RHOQ
       VK=RHOQ/DBLE(NPOWG)
       DO 20 IP=0,NPOWG
       RHO=VK*DBLE(IP)
       XI(1:NBDM)=XIP+RHO*SLOP
       CALL SHAPEF(NDIM,NODE,XI,CK,XP,RI,R2,CDL,SHAP)
       IF(IP.NE.0) GOTO 10
       CALL DSHAPE(NDIM,NBDM,NODE,XI,CK,COSN,FJCB,CDL,GCD)   !dX/dXI
       RI=MATMUL(GCD,SLOP); COEFG(0)=DOT_PRODUCT(RI,RI)  ! Eq.(3-6-58)
       GOTO 20
 10    COEFG(IP)=((DSQRT(R2)/RHO)**2-COEFG(0))/RHO     ! Eq.(3-6-60)
       RMAT(IP,1)=1.D0        ! Eq.(3-6-42)
       DO JP=2,NPOWG; RMAT(IP,JP)=RMAT(IP,JP-1)*RHO; ENDDO
 20    CONTINUE
       CALL INVSOLVR(NPOWG,NPOWG,RMAT,NPOWG,-1)       ! INVERSE MATRIX [R] IN Eq.(3-6-59)
       COEFG(1:NPOWG)=MATMUL(RMAT,COEFG(1:NPOWG))     ! SOLVING Eq.(3-6-59) FOR {G}
      ENDIF
      END

      SUBROUTINE COEF_B(NDIM,NBDM,NODE,BETA,NPOWG,NPOWF,NF,CDL,CK,      &
     &                  XP,XIP,XIQ,SLOP,RHOQ,RI,COSN,GCD,COEFG,COEFB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XP(NDIM),CK(3,*),DRDX(NDIM),CDL(*),COSN(NDIM),          &
     &          GCD(NDIM,*),XIP(NBDM),XIQ(NBDM),XI(NBDM),SLOP(NBDM),    &
     &          RI(NDIM),X(NDIM),RMAT(NPOWF,NPOWF),COEFG(0:NPOWG),      &
     &          COEFB(0:11,NF),SHAP(NODE),FQ(NF),A(NDIM)
      VK=RHOQ/DBLE(NPOWF)
      DO 20 IP=0,NPOWF
      RHO=VK*DBLE(IP)
      XI=XIP+RHO*SLOP
      ROBAR=COEFG(0); DO M=1,NPOWG; ROBAR=ROBAR+COEFG(M)*RHO**M; ENDDO
      ROBAR=DSQRT(ROBAR)       ! Eq.(3-6-28)
      CALL SHAPEF(NDIM,NODE,XI,CK,XP,RI,R2,CDL,SHAP)
      R=DSQRT(R2); X=XP+RI
      CALL DSHAPE(NDIM,NBDM,NODE,XI,CK,COSN,FJCB,CDL,GCD)
      IF(RHO.GT.1.0D-10)THEN    
       DRDX=RI/R
      ELSE
       DO 10 I=1,NDIM; A(I)=0.D0
       DO 10 J=1,NBDM            
       A(I)=A(I)+GCD(I,J)*SLOP(J)
10     CONTINUE                  
       GM=DSQRT(DOT_PRODUCT(A,A))
       DRDX=A/GM                 ! Eq.(3-6-74)
      ENDIF 
      DRDN=DOT_PRODUCT(COSN,DRDX)
      CALL F_BAR(NDIM,NBDM,DRDX,COSN,R,DRDN,XI,SHAP,XP,X,NF,FQ)
      COEFB(IP,:)=FQ*FJCB/ROBAR**BETA
      IF(IP.EQ.0) GOTO 20
      COEFB(IP,:)=(COEFB(IP,:)-COEFB(0,:))/RHO
      RMAT(IP,1)=1.D0
      DO JP=2,NPOWF; RMAT(IP,JP)=RMAT(IP,JP-1)*RHO; ENDDO
 20   CONTINUE
      CALL INVSOLVR(NPOWF,NPOWF,RMAT,NPOWF,-1)
      COEFB(1:NPOWF,:)=MATMUL(RMAT,COEFB(1:NPOWF,:))    ! Bk IN Eq.(3-6-62)
      END

      SUBROUTINE GAUSSV(NGAUS,GP,GW)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION GP(NGAUS),GW(NGAUS)
  !   NGAUS:  The number of Gauss integration points (from 2 to 10)
  !   GP:     The coordinates of Gauss points
  !   GW:     The weighting factors of Gauss points
      SELECT CASE(NGAUS)
      CASE(2); GP(1)=-0.57735026918962576451D0; GW(1)=1.0000000000000D0
      CASE(3); GP(1)=-0.77459666924148337704D0; GP(2)=0.0000000000000D0
       GW(1)=0.555555555555555555556D0; GW(2)=0.888888888888888888889D0
      CASE(4)
       GP(1)=-0.86113631159405257522D0; GP(2)=-0.33998104358485626480D0
       GW(1)=0.34785484513745385737D0;  GW(2)=0.65214515486254614263D0
      CASE(5)
       GP(1)=-0.90617984593866399280D0; GP(2)=-0.53846931010568309104D0
       GP(3)=0.00000000000000000000D0;  GW(1)=0.23692688505618908751D0
       GW(2)=0.47862867049936646804D0;  GW(3)=0.56888888888888888889D0
      CASE(6)
       GP(1)=-0.93246951420315202781D0; GP(2)=-0.66120938646626451366D0
       GP(3)=-0.23861918608319690863D0; GW(1)=0.17132449237917034504D0
       GW(2)=0.36076157304813860757D0;  GW(3)=0.46791393457269104739D0
      CASE(7)
       GP(1)=-0.94910791234275852453D0; GP(2)=-0.74153118559939443986D0
       GP(3)=-0.40584515137739716691D0; GP(4)=0.00000000000000000000D0
       GW(1)=0.12948496616886969327D0;  GW(2)=0.27970539148927666790D0
       GW(3)=0.38183005050511894495D0;  GW(4)=0.41795918367346938776D0
      CASE(8)
       GP(1)=-0.96028985649753623168D0; GP(2)=-0.79666647741362673959D0
       GP(3)=-0.52553240991632898582D0; GP(4)=-0.18343464249564980494D0
       GW(1)=0.10122853629037625915D0;  GW(2)=0.22238103445337447054D0
       GW(3)=0.31370664587788728734D0;  GW(4)=0.36268378337836198297D0
      CASE(9)
       GP(1)=-0.96816023950762608984D0; GP(2)=-0.83603110732663579430D0
       GP(3)=-0.61337143270059039731D0; GP(4)=-0.32425342340380892904D0
       GP(5)=0.00000000000000000000D0;  GW(1)=0.08127438836157441197D0
       GW(2)=0.18064816069485740406D0;  GW(3)=0.26061069640293546232D0
       GW(4)=0.31234707704000284007D0;  GW(5)=0.33023935500125976317D0
      CASE(10)
       GP(1)=-0.97390652851717172008D0; GP(2)=-0.86506336668898451073D0
       GP(3)=-0.67940956829902440623D0; GP(4)=-0.43339539412924719080D0
       GP(5)=-0.14887433898163121089D0; GW(1)=0.06667134430868813759D0
       GW(2)=0.14945134915058059315D0;  GW(3)=0.21908636251598204400D0
       GW(4)=0.26926671930999635509D0;  GW(5)=0.29552422471475287017D0
      END SELECT
      KGAUS=NGAUS/2; DO IGASH=1,KGAUS; JGASH=NGAUS+1-IGASH
      GP(JGASH)=-GP(IGASH); GW(JGASH)=GW(IGASH); ENDDO
      END 

      SUBROUTINE INVSOLVR(NROW,NCOL,A,N,INDIC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NROW,NCOL),W(N),IROW(N)
      EPS=1.D-8             ! The tolerance of the minimum pivot
      DO 40 K=1,N; KM1=K-1; PIVOT=0.D0
      DO 20 I=1,N; IF(K.EQ.1) GOTO 10
      DO ISCAN=1,KM1; IF(I.EQ.IROW(ISCAN)) GOTO 20; ENDDO
  10  IF(DABS(A(I,K)).LE.DABS(PIVOT)) GOTO 20
      PIVOT=A(I,K); IROW(K)=I
  20  CONTINUE
      IF(DABS(PIVOT).GT.EPS) GOTO 30; STOP 9999
  30  IROWK=IROW(K)
      A(IROWK,1:NCOL)=A(IROWK,1:NCOL)/PIVOT; A(IROWK,K)=1.D0/PIVOT
      DO 40 I=1,N; AIK=A(I,K); IF(I.EQ.IROWK) GOTO 40
      A(I,1:NCOL)=A(I,1:NCOL)-AIK*A(IROWK,1:NCOL); A(I,K)=-AIK/PIVOT
  40  CONTINUE
      IF(INDIC.LT.0) GOTO 60
      DO 50 IX=N+1,NCOL; W(1:N)=A(IROW(1:N),IX); DO 50 I=1,N
  50  A(I,IX)=W(I)
      IF(INDIC.GT.0) RETURN
  60  DO 70 J=1,N; W(1:N)=A(IROW(1:N),J); DO 70 I=1,N
  70  A(I,J)=W(I)
      DO 80 I=1,N; W(IROW(1:N))=A(I,1:N); DO 80 J=1,N
  80  A(I,J)=W(J)
      END

      SUBROUTINE SHAPEF(NDIM,NODE,X,CK,XP,RI,RQ2,C,SP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SP(NODE),CK(3,*),XP(*),RI(*),C(*),X(*)
      IF(NODE.GT.3) GOTO 4
  !                  2-noded line element
      SP(1)=0.5*(1.-X(1)); SP(2)=0.5*(1.+X(1))
      IF(NODE.EQ.2) GOTO 50
  !                  3-noded line element
      SP(1)=-X(1)*SP(1); SP(2)=X(1)*SP(2); SP(3)=1.-X(1)*X(1)
      GOTO 50
  !                  4-noded quadrilateral element
  4   DO I=1,4
       SP(I)=0.25*(1.+C(2*I-1)*X(1))*(1.+C(2*I)*X(2))
      ENDDO
      IF(NODE.EQ.8) THEN
  !                  8 noded-element (square element)
       DO 15 I=1,4; L=2*I-1; SP(I)=SP(I)*(C(L)*X(1)+C(L+1)*X(2)-1.D0)
       WL=C(L+8)*X(1)+C(L+9)*X(2)
  15   SP(I+4)=.5D0*(WL+1.D0)*(1.D0-(C(L+8)*X(2))**2-(C(L+9)*X(1))**2)
      ELSEIF(NODE.EQ.9) THEN
  !                  9 noded-element (square element)
       DO 20 I=1,4; L=2*I-1; SP(I)=SP(I)*C(L)*X(1)*C(L+1)*X(2)
       WL=C(L+8)*X(1)+C(L+9)*X(2)
  20   SP(I+4)=0.5D0*WL*(WL+1.D0)*(1.D0-(C(L+8)*X(2))**2-               &
     &                  (C(L+9)*X(1))**2)
       SP(9)=(1.D0-X(1)*X(1))*(1.D0-X(2)*X(2))
      ENDIF
  !                  Calculate r and its vector components
  50  RQ2=0.; DO 70 I=1,NDIM; RI(I)=-XP(I)
      DO ID=1,NODE; RI(I)=RI(I)+SP(ID)*CK(I,ID); ENDDO
  70  RQ2=RQ2+RI(I)*RI(I)
      END

      SUBROUTINE DSHAPE(NDIM,NBDM,NODE,X,CK,COSN,FJCB,C,GD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),CK(3,*),DN(2,NODE),GD(3,*),COSN(*),C(*),GR(3)
      IF(NODE.GT.3) GOTO 5
      DN(1,1)=-0.5; DN(1,2)=0.5        ! 2-noded line element
      IF(NODE.EQ.2) GOTO 30
      DN(1,1)=-0.5*(1.-2.*X(1))        ! 3-noded line element
      DN(1,2)=0.5*(1.+2.*X(1)); DN(1,3)=-2.*X(1); GOTO 30
   5  DO 10 I=1,4; I0=2*(I-1)          ! 4-noded quadr. element
      DN(1,I)=0.25*C(I0+1)*(1.+C(I0+2)*X(2))
  10  DN(2,I)=0.25*C(I0+2)*(1.+C(I0+1)*X(1))
      IF(NODE.EQ.8) THEN
!                  8 noded-element (square element)
       DO 15 I=1,4; L=2*I-1; S=C(L)*X(1)+C(L+1)*X(2)-1.D0
       DN(1,I)=DN(1,I)*S+0.25D0*(1.D0+C(L)*X(1))*(1.D0+C(L+1)*X(2))*C(L)
       DN(2,I)=DN(2,I)*S+0.25D0*(1.D0+C(L)*X(1))*(1.D0+C(L+1)*X(2))*    &
     &                 C(L+1)
       S=1.D0+C(L+8)*X(1)+C(L+9)*X(2)      
       T=1.D0-(C(L+8)*X(2))**2-(C(L+9)*X(1))**2
       DN(1,I+4)=0.5D0*C(L+8)*T-C(L+9)*C(L+9)*X(1)*S
  15   DN(2,I+4)=0.5D0*C(L+9)*T-C(L+8)*C(L+8)*X(2)*S
      ELSEIF(NODE.EQ.9) THEN
  !                 9 noded-element (square element)
       DO 20 I=1,4; L=2*I-1
       DN(1,I)=DN(1,I)*C(L+1)*X(2)*(1.+2.D0*C(L)*X(1))
       DN(2,I)=DN(2,I)*C(L)*X(1)*(1.+2.D0*C(L+1)*X(2))
       S=C(L+8)*X(1)+C(L+9)*X(2)      
       SS=S*(1.D0+S)
       T=(1.D0-(C(L+8)*X(2))**2-(C(L+9)*X(1))**2)*(1.D0+2.D0*S)
       XI2=C(L+8)*C(L+8); ET2=C(L+9)*C(L+9)
       DN(1,I+4)=0.5D0*(C(L+8)*T-2.D0*ET2*X(1)*SS)
  20   DN(2,I+4)=0.5D0*(C(L+9)*T-2.D0*XI2*X(2)*SS)
       DN(1,9)=-2.D0*X(1)*(1.D0-X(2)*X(2))
       DN(2,9)=-2.D0*X(2)*(1.D0-X(1)*X(1))
      ENDIF
  30  DO 50 I=1,NDIM; DO 50 J=1,NBDM; GD(I,J)=0.; DO 50 ID=1,NODE
50    GD(I,J)=GD(I,J)+DN(J,ID)*CK(I,ID)
      IF(NDIM.EQ.NBDM) THEN
       GOTO (51,52,53),NDIM
  51   FJCB=DABS(GD(1,1))
       RETURN          ! For internal cell integral
  52   FJCB=GD(1,1)*GD(2,2)-GD(1,2)*GD(2,1)
       RETURN
  53   FJCB=GD(1,1)*(GD(2,2)*GD(3,3)-GD(3,2)*GD(2,3))                   &
     &     -GD(1,2)*(GD(2,1)*GD(3,3)-GD(3,1)*GD(2,3))+                  &
     &      GD(1,3)*(GD(2,1)*GD(3,2)-GD(2,2)*GD(3,1))
       RETURN
      ENDIF
      IF(NODE.GT.3) GOTO 60
      GR(1)=GD(2,1); GR(2)=-GD(1,1)             ! For 2D normals
      FJCB=DSQRT(DOT_PRODUCT(GD(1:NDIM,1),GD(1:NDIM,1))) ! LINE JACOBIAN
      GOTO 70
  60  GR(1)=GD(2,1)*GD(3,2)-GD(3,1)*GD(2,2)     ! For 3D normals
      GR(2)=GD(3,1)*GD(1,2)-GD(1,1)*GD(3,2)
      GR(3)=GD(1,1)*GD(2,2)-GD(2,1)*GD(1,2)
      FJCB=DSQRT(GR(1)*GR(1)+GR(2)*GR(2)+GR(3)*GR(3)) ! 3D JACOBIAN
  70  COSN(1:NDIM)=GR(1:NDIM)/FJCB              ! 2D and 3D Normal
      END
      
      SUBROUTINE F_BAR(NDIM,NBDM,DRDX,COSN,R,DRDN,XI,SHAP,XP,XQ,NF,FB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XP(NDIM),X(NDIM),XI(NBDM),DRDX(NDIM),COSN(NDIM),SHAP(*),&
     &          FB(NF),XQ(NDIM)
      COMMON/RIM_COEF/PI,DLT(3,3),CNU

!     COMPUTING Gij (The first example (2D))
!      FB(1)=-(1.-2.*DRDX(1)*DRDX(1))/(2.*PI)
!      FB(2)=2.*DRDX(1)*DRDX(2)/(2.*PI)
!      FB(3)=-(1.-2.*DRDX(2)*DRDX(2))/(2.*PI)

!     CYLINDER SURFACE ELEMENT (The second example (3D))
      FB(1)=-(3.*DRDX(3)*DRDN-COSN(3))/(4.*PI)    ! GUIG 4.2


!        IF(R.GT.1.D-18) FB(1)=R*DLOG(R)
!        IF(R.LT.1.D-18) FB(1)=0.


!	  FB(1)=-DRDX(1)      ! Karami 4.3

!        FB(1)=DRDX(1)    
!        FB(2)=DRDX(2)    

!	  FB=DRDX

!!      FORM Tij OF 2D
!      CON=-1./(4.*PI*(1.-CNU)); PR2=1.-2.*CNU; N=0
!      DO I=1,NDIM; DO J=1,NDIM; N=N+1
!       FB(N)=CON*(PR2*(COSN(I)*DRDX(J)-COSN(J)*DRDX(I))+               &
!     &           (PR2*DLT(I,J)+2*DRDX(I)*DRDX(J))*DRDN)
!      ENDDO; ENDDO
      
!     FORM Tij OF 3D
     ! NBDM=NDIM-1; PR1=1.-CNU; PR2=1.-2.*CNU
     ! CON=-1./(4.*NBDM*PI*PR1)
     ! N=0
     ! DO I=1,NDIM; DO J=1,NDIM
     !  N=N+1
     !  FB(N)=CON*(DRDN*(PR2*DLT(I,J)+NDIM*DRDX(I)*DRDX(J))-             &
     !& PR2*(DRDX(I)*COSN(J)-DRDX(J)*COSN(I)))
     ! ENDDO; ENDDO

!      FB(1)=DRDX(2)    ! Karami 4.1
!     FB(1)=1.D0  ! GUIG 4.1


!      FB(1)=DRDX(1)  ! KARAMI3

!       FB(1)=(1-2.*DRDX(1)**2)
!       FB(2)=(1-2.*DRDX(2)**2)

!     COMPUTING Gij
!      FB(1)=(1.-2.*DRDX(1)*DRDX(1))/(2.*PI)
!      FB(2)=-2.*DRDX(1)*DRDX(2)/(2.*PI)
!      FB(3)=(1.-2.*DRDX(2)*DRDX(2))/(2.*PI)

!     COMPUTING Uijk
!      CON=1./(4.*PI*(1.-CNU)); PR2=1.-2.*CNU
!      PHI2=SHAP(3)        
!      FB(1)=CON*PHI2*(PR2*DRDX(1)+2.*DRDX(1)**3)
!      FB(2)=CON*PHI2*(PR2*DRDX(1)+2.*DRDX(1)*DRDX(2)**2)
!      FB(3)=CON*PHI2*(PR2*DRDX(2)+2.*DRDX(2)**3)

!     COMPUTING Tijk 
!      E=1.; G=E/(2.*(1.+CNU)); CON=G/(2.*PI*(1.-CNU))
!      C=1.-2.*CNU; D=1.-4.*CNU; K=2
!      M=0; DO 10 I=1,NDIM; DO 10 J=I,NDIM; M=M+1
!      TEM=2.*DRDN*(C*DLT(I,J)*DRDX(K)+CNU*(DLT(I,K)*DRDX(J)+DLT(J,K)*DRDX(I))-   &
!     & 4.*DRDX(I)*DRDX(J)*DRDX(K))+2.*CNU*(COSN(I)*DRDX(J)*DRDX(K)+COSN(J)*DRDX(I)*       &
!     & DRDX(K))+C*(2.*COSN(K)*DRDX(I)*DRDX(J)+COSN(J)*DLT(I,K)+COSN(I)*          &
!     & DLT(J,K))-D*COSN(K)*DLT(I,J)
!      FB(M)=CON*TEM*SHAP(M)       
!  10  CONTINUE

      END
