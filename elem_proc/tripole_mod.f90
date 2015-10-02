module tripole_mod
    use mod_func!!spfunc is used 
    use extend_mesh
    implicit none
    !real(8) :: xc,yc,zc

    integer :: NOSAMP
    real(8) :: XYNOD(3,50),DXYNOD(6,50),SAMNOD(50,0:8)
contains
    SUBROUTINE TRIPOL(NODNUM,NCNE,XYZT,DXYZT)

    IMPLICIT NONE
    INTEGER NCNE,I,J,LK,LJ,LI,IXRHO,IYRHO,NOTRI,NODNUM
    REAL*8 XIQ(4),WIQ(4)
    REAL*8 SF(8),DSF(2,8),XJ(3,3)
    REAL*8 SITRI(3,3),ETATRI(3,3)
    REAL*8 XYZT(3,8),DXYZT(3,8)
    REAL*8 X,Y,Z,DX,DY,DZ,SI,ETA,XRHO,YRHO,AG,&
      & DET,DET1,DET2,DET3,DUM,DUM1,DUM2,DUM3 
    REAL(8) :: XYZC(3)
    
!C
!C
!C
!!C** matrix XIQ store the Gauss-Legendre sampling points (4) for the 
!!C   quadrilateral element 
!C
      DATA XIQ/-0.861136311594053D0,-0.339981043584856D0,&
              & 0.339981043584856D0, 0.861136311594053D0/
!C
!!C** matrix WIQ store the Gauss-Legendre weighting fa!Ctors for the 
!!C   quadrilateral element 
!C
      DATA WIQ/0.347854845137454D0,0.652145154862546D0,&
            & 0.652145154862546D0,0.347854845137454D0/
!C
!C
!!C** identity the local nodal position of INODE IN IELEM
!C
        IF(NCNE.EQ.8) THEN
          IF(NODNUM.EQ.1.OR.NODNUM.EQ.3.OR.&
         & NODNUM.EQ.5.OR.NODNUM.EQ.7) THEN
           NOTRI=2
          ELSE
           NOTRI=3
          ENDIF
        ELSE
          IF(NODNUM.LE.3)THEN
            NOTRI=1
          ELSE
            NOTRI=2
          ENDIF
        ENDIF
!C
!!C** identity the local coordinate(SI,ETA) of the trianges
!C
      IF(NCNE.EQ.8) THEN
!C
!!C** quadrilateral element
!C
      IF(NODNUM.EQ.1) THEN
!C
      SITRI (1,1)=-1.0D0
      SITRI (1,2)= 1.0D0
      SITRI (1,3)=-1.0D0
!C
      SITRI (2,1)=-1.0D0
      SITRI (2,2)= 1.0D0
      SITRI (2,3)= 1.0D0
!C
      SITRI (3,1)= 0.0D0
      SITRI (3,2)= 0.0D0
      SITRI (3,3)= 0.0D0
!C
      ETATRI(1,1)=-1.0D0
      ETATRI(1,2)= 1.0D0
      ETATRI(1,3)= 1.0D0
!C
      ETATRI(2,1)=-1.0D0
      ETATRI(2,2)=-1.0D0
      ETATRI(2,3)= 1.0D0
!C
      ETATRI(3,1)= 0.0D0
      ETATRI(3,2)= 0.0D0
      ETATRI(3,3)= 0.0D0
!C
      ELSE IF(NODNUM.EQ.3) THEN
!C
      SITRI (1,1)= 1.0D0
      SITRI (1,2)=-1.0D0
      SITRI (1,3)=-1.0D0
!C
      SITRI (2,1)= 1.0D0
      SITRI (2,2)= 1.0D0
      SITRI (2,3)=-1.0D0
!C
      SITRI (3,1)= 0.0D0
      SITRI (3,2)= 0.0D0
      SITRI (3,3)= 0.0D0
!C
      ETATRI(1,1)=-1.0D0
      ETATRI(1,2)= 1.0D0
      ETATRI(1,3)=-1.0D0
!C
      ETATRI(2,1)=-1.0D0
      ETATRI(2,2)= 1.0D0
      ETATRI(2,3)= 1.0D0
!C
      ETATRI(3,1)= 0.0D0
      ETATRI(3,2)= 0.0D0
      ETATRI(3,3)= 0.0D0
!C
      ELSE IF(NODNUM.EQ.5) THEN
!C
      SITRI (1,1)= 1.0D0
      SITRI (1,2)=-1.0D0
      SITRI (1,3)=-1.0D0
!C
      SITRI (2,1)= 1.0D0
      SITRI (2,2)=-1.0D0
      SITRI (2,3)= 1.0D0
!C
      SITRI (3,1)= 0.0D0
      SITRI (3,2)= 0.0D0
      SITRI (3,3)= 0.0D0
!C
      ETATRI(1,1)= 1.0D0
      ETATRI(1,2)= 1.0D0
      ETATRI(1,3)=-1.0D0
!C
      ETATRI(2,1)= 1.0D0
      ETATRI(2,2)=-1.0D0
      ETATRI(2,3)=-1.0D0
!C
      ETATRI(3,1)= 0.0D0
      ETATRI(3,2)= 0.0D0
      ETATRI(3,3)= 0.0D0
!C
      ELSE IF(NODNUM.EQ.7) THEN
!C
      SITRI (1,1)=-1.0D0
      SITRI (1,2)=-1.0D0
      SITRI (1,3)= 1.0D0
!C
      SITRI (2,1)=-1.0D0
      SITRI (2,2)= 1.0D0
      SITRI (2,3)= 1.0D0
!C
      SITRI (3,1)= 0.0D0
      SITRI (3,2)= 0.0D0
      SITRI (3,3)= 0.0D0
!C
      ETATRI(1,1)= 1.0D0
      ETATRI(1,2)=-1.0D0
      ETATRI(1,3)=-1.0D0
!C
      ETATRI(2,1)= 1.0D0
      ETATRI(2,2)=-1.0D0
      ETATRI(2,3)= 1.0D0
!C
      ETATRI(3,1)= 0.0D0
      ETATRI(3,2)= 0.0D0
      ETATRI(3,3)= 0.0D0
!C
      ELSE IF(NODNUM.EQ.2) THEN
!C
      SITRI (1,1)= 0.0D0
      SITRI (1,2)=-1.0D0
      SITRI (1,3)=-1.0D0
!C
      SITRI (2,1)= 0.0D0
      SITRI (2,2)= 1.0D0
      SITRI (2,3)=-1.0D0
!C
      SITRI (3,1)= 0.0D0
      SITRI (3,2)= 1.0D0
      SITRI (3,3)= 1.0D0
!C
      ETATRI(1,1)=-1.0D0
      ETATRI(1,2)= 1.0D0
      ETATRI(1,3)=-1.0D0
!C
      ETATRI(2,1)=-1.0D0
      ETATRI(2,2)= 1.0D0
      ETATRI(2,3)= 1.0D0
!C
      ETATRI(3,1)=-1.0D0
      ETATRI(3,2)=-1.0D0
      ETATRI(3,3)= 1.0D0
!C
      ELSE IF(NODNUM.EQ.4) THEN
!C
      SITRI (1,1)= 1.0D0
      SITRI (1,2)= 1.0D0
      SITRI (1,3)=-1.0D0
!C
      SITRI (2,1)= 1.0D0
      SITRI (2,2)=-1.0D0
      SITRI (2,3)=-1.0D0
!C
      SITRI (3,1)= 1.0D0
      SITRI (3,2)=-1.0D0
      SITRI (3,3)= 1.0D0
!C
      ETATRI(1,1)= 0.0D0
      ETATRI(1,2)= 1.0D0
      ETATRI(1,3)= 1.0D0
!C
      ETATRI(2,1)= 0.0D0
      ETATRI(2,2)= 1.0D0
      ETATRI(2,3)=-1.0D0
!C
      ETATRI(3,1)= 0.0D0
      ETATRI(3,2)=-1.0D0
      ETATRI(3,3)=-1.0D0
!C
      ELSE IF(NODNUM.EQ.6) THEN
!C
      SITRI (1,1)= 0.0D0
      SITRI (1,2)=-1.0D0
      SITRI (1,3)=-1.0D0
!C
      SITRI (2,1)= 0.0D0
      SITRI (2,2)=-1.0D0
      SITRI (2,3)= 1.0D0
!C
      SITRI (3,1)= 0.0D0
      SITRI (3,2)= 1.0D0
      SITRI (3,3)= 1.0D0
!C
      ETATRI(1,1)= 1.0D0
      ETATRI(1,2)= 1.0D0
      ETATRI(1,3)=-1.0D0
!C
      ETATRI(2,1)= 1.0D0
      ETATRI(2,2)=-1.0D0
      ETATRI(2,3)=-1.0D0
!C
      ETATRI(3,1)= 1.0D0
      ETATRI(3,2)=-1.0D0
      ETATRI(3,3)= 1.0D0
!C
      ELSE IF(NODNUM.EQ.8) THEN
!C
      SITRI (1,1)=-1.0D0
      SITRI (1,2)= 1.0D0
      SITRI (1,3)=-1.0D0
!C
      SITRI (2,1)=-1.0D0
      SITRI (2,2)= 1.0D0
      SITRI (2,3)= 1.0D0
!C
      SITRI (3,1)=-1.0D0
      SITRI (3,2)=-1.0D0
      SITRI (3,3)= 1.0D0
!C
      ETATRI(1,1)= 0.0D0
      ETATRI(1,2)= 1.0D0
      ETATRI(1,3)= 1.0D0
!C
      ETATRI(2,1)= 0.0D0
      ETATRI(2,2)=-1.0D0
      ETATRI(2,3)= 1.0D0
!C
      ETATRI(3,1)= 0.0D0
      ETATRI(3,2)=-1.0D0
      ETATRI(3,3)=-1.0D0
!C
      ENDIF
!C
      ELSE
!C
!!C** triangular element ***************
!!C         
      IF(NODNUM.EQ.1) THEN
!C
      SITRI (1,1)= 0.0D0
      SITRI (1,2)= 1.0D0
      SITRI (1,3)= 0.0D0
!C
      SITRI (2,1)= 0.0D0
      SITRI (2,2)= 0.0D0
      SITRI (2,3)= 0.0D0
!C
      SITRI (3,1)= 0.0D0
      SITRI (3,2)= 0.0D0
      SITRI (3,3)= 0.0D0
!C
      ETATRI(1,1)= 0.0D0
      ETATRI(1,2)= 0.0D0
      ETATRI(1,3)= 1.0D0
!C
      ETATRI(2,1)= 0.0D0
      ETATRI(2,2)= 0.0D0
      ETATRI(2,3)= 0.0D0
!C
      ETATRI(3,1)= 0.0D0
      ETATRI(3,2)= 0.0D0
      ETATRI(3,3)= 0.0D0
!C
      ELSE IF(NODNUM.EQ.2) THEN
!C
      SITRI (1,1)= 1.0D0
      SITRI (1,2)= 0.0D0
      SITRI (1,3)= 0.0D0
!C
      SITRI (2,1)= 0.0D0
      SITRI (2,2)= 0.0D0
      SITRI (2,3)= 0.0D0
!C
      SITRI (3,1)= 0.0D0
      SITRI (3,2)= 0.0D0
      SITRI (3,3)= 0.0D0
!C
      ETATRI(1,1)= 0.0D0
      ETATRI(1,2)= 1.0D0
      ETATRI(1,3)= 0.0D0
!C
      ETATRI(2,1)= 0.0D0
      ETATRI(2,2)= 0.0D0
      ETATRI(2,3)= 0.0D0
!C
      ETATRI(3,1)= 0.0D0
      ETATRI(3,2)= 0.0D0
      ETATRI(3,3)= 0.0D0
!C
      ELSE IF(NODNUM.EQ.3) THEN
!C
      SITRI (1,1)= 0.0D0
      SITRI (1,2)= 0.0D0
      SITRI (1,3)= 1.0D0
!C
      SITRI (2,1)= 0.0D0
      SITRI (2,2)= 0.0D0
      SITRI (2,3)= 0.0D0
!C
      SITRI (3,1)= 0.0D0
      SITRI (3,2)= 0.0D0
      SITRI (3,3)= 0.0D0
!C
      ETATRI(1,1)= 1.0D0
      ETATRI(1,2)= 0.0D0
      ETATRI(1,3)= 0.0D0
!C
      ETATRI(2,1)= 0.0D0
      ETATRI(2,2)= 0.0D0
      ETATRI(2,3)= 0.0D0
!C
      ETATRI(3,1)= 0.0D0
      ETATRI(3,2)= 0.0D0
      ETATRI(3,3)= 0.0D0
!C
      ELSE IF(NODNUM.EQ.4) THEN
!C
      SITRI (1,1)= 0.5D0
      SITRI (1,2)= 0.0D0
      SITRI (1,3)= 0.0D0
!C
      SITRI (2,1)= 0.5D0
      SITRI (2,2)= 1.0D0
      SITRI (2,3)= 0.0D0
!C
      SITRI (3,1)= 0.0D0
      SITRI (3,2)= 0.0D0
      SITRI (3,3)= 0.0D0
!C
      ETATRI(1,1)= 0.0D0
      ETATRI(1,2)= 1.0D0
      ETATRI(1,3)= 0.0D0
!C
      ETATRI(2,1)= 0.0D0
      ETATRI(2,2)= 0.0D0
      ETATRI(2,3)= 1.0D0
!C
      ETATRI(3,1)= 0.0D0
      ETATRI(3,2)= 0.0D0
      ETATRI(3,3)= 0.0D0
!C
      ELSE IF(NODNUM.EQ.5) THEN
!C
      SITRI (1,1)= 0.5D0
      SITRI (1,2)= 0.0D0
      SITRI (1,3)= 1.0D0
!C
      SITRI (2,1)= 0.5D0
      SITRI (2,2)= 0.0D0
      SITRI (2,3)= 0.0D0
!C
      SITRI (3,1)= 0.0D0
      SITRI (3,2)= 0.0D0
      SITRI (3,3)= 0.0D0
!C
      ETATRI(1,1)= 0.5D0
      ETATRI(1,2)= 0.0D0
      ETATRI(1,3)= 0.0D0
!C
      ETATRI(2,1)= 0.5D0
      ETATRI(2,2)= 1.0D0
      ETATRI(2,3)= 0.0D0
!C
      ETATRI(3,1)= 0.0D0
      ETATRI(3,2)= 0.0D0
      ETATRI(3,3)= 0.0D0
!C
      ELSE IF(NODNUM.EQ.6) THEN
!C
      SITRI (1,1)= 0.0D0
      SITRI (1,2)= 0.0D0
      SITRI (1,3)= 1.0D0
!C
      SITRI (2,1)= 0.0D0
      SITRI (2,2)= 1.0D0
      SITRI (2,3)= 0.0D0
!C
      SITRI (3,1)= 0.0D0
      SITRI (3,2)= 0.0D0
      SITRI (3,3)= 0.0D0
!C
      ETATRI(1,1)= 0.5D0
      ETATRI(1,2)= 0.0D0
      ETATRI(1,3)= 0.0D0
!C
      ETATRI(2,1)= 0.5D0
      ETATRI(2,2)= 0.0D0
      ETATRI(2,3)= 1.0D0
!C
      ETATRI(3,1)= 0.0D0
      ETATRI(3,2)= 0.0D0
      ETATRI(3,3)= 0.0D0
!C
      ENDIF
!C
      ENDIF
!C
      NOSAMP=0
      DO 100 I=1,NOTRI
      IF(NCNE.EQ.8) THEN
      AG=2.0D0
         IF(NOTRI.EQ.3) THEN
            IF(I.EQ.1.OR.I.EQ.3) AG=1.0D0
         ENDIF
      ELSE
         AG=0.5D0
            IF(NOTRI.EQ.2) AG=0.25D0
      ENDIF
      DO 110 IXRHO=1,4
      DO 120 IYRHO=1,4
      NOSAMP=NOSAMP+1
!C
!!C**  calculate the shape function at the sampling points(XRHO,YRHO)
!C
      XRHO = XIQ(IXRHO)
      YRHO = XIQ(IYRHO)
!C
      DUM1=0.5D0*(1.0D0-XRHO)
      DUM2=0.25D0*(1.0D0+XRHO)*(1.0D0-YRHO)
      DUM3=0.25D0*(1.0D0+XRHO)*(1.0D0+YRHO)
!C
      SI =DUM1*SITRI (I,1)+DUM2*SITRI (I,2)+DUM3*SITRI (I,3)
      ETA=DUM1*ETATRI(I,1)+DUM2*ETATRI(I,2)+DUM3*ETATRI(I,3)
!C
!C
      IF(NCNE.EQ.8) THEN
    CALL SPFUNC8(SI,ETA,SF,DSF)
      ELSE
    CALL SPFUNC6(SI,ETA,SF,DSF)
      ENDIF
!C
!!C** evaluate the Jacobian matrix at (SI,ETA)
!C
        DO 130 LI=1,2
        DO 130 LJ=1,3
        DUM=0.0D0
        DO 140 LK=1,NCNE
        DUM=DUM+DSF(LI,LK)*XYZT(LJ,LK)
140     CONTINUE
130     XJ(LI,LJ)=DUM
!C
!!C** compute the determinant of the Jacobian maxtix at (XRHO,YRHO)
!C
        DET1=XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2) 
        DET2=XJ(1,1)*XJ(2,3)-XJ(1,3)*XJ(2,1) 
        DET3=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1) 
        DET=DSQRT(DET1*DET1+DET2*DET2+DET3*DET3)
!C
!C
        X=0.0D0
        Y=0.0D0
        Z=0.0D0
!C
        DX=0.0D0
        DY=0.0D0
        DZ=0.0D0
!C
        DO 150 LK=1, NCNE
        X=X+SF(LK)*XYZT(1,LK)
        Y=Y+SF(LK)*XYZT(2,LK)
        Z=Z+SF(LK)*XYZT(3,LK)
        DX=DX+SF(LK)*DXYZT(1,LK)
        DY=DY+SF(LK)*DXYZT(2,LK)
        DZ=DZ+SF(LK)*DXYZT(3,LK)
150     CONTINUE
!C
        XYNOD(1,NOSAMP)=X
        XYNOD(2,NOSAMP)=Y
        XYNOD(3,NOSAMP)=Z
        DXYNOD(1,NOSAMP)=DX
        DXYNOD(2,NOSAMP)=DY
        DXYNOD(3,NOSAMP)=DZ 
        DXYNOD(4,NOSAMP)=(Y-YC)*DZ-(Z-ZC)*DY
        DXYNOD(5,NOSAMP)=(Z-ZC)*DX-(X-XC)*DZ
        DXYNOD(6,NOSAMP)=(X-XC)*DY-(Y-YC)*DX
!C
        SAMNOD(NOSAMP,0)=0.25D0*AG*(1.0+XRHO)*WIQ(IXRHO)*WIQ(IYRHO)*DET
!C
        DO 160 J=1,NCNE
        SAMNOD(NOSAMP,J)=SF(J)*SAMNOD(NOSAMP,0)
160     CONTINUE
!C
120     CONTINUE
110     CONTINUE
100     CONTINUE
!C
!C
        RETURN
        END SUBROUTINE TRIPOL
end module
