module hi_shape_function

contains

      function new_norm2(r)
        real(8) :: r(2),new_norm2

        new_norm2 = dsqrt(r(1)*r(1)+r(2)*r(2))

        return

      end function

      function new_norm3(r)
        real(8) :: r(3),new_norm3

        new_norm3 = dsqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))

        return

      end function

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
      
      end SUBROUTINE GAUSSV

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

SUBROUTINE SHAPEF(num_dim,elem_nd_count,C,cnr_glb_mtx,Q_lcl,P_glb,RI,SP)
    ! P is start point
    ! Q is end point
    ! C is corner_local_matrix

    ! given Q_local, P_global,and information about the element
    ! like elem_nd_count,dimension and corner node global and local position

    ! we output RI,SP

    ! RI is the vector from P to Q in global coordinate, R2 is the square of distance
    ! SP is the shape func based on Q_local
    
    implicit none

    integer,intent(in) :: num_dim,elem_nd_count
    real(8),intent(in) :: cnr_glb_mtx(3,*),P_glb(*),C(*),Q_lcl(*)

    real(8),intent(out) :: SP(elem_nd_count),RI(*)
    integer :: i,l
    real(8) :: wl

    IF(elem_nd_count.GT.3) GOTO 4
!                  2-noded line element
    SP(1)=0.5*(1.-Q_lcl(1)); SP(2)=0.5*(1.+Q_lcl(1))
    IF(elem_nd_count.EQ.2) GOTO 50
!                  3-noded line element
    SP(1)=-Q_lcl(1)*SP(1); SP(2)=Q_lcl(1)*SP(2); SP(3)=1.-Q_lcl(1)*Q_lcl(1)
    GOTO 50
!
!                  4-noded quadrilateral element
4   DO I=1,4
     SP(I)=0.25*(1.+C(2*I-1)*Q_lcl(1))*(1.+C(2*I)*Q_lcl(2))
    ENDDO
    IF(elem_nd_count.EQ.8) THEN
!                  8 noded-element (square element)
         DO 15 I=1,4; L=2*I-1; SP(I)=SP(I)*(C(L)*Q_lcl(1)+C(L+1)*Q_lcl(2)-1.D0)
         WL=C(L+8)*Q_lcl(1)+C(L+9)*Q_lcl(2)
15       SP(I+4)=.5D0*(WL+1.D0)*(1.D0-(C(L+8)*Q_lcl(2))**2-(C(L+9)*Q_lcl(1))**2)
    ELSEIF(elem_nd_count.EQ.9) THEN
!                  9 noded-element (square element)
     DO 20 I=1,4; L=2*I-1; SP(I)=SP(I)*C(L)*Q_lcl(1)*C(L+1)*Q_lcl(2)
     WL=C(L+8)*Q_lcl(1)+C(L+9)*Q_lcl(2)
20   SP(I+4)=0.5D0*WL*(WL+1.D0)*(1.D0-(C(L+8)*Q_lcl(2))**2-               &
   &                  (C(L+9)*Q_lcl(1))**2)
     SP(9)=(1.D0-Q_lcl(1)*Q_lcl(1))*(1.D0-Q_lcl(2)*Q_lcl(2))
    ENDIF

    !=====================================================================================
!                  Calculate r and its vector components

50  DO I=1,num_dim
        RI(I)=-P_glb(I)+dot_product(SP(1:elem_nd_count),cnr_glb_mtx(I,1:elem_nd_count))
    end do
    ! get Q_global using Q_local, then calculate Ri is vector from P_global to Q_global

    END SUBROUTINE SHAPEF


    !=============================================================================

SUBROUTINE DSHAPE(NDIM,NODE,C,CK,X,COSN,FJCB,GD) 

    ! C is cnr_lcl_mtx,CK is cnr_glb_mtx
    ! X is local node
    ! FJCB is ...........
    ! GD is ..........
    
    IMPLICIT REAL*8 (A-H,O-Z)
    DIMENSION X(*),CK(3,*),DN(2,NODE),GD(3,*),COSN(*),C(*),GR(3)

    NBDM = NDIM - 1

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
    !========================================================

30  DO 50 I=1,NDIM; DO 50 J=1,NBDM; GD(I,J)=0.; DO 50 ID=1,NODE
50    GD(I,J)=GD(I,J)+DN(J,ID)*CK(I,ID) ! partial shape over partial ksi/eta

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

    END SUBROUTINE DSHAPE

end module hi_shape_function