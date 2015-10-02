    subroutine comp_coef_GH(npowg,npw,COEFG,COEFH)
        ! set hi_beta to shared variable
        implicit none

        integer,intent(in) ::  NPOWG,npw

        real(8),intent(in) :: COEFG(0:NPOWG)
        real(8),intent(out) :: COEFH(0:NPW)
        real(8) :: COEFC(0:NPW)
        integer :: i,n,ip,jp
   
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
    end subroutine
    
    SUBROUTINE compute_coeff_G(NDIM,NBDM,NODE,NPOWG,XP,XIP,XIQ,         &
     &                  COEFG)
      ! THIS ROUTINE DETERMINES COEFFICIENTS Gm USING Eqs.(3-6-58)~(3-6-60)
        implicit none 


        integer ::  NDIM,nBDM,NODE,NPOWG
        real(8) :: drdn,xbar2,xxbar,vk,rhoq,rho,r2,fjcb


        real(8) ::XP(NDIM),XIP(NBDM),XIQ(NBDM),COSN(NDIM),        &
            &          GCD(NDIM,NBDM),XI(NDIM),SLOP(NBDM),RI(NDIM),     &
            &          SHAP(NODE),RMAT(NPOWG,NPOWG),COEFG(0:NPOWG)

        integer :: ip,jp


      IF(NODE.EQ.2) THEN         ! Eq.(24)
!        COEFG(0)=0.25*((cnr_glb_mtx(1,1)-cnr_glb_mtx(1,2))**2+(cnr_glb_mtx(2,1)-cnr_glb_mtx(2,2))**2)  
!       ELSEIF(NODE.EQ.3) THEN             ! FORM Eq.(3-6-33)
!        RI=cnr_glb_mtx(:,1)+cnr_glb_mtx(:,2)-2.*cnr_glb_mtx(:,3)     ! Xbar
!        XI=cnr_glb_mtx(:,1)-cnr_glb_mtx(:,2)                ! X1i-X2i 
!        XBAR2=RI(1)*RI(1)+RI(2)*RI(2)     ! Xbar^2
!        XXBAR=XI(1)*RI(1)+XI(2)*RI(2)
!        COEFG(0)=0.25*(XI(1)*XI(1)+XI(2)*XI(2))-XXBAR*XIP(1)             &
!      &         +XBAR2*XIP(1)*XIP(1)
!        COEFG(1)=(-0.5*XXBAR+XBAR2*XIP(1))*XIQ(1)
!        COEFG(2)=0.25*XBAR2
        ELSE      ! QUADRILATERAL ELEMENTS
            SLOP=XIQ-XIP
            RHOQ=DSQRT(DOT_PRODUCT(SLOP,SLOP))
            SLOP=SLOP/RHOQ
            VK=RHOQ/DBLE(NPOWG)
            DO 20 IP=0,NPOWG
            RHO=VK*DBLE(IP)
            XI(1:NBDM)=XIP+RHO*SLOP
            CALL SHAPEF(NDIM,NODE,cnr_lcl_mtx,cnr_glb_mtx,XI,XP,RI,SHAP)
            IF(IP.NE.0) GOTO 10
            CALL DSHAPE(NDIM,NODE,cnr_lcl_mtx,cnr_glb_mtx,XI,COSN,FJCB,GCD)   !dX/dXI
            RI=MATMUL(GCD,SLOP); COEFG(0)=DOT_PRODUCT(RI,RI)  ! Eq.(3-6-58)
            GOTO 20
!  10    COEFG(IP)=((DSQRT(R2)/RHO)**2-COEFG(0))/RHO     ! Eq.(3-6-60)
 10    COEFG(IP)=((norm2(RI)/RHO)**2-COEFG(0))/RHO     ! Eq.(3-6-60)

       RMAT(IP,1)=1.D0        ! Eq.(3-6-42)
       DO JP=2,NPOWG; RMAT(IP,JP)=RMAT(IP,JP-1)*RHO; ENDDO
 20    CONTINUE
       
       CALL INVSOLVR(NPOWG,NPOWG,RMAT,NPOWG,-1)       ! INVERSE MATRIX [R] IN Eq.(3-6-59)
       COEFG(1:NPOWG)=MATMUL(RMAT,COEFG(1:NPOWG))     ! SOLVING Eq.(3-6-59) FOR {G}
        ENDIF
    end subroutine compute_coeff_g
