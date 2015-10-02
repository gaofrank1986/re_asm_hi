
     
    subroutine compute_coeff_gh(ndim,nbdm,npw,node,npowg,xp,xip,xiq,  &
     &                    coefg,coefh)
        implicit none

        integer,intent(in) ::  ndim,nbdm,node,npw,npowg
        real(8),intent(in) :: xp(ndim),xip(nbdm),xiq(nbdm)
        real(8),intent(out) :: COEFG(0:NPOWG),COEFH(0:NPW)
        real(8) :: drdn

        real(8) ::COSN(NDIM),GCD(NDIM,NBDM),XI(NDIM),RI(NDIM),SHAP(NODE),     &
        &         COEFC(0:NPW)

        integer :: i,n,ip,jp
       
        ! DETERMINE COEFFICIENTS Gn
       
        !COEFG=0.D0
       
        call compute_coeff_g(ndim,nbdm,node,npowg,xp,xip,xiq,    &
        &           coefg)
        
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
   
    end subroutine
    subroutine compute_coeff_g(ndim,nbdm,node,npowg,xp,xip,xiq,         &
     &                  coefg)
      ! THIS ROUTINE DETERMINES COEFFICIENTS Gm USING Eqs.(3-6-58)~(3-6-60)
        implicit none 

        integer,intent(in) ::  ndim,nbdm,node,npowg
        real(8),intent(in) :: xp(ndim),xip(nbdm),xiq(nbdm)
        real(8) :: drdn,xbar2,xxbar,vk,rhoq,rho,r2,fjcb

        real(8) ::     cosn(ndim),        &
            &          gcd(ndim,nbdm),xi(ndim),slop(nbdm),ri(ndim),     &
            &          shap(node),rmat(npowg,npowg),coefg(0:npowg)

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
            slop=xiq-xip
            !rhoq=dsqrt(dot_product(slop,slop))
            rhoq = norm2(slop)
            slop=slop/rhoq
            vk=rhoq/dble(npowg)
            do 20 ip=0,npowg
            rho=vk*dble(ip)
            xi(1:nbdm)=xip+rho*slop
            call shapef(ndim,node,cnr_lcl_mtx,cnr_glb_mtx,xi,xp,ri,shap)
            if(ip.ne.0) goto 10
            call dshape(ndim,node,cnr_lcl_mtx,cnr_glb_mtx,xi,cosn,fjcb,gcd)   !dx/dxi
            ri=matmul(gcd,slop); coefg(0)=dot_product(ri,ri)  ! eq.(3-6-58)
            goto 20
!  10    coefg(ip)=((dsqrt(r2)/rho)**2-coefg(0))/rho     ! eq.(3-6-60)

 10    coefg(ip)=((norm2(ri)/rho)**2-coefg(0))/rho     ! eq.(3-6-60)

       rmat(ip,1)=1.d0        ! eq.(3-6-42)
       do jp=2,npowg; rmat(ip,jp)=rmat(ip,jp-1)*rho; enddo
 20    continue
       
       call invsolvr(npowg,npowg,rmat,npowg,-1)       ! inverse matrix [r] in eq.(3-6-59)
       coefg(1:npowg)=matmul(rmat,coefg(1:npowg))     ! solving eq.(3-6-59) for {g}
        endif
    end subroutine compute_coeff_g
