!module hi_integrand

!contains
    !subroutine f_bar(ndim,nbdm,drdx,cosn,r,drdn,xi,shap,xp,xq,nf,fb)
        !nidm,nbdm,nf : dimension,dimension for mapped space,number integrand
        !               component 
        !drdx,cosn,r,dr,dn
        ! xi          : src local position
        ! xp          : src point glb position
        ! xq          : integration point glb position
        !implicit none
        !integer,intent(in) :: ndim,nbdm,nf
        !real(8),intent(in) :: xp(ndim),xi(nbdm),drdx(ndim),cosn(ndim),shap(*),xq(ndim)
        
    subroutine f_integrand(ndim,nf,cosn,drdx,drdn,shap,fb)
        implicit none
        integer,intent(in) :: ndim,nf
        real(8),intent(in) :: drdx(ndim),drdn,cosn(ndim),shap(*)
        real(8),intent(out) :: fb(nf)       

        real(8) :: x(ndim),tmp
        integer :: id
      
        tmp =(3.*DRDX(3)*DRDN-COSN(3))/(4.*hi_PI)    ! GUIG 4.2
        !attention neg sign removed in tmp to match teng's result
        do id =1,8
            FB(id) = tmp*shap(id)
        end do
    end subroutine

!end module 
