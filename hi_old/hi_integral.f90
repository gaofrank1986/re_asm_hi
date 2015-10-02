module hi_SIEPPEM

    use hi_const
    use hi_shape_function
    use hi_target_func
    
    implicit none 
    
    integer,protected    ::  num_dim,num_node,num_nrml,num_elem
    integer,protected    ::  elem_nd_count,num_target_func
    real(8),protected    ::  hi_beta 


    real(8),allocatable,private ::  node_matrix(:,:),normal_matrix(:,:),src_local_list(:,:)
    !node_matrix(1:3,node_id),normal_matrix(1:3,nrml_id)
    integer,allocatable,private ::  elem_matrix(:,:),src_flag(:)
    !elem_matrix(1:8,elem_id)
    real(8),allocatable,private ::  full_mesh_matrix(:,:,:)
    !full_mesh_matrix(1:3,1:8,elem_id)

    
    integer,private :: model_readed_flag = 0 ! 0 for not readed
    integer,private :: external_src_ctr_flag = 0! if or not use external src ctr input

    integer,private,parameter :: NPW = 4
    real(8),private,allocatable :: value_list(:,:)
    integer,private :: n_pwr_g = -1
    real(8),allocatable,private :: cnr_glb_mtx(:,:) !corner_global_matrix
    real(8),private :: src_lcl_preset(2),src_glb_preset(3)
   
    real(8),private ::  src_glb(3),src_ctr_glb(3)
    !src coordinate in global,src center point in global
    
contains

    subroutine eval_singular_elem(this_elem_id,num_edge,result,src_preset_flag)
        !  removed input variable ndsid.......not used---July 25th
        !  add result variable initialisation-----July 26th
        !  Add src_preset function----------August 1st
        !  what about output shape,input src_ctr_glb 
        implicit none

        integer,intent(in) :: num_edge,this_elem_id,src_preset_flag
        real(8),intent(out) :: result(num_target_func)

        real(8) :: GPL(IABS(NGL)),GWL(IABS(NGL)) ! gaussian points and weights     

        real(8) :: src_lcl(num_dim-1),pt_intg(num_dim-1),pt_intg_tmp(num_dim - 1)  
        ! src point, integration point , temporary integration point     

        real(8) :: end_nodes(2,2),ri(3),SF_src(elem_nd_count),RINT(num_target_func)
        ! end nodes recorder, shape function
        ! rho integration result

        real(8) :: COEF_G(0:n_pwr_g),COEF_H(0:NPW)

        !====================================================================
        integer :: id,ie,tmp,i_edge,ks,src_identifier,igl,num_converge

        integer :: unfixed_cmp,fixed_cmp

        real(8) :: vlc2,rho_q,comt
        real(8) :: fk,wfa,sign_value,ri_tmp,full_step_size,sum,step_size,drdnp
        real(8) :: diff_1,diff_2

        integer :: debug_flag,debug_file_id
        
        debug_file_id = 109
        debug_flag = 0

 
        print *,"============in eval_SINGULAR_ELEM=============="

        result = 0.
 
        ! cnr_glb_mtx is allocated when read in
        cnr_glb_mtx = full_mesh_matrix(1:num_dim,1:elem_nd_count,this_elem_id)
       
        ! need set up src_lcl and src_glb
        if (src_preset_flag .eq. 0) then
            src_identifier = src_flag(this_elem_id)
            
            if (src_identifier < 0) then 
                ! use local src info if identifier less than zero
                src_lcl = src_local_list(:,this_elem_id) !!!!!====== attention how src_lcl_list is read
                call SHAPEF(num_dim,elem_nd_count,cnr_lcl_mtx,cnr_glb_mtx,src_lcl, &
                            & src_ctr_glb,ri,SF_src)
              
                ! shape function based src_lcl
                ! return ri and SF_src

                write(110,*)' source pt =',src_lcl
                write(110,*)' global pt =',src_glb
            else    
                ! src_identifier is 0 will not be called for eval_singular-----
                src_glb(1:num_dim)=cnr_glb_mtx(1:num_dim,src_identifier)
                ! attention when > 0, here use node index in a element, not global node id 
                ri = src_glb - src_ctr_glb   
                if (num_dim == 2) then
                    src_lcl(1)=(-1)**src_identifier+src_identifier/3
                else ! num_dim = 3
                    src_lcl = cnr_lcl_mtx( 2*src_identifier-1 : 2*src_identifier )
                end if
                write(110,*)' source pt =',src_lcl
                write(110,*)' global pt =',src_glb

            end if 
        else
            src_lcl = src_lcl_preset
            src_glb = src_glb_preset
            ri = src_glb - src_ctr_glb 
            write(110,*)' source pt =',src_lcl
            write(110,*)' global pt =',src_glb           
               
        end if
       
        


        !===================FINISH prepare===============================

        if (debug_flag .eq. 1) then
            write(110 ,*) "========output glb mtx in line 123,hi_integral,order switched=========="
            do id = 1,8
                write (110,*) 'node',id,'is',cnr_glb_mtx(1:3,id)
            end do
            write (debug_file_id,*) '================= corner_global matrix ==============='
            write (debug_file_id,*) cnr_glb_mtx
            write (debug_file_id,*) '================= src point global pos ==============='
            write (debug_file_id,*) src_glb
            write (debug_file_id,*) '================= src point local  pos,will changed if ... ==============='
            write (debug_file_id,*) src_lcl
            write (debug_file_id,*) '================= src center globl pos ==============='
            write (debug_file_id,*) src_ctr_glb
            write (debug_file_id,*) '================= ri ==============='
            write (debug_file_id,*) ri

        endif
 

        if (num_dim == 2) then

            !   Evaluate 2D SINGULAR INTEGRALS

!             DO ISUB=1,2 ! two division

!                 pt_intg(1)=(-1.D0)**ISUB

!                 IF(DABS(pt_intg(1)-src_lcl(1)).LT.TOL) then 
!                     print *,"skip this sub!!!!"
!                 end if

!                 !   SET UP COEFFICIENTS Gn ADN Hn
!                 CALL COEFS_GH(NDIM,NBDM,NODE,hi_beta,n_pwr_g,cnr_lcl_mtx,cnr_glb_mtx,XP,XIP,XIQ,        &
!                 &              NPW,COEFG,COEFH,NF)
!                 CALL INT_RHO(NDIM,NBDM,NODE,hi_beta,NF,n_pwr_g,NPW,cnr_lcl_mtx,cnr_glb_mtx,src_ctr_glb,         &
!                 &               XP,XIP,XIQ,SF_src,COEFG,COEFH,RINT)
!                 V1E=V1E+RINT   
!             end do 
      
        else 
            !-----------------------------------------------------------------------
            call GAUSSV(iabs(NGL),GPL,GWL)
!             if (debug_flag .eq. 1) then
!                 write (debug_file_id,*) '================= GPL ==============='
!                 write (debug_file_id,*) GPL
!                 write (debug_file_id,*) '================= GWL ==============='
!                 write (debug_file_id,*) GWL
!             endif    

            WFA=DSQRT(hi_beta*2.D0/3.D0+0.4D0)*DLOG(DABS(TOLGP)/2.D0)   
            FK=3.D0/8.D0*(-10.D0*DBLE(iabs(NGL))/WFA-1.D0)**(4.D0/3.D0)

!             if (debug_flag .eq. 1) then
!                 write (debug_file_id,*) '================= wfa ==============='
!                 write (debug_file_id,*) wfa
!                 write (debug_file_id,*) '================= fk ==============='
!                 write (debug_file_id,*) fk

!             endif
      
            do i_edge = 1,num_edge ! ITERATE through each edge

                KS=KSB(i_edge)
                IF(DABS(src_lcl(IABS(KS))-DBLE(KS)/DABS(DBLE(KS))).LT.TOL) then
                    print *,"Current edge iteration skipped! elem_id = ",this_elem_id," edge =",i_edge
                    if (debug_flag .eq. 1) then
                        write (debug_file_id,*) "Current edge iteration skipped! elem_id = ",this_elem_id," edge =",i_edge
                    endif
                    !------------------------------------------
                    goto 100
                end if
                
                do id = 1,2
                    tmp = node_grp_by_edge(3*(i_edge - 1) + ID)! determine which group of node to used
                    end_nodes(1:2,ID)=cnr_lcl_mtx(2*tmp-1 : 2*tmp) !get local node from corner table
                end do

                !==================================================
                ! for first and third edge, vertical component unchanged,1st comp change
                ! for second and forth edge, horizontal component unchanged,2rd comp change
                !       4----7----3
                !       |         |
                !       8    9    6
                !       |         |
                !       1----5----2

                !         7     6     5

                !         8           4

                !         1     2     3
               

                unfixed_cmp=1
                if (i_edge/2*2 .EQ. i_edge) unfixed_cmp=2 
                ! i_edge can is mutiple of 2
                fixed_cmp=3-unfixed_cmp
                !==================================================
                
                sign_value = DSIGN(1.D0,end_nodes(unfixed_cmp,2)-end_nodes(unfixed_cmp,1))
                !dsign(a,b) a time sign of b,end_node(:,id)
                ! sign_value is sign of second node minus first node,actually shows the direction

                VLc2=(end_nodes(fixed_cmp,2)-src_lcl(fixed_cmp))**2

!                  if (debug_flag .eq. 1) then
!                     write (debug_file_id,*) '======================================='
!                     write (debug_file_id,*) "current element",this_elem_id
!                     write (debug_file_id,*) "i_edge",i_edge
!                     write (debug_file_id,*) 'unfixed_cmp',unfixed_cmp
!                     write (debug_file_id,*) 'vlc2',vlc2
!                     write (debug_file_id,*) '========== start num_converge loop ============================='
!                 end if

                pt_intg_tmp=end_nodes(:,1) 
                !tmp integration point
                pt_intg(fixed_cmp)=pt_intg_tmp(fixed_cmp)
                ! the other xiq component is not initialised

                do num_converge=1,500 
                !control the maxium step to go thru one edge
                !also controled by step-size, example finished in less than 10 step
                !pt_intg(unfixed_cmp) is updated each time

                   
                    diff_1 = src_lcl(unfixed_cmp)-pt_intg_tmp(unfixed_cmp)!between tmp src and end node
                    diff_2 = end_nodes(unfixed_cmp,2)-pt_intg_tmp(unfixed_cmp)!between two nodes

                    if (sign_value*diff_2 < 1.D-8) then
                        print *, "if pt_intg_tmp coordinate out of range,exit num_converge loop"
                        !if pt_intg_tmp goes out of the edge,then stop
                        goto 100
                    end if
                    
                    ri_tmp = new_norm2(pt_intg_tmp-src_lcl) 
                    ! this is r from temporary integration point to local src point

                    if (0..GE.sign_value*diff_1) then
                        full_step_size = FK*ri_tmp ! fk is a factor?
                    else 
                        SUM=VLc2*(1.-FK**2)+(diff_1)**2
                        IF(SUM < 0.D0) THEN
                            full_step_size=sign_value*diff_1
                        ELSE
                            full_step_size=FK*(FK*sign_value*(diff_1)-DSQRT(SUM))/(FK*FK-1.)
                        ENDIF
                    endif

                
                    IF (sign_value*sign_value*full_step_size+1.D-8 > sign_value*diff_2)  full_step_size=sign_value*(diff_2)

                    step_size = 0.5D0*sign_value*full_step_size
                    !step size control,step size is half the full step size,also consider sign
                   
                    print *,"step_size = ",step_size
                    !==================================================================
                    !--- compute integral below

                    do IGL = 1,iabs(NGL) ! gaussian sampling points
                        ! this method change the double integral to line integral
                        
                        pt_intg(unfixed_cmp)=pt_intg_tmp(unfixed_cmp)+step_size*(1.D0+GPL(IGL))
                        ! update integration point position
                        
                        RHO_Q=new_norm2(pt_intg-src_lcl)
                        ! recaluate rho_q
                        
                        DRDNP=DABS(pt_intg(fixed_cmp)-src_lcl(fixed_cmp))/RHO_Q !sin(theta)

                        call compute_coeff_GH(num_dim,num_dim - 1,elem_nd_count,n_pwr_g,src_glb &
                                            & ,src_lcl,pt_intg,COEF_G,COEF_H)
                        !             get coef_g and coef_h
                        !             why not inside integrate_Rho?????

!                         if (debug_flag .eq. 1) then
!                             write (debug_file_id,*) "================coef_g=========="
!                              write (debug_file_id,*) coef_g
!                             write (debug_file_id,*) '================coef_h==========='
!                             write (debug_file_id,*) coef_h
!                         end if

                        call integrate_rho(n_pwr_g,src_lcl,pt_intg,coef_g,coef_h,RINT)

                        result = result + (DABS(step_size)*GWL(IGL)*DRDNP/RHO_Q)*RINT
                        ! Equation (3-6-50)

                    end do ! igl =1,iabs(ngl)
                    
                    pt_intg_tmp(unfixed_cmp)=pt_intg_tmp(unfixed_cmp)+sign_value*full_step_size
!                       if (debug_flag .eq. 1) then
!                             write (debug_file_id,*) "================result=========="
!                              write (debug_file_id,*) result
!                             write (debug_file_id,*) '================pt_intg_tmp==========='
!                             write (debug_file_id,*) pt_intg_tmp
!                         end if
                end do ! num_converge = 1,500

                if (debug_flag .eq. 1) then
                            write (debug_file_id,*) "================end num_converge loop =========="
                end if

100          end do ! i_edge = 1,num_edge
        end if



        
!             write (110,*) "==========final result========="
            call swap_result(result)
!             write (110,*) result
!         end if


    end subroutine



    SUBROUTINE integrate_RHO(n_pwr_g,src_lcl,pt_intg,coef_g,coef_h,result)      
      ! changed cnr_glb_mtx to private variable shared in module
        implicit none 

        real(8),intent(in)  :: src_lcl(num_dim-1),pt_intg(num_dim-1)
        integer,intent(in)  :: n_pwr_g
        real(8),intent(out) :: result(num_target_func)

        real(8)  :: COSN(num_dim),RI(num_dim),GCD(num_dim,num_dim-1)
        real(8)  :: COEF_G(0:n_pwr_g),COEF_H(0:NPW),COEF_B(0:11,num_target_func)

        integer :: k ,beta1,nbdm,npowf

        real(8) :: SLOP(num_dim-1),rho_q,E_k,pw,nbeta

        integer :: n_pwr_k

        !=============================================

        result=0.D0 ! output initialization 
        SLOP = pt_intg - src_lcl
        rho_q = new_norm2(SLOP)
        SLOP = SLOP/rho_q ! normalized vector, cos(theta),sin(theta)

        !!! -------------Compute n_pwr_k
        n_pwr_k = int(3+2.1214*RHO_Q)    ! NPOWF IS FROM 3 TO 9
        !order of power expansion, for parameter K in equation (3-6-62)
        if (n_pwr_k.LT.(hi_beta-num_dim+1)) then
            n_pwr_k = hi_beta-num_dim+1
        end if       

        !!! - ----------End computing pwr_k

        call compute_coeff_B(elem_nd_count,n_pwr_g,n_pwr_k, &
                & src_glb,src_lcl,pt_intg,SLOP,RHO_Q,RI,COSN,GCD,COEF_G,COEF_B)     

        ! Case 1 for Ek, 0 <= k <= lamda - 3        

        do k =0 ,int(hi_beta)-num_dim

            PW= hi_beta - k - (num_dim - 1) ! the power coefficient of rho_q
            E_k = (1.D0/RHO_Q**PW-COEF_H(INT(PW)))/(-PW)
            result = result + E_k*COEF_B(K,:)
!              print *,"case 1"
             !print *,"E_k = ",E_k

        end do

        ! Case 2 for Ek, k = lamda - 2
        k = hi_beta - (num_dim - 1)
        IF (k.GE.0) then
            E_k = (DLOG(rho_q) - DLOG(COEF_H(0)))
            result = result + COEF_B(k,:)*E_k
        end if 

        ! Case 3 for Ek,  lamda - 2 <= k <= n_pwr_k
        do k = int(hi_beta-(num_dim-1)+1),n_pwr_k

            PW= hi_beta - k - (num_dim - 1) ! the power coefficient of rho_q

            !PW = k - hi_beta +(num_dim - 1) ! Equ (3-6-64) Ek, power k+2 - lamda
            E_k = rho_q**(-PW)/(-PW)
            result = result + E_k*COEF_B(k,:)
        end do


    END SUBROUTINE 



    subroutine compute_coeff_B(NODE,NPOWG,NPOWF,XP,XIP,XIQ,SLOP,RHOQ,RI,COSN,GCD,COEFG,COEFB)

    !cnr_glb_mtx : corner global matrix changed to shared private variable in module
    ! call c_COEF_B(num_dim,num_dim-1,elem_nd_count,beta,n_pwr_g,n_pwr_k,cnr_glb_mtx,      &
    !    &                  src_glb,src_lcl,pt_intg,SLOP,RHO_Q,RI,COSN,GCD,COEF_G,COEF_B)
      
        implicit none

        integer ::  NODE,NPOWG,NPOWF
       
      
        real(8),intent(in) :: XP(num_dim),XIP(num_dim - 1),XIQ(num_dim - 1),SLOP(num_dim - 1), &
                    & COSN(num_dim)


        real(8) :: DRDX(num_dim),XI(num_dim - 1),X(num_dim),RMAT(NPOWF,NPOWF), &
                    &   COEFG(0:NPOWG),COEFB(0:11,num_target_func),SF_iter(NODE),&
                    &   FQ(num_target_func),A(num_dim), GCD(num_dim,*),RI(num_dim)
      
        integer :: i,j,ip,m,jp
        integer :: NDIM,NBDM
        real(8) :: vk,robar
        real(8) :: rhoq,rho,r2,r,fjcb,gm,drdn

        NDIM = num_dim
        NBDM = num_dim - 1


      
        
        VK=RHOQ/DBLE(NPOWF) ! divide rho_q to npowf parts

        DO 20 IP=0,NPOWF 

            RHO=VK*DBLE(IP) 
            XI=XIP+RHO*SLOP !!!!!!!!!!!!!!!!! xi updated here!!!!!!!!!!!!!!!!!!
            ROBAR=COEFG(0)

            DO M=1,NPOWG; ROBAR=ROBAR+COEFG(M)*RHO**M; ENDDO

            ROBAR=DSQRT(ROBAR)       ! Eq.(3-6-28)

            CALL SHAPEF(NDIM,NODE,cnr_lcl_mtx,cnr_glb_mtx,XI,XP,RI,SF_iter)
            !shape func based on xi
            !ri give distance vector from xp(it is src_glb) to xi (update along rho)
            R = new_norm3(RI)

            X=XP+RI
          
            CALL DSHAPE(NDIM,NODE,cnr_lcl_mtx,cnr_glb_mtx,XI,COSN,FJCB,GCD)
            !SUBROUTINE DSHAPE(NDIM,NODE,C,CK,X,COSN,FJCB,GD)
            ! GCD gives the normal n vector
            ! dshape based on xi

            IF(RHO.GT.1.0D-10)THEN    
                DRDX=RI/R
            ELSE
                A = 0.D0
    !             DO 10 I=1,NDIM; A(I)=0.D0
    !             DO 10 J=1,NBDM        
                forall (i = 1:ndim)    
                    A(i)=A(i)+dot_product(GCD(i,1:nbdm),SLOP(1:nbdm))
                end forall                  
            
                GM=DSQRT(DOT_PRODUCT(A,A))
                DRDX=A/GM                 ! Eq.(3-6-74)
            endif 

            ! dr/dx is defined above
            
            DRDN = DOT_PRODUCT(COSN,DRDX)  !!!!  notice dr/dn is defined  here
            
            CALL F_BAR(NDIM,NBDM,DRDX,COSN,R,DRDN,XI,SF_iter,XP,X,num_target_func,FQ)
            !!SUBROUTINE F_BAR(NDIM,NBDM,DRDX,COSN,R,DRDN,XI,SHAP,XP,XQ,NF,FB)

            COEFB(IP,:) = FQ*FJCB/ROBAR**hi_beta

     
            IF(IP.EQ.0) GOTO 20
            COEFB(IP,:)=(COEFB(IP,:)-COEFB(0,:))/RHO

            RMAT(IP,1)=1.D0
            DO JP=2,NPOWF; RMAT(IP,JP)=RMAT(IP,JP-1)*RHO; ENDDO

20      CONTINUE
        


        CALL INVSOLVR(NPOWF,NPOWF,RMAT,NPOWF,-1)
   

        forall (i = 1:npowf,j=1:num_target_func)

        coefb(i,j) = DOT_PRODUCT(rmat(i,:),coefb(1:npowf,j))

        end forall
        !COEFB(1:NPOWF,:)=MATMUL(RMAT,COEFB(1:NPOWF,:))    ! Bk IN Eq.(3-6-62)
        ! coefb has size 12xnum_target_func,but only from index 1 to 7 has useful info

    end subroutine compute_coeff_B



    subroutine compute_coeff_GH(NDIM,NBDM,NODE,NPOWG,XP,XIP,XIQ,  &
     &                    COEFG,COEFH)
        ! set hi_beta to shared variable
        implicit none

        integer ::  NDIM,nBDM,NODE,NPOWG
        real(8) :: drdn

        real(8) ::XP(NDIM),XIP(NBDM),XIQ(NBDM),COSN(NDIM),        &
        &          GCD(NDIM,NBDM),XI(NDIM),RI(NDIM),SHAP(NODE),     &
        &          COEFG(0:NPOWG),COEFC(0:NPW),COEFH(0:NPW)

        integer :: i,n,ip,jp
       
        ! DETERMINE COEFFICIENTS Gn
       
        COEFG=0.D0
       
        CALL compute_coeff_G(NDIM,NBDM,NODE,NPOWG,XP,XIP,XIQ, DRDN,COSN,    &
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
   
    end subroutine

     


    SUBROUTINE compute_coeff_G(NDIM,NBDM,NODE,NPOWG,XP,XIP,XIQ,         &
     &                  DRDN,COSN,GCD,RI,SHAP,COEFG)
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
 10    COEFG(IP)=((new_norm3(RI)/RHO)**2-COEFG(0))/RHO     ! Eq.(3-6-60)

       RMAT(IP,1)=1.D0        ! Eq.(3-6-42)
       DO JP=2,NPOWG; RMAT(IP,JP)=RMAT(IP,JP-1)*RHO; ENDDO
 20    CONTINUE
       
       CALL INVSOLVR(NPOWG,NPOWG,RMAT,NPOWG,-1)       ! INVERSE MATRIX [R] IN Eq.(3-6-59)
       COEFG(1:NPOWG)=MATMUL(RMAT,COEFG(1:NPOWG))     ! SOLVING Eq.(3-6-59) FOR {G}
        ENDIF
    end subroutine compute_coeff_g



    subroutine read_model_from_WAVDUT()
        
        use MVAR_MOD

        implicit none

        integer :: i,j,ie,id,tmp,tmp1

        print *,"------------Start Reading Model from WAVDUT-------------"

        num_dim = 3
        num_node = NNODE
        num_elem = NELEM
        num_nrml = NNODED
        elem_nd_count = 8 !NCN(IELEM)!! to be changed
        hi_beta = 3.
        num_target_func = 8

        allocate(node_matrix(num_dim,num_node))
        allocate(normal_matrix(num_dim,num_nrml))

        allocate(elem_matrix(elem_nd_count,num_elem))
        allocate(value_list(num_elem,num_target_func))
        allocate(cnr_glb_mtx(num_dim,elem_nd_count))

        if (num_dim == 2) ngl = 1

        forall (i = 1:num_dim,j=1:num_node)
            node_matrix(i,j) = XYZ(i,j)
        end forall

        forall (i = 1:num_dim,j=1:num_nrml)
            normal_matrix(i,j) = DXYZ(i,j)
        end forall

        forall(i = 1:num_elem,j = 1:elem_nd_count)
            elem_matrix(j,i) = NCON(i,j)!triangle
            ! there is index order problem
!              elem_mtx_nrml(i,j) = NCOND(i,j)
        end forall

        ! for first and third edge, vertical component unchanged,1st comp change
                ! for second and forth edge, horizontal component unchanged,2rd comp change
                !       4----7----3
                !       |         |
                !       8    9    6
                !       |         |  new position
                !       1----5----2

                !         7     6     5

                !         8           4 old position

                !         1     2     3
        !switch order
        print *,"Transfering elem matrix from old order to new order"
            do i = 1,num_elem
                if (elem_nd_count .eq. 8) then
!                     tmp = elem_matrix(2,i)
!                     elem_matrix(2,i) = elem_matrix(5,i)
!                     tmp1 = elem_matrix(3,i)
!                     elem_matrix(3,i) = tmp
!                     tmp = elem_matrix(4,i)
!                     elem_matrix(4,i) = elem_matrix(6,i)
!                     elem_matrix(5,i) = tmp1
!                     elem_matrix(6,i) = elem_matrix(7,i)
!                     elem_matrix(7,i) = tmp



                endif
                if (NCN(i) .eq. 8) then
                tmp = elem_matrix(2,i)
                elem_matrix(2,i) = elem_matrix(3,i)
                elem_matrix(3,i) = elem_matrix(5,i)
                elem_matrix(5,i) = tmp
                tmp = elem_matrix(4,i)
                elem_matrix(4,i) = elem_matrix(7,i)
                elem_matrix(7,i) = elem_matrix(6,i)
                elem_matrix(6,i) = tmp
            endif
            end do
        !-------------------- Data Manipulation------------------

        allocate(full_mesh_matrix(num_dim,elem_nd_count,num_elem))
            
        forall (ie = 1:num_elem,id = 1:elem_nd_count)

                full_mesh_matrix(1:num_dim,id,ie)=node_matrix(1:num_dim,elem_matrix(id,ie))
                ! reorganize nodes coordinate by element node order
        end forall

        !------------- Initialization-------------------------
        model_readed_flag = 0
        value_list = 0

        print *,"------------Stop Reading Model from WAVDUT-------------"
        !print *, node_matrix
!         print *,"------------elem matrix-------------"

        !print *, elem_matrix



    end subroutine read_model_from_WAVDUT

    subroutine swap_result(result)
        implicit none
        real(8) :: result(*),tmp,tmp1

                !       4----7----3
                !       |         |
                !       8    9    6
                !       |         |
                !       1----5----2

                !         7     6     5

                !         8           4

                !         1     2     3

        tmp = result(2)
        result(2) = result(5)
        tmp1 = result(3)
        result(3) = tmp
        tmp = result(4)
        result(4) = result(6)
        result(5) = tmp1
        result(6) = result(7)
        result(7) = tmp
    end subroutine


    subroutine read_model_from_DAT()

        implicit none
        integer :: ip,ie,tmp,i,id        

        if (model_readed_flag == 0) then
            print *,"------------Start Reading Model-------------"
            OPEN(5,FILE='SIEPPEM.DAT',STATUS='OLD')

            read (5,*) num_dim,num_node,num_elem,elem_nd_count,hi_beta,num_target_func
            ! dimesnion
            ! number of node
            ! number of element
            ! number of node per element
            ! beta is the power of r in target equation
            ! number of target func components

            allocate(node_matrix(num_dim,num_node))
            allocate(elem_matrix(elem_nd_count,num_elem))
            allocate(src_flag(num_elem))
            allocate(src_local_list(2,num_elem))
            allocate(value_list(num_elem,num_target_func))

            if (num_dim == 2) ngl = 1
    
         !    Input nodal coordinates and element connectivity
          
            DO IP = 1,num_node
                READ(5,*) tmp,(node_matrix(I,tmp),I=1,num_dim)                 ! Card set 2
            end do  


            DO IE = 1,num_elem
                READ(5,*) tmp,(elem_matrix(ID,tmp),ID=1,elem_nd_count),src_flag(tmp)    ! Card set 3
            end do
            !==========================
            !====src_flag
            ! if = 0 src not on elem
            ! if > 0 src is given in global coordinate, use node with id (src_flag)
            ! if < 0 src is given in local coordinate, use local src list given in card set 4

          
            READ(5,*) (src_glb(i),i=1,num_dim)                       ! Card set 4  
             ! read src x,y,z coord
             ! there seems a error, src_glb should be an array of global coordinate
             ! since src_glb cannot remain unchanged for different element
        
            DO IE=1,num_elem
                if (src_flag(ie) < 0) then
                    read(5,*) (src_local_list(i,ie),i=1,num_dim-1)
                end if
                ! src local position given in elements input order
            end do

            close(5)

            print *,"------------Finish Reading Model-------------"
            !------------Some initialisation of data

            allocate(full_mesh_matrix(num_dim,elem_nd_count,num_elem))
            
            forall (ie = 1:num_elem,id = 1:elem_nd_count)

                    full_mesh_matrix(1:num_dim,id,ie)=node_matrix(1:num_dim,elem_matrix(id,ie))
                    ! reorganize nodes coordinate by element node order
            end forall

            model_readed_flag = 0
            value_list = 0

        else
            print *,"--Attention! Reading process skipped,model already loaded---"
        end if

    end subroutine read_model_from_DAT

     subroutine set_npwg(number)
        
        implicit none

        integer,intent(in) :: number

        n_pwr_g = number
        print *,"n_power_g is set to ",number

    end subroutine

    subroutine set_src_preset(ksi,eta,glb,ctr_glb)
        
        implicit none

        real(8),intent(in) :: ksi,eta,glb(3),ctr_glb(3)

        src_lcl_preset(1) = ksi
        src_lcl_preset(2) = eta
        src_glb_preset = glb
        src_ctr_glb = ctr_glb
        print *,"src preseted as",src_lcl_preset
        print *,"src preseted as",src_glb_preset
        print *,"src_ctr_glb preseted as",src_ctr_glb


    end subroutine

    subroutine debug_test()
        
        implicit none

        print *,"elem matrix",elem_matrix

    end subroutine
    


   
!     subroutine RIM_ELEMS()

!         implicit none 

!         integer  :: num_edge,ie,id

! !         EXTERNAL INT_ELEM

!         IF(num_dim == 2) n_pwr_g = (elem_nd_count/3)*2               ! 0,  2
!         IF(num_dim == 3) n_pwr_g = elem_nd_count/2+(elem_nd_count/9)*2        ! 2,  4,  6 
!         ! refer to  Equ. 3-6-56 for parameter m
     
!         num_edge = 2 * (num_dim - 1 ) ! 4 -----how many edges

!         ! NDSID=2+elem_nd_count/8 !3 !removed not used=== July 25th====

!         src_ctr_glb = 0. !src center global, define the center of src for calculating r

!         do ie = 1,num_elem ! can introduce parallel here!!!!!!!!!!!!!!!!!
!             If (src_flag(ie) == 0) THEN    
!                 ! EVALUATE INTEGRAL OVER REGULAR ELEMENT
!                 !CALL ADAPTINT_ELEM(ie,src_ctr_glb,cnr_lcl_mtx,value_list(ie),GPR,GWR,GPL,GWL,INT_ELEM)
!                 print *," need evaluate integral over element,src not on element"            
!             else     
!                 CALL eval_SINGULAR_ELEM(ie,num_edge,value_list(ie,:),0)
!             end if 
!         end do

!         print *,"==========final result========="
!         print *,sum(value_list)
!         print *,"==============================="
!         print *,"end of RIM_ELEMS"
!     end subroutine

!     subroutine initialise_ELEMS(x0,y0,z0)

!         implicit none 

!         real,intent(in) :: x0,y0,z0
!         integer  :: num_edge,ie,id

!         !         EXTERNAL INT_ELEM

!         !-------------------------------------------

!         IF(num_dim == 2) n_pwr_g = (elem_nd_count/3)*2               ! 0,  2
!         IF(num_dim == 3) n_pwr_g = elem_nd_count/2+(elem_nd_count/9)*2        ! 2,  4,  6 
!         ! refer to  Equ. 3-6-56 for parameter m
     
!         num_edge = 2 * (num_dim - 1 ) ! 4 -----how many edges

!         external_src_ctr_flag = 1
        
!         if (external_src_ctr_flag .eq. 0) then
!             src_ctr_glb = 0. !src center global, define the center of src for calculating r
!         else
!             src_ctr_glb(1) = x0
!             src_ctr_glb(2) = y0
!             src_ctr_glb(3) = z0
!         end if

!         !         do ie = 1,num_elem ! can introduce parallel here!!!!!!!!!!!!!!!!!
!         !             If (src_flag(ie) == 0) THEN    
!         !                 ! EVALUATE INTEGRAL OVER REGULAR ELEMENT
!         !                 !CALL ADAPTINT_ELEM(ie,src_ctr_glb,cnr_lcl_mtx,value_list(ie),GPR,GWR,GPL,GWL,INT_ELEM)
!         !                 print *," need evaluate integral over element,src not on element"            
!         !             else     
!         !                 CALL eval_SINGULAR_ELEM(ie,num_edge,value_list(ie,:))
!         !             end if 
!         !         end do

!         !         print *,"==========final result========="
!         !         print *,sum(value_list)
!         !         print *,"==============================="
!         !         print *,"end of RIM_ELEMS"
!     end subroutine
end module
      


