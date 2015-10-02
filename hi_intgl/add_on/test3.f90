        subroutine compute_coeff_B(NDIM,NF,NODE,NPOWG,NPOWF,XP,XIP,&
                                & XIQ,SLOP,RHOQ &
                                &,COEFG,COEFB)
      
        implicit none

        integer,intent(in) ::  NODE,NPOWG,NPOWF,NDIM,NF
       
      
        real(8),intent(in) :: XP(NDIM),XIP(NDIM - 1)&
                    & ,XIQ(NDIM - 1),SLOP(NDIM - 1), &
                                        & COEFG(0:NPOWG)
        real(8):: ri(NDIM), COSN(NDIM), GCD(3,node)

        real(8),intent(in) :: rhoq
        real(8),intent(out) :: COEFB(0:11,NF)

        real(8) :: DRDX(NDIM),XI(NDIM - 1),X(NDIM)&
                    &,RMAT(NPOWF,NPOWF), &
                    &   SF_iter(NODE),&
                    &   FQ(NF),A(NDIM)      
        integer :: i,j,ip,m,jp
        integer :: NBDM
        real(8) :: vk,fjcb,robar
        real(8) :: rho,r2,r,gm,drdn

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
            R = norm2(RI)

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
            
            !CALL F_BAR(NDIM,NBDM,DRDX,COSN,R,DRDN,XI,SF_iter,XP,X,NF,FQ)
            call f_integrand(ndim,nf,cosn,drdx,drdn,sf_iter,fq)
            COEFB(IP,:) = FQ*FJCB/ROBAR**hi_beta

     
            IF(IP.EQ.0) GOTO 20
            COEFB(IP,:)=(COEFB(IP,:)-COEFB(0,:))/RHO

            RMAT(IP,1)=1.D0
            DO JP=2,NPOWF; RMAT(IP,JP)=RMAT(IP,JP-1)*RHO; ENDDO

20      CONTINUE
        


        CALL INVSOLVR(NPOWF,NPOWF,RMAT,NPOWF,-1)
   

        forall (i = 1:npowf,j=1:NF)

        coefb(i,j) = DOT_PRODUCT(rmat(i,:),coefb(1:npowf,j))

        end forall
        !COEFB(1:NPOWF,:)=MATMUL(RMAT,COEFB(1:NPOWF,:))    ! Bk IN Eq.(3-6-62)
        ! coefb has size 12xnum_intgd,but only from index 1 to 7 has useful info
      
        

     end subroutine compute_coeff_B
