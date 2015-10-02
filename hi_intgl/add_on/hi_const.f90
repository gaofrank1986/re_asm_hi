!module hi_const
    integer,parameter,private   :: ngr = -10
    integer,private             :: ngl = -10 ! how many gaussian point on 

    real(8),parameter,private   :: dlt(3,3) = reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))
    ! d-matrix--unix matix
    
    real(8),parameter,private   :: hi_pi  = 4.d0*datan(1.d0)
    real(8),parameter,private   :: cnu = 0.3
    real(8),parameter,private   :: TOLGP=1.D-14


    real(8),parameter,private :: cnr_lcl_mtx(18) =(/-1.,-1.,1.,-1.,1.,1.,-1.,1.,0.,-1.,1.,0.,0.,1.,          &
        &         -1.,0.,0.,0./) ! 9-point element

    integer,parameter,private :: node_grp_by_edge(12) = (/1,2,5, 2,3,6, 3,4,7, 4,1,8/)

    real(8),parameter,private :: tol=1.d-5 !tolerence?

    integer,parameter,private :: ksb(4)=(/-2,1,2,-1/) !some matrix,delcared
    ! IF(DABS(XIP(IABS(KS))-DBLE(KS)/DABS(DBLE(KS))).LT.TOL) GOTO 80  
    ! edge 1 , xip(2) != -1
    ! edge 2 , xip(1) != 1
    ! edge 3 , xip(2) != 1
    ! edge 4 , xip(1) != -1
    ! why?

!end module hi_const
