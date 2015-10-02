module hi_intg 

    use hi_mod_funcs
    
    implicit none 
    
    include './add_on/hi_const.f90'    
    integer,protected    ::  num_dim,num_node,num_nrml,num_elem
    integer,protected    ::  elem_nd_count,num_intgd
    real(8),protected    ::  hi_beta 


    real(8),allocatable,private ::  node_matrix(:,:),normal_matrix(:,:),src_local_list(:,:)
    !node_matrix(1:3,node_id),normal_matrix(1:3,nrml_id)
    integer,allocatable,private ::  elem_matrix(:,:),src_flag(:)
    !elem_matrix(1:8,elem_id)
    real(8),allocatable,private ::  full_mesh_matrix(:,:,:)
    !full_mesh_matrix(1:3,1:8,elem_id)

    
    integer,private :: model_readed_flag = 0 ! 0 for not readed
    !integer,private :: external_src_ctr_flag = 0! if or not use external src ctr input

    integer,private,parameter :: NPW = 4
    real(8),private,allocatable :: value_list(:,:)
    integer,private :: n_pwr_g = -1
    real(8),allocatable,private :: cnr_glb_mtx(:,:) !corner_global_matrix
    real(8),private :: src_lcl_preset(2),src_glb_preset(3)
   
    real(8),private ::  src_glb(3),src_ctr_glb(3)
    !src coordinate in global,src center point in global
    !private :: read_model_from_WAVDUT
contains
    include './add_on/test6.f90'
    include './add_on/test3.f90'
    include './add_on/test2.f90'
    include './add_on/test.f90'
    include './add_on/hi_integrand.f90'        
    include './add_on/run_thru_elems.f90'    
    
    subroutine get_node_matrix(nd,ex_node_matrix)
        implicit none
        integer ::nd
        real(8),intent(out) :: ex_node_matrix(3,nd)
        
        ex_node_matrix = node_matrix
    end subroutine 
    

    subroutine read_model_from_DAT()

        implicit none

        integer :: ip,ie,tmp,i,id        

        if (model_readed_flag == 0) then
            print *,"------------Start Reading Model-------------"
            OPEN(5,FILE='SIEPPEM.DAT',STATUS='OLD')

            read (5,*) num_dim,num_node,num_elem,elem_nd_count,hi_beta,num_intgd
            ! number of node per element
            ! beta is the power of r in target equation
            ! number of target func components
            allocate(node_matrix(num_dim,num_node))
            allocate(elem_matrix(elem_nd_count,num_elem))
            allocate(src_flag(num_elem))
            allocate(src_local_list(2,num_elem))
            allocate(value_list(num_elem,num_intgd))

            if (num_dim == 2) ngl = 1
    
         !    Input nodal coordinates and element connectivity
          
            do ip = 1,num_node
                read(5,*) tmp,(node_matrix(i,tmp),i=1,num_dim)                 ! card set 2
            end do  

            do ie = 1,num_elem
                read(5,*) tmp,(elem_matrix(id,tmp),id=1,elem_nd_count),src_flag(tmp)    ! card set 3
            end do
            !====src_flag
            ! if = 0 not valid
            ! if > 0 src is given in global coordinate, use node with id (src_flag)
            ! if < 0 src is given in local coordinate, use local src list given in card set 4

            read(5,*) (src_glb(i),i=1,num_dim)                       ! card set 4  
             ! read src x,y,z coord
             ! there seems a error, src_glb should be an array of global coordinate
             ! since src_glb cannot remain unchanged for different element
        
            do ie=1,num_elem
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

            model_readed_flag = 1 
            value_list = 0

        else
            print *,"--Attention! Reading process skipped,model already loaded---"
        end if

    end subroutine read_model_from_DAT

    subroutine read_model_from_WAVDUT()
        !use mesh
         USE MVAR_MOD
        implicit none

        integer :: i,j,ie,id,tmp,tmp1

        print *,"------------Start Reading Model from WAVDUT-------------"

        num_dim = 3
        num_node = NNODE
        num_elem = NELEM
        num_nrml = NNODED
        elem_nd_count = 8 !NCN(IELEM)!! to be changed
        hi_beta = 3.
        num_intgd = 8

        allocate(node_matrix(num_dim,num_node))
        allocate(normal_matrix(num_dim,num_nrml))

        allocate(elem_matrix(elem_nd_count,num_elem))
        allocate(value_list(num_elem,num_intgd))
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

        !switch order
        print *,"Transfering elem matrix from old order to new order"
            do i = 1,num_elem
                if (elem_nd_count .eq. 8) then
                    print *,"god knows why@-@"
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
        model_readed_flag = 1
        value_list = 0

        print *,"------------Stop Reading Model from WAVDUT-------------"
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


    
end module
      


