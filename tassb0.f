C  TASSB0
C *********************************************************
C *                                                       *
C * Calculate the element contribution and assembly the   *
C * coefficients of the corresponding system of equationn *
C *                                                       *
C *********************************************************
C
        SUBROUTINE TASSB0
C  
	  USE MVAR_MOD
	  USE PVAR_MOD
	  USE MFUNC_mod
        USE SEBSM_MOD
        !use elem_intgl
        
        IMPLICIT   NONE  
	  INTEGER  INODE,IELEM,J,JNODE,IND,INDD,IP,
	1	       I,II,IS,JNCON,KNCON,L
	  REAL*8  XP,YP,ZP,XSB,YSB,ZSB,R
	  REAL*8  DX,DY,DZ,NX,NY,NZ
	  REAL*8 RSN(4,4),EX(4),EY(4)
	  
        REAL*8  BMATRIX(4,8),AMATRIX(4,8),BMAT(4)

	  REAL*8  S_ANGLE
!	  REAL*8  POXY 
	  REAL*8  DSIGN
!
        REAL*8  CELE31(4),CELE32(4),CELE33(4),AL1(4)
	  REAL*8  DPOX,DPOY,DPOZ,DPDN,PHI


        DATA RSN /1.,  1.,  1.,  1., 
     1            1., -1.,  1., -1.,
     2            1.,  1., -1., -1., 
     3            1., -1., -1.,  1./ 
!
	  DATA EX /  1.0d0,  1.0d0, -1.0d0, -1.0d0/                       
        DATA EY /  1.0d0, -1.0d0, -1.0d0,  1.0d0/
!
!  ----------------------------------------------------
!
	  WRITE(10, *)   ' IN TASSB0 '
	  DSDT(:)=0.0
!	  
!  ----------------------------------------------------
!	  	  
        DO 50 INODE=1, NNODE 
        L=0
        DO 40 IELEM=1,  NELEM
        DO 30 J=1,      NCN(IELEM)
        IF(INODE.EQ.NCON(IELEM,J)) THEN
        L=L+1
        NODELE(INODE,L)=IELEM
        NODELJ(INODE,L)=J
        ENDIF
30      CONTINUE
40      CONTINUE
        NODNOE(INODE)=L
!                          
        NODQUA(INODE)=0
        IF( NSYS .GE. 2) THEN
          IF( DABS(XYZ(2,INODE)).LT.1.0E-06 ) THEN
          NODQUA(INODE)=2
          END IF
        END IF
!
        IF( NSYS .EQ. 4) THEN
          IF( DABS(XYZ(1,INODE)).LT.1.0E-06.AND.
     1        DABS(XYZ(2,INODE)).LT.1.0E-06) THEN
           NODQUA(INODE)=5
          ELSE IF( DABS(XYZ(1,INODE)).LT.1.0E-06 ) THEN
           NODQUA(INODE)=4
          ENDIF
        END IF
!
50      CONTINUE
!                
! ***************************************

        DO  INODE=1, NNODE 
          write (14,51) inode,nodele(inode,1:5)
        end do
51    FORMAT(6I5)
!
        line_sum=0.0d0
        AMATA(:,:,:)=(0.0D0,0.0D0)
        BMATA(:,:)=(0.0D0,0.0D0)

!
! =======================================================================
! 
 	  WRITE(9, *) '   INODE      XP       YP       ZP       S_ANGLE'
        WRITE(10,*) '  INODE      ANGLE         A3   ',
     1               '       C31           C32           C33'
!
        DO  500   INODE=1,  NNF   ! Source point is on the free surface
!
	   XP=XYZ(1,INODE)
         YP=XYZ(2,INODE)
         ZP=XYZ(3,INODE) 
!
  	   CELE31(:)=0.0d0
	   CELE32(:)=0.0d0
	   CELE33(:)=0.0d0
	   AL1(:)=0.0d0	  
	  
!      
     	  CALL SOLIDANGLE(INODE,NNODE,NELEM,NCN,NCON,NODQUA,
     1                        H,XYZ,DXYZE,S_ANGLE)    
     
         S_ANGLE=1.0d0-S_ANGLE
	   WRITE(9,102)  INODE, XP, YP, ZP, S_ANGLE
      	 WRITE(*,102)  INODE, XP, YP, ZP, S_ANGLE
	   WRITE(101,102)  INODE, XP, YP, ZP, S_ANGLE
        
102	  FORMAT(I6,3F12.4,F15.6) 
!
	   ANGLE(INODE)=S_ANGLE
!
        DO   IP=1,  NSYS 
	    AMATA(INODE,INODE,IP)= ANGLE(INODE)
        ENDDO
!
!  ---------------------------
!  Integration on the free surface
!
        DO  300   IELEM=1,  NELEMF

        WRITE(101,*)
        WRITE(101,*) ' IELEM=',IELEM
        WRITE(101,*)  ' XP,YP,ZP=',XP,YP,ZP

         II=0   
         DO I=1, NODNOE(INODE)
          IF(IELEM .EQ. NODELE(INODE,I)) THEN
          II=II+1
          ENDIF
         ENDDO
C
C Using SING_ELE1 if the source point is in the element or its mirror
C     elements about any symmetrical axis, otherwise using NORM_ELE1
C

        IF (II .EQ. 0)   THEN 
         CALL NORM_ELE1(IELEM,XP,YP,ZP,AMATRIX,BMATRIX)
!         write(101,*) ' After NORM_ELE1'
        ELSE IF (II .NE. 0)   THEN 
          CALL SING_ELE1(INODE,IELEM,NODQUA(INODE),XP,YP,ZP,
     1                   AMATRIX,BMATRIX)
!         write(101,*) ' After SING_ELE1'
        END IF 

!       write(101,*) ' BMATRIX=',BMATRIX(1,:)
!       write(101,*) ' AMATRIX=',AMATRIX(1,:)

!
! ----------------------------
!
!
        DO  110   J=1,  NCN(IELEM) 
         JNCON=NCON(IELEM,J)
         KNCON=NCOND(IELEM,J)  
        DO  110   IP=1, NSYS          
          XSB=EX(IP)*XYZ(1,JNCON)
	    YSB=EY(IP)*XYZ(2,JNCON)
	    ZSB=       XYZ(3,JNCON)
          NX=EX(IP)*DXYZ(1,KNCON)
	    NY=EY(IP)*DXYZ(2,KNCON)
	    NZ=       DXYZ(3,KNCON)
	    
          CALL DINP(XSB,YSB,ZSB,DPOX,DPOY,DPOZ)       
           DPDN=DPOX*NX+DPOY*NY+DPOZ*NZ
            	 
         DO   IS=1, NSYS    
           AMATA(INODE,JNCON,IP)=AMATA(INODE,JNCON,IP)+
     1                           RSN(IS,IP)*BMATRIX(IS,J)		
        
           BMATA(INODE,IP)=BMATA(INODE,IP)+RSN(IS,IP)*AMATRIX(IS,J)
     1                        *POXY(XSB,YSB,ZSB)
         ENDDO        
!
! -------------------------
!
          CALL DINP0(0,XSB,YSB,ZSB,PHI,DPOX,DPOY,DPOZ)       
           DPDN=DPOX*NX+DPOY*NY+DPOZ*NZ
         DO    IS=1, NSYS             
           AL1(IP)=AL1(IP)-RSN(IS,IP)*BMATRIX(IS,J)*DPDN	
           AL1(IP)=AL1(IP)+RSN(IS,IP)*AMATRIX(IS,J)*PHI
         ENDDO

!
          CALL DINP0(1,XSB,YSB,ZSB,PHI,DPOX,DPOY,DPOZ)       
           DPDN=DPOX*NX+DPOY*NY+DPOZ*NZ
         DO    IS=1, NSYS             
           CELE31(IP)=CELE31(IP)-RSN(IS,IP)*BMATRIX(IS,J)*DPDN	
           CELE31(IP)=CELE31(IP)+RSN(IS,IP)*AMATRIX(IS,J)*PHI
         ENDDO
!
          CALL DINP0(2,XSB,YSB,ZSB,PHI,DPOX,DPOY,DPOZ)       
           DPDN=DPOX*NX+DPOY*NY+DPOZ*NZ
         DO    IS=1, NSYS             
           CELE32(IP)=CELE32(IP)-RSN(IS,IP)*BMATRIX(IS,J)*DPDN	
           CELE32(IP)=CELE32(IP)+RSN(IS,IP)*AMATRIX(IS,J)*PHI
         ENDDO
!
          CALL DINP0(3,XSB,YSB,ZSB,PHI,DPOX,DPOY,DPOZ)       
           DPDN=DPOX*NX+DPOY*NY+DPOZ*NZ
         DO    IS=1, NSYS             
           CELE33(IP)=CELE33(IP)-RSN(IS,IP)*BMATRIX(IS,J)*DPDN	
           CELE33(IP)=CELE33(IP)+RSN(IS,IP)*AMATRIX(IS,J)*PHI
         ENDDO

!
! -------------------------
!
110     CONTINUE        
!
          write(101,*) '  AL1=', AL1(1)
       


300     CONTINUE  
!	  
!  --------------------------
!  Integration on the body surface

        DO  400   IELEM=NELEMF+1,  NELEM
        WRITE(101,*)
        WRITE(101,*)  ' IELEM=',IELEM
        WRITE(101,*)  ' XP,YP,ZP=',XP,YP,ZP
         
         II=0
         DO  I=1, NODNOE(INODE)
         IF(IELEM .EQ. NODELE(INODE,I)) THEN
         II=II+1
         ENDIF
         ENDDO
!
        IF (II .EQ. 0)   THEN 
         CALL NORM_ELE1(IELEM,XP,YP,ZP,AMATRIX,BMATRIX)
!         write(101,*) ' After NORM_ELE1'
        ELSE IF (II .NE. 0)   THEN 
          CALL SING_ELE1(INODE,IELEM,NODQUA(INODE),XP,YP,ZP,
     1                   AMATRIX,BMATRIX)
!         write(101,*) ' After SING_ELE1'
        END IF                
! 
!       write(101,*) ' BMATRIX=',BMATRIX(1,:)
!       write(101,*) ' AMATRIX=',AMATRIX(1,:)
!
        DO  320   J=1,  NCN(IELEM) 
          JNCON=NCON(IELEM,J)
          KNCON=NCOND(IELEM,J)
         DO  320   IP=1, NSYS 
          XSB=EX(IP)*XYZ(1,JNCON)
	    YSB=EY(IP)*XYZ(2,JNCON)
	    ZSB=       XYZ(3,JNCON)
   
          NX=EX(IP)*DXYZ(1,KNCON)
	    NY=EY(IP)*DXYZ(2,KNCON)
	    NZ=       DXYZ(3,KNCON)

          CALL DINP(XSB,YSB,ZSB,DPOX,DPOY,DPOZ)       
          DPDN=DPOX*NX+DPOY*NY+DPOZ*NZ 
 
         DO  IS=1, NSYS    
	   IF(JNCON .GT. NNF)  THEN
          AMATA(INODE,JNCON,IP)=AMATA(INODE,JNCON,IP)-
     1                          RSN(IS,IP)*AMATRIX(IS,J)
         ELSE
          PHI=POXY(XSB,YSB,ZSB)
 	    BMATA(INODE,IP)=BMATA(INODE,IP)+RSN(IS,IP)*AMATRIX(IS,J)*PHI   !  * ******
         ENDIF

 	    BMATA(INODE,IP)=BMATA(INODE,IP)-RSN(IS,IP)*BMATRIX(IS,J)*DPDN  !  * ******
        ENDDO
!
! -------------------------
!
          CALL DINP0(0,XSB,YSB,ZSB,PHI,DPOX,DPOY,DPOZ)       
           DPDN=DPOX*NX+DPOY*NY+DPOZ*NZ
         DO    IS=1, NSYS             
           AL1(IP)=AL1(IP)-RSN(IS,IP)*BMATRIX(IS,J)*DPDN	
           AL1(IP)=AL1(IP)+RSN(IS,IP)*AMATRIX(IS,J)*PHI
         ENDDO
!
          CALL DINP0(1,XSB,YSB,ZSB,PHI,DPOX,DPOY,DPOZ)       
           DPDN=DPOX*NX+DPOY*NY+DPOZ*NZ
         DO    IS=1, NSYS             
           CELE31(IP)=CELE31(IP)-RSN(IS,IP)*BMATRIX(IS,J)*DPDN	
           CELE31(IP)=CELE31(IP)+RSN(IS,IP)*AMATRIX(IS,J)*PHI
         ENDDO
!
          CALL DINP0(2,XSB,YSB,ZSB,PHI,DPOX,DPOY,DPOZ)       
           DPDN=DPOX*NX+DPOY*NY+DPOZ*NZ
         DO    IS=1, NSYS             
           CELE32(IP)=CELE32(IP)-RSN(IS,IP)*BMATRIX(IS,J)*DPDN	
           CELE32(IP)=CELE32(IP)+RSN(IS,IP)*AMATRIX(IS,J)*PHI
         ENDDO
!
          CALL DINP0(3,XSB,YSB,ZSB,PHI,DPOX,DPOY,DPOZ)       
           DPDN=DPOX*NX+DPOY*NY+DPOZ*NZ
         DO    IS=1, NSYS             
           CELE33(IP)=CELE33(IP)-RSN(IS,IP)*BMATRIX(IS,J)*DPDN	
           CELE33(IP)=CELE33(IP)+RSN(IS,IP)*AMATRIX(IS,J)*PHI
         ENDDO
!
! -------------------------
!       
320      CONTINUE
!
          write(101,*) '  AL1=', AL1(1)

400     CONTINUE
!
!  --------------------------
!
        FrA3(INODE)=AL1(1)
        FrC31(INODE)=CELE31(1)
        FrC32(INODE)=CELE32(1)
        FrC33(INODE)=CELE33(1)
   
         WRITE(10,620) INODE,ANGLE(INODE),FrA3(INODE),
     1                 FrC31(INODE),FrC32(INODE),FrC33(INODE)

         PHI=POXY(XP,YP,ZP)
         CALL DINP(XP,YP,ZP,DPOX,DPOY,DPOZ)       
         BMATA(INODE,1)=BMATA(INODE,1)-FrA3(INODE)*PHI-
     1                 FrC31(INODE)*DPOX-FrC32(INODE)*DPOY
!        
500     CONTINUE
!
! =======================================================================
!    Source point is on the body surface
!
        DO  1000   INODE=NNF+1, NNODE   
!

	   XP=XYZ(1,INODE)
         YP=XYZ(2,INODE)
         ZP=XYZ(3,INODE) 
!	        
     	   CALL SOLIDANGLE(INODE,NNODE,NELEM,NCN,NCON,NODQUA,
     1                    H,XYZ,DXYZE,S_ANGLE) 
     
         S_ANGLE=1.0d0-S_ANGLE
    
	   WRITE(9,102)  INODE, XP, YP, ZP, S_ANGLE
      	 WRITE(*,102)  INODE, XP, YP, ZP, S_ANGLE
!
	   ANGLE(INODE)=S_ANGLE
!
        DO    IP=1,  NSYS 
	    AMATA(INODE,INODE,IP)= ANGLE(INODE)
        ENDDO        
! ------------------------------
! Intergration on the free surface
! 
        DO  800   IELEM=1,  NELEMF
        II=0   
C
C Using TNORWP0 to integrate 
C
         CALL NORM_ELE0(IELEM,XP,YP,ZP,AMATRIX,BMATRIX)

!  ------
!
        DO  710   J=1,  NCN(IELEM) 
         JNCON=NCON(IELEM,J)
         KNCON=NCOND(IELEM,J)
        DO  710   IP=1, NSYS  
          XSB=EX(IP)*XYZ(1,JNCON)
	    YSB=EY(IP)*XYZ(2,JNCON)
	    ZSB=       XYZ(3,JNCON)
   
          NX=EX(IP)*DXYZ(1,KNCON)
	    NY=EY(IP)*DXYZ(2,KNCON)
	    NZ=       DXYZ(3,KNCON)
          CALL DINP(XSB,YSB,ZSB,DPOX,DPOY,DPOZ)       
          DPDN=DPOX*NX+DPOY*NY+DPOZ*NZ 
      	 
         DO  710   IS=1, NSYS    
           AMATA(INODE,JNCON,IP)=AMATA(INODE,JNCON,IP)+
     1                           RSN(IS,IP)*BMATRIX(IS,J)
	     BMATA(INODE,IP)=BMATA(INODE,IP)+RSN(IS,IP)*AMATRIX(IS,J)*
     1               POXY(XSB,YSB,ZSB)	 	 
710     CONTINUE  
!

800     CONTINUE
	  
!  --------------------------
! Intergration on the body surface    
! 
        DO  900   IELEM=1+NELEMF, NELEM
        II=0   
C
C Using TSING if the source point is in the element or its mirror
C     elements about any symmetrical axis, otherwise using TINBOD
C
         DO  I=1, NODNOE(INODE)
         IF(IELEM .EQ. NODELE(INODE,I)) THEN
         II=II+1
         ENDIF
         ENDDO
!
        IF (II .EQ. 0)   THEN 
         CALL NORM_ELE0(IELEM,XP,YP,ZP,AMATRIX,BMATRIX)
        ELSE IF (II .NE. 0)   THEN 
         CALL SING_ELE0(INODE,IELEM,NODQUA(INODE),XP,YP,ZP,
     1                   AMATRIX,BMATRIX)
        END IF                
!

         DO 820 J=1,  NCN(IELEM) 
          JNCON=NCON(IELEM,J)         
          KNCON=NCOND(IELEM,J)
         DO 820 IP=1, NSYS                   
           XSB=EX(IP)*XYZ(1,JNCON)
	     YSB=EY(IP)*XYZ(2,JNCON)
	     ZSB=       XYZ(3,JNCON)
           NX=EX(IP)*DXYZ(1,KNCON)
	     NY=EY(IP)*DXYZ(2,KNCON)
	     NZ=       DXYZ(3,KNCON)
          CALL DINP(XSB,YSB,ZSB,DPOX,DPOY,DPOZ)       
          DPDN=DPOX*NX+DPOY*NY+DPOZ*NZ 
          
           DO  IS=1, NSYS    
	     IF(JNCON .GT. NNF)  THEN
           AMATA(INODE,JNCON,IP)=AMATA(INODE,JNCON,IP)-
     1                          RSN(IS,IP)*AMATRIX(IS,J)
	     ELSE	
	      PHI=POXY(XSB,YSB,ZSB)     	    
	      BMATA(INODE,IP)=BMATA(INODE,IP)+RSN(IS,IP)*AMATRIX(IS,J)*PHI  ! *******
	     ENDIF
            BMATA(INODE,IP)=BMATA(INODE,IP)-
     1                      RSN(IS,IP)*BMATRIX(IS,J)*DPDN	
           ENDDO
820     CONTINUE
!
900     CONTINUE
!
1000     CONTINUE
!
! =============================================

        IF( NSYS .EQ. 2) THEN
!
	   DO INODE=1, NNF
	    IF(NODQUA(INODE) .EQ. 2) THEN
	     AMATA(INODE,INODE,2)=1.0E20	 
		ENDIF
	   ENDDO
!
	  ELSE IF( NSYS .EQ. 4) THEN
!
	   DO INODE=1, NNF
	    IF(NODQUA(INODE) .EQ. 2) THEN
	     AMATA(INODE,INODE,2)=1.0E20
	     AMATA(INODE,INODE,4)=1.0E20	    
	    ELSE IF(NODQUA(INODE) .EQ. 4) THEN
	     AMATA(INODE,INODE,3)=1.0E20
	     AMATA(INODE,INODE,2)=1.0E20
		  ELSE IF(NODQUA(INODE) .EQ. 5) THEN
	     AMATA(INODE,INODE,2)=1.0E20
	     AMATA(INODE,INODE,3)=1.0E20	    
	     AMATA(INODE,INODE,4)=1.0E20	    
		ENDIF
	   ENDDO
!
	  ENDIF
!
            do i = 1,nnode
                 do j = 1,nnode
                   write(400,*) amata(i,j,1:nsys)
            end do;end do
            do i = 1,nnode
            write(401,*) bmata(i,1:nsys)
            end do
! =============================================
!
	 DO IP=1, NSYS
        WRITE(101, *) '  IP=',IP
	  WRITE(101, *) '    INODE=',nNODE,'      AMATA' 
        DO INODE=1,  NNODE
         WRITE(101, *) '  INODE=',INODE
         DO IND= 1,  NNODE 
         WRITE(101, 620) IND,XYZ(1,IND),XYZ(2,IND),XYZ(3,IND),
     1                    AMATA(INODE,IND,IP)
         ENDDO      
        ENDDO
       ENDDO
!
!
        WRITE(102, *) '  =========== Before RLUDCMP =============='
	 DO IP=1, NSYS
        WRITE(102, *) '  IP=',IP
	  WRITE(102, *) '    INODE          BMATA' 
         DO IND= 1,  NNODE 
         WRITE(102, 620) IND,XYZ(1,IND),XYZ(2,IND),XYZ(3,IND),
     1                    BMATA(IND,IP) 
         ENDDO      
       ENDDO
!       
!
!
        WRITE(102, *) 
        WRITE(102, *)
        WRITE(102, *) '  =========== After RLUDCMP =============='
	   DO IP=1, NSYS
           WRITE(6, *) '  IP=',IP,'    Before RLUDCMP'
	     CALL RLUDCMP(IP,AMATA,NNODE,NNODE,NSYS,INDX,DSIGN)  
	   ENDDO
!
         DO IS=1, NSYS   
           CALL RLUBKSB(IS,AMATA,NNODE,NNODE,1,NSYS,1,INDX,BMATA)
         ENDDO
!         
!
	 DO IP=1, NSYS
        WRITE(102, *) '  IP=',IP
	  WRITE(102, *) '    INODE          BMATA' 
         DO IND= 1,  NNODE 
         WRITE(102, 620) IND,XYZ(1,IND),XYZ(2,IND),XYZ(3,IND),
     1                    BMATA(IND,IP) 
         ENDDO      
       ENDDO
!
! =======================================================                 
! ** output the results
!      
        DO 1500 IP=1, NSYS 
        DO 1360 IND=1, NNODE 
        BMAT(IP)=(0.0D0, 0.0D0)
        DO 1350 IS=1, NSYS
1350     BMAT(IP)=BMAT(IP)+BMATA(IND,IS)*RSN(IP,IS)
        BMAT(IP)=BMAT(IP)/NSYS
        UNKN(IND,IP)=BMAT(IP)
1360     CONTINUE
!
1500     CONTINUE


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!用于输出Ax=B中的x
!       IF(TIMERK .GT. 1.0) THEN
!  
         WRITE(9, *)
	   WRITE(9, *)  '  Direvative of potential on the upper surface'
         WRITE(9, *)  '    INODE     XP      YP     ZP     Unkn    DPDZ'
           
        DO  INODE=1, NNF 
        XP=XYZ(1,INODE)
        YP=XYZ(2,INODE)
        ZP=XYZ(3,INODE) 
        CALL DINP(XP,YP,ZP,DPOX,DPOY,DPOZ)       

        WRITE(9, 620)  INODE,XP,YP,ZP,unkn(INODE,1),DPOZ
        ENDDO
        
         WRITE(9, *)
	   WRITE(9, *)  '  Potential on the body surface'
         WRITE(9, *)  '    INODE     XP      YP     ZP     Unkn    POXY'

        DO  INODE=1+NNF, NNODE
        XP=XYZ(1,INODE)
        YP=XYZ(2,INODE)
        ZP=XYZ(3,INODE) 
        WRITE(9, 620)  INODE,XP,YP,ZP,unkn(INODE,1),POXY(XP,YP,ZP)
        ENDDO

        do i = 1,NELEMF
          write (200,999) i,ncon(i,1:8)
C           x = 0.5*(XYZ(1,NCON(2))+XYZ(1,NCON(6)))
C           Y = 0.5*(XYZ(2,NCON(2))+XYZ(2,NCON(6)))
        end DO

         do i = NELEMF+1,NELEM
          write (201,999) i,ncon(i,1:8)
C           x = 0.5*(XYZ(1,NCON(2))+XYZ(1,NCON(6)))
C           Y = 0.5*(XYZ(2,NCON(2))+XYZ(2,NCON(6)))
        end DO


c 610     FORMAT(10x,'NNODE=',I6,/,2X,'IND',10X,'AMATA(IND,IND,1)')
!    
 999     format(9i4)
 610   FORMAT(10x,'NNODE     BMATA(IND,1)     T=',F14.6)
 620   FORMAT(1X,I4,6(1X,F13.6))      
 630   FORMAT(10x,'NNODE     UNKN(IND,1)     T=',F14.6)     
 640   FORMAT(1X,I4,2(2X,F13.6,2X,F13.6))                     
	            
                      
      RETURN
      END


