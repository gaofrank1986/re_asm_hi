!  
!  G=-1/4/pi*(1/r+1/r1) 
!
!
        SUBROUTINE TGRN(H,X,X0,Y,Y0,Z,Z0,XHF) 
! 
        IMPLICIT NONE
		INTEGER I

	    REAL*8,INTENT(IN)::  H,X,Y,Z,X0,Y0,Z0
		REAL*8,INTENT(OUT):: XHF(4)
	    REAL*8 DX,DY,DZ,DZ1
		REAL*8 RXY2,DZ02,DZ12,SR2,ST2,SR,SR1,PI,PI4
!
	    DATA PI/3.14159265358979/ 
!C 
!C H: water depth, negtive for infinity water depth 
!C
        PI4 = 4.0D0*PI

        DX =X-X0
        DY =Y-Y0
        DZ =Z- Z0
        DZ1=Z+ Z0 +2.0*H


		RXY2=DX*DX+DY*DY
	    DZ02=DZ*DZ
		DZ12=DZ1*DZ1
!
	   IF(H .GT. 0) THEN
        SR2=RXY2+DZ02                      
        ST2=RXY2+DZ12           
        SR =DSQRT(SR2) 
        SR1=DSQRT(ST2)   
! 
        XHF(1)= 1.D0/SR +  1.D0/SR1 
        XHF(2)=-DX/SR**3 -  DX /SR1**3 
        XHF(3)=-DY/SR**3 -  DY /SR1**3 
        XHF(4)=-DZ/SR**3 -  DZ1/SR1**3    
! 
	   ELSE
        SR2=RXY2+DZ02                      
        SR =DSQRT(SR2) 
! 
        XHF(1)= 1.D0/SR  
        XHF(2)=-DX/SR**3  
        XHF(3)=-DY/SR**3  
        XHF(4)=-DZ/SR**3     
! 
	   ENDIF
!
       DO 120   I=1,  4
         XHF(I)=-XHF(I)/PI4 
120	   CONTINUE
!
        RETURN 
        END 
		    
		    
!  
! ==================================================================
!
!  G=-1/4/pi*(1/r+1/r1) 
! XF(1) : dG(x,x0)/dz0
! XF(2) :d^2G/dx/dz0;    XF(3):d^2G/dy/dz0；   XF(4):d^2G/dz/dz0
!
! ==================================================================
!
        SUBROUTINE DTGRN(H,X,X0,Y,Y0,Z,Z0,XHF) 
! 
        IMPLICIT NONE
		INTEGER I

	    REAL*8,INTENT(IN)::  H,X,Y,Z,X0,Y0,Z0
		REAL*8,INTENT(OUT):: XHF(4)
	    REAL*8 DX,DY,DZ,DZ1
		REAL*8 RXY2,DZ02,DZ12,SR,SR2,SR3,SR5,ST,ST2,ST3,ST5
		REAL*8 PI,PI4
!
	    DATA PI/3.14159265358979/ 
! 
! H: water depth, negtive for infinity water depth 
!
        PI4 = 4.0D0*PI

        DX=X-X0
        DY=Y-Y0
        DZ =Z- Z0
        DZ1=Z+ Z0 +2.0*H

		RXY2=DX*DX+DY*DY
	    DZ02=DZ*DZ
		DZ12=DZ1*DZ1
!
	   IF(H .GT. 0) THEN
        SR2=RXY2+DZ02                      
        SR =DSQRT(SR2) 
        SR3=SR*SR2
        SR5=SR3*SR2

        ST2=RXY2+DZ12           
        ST=DSQRT(ST2)   
        ST3=ST*ST2
        ST5=ST3*ST2
! 
        XHF(1)=DZ/SR3 -DZ1/ST3 
        XHF(2)=3.0d0*DX*(-DZ/SR5 + DZ1/ST5)
        XHF(3)=3.0d0*DY*(-DZ/SR5 + DZ1/ST5) 
        XHF(4)=-3.0d0*DZ02/SR5+1.0d0/SR3+3.0d0*DZ12/ST5-1.0d0/ST3    
! 
	   ELSE
        SR2=RXY2+DZ02                      
        SR =DSQRT(SR2) 
        SR3=SR*SR2
        SR5=SR3*SR2
! 
        XHF(1)=DZ/SR3 
        XHF(2)=-3.0d0*DX*DZ/SR5 
        XHF(3)=-3.0d0*DY*DZ/SR5  
        XHF(4)=-3.0d0*DZ02/SR5+1.0d0/SR3 
! 
	   ENDIF
!
       DO 120   I=1,  4
         XHF(I)=-XHF(I)/PI4 
120	   CONTINUE
!
        RETURN 
        END 
		           


!  
! ==================================================================
!
!  G=-1/4/pi*1/r 
! XF(1) : dG(x,x0)/dz0
! XF(2) :d^2G/dx/dz0;    XF(3):d^2G/dy/dz0；   XF(4):d^2G/dz/dz0
!
! ==================================================================
!
        SUBROUTINE DTGRN0(H,X,X0,Y,Y0,Z,Z0,XHF) 
! 
        IMPLICIT NONE
		INTEGER I

	    REAL*8,INTENT(IN)::  H,X,Y,Z,X0,Y0,Z0
		REAL*8,INTENT(OUT):: XHF(4)
	    REAL*8 DX,DY,DZ
		REAL*8 RXY2,DZ02,SR,SR2,SR3,SR5
		REAL*8 PI,PI4
!
	    DATA PI/3.14159265358979/ 
! 
! H: water depth, negtive for infinity water depth 
!
        DX=X-X0
        DY=Y-Y0
        DZ=Z-Z0

		RXY2=DX*DX+DY*DY
	    DZ02=DZ*DZ
!
        SR2=RXY2+DZ02                      
        SR =DSQRT(SR2) 
        SR3=SR*SR2
        SR5=SR3*SR2
! 
        XHF(1)=DZ/SR3 
        XHF(2)=-3.0d0*DX*DZ/SR5 
        XHF(3)=-3.0d0*DY*DZ/SR5  
        XHF(4)=-3.0d0*DZ02/SR5+1.0d0/SR3     
! 
        PI4 = 4.0D0*PI
!
       DO 120   I=1,  4
         XHF(I)=-XHF(I)/PI4 
120	   CONTINUE
!
        RETURN 
        END 
		           
!  
!
!===================================================
!  G=-1/4/pi*(1/r1) 
! XF(1) : dG(x,x0)/dz0
! XF(2) :d^2G/dx/dz0;     XF(3):d^2G/dy/dz0；      XF(4):d^2G/dz/dz0
!
!
        SUBROUTINE DTGRN1(H,X,X0,Y,Y0,Z,Z0,XHF) 
! 
        IMPLICIT NONE
		INTEGER I

	    REAL*8,INTENT(IN)::  H,X,Y,Z,X0,Y0,Z0
		REAL*8,INTENT(OUT):: XHF(4)
	    REAL*8 DX,DY,DZ,DZ1
		REAL*8 RXY2,DZ02,DZ12,SR2,ST2,SR,SR3,SR5,ST,ST3,ST5
		REAL*8 PI,PI4
!
	    DATA PI/3.14159265358979/ 
!C 
!C H: water depth, negtive for infinity water depth 
!C
        PI4 = 4.0D0*PI

        DX=X-X0
        DY=Y-Y0
        DZ1 =Z+ Z0 +2.0*H


		RXY2=DX*DX+DY*DY
		DZ12=DZ1*DZ1
!
	   IF(H .GT. 0) THEN

        ST2=RXY2+DZ12           
        ST=DSQRT(ST2)   
        ST3=ST*ST2
        ST5=ST3*ST2
! 
        XHF(1)= -DZ1/ST3 
        XHF(2)=3.0d0*DX*DZ1/ST5
        XHF(3)=3.0d0*DY*DZ1/ST5 
        XHF(4)=3.0d0*DZ12/ST5-1.0d0/ST3    
! 
	   ELSE
! 
        XHF(:)=0.0d0 
! 
	   ENDIF
!
       DO 120   I=1,  4
         XHF(I)=-XHF(I)/PI4 
120	   CONTINUE
!
        RETURN 
        END 
		     

!===================================================
!   G=-1/4/pi*(1/r+1/r1)      ！！没有写完整
! XF(1) : G(x,x0)
! XF(2) :dG/dx;     XF(3):dG/dy；      XF(4):dG/dz
! XF(5):dG/dx/dx0;  XF(6):dG/dx/dy0;   XF(7):dG/dx/dz0
! XF(8):dG/dy/dx0;  XF(9):dG/dy/dy0;   XF(10):dG/dy/dz0
! XF(11):dG/dz/dx0; XF(12):dG/dz/dy0;  XF(13):dG/dz/dz0
!

      SUBROUTINE GREEN1(H,X,X0,Y,Y0,Z,Z0,XFF) 
      IMPLICIT NONE

  	  REAL*8,INTENT(IN)::  H,X,X0,Y,Y0,Z,Z0
	  REAL*8,INTENT(OUT):: XFF(13)

	  REAL*8 DX,DY,DZ,DZ1,RXY
	  REAL*8 DX2,DY2,DZ02,DZ12
	  REAL*8 RXY2,SR2,ST2,R,R1
	  REAL*8 PI,PI4

!
	  DATA PI/3.14159265358979/ 
! 
! H: water depth, negtive for infinity water depth 
! 
        PI4 = 4.*PI


	   IF(H .GT. 0) THEN

      DX = X-X0
      DY = Y-Y0
      DZ = Z-Z0
      DZ1= Z+Z0+2*H

      DX2 = DX*DX
      DY2 = DY*DY
      RXY2= DX2 + DY2
      RXY = DSQRT(RXY2)

      DZ02= DZ*DZ
      DZ12= DZ1*DZ1
      SR2=RXY2+DZ02                     
      ST2=RXY2+DZ12          
  
      R =DSQRT(SR2)
      R1=DSQRT(ST2)  
        
      XFF(1) =(-1./R -1./R1)/PI4
	  XFF(2) =(DX/R**3+ DX/R1**3)/PI4
	  XFF(3) =(DY/R**3+ DY /R1**3)/PI4
	  XFF(4) =(DZ/R**3+ DZ1/R1**3)/PI4   

	  XFF(5) =((3.*DX2-SR2)/R**5+(3.*DX2-ST2)/R1**5)/PI4
   	  XFF(6) =(3.*DX*DY*(1./R**5+1./R1**5))/PI4
      XFf(7) =(3.*DX*(DZ/R**5))/PI4                       !-DZ1/R1**5))/PI4

      XFF(8) =(3.*DX*DY*(1./R**5+1./R1**5))/PI4
      XFF(9)=((3.*DY2-SR2)/R**5+(3.*DY2-ST2) /R1**5)/PI4
      XFf(10)=(3.*DY*(DZ/R**5))/PI4                       !-DZ1/R1**5))/PI4

      XFF(11)=(3.*DX*(DZ/R**5+DZ1/R1**5))/PI4
      XFF(12)=(3.*DY*(DZ/R**5+DZ1/R1**5))/PI4
      XFf(13)=((3.*DZ02-SR2)/R**5)/PI4                    !-(3.*DZ12-ST2) /R1**5)/PI4
    
      ELSE

     DZ  = Z-Z0
  
     DX2 = DX*DX
     DY2 = DY*DY
     RXY2= DX2 + DY2
     RXY = DSQRT(RXY2)

     DZ02= DZ*DZ
     SR2=RXY2+DZ02                     
  
     R =DSQRT(SR2)
        
     XFF(1) =-1./R /PI4
	 XFF(2) =DX/R**3/PI4
	 XFF(3) =DY/R**3/PI4
	 XFF(4) =DZ/R**3/PI4   

 	 XFF(5) =(3.*DX2-SR2)/R**5/PI4
     XFF(6) =3.*DX*DY/R**5/PI4
     XFf(7) =3.*DX*DZ/R**5/PI4                       !-DZ1/R1**5))/PI4

     XFF(8) =3.*DX*DY/R**5/PI4
     XFF(9) =(3.*DY2-SR2)/R**5/PI4
     XFf(10)=3.*DY*DZ/R**5/PI4                       !-DZ1/R1**5))/PI4

     XFF(11)=3.*DX*DZ/R**5/PI4
     XFF(12)=3.*DY*DZ/R**5/PI4
     XFf(13)=(3.*DZ02-SR2)/R**5/PI4                 !-(3.*DZ12-ST2) /R1**5)/PI4

      ENDIF
  

     
	RETURN
   END 