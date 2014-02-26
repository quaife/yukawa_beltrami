  PROGRAM testcon
! ------------------------------------------------------
! This is a test program for the Fortran 90 
! module Conical.
! ------------------------------------------------------
! 1)First, we compare the values obtained from 
!   the routine conic with 25 pre-computed values of the 
!   functions P^m_(-1/2+i tau) (x).
!   Relative errors are shown.
! 2)Second, we test the three-term recurrence relations:
!
!   For -1<x<1:  
!      pmu+1=-2*mu*x*pmu/sqrt(1-x*x)+(tau*tau+1/4+&
!             mu*(mu-1))*pmu-1
!   For x>1:
!      pmu+1=2*mu*x*pmu/sqrt(x*x-1)-(tau*tau+1/4+&
!             mu*(mu-1))*pmu-1 
!   
!   for several values of (m,tau,x).
! ------------------------------------------------------ 
  USE Conical
  USE someconstants
  IMPLICIT NONE
  REAL(r8) x,tau
  REAL(r8) x1(25),tau1(25),pcon(25)
  REAL(r8) pmtau,pmtau1,pmtau2,pmtau3,xargu,pmtest,ermax,err1
  INTEGER  mum(25),i,j,mu,mu1,mu2,mu3,ierr1,ierr2,ierr3
  DATA x1/-0.9_r8,0.0_r8,0.9_r8,1.1_r8,3.0_r8,5.0_r8,6.0_r8,7.0_r8,8.0_r8,&
          9.0_r8,10.0_r8,11.0_r8,13.0_r8,15.0_r8,17.0_r8,19.0_r8,21.0_r8,&
          23.0_r8,25_r8,30.0_r8,40.0_r8,50.0_r8,70.0_r8,90.0_r8,100.0_r8/
  DATA tau1/0.1_r8,0.5_r8,1.0_r8,3.0_r8,7.0_r8,9.0_r8,10.0_r8,11.0_r8,13.0_r8,&
           15.0_r8,17.0_r8,19.0_r8,21.0_r8,23.0_r8,25_r8,30.0_r8,35.0_r8,40.0_r8,&
           50.0_r8,55.0_r8,60.0_r8,70.0_r8,80.0_r8,90.0_r8,100.0_r8/
  DATA mum/0,1,2,5,7,10,15,20,22,25,30,33,35,40,45,50,55,60,65,70,75,80,85,&
           90,100/
  DATA pcon/1.8900720543540235_r8,0.5807488453711801_r8,&
            0.10918150796576438_r8,3.766101638152348_r8,&
           -297122.06012454815_r8,-1.8418257208922422e+9_r8,&
           -4.1221857032656785e+15_r8,-2.337706510495258e+22_r8,&
           -9.205609694695491e+25_r8,1.9573607117438758e+31_r8,&
           -3.517059633224470e+39_r8,-1.2857020029215802e+45_r8,&
           -1.7684756798987016e+49_r8,-1.9584930819993896e+58_r8,&
           -3.464798870971219e+67_r8,5.197774400454444e+77_r8,&
           -6.349599741729241e+88_r8,6.1958997557824016e+100_r8,&
           -5.7767267898372044e+113_r8,8.149772408354674e+125_r8,&
            2.178951636437889e+137_r8,-1.8017434379090414e+151_r8,&
            1.6655668080406992e+164_r8,9.350636511547676e+178_r8,&
            1.7249272890028403e+203_r8/
  open(unit=1,file='testf.dat',status='unknown') 
! ------------------------------------------------
! We check the values of the functions
! ------------------------------------------------
  write(1,*)'ERR means relative error'
  write(1,*)'***************************************************'
  write(1,*)'Test of the 25 pre-computed values of the functions'
  write(1,*)'***************************************************'
  write(1,30)'x','tau','m','ERR(Pconic)'
  DO i=1,25
    x=x1(i)
    tau=tau1(i)
    mu=mum(i) 
    ierr1=0
    CALL conic(x,mu,tau,pmtau,ierr1)
    err1=abs(1.0_r8-pmtau/pcon(i))
    write(1,31)x,tau,mu,err1
  ENDDO
 ! ---------------------------------------------------------
 !  Test of the three-term recurrence relations
 ! ---------------------------------------------------------
 ! We test the relations for 
 !   a)  -1<x<1,  0<tau<=100, 0<=m<=40
 !   b)     x>1,  0<tau<=100, 0<=m<=100       
 ! ---------------------------------------------------------
  write(1,*)' '
  write(1,*)'----------------------------------------'
  write(1,*)'Test of the TTRRs for P^m_(-1/2+itau)(x)' 
  write(1,*)'----------------------------------------'
  ermax=0.0_r8
  DO i=1,100,1
    x=-0.9_r8+i*1.0_r8
    mu=i
    DO j=0,99
      tau=0.1_r8+j  
      ierr1=0
      ierr2=0
      ierr3=0
      mu1=mu-1
      mu2=mu
      mu3=mu+1
      CALL conic(x,mu1,tau,pmtau1,ierr1)
      CALL conic(x,mu2,tau,pmtau2,ierr2)
      CALL conic(x,mu3,tau,pmtau3,ierr3)
      IF ((ierr1==0).AND.(ierr2==0).AND.(ierr3==0)) THEN
        IF (x<1.0_r8) THEN
          xargu=1.0_r8-x*x 
          pmtest=-2.0_r8*mu*x*pmtau2/sqrt(xargu)+(tau*tau+0.25_r8+&
                 mu*mu1)*pmtau1
        ELSE
          xargu=x*x-1.0_r8 
          pmtest=2.0_r8*mu*x*pmtau2/sqrt(xargu)-(tau*tau+0.25_r8+&
                 mu*mu1)*pmtau1
        ENDIF
        err1=abs(1.0_r8-pmtest/pmtau3)
        IF (err1>ermax) THEN
          ermax=err1
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  write(1,*)' '
  write(1,*)'Maximum value of the TTRR check =',ermax
 30 format(6x,a1,13x,a3,7x,a8,3x,a12)
 31 format(d16.9,1x,d16.9,1x,i4,1x,d16.9)    
  END PROGRAM testcon

     


