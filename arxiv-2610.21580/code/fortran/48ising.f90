!----------------------!
 module systemvariables
 implicit none

 integer :: lx 
 integer :: ly                       ! system length
 integer :: n                        ! number of spins (n=lx*ly)
 integer :: nn1                      ! number of spins (nn1=lx*ly)
 integer :: nb                       ! number of bonds (nb=2*n)
 real(8) :: pflip(-1:1,-5:5)         ! flip probabilities
 integer, allocatable :: bsites(:,:)
 integer, allocatable :: spin(:)     ! spin array
 real(8) :: j1
 real(8) :: j2
 real(8) :: h
 real(8) :: t
 end module systemvariables
!--------------------------!

!-------------------------------------------------------------------!
!Metropolis-algorithm simulation of the two-dimensional Ising model.!
!Reads the following from a file 'ising.in':                        !
! l,t,h = system length, temperature (T/J), field (h.J)             !
! bins, binsteps = number of bins, MC steps per bin                 !
!Writes these bin averages to the file 'bindata.dat':               ! 
! <E>/N, <E**2>/N**2, <|M|>/N, <M**2>/N**2                          !
!These results should be processed by the program average.f90.      !
!A single random numbers seeds (3) are read from a file 'seed.in'.
!-------------------------------------------------------------------!
 program ising2d
!----------------------------------!
 use systemvariables; implicit none

 integer :: i,j,bins,relaxsteps,binsteps,seed

 open(10,file='ising.in',status='old')
 read(10,*)lx,ly,t,h,j2,j1
 read(10,*)bins,relaxsteps,binsteps
 close(10)

 call initialize()
 call makelattice()
 do i=1,relaxsteps
    call mcstep
 end do
 do j=1,bins
    call resetdatasums
    do i=1,binsteps
       call mcstep
	   call measure()
       !if (mod(binsteps,10)==0) call measure()       
    end do
    call writebindata(binsteps)
 end do
 call test
 call writeaverage
 call values
 deallocate(spin)
 deallocate(bsites)   
!-------------------!
 end program ising2d
!-------------------!
!------------------------!
 subroutine makelattice()
!------------------------!
 use systemvariables; implicit none

 integer :: x1,x2,y1,y2,s,s1,s2

 nn1=lx*ly
 n=nn1
 nb=nn1*8
 allocate(bsites(2,nb))
 
 do y1=0,ly-1
 do x1=0,lx-1
    s=1+x1+y1*lx
    bsites(1,4*s-3)=4*s-3
    bsites(2,4*s-3)=4*s-2
    bsites(1,4*s-2)=4*s-2
    bsites(2,4*s-2)=4*s-1
    bsites(1,4*s-1)=4*s-1
    bsites(2,4*s-1)=4*s
    bsites(1,4*s)=4*s
    bsites(2,4*s)=4*s-3
    s1=1+mod(x1+1,lx)+y1*lx
    s2=1+mod(x1-1+lx,lx)+y1*lx
    bsites(1,s+4*n)=4*s-2
    bsites(2,s+4*n)=4*s1
    bsites(1,s+5*n)=4*s
    bsites(2,s+5*n)=4*s2-2
    s1=1+x1+mod(y1+1,ly)*lx
    s2=1+x1+mod(y1+ly-1,ly)*lx
    bsites(1,s+6*n)=4*s-1
    bsites(2,s+6*n)=4*s1-3   
    bsites(1,s+7*n)=4*s-3
    bsites(2,s+7*n)=4*s2-1     
 enddo
 enddo

 end subroutine makelattice
!--------------------------!

!--------------------------------------------!
!Carries out one Monte Carlo step, defined as! 
!n flip attempts of randomly selected spins. !
!--------------------------------------------!
 subroutine mcstep
!----------------------------------!
 use systemvariables; implicit none

 integer :: i,s,x,y,s1,s2,s3,s4,s0
 real(8), external :: ran

 do i=1,n
    s=int(ran()*n)
    x=mod(s,lx); y=s/lx
    s1=1+x+mod(y+1,ly)*lx
    s2=1+mod(x+1,lx)+y*lx
    s3=1+x+mod(y+ly-1,ly)*lx
    s4=1+mod(x-1+lx,lx)+y*lx
    s0=1+x+y*lx
    if (ran()<pflip(spin(4*s0-1),spin(4*s0-2)+spin(4*s0)+3*spin(4*s1-3))) then
        spin(4*s0-1)=-spin(4*s0-1)
    endif
    if (ran()<pflip(spin(4*s0-2), spin(4*s0-1)+spin(4*s0-3)+3*spin(4*s2))) then
		spin(4*s0-2)=-spin(4*s0-2)
	endif
    if (ran()<pflip(spin(4*s0-3), spin(4*s0-2)+spin(4*s0)+3*spin(4*s3-1))) then
		spin(4*s0-3)=-spin(4*s0-3)
	endif
    if (ran()<pflip(spin(4*s0), spin(4*s0-1)+spin(4*s0-3)+3*spin(4*s4-2))) then 
		spin(4*s0)=-spin(4*s0)
	endif
end do

 end subroutine mcstep
!---------------------!

!--------------------------------------------------------!
!Measures the energy (enrg1), the squared energy (enrg2),! 
!the absolute value of the magnetization (magn1), and the!
!squared magnetization (magn2).                          !
!--------------------------------------------------------!
 subroutine measure()
!----------------------------------!
 use systemvariables; implicit none

 real(8) :: enrg1,enrg2,magn1,magn2,magn3
 common/measurments/enrg1,enrg2,magn1,magn2,magn3

 integer :: s,x,y,e,m,s1,s2,s3,s4
 e=0
 m=0
 do s=1,n
    x=mod(s-1,lx); y=s/ly
    s1=1+x+mod(y+1,ly)*lx
    s2=1+mod(x+1,lx)+y*lx
    s3=1+x+mod(y+ly-1,ly)*lx
    s4=1+mod(x-1+lx,lx)+y*lx
    e=e+(spin(4*s-1)*(j1*spin(4*s-2)+j1*spin(4*s)+j2*spin(4*s1-3))&
            +spin(4*s-2)*(j1*spin(4*s-1)+j1*spin(4*s-3)+j2*spin(4*s2))&
            +spin(4*s-3)*(j1*spin(4*s-2)+j1*spin(4*s)+j2*spin(4*s3-1))&
            +spin(4*s)*(j1*spin(4*s-3)+j1*spin(4*s-1)+j2*spin(4*s4-2)))/2
    m=m+((-1)**(x+y))*(-spin(4*s)+spin(4*s-1)-spin(4*s-2)+spin(4*s-3))
 enddo
 enrg1=enrg1+dble(e-m*h)
 enrg2=enrg2+dble(e-m*h)**2
 magn1=magn1+dble(abs(m))
 magn2=magn2+dble(m)**2
 magn3=magn3+dble(m)**4
 end subroutine measure
!----------------------!

!-------------------------------------------------!
!Writes accumulated data to the file 'bindata.dat'!
!-------------------------------------------------!
 subroutine writebindata(steps)
!------------------------------!
 use systemvariables; implicit none

 real(8) :: enrg1,enrg2,magn1,magn2,magn3
 common/measurments/enrg1,enrg2,magn1,magn2,magn3

 integer :: steps
 
 open(1,file='bindata.dat',status='unknown',position='append')
 enrg1=enrg1/(dble(steps)*dble(n*4))
 enrg2=enrg2/(dble(steps)*dble(n*4)**2)
 magn1=magn1/(dble(steps)*dble(n*4))
 magn2=magn2/(dble(steps)*dble(n*4)**2)
 magn3=magn3/(dble(steps)*dble(n*4)**4)
 write(1,1)enrg1,enrg2,magn1,magn2,magn3
 1 format(5f18.12)
 close(1)

 end subroutine writebindata
!---------------------------!

!--------------------------------------------!
!Sets the data accumulation variables to zero!
!--------------------------------------------!
 subroutine resetdatasums
!------------------------!
 implicit none

 real(8) :: enrg1,enrg2,magn1,magn2,magn3
 common/measurments/enrg1,enrg2,magn1,magn2,magn3

 enrg1=0.d0
 enrg2=0.d0
 magn1=0.d0
 magn2=0.d0
 magn3=0.d0

 end subroutine resetdatasums
!----------------------------!

!------------------------------------------------------------------!
!Various initializations: system size n=l*l, flipping probabilities!
!pflip, random number generator, initial random spin state.        !
!------------------------------------------------------------------!
 subroutine initialize()
!--------------------------!
 use systemvariables
 implicit none

 integer :: i,j,ns,nn

 integer, allocatable :: seed(:)
 real(8), external :: ran

 n=lx*ly
 do i=-2,2
    pflip(-1,i+3)=exp(-2.*(i*j1+j2+h)/t)
    pflip(-1,i-3)=exp(-2.*(i*j1-j2+h)/t)
    pflip(+1,i+3)=exp(2.*(i*j1+j2+h)/t)
    pflip(+1,i-3)=exp(2.*(i*j1-j2+h)/t)
 enddo
 
 call initran(1)
 nn=4*n
 allocate (spin(1:nn))
 do i=1,nn
    spin(i)=2*int(2.d0*ran())-1
 enddo

 end subroutine initialize
!-------------------------!

!----------------------!
 real(8) function ran()
!----------------------------------------------!
! 64-bit congruental generator                 !
! iran64=oran64*2862933555777941757+1013904243 !
!----------------------------------------------!
 implicit none

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64

 ran64=ran64*mul64+add64
 ran=0.5d0+dmu64*dble(ran64)

 end function ran
!----------------!

!---------------------!
 subroutine initran(w)
!---------------------!
 implicit none

 integer(8) :: irmax
 integer(4) :: w,nb,b

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64
      
 irmax=2_8**31
 irmax=2*(irmax**2-1)+1
 mul64=2862933555777941757_8
 add64=1013904243
 dmu64=0.5d0/dble(irmax)

 open(10,file='seed.in',status='old')
 read(10,*)ran64
 close(10)
 if (w.ne.0) then
    open(10,file='seed.in',status='unknown')
    write(10,*)abs((ran64*mul64)/5+5265361)
    close(10)
 endif

 end subroutine initran
!----------------------!

!----measure average----!
subroutine test
    implicit none
    integer::i,j,num,end_num
    real,ALLOCATABLE,dimension(:,:)::list
    real::list_sum1,list_sum2,list_sum3,list_sum4,list_sum5
    common/tests/list_sum1,list_sum2,list_sum3,list_sum4,list_sum5,num
    character(len=4)::temp
    open(12,file='bindata.dat',status='old',action='read')
    num=0
    conuter_number:&
    do
        read(12,*,iostat=end_num) temp
        if(end_num<0)then  !!确定文件末
            exit
        else
            num=num+1
            if(temp=='    ')then  !!去除空白行计数
                num=num-1
            endif
        endif
    enddo conuter_number
    rewind(12)
    allocate(list(num,5))
    do i=1,num
        read(12,*) (list(i,j),j=1,5)    
    enddo
    close(12)
    list_sum1=0.
    list_sum2=0.
    list_sum3=0.
    list_sum4=0.
    list_sum5=0.
    list_sum1=sum(list(:,1))
    list_sum2=sum(list(:,2))
    list_sum3=sum(list(:,3))
    list_sum4=sum(list(:,4))
    list_sum5=sum(list(:,5))
    DEALLOCATE(list)

end subroutine test

!-------------------------------------------------!
!Writes accumulated data to the file 'test.dat'!
!-------------------------------------------------!
 subroutine writeaverage
!------------------------------!
 use systemvariables; implicit none

 integer :: num
 real :: list_sum1,list_sum2,list_sum3,list_sum4,list_sum5,average1,average2,average3,average4,average5 
 common/tests/list_sum1,list_sum2,list_sum3,list_sum4,list_sum5,num
 common/thevalues/average1,average2,average3,average4,average5
 integer :: j

 open(1,file='test.dat',status='unknown',position='append')
 average1=list_sum1/real(num)
 average2=list_sum2/real(num)
 average3=list_sum3/real(num)
 average4=list_sum4/real(num)
 average5=list_sum5/real(num)
 write(1,1)average1,average2,average3,average4,average5 
 1 format(5f18.12)
 close(1)

 end subroutine writeaverage
!----------------------------------------------------------------------------------!
!Writes the value of heat capacity,magnetic susceptibility to the file 'values.dat'!
!----------------------------------------------------------------------------------!
 subroutine values
!------------------------------!
 use systemvariables; implicit none

 real :: average1,average2,average3,average4,average5
 common/thevalues/average1,average2,average3,average4,average5
 real :: capacity,susceptibility,Binder_cumulan
 open(1,file='values.dat',status='unknown',position='append')
 capacity=(average2-(average1)**2)/(t**2)
 susceptibility=(average4-average3**2)*dble(n)/t
 Binder_cumulan=(1-average5/((average4**2)*3))
 write(1,1)capacity,susceptibility,Binder_cumulan
 1 format(4f18.12)
 close(1)
 end subroutine values
