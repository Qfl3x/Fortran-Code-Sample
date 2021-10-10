program deform_energy
  
  implicit none

  integer,parameter :: dp=kind(1.6) !Used to explicitly fix the precision
                                    !throughout the program, avoiding address issues.
  integer :: i
  real(kind=dp),dimension(11) :: w_arr11,x_arr11
  real(kind=dp),dimension(10) :: w_arr10,x_arr10

  open(1,file='gauss11.dat',status='old')  !n=11 Gaussian Quadrature (Legendre)
  open(2,file='gauss10.dat',status='old')   !n=10 Gaussian Quadrature (Legendre)
  do i=1,11
     read(1,*) w_arr11(i),x_arr11(i)
  end do
  do i=1,10
     read(2,*) w_arr10(i),x_arr10(i)
  end do
  
  close(1)
  close(2)
  call c_var(236,92,w_arr11,x_arr11,w_arr10,x_arr10,"U236.dat")
  call c_var(240,94,w_arr11,x_arr11,w_arr10,x_arr10,"Pu240.dat")


end program deform_energy

real(kind=kind(1.6)) function sum_func(arr,init,step)  !Function for summing over the odd/even elements of an array
  implicit none
  integer, parameter :: dp = kind(1.6)
  integer :: init,step,i
  real(kind=dp),dimension(501) :: arr
  sum_func=0

  do i=init,500,step
     sum_func=sum_func+arr(i)
  end do
end function sum_func

function gauss_quad(func,typ,params,w_arr,x_arr,w_arr10,x_arr10)   !Function calculating the integral using Legendre-Gauss
  !func: Integrand.
  !typ: integer parameter to differentiate between different integrals: (Because of different number of parameters)
  !      - 0 for B_s
  !      - 1 for B_c
  !      - 2 for Phi
  !params: the extra parameters to be input to 'func'
  !w_arr,x_arr: weights and nodes lists (resp.) for the n=11 Legendre-Gauss
  !w_arr10,x_arr10: weights and nodes lists (resp.) for the n=10 Legendre-Gauss
  !Note: w_arr/10 x_arr/10 are sometimes put in as dummy variables when not needed, these are only needed
  ! explicitly for E_c where Phi needs the n=10 parameters (as explained previously)
  !=== func syntax ===!
  ! func(x,[params],[w_arr,x_arr]) x:=> looped parameter
  ! [params]:=> extra paremeters for func
  ! [w_arr,x_arr]:=> w_arr and x_arr arrays passed to next function for B_c/E_c
  implicit none
  integer, parameter :: dp=kind(1.6)
  real(kind=dp), external :: func
  real(kind=dp),dimension(10) :: params
  real(kind=dp),dimension(10),intent(in) :: w_arr10,x_arr10
  integer :: i,typ
  real(kind=dp),dimension(11),intent(in) :: w_arr,x_arr
  real(kind=dp) :: gauss_quad
  !typ: 0 -> B_s ; 1 -> B_c 2-> Phi
  if (typ .EQ. 0) then ! B_s Only w_arr,x_arr are relevant.
     do i=1,11
        gauss_quad=gauss_quad+w_arr(i)*func(x_arr(i),params(1),params(2))
     end do
  elseif (typ .EQ. 1) then !B_c/E_c ALL 4 ARRAYS ARE RELEVANT
     do i=1,11
        gauss_quad=gauss_quad+w_arr(i)*func(x_arr(i),params(1),params(2), &
             params(3),params(4),w_arr10,x_arr10)
     end do
  elseif (typ .EQ. 2) then !Phi only w_arr10,x_arr10 are relevant
     do i=1,10
        gauss_quad=gauss_quad+w_arr10(i)*func(x_arr10(i),params(1),params(2), &
             params(3),params(4),params(5),params(6))
     end do
  else
     do i=1,11
        gauss_quad=gauss_quad+w_arr(i)*func(x_arr(i))
     end do

  end if
end function gauss_quad



real(kind=kind(1.6)) function F(a,b) !First kind elliptical integral
  implicit none

  integer, parameter :: dp = kind(1.6)
  real(kind=dp) :: a,b
  real(kind=dp),parameter :: pi=4.0*atan(1.0)
  real(kind=dp) :: linspace_seq=pi/1000
  real(kind=dp),dimension(501)::F_arr,phi_arr
  real(kind=dp),dimension(10) :: F_in
  real(kind=dp) :: sum_func,k

  integer :: i
  phi_arr(1)=0.0
  a=abs(a)
  b=abs(b)
  k=1-b/a
  do i=2,501
     phi_arr(i)=phi_arr(i-1)+linspace_seq
  end do
  phi_arr(501)=pi/2

  F_arr=(sqrt(1-(k)*(sin(phi_arr))**2.))**(-1.0)

  F=F_arr(1)+F_arr(501)+4.0*sum_func(F_arr,2,2)+2.0*sum_func(F_arr,3,2)
  F=F*(linspace_seq/3.)/sqrt(a) 


end function F


real(kind=kind(1.6)) function E(a,b) !Second kind elliptical integral
  implicit none
  integer, parameter :: dp = kind(1.6)
  real(kind=dp) :: a,b
  real(kind=dp),parameter :: pi=4.0*atan(1.0)
  real(kind=dp) :: linspace_seq=pi/1000
  real(kind=dp),dimension(501)::E_arr,phi_arr,E_in
  integer :: i
  real(kind=dp) :: sum_func,k
  phi_arr(1)=0.0
  a=abs(a)
  b=abs(b)

  do i=2,501
     phi_arr(i)=phi_arr(i-1)+linspace_seq
  end do
  phi_arr(501)=pi/2.

  E_arr=sqrt(a-(a-b)*(sin(phi_arr))**2.)

  E=E_arr(1)+E_arr(501)+4.*sum_func(E_arr,2,2)+2.*sum_func(E_arr,3,2)
  E=E*((linspace_seq)/(3)) 
end function E

real(kind=(kind(1.6))) function Phi_in(x,u,v,c,R0,A,B)
  !This function represents the integrand of Phi. This function DEPENDS on F and E
  integer, parameter :: dp = kind(1.6)
  real(kind=dp) :: u,v,c,R0,A,B,F,E,sum_func,der,x,vp
  real(kind=dp),dimension(10) :: Phi_arr,der_arr
  real(kind=dp) ::linspace_seq,aa,bb



  der=(1-x**2.)*(2*B*x)-2*x*(A+B*x**2.)
  vp=sqrt((1-x**2.)*(A+B*x**2.))
  Phi_in=vp**2.-x**2.-(x-u)**2. - (x-u)*der

  aa=((x-u)**2.+(vp-v)**2.)
  bb=((x-u)**2.+(vp+v)**2.)
  Phi_in=Phi_in*F(aa,bb) + E(aa,bb)

end function Phi_in


real(kind=kind(1.6)) function E_c_in(x,A,B,c,R0,w_arr,x_arr)
  !Integrand of E_c, DEPENDS on Phi_in (subsequently F and E)
  implicit none
  integer, parameter :: dp = kind(1.6)
  integer :: At
  real(kind=dp),external :: gauss_quad,Phi_in
  real(kind=dp) :: x,v,A,B,c,R0,der
  real(kind=dp),dimension(10) :: params
  real(kind=dp),dimension(10),intent(in) :: w_arr,x_arr

  der=B*x*(1-x**2.)-x*(A+B*x**2.)
  v=(1-x**2.)*(A+B*(x**2.))
  params(1)=x
  params(2)=v
  params(3)=c
  params(4)=R0
  params(5)=A
  params(6)=B


  E_c_in=gauss_quad(Phi_in,2,params,[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.], &
       [1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],w_arr,x_arr)*(v**2.-x*der)
  !write(*,*) "E_c", E_c_in
end function E_c_in

real(kind=kind(1.6)) function B_c(At,c,w_arr11,x_arr11,w_arr10,x_arr10)
  !B_c calculator, DEPENDS on E_c and its dependencies.
  implicit none
  integer,parameter :: dp = kind(1.6)
  real(kind=dp) :: c,A,B,R0,v,der
  real(kind=dp),dimension(10) :: params
  integer :: At
  real(kind=dp),external :: gauss_quad,E_c_in
  real(kind=dp),dimension(11),intent(in) :: w_arr11,x_arr11
  real(kind=dp),dimension(10),intent(in) :: w_arr10,x_arr10
  real(kind=dp) :: pi=4.0*atan(1.0)


  R0=0.53E-10*(At**(1./3.))
  A=c**(-3)-0.1*(c-1)
  B=0.5*(c-1)


  !print*, c,At
  !print*, x_arr10

  params(1)=A
  params(2)=B
  params(3)=c
  params(4)=R0
  B_c=gauss_quad(E_c_in,1,params,w_arr11,x_arr11,w_arr10,x_arr10)
  B_c=B_c*(3.*(c**5.))/(5.*8.*pi )
  print*,"B_c",B_c
end function B_c

real(kind=kind(1.6)) function E_s_in(x,A,B)
  implicit none
  integer,parameter :: dp=kind(1.6)
  real(kind=dp) ::x,A,B,v,der

  v=sqrt((1-x**2.)*(A+B*(x**2.)))

  der=B*x*(1-x**2.)-x*(A+B*(x**2.))
  E_s_in=(v**2. + der**2.)**(0.5)
return

end function E_s_in

real(kind=kind(1.6)) function B_s(At,c,w_arr,x_arr)
  integer, parameter :: dp = kind(1.6)
  real(kind=dp):: pi=4.0*atan(1.0)
  real(kind=dp) :: c,A,B
  real(kind=dp) :: sum_func
  real(kind=dp),dimension(11),intent(in) :: w_arr,x_arr
  real(kind=dp),dimension(10) :: params
  real(kind=dp),external :: gauss_quad,E_s_in

  integer :: At
  A=c**(-3.)-0.1*(c-1)
  B=0.5*(c-1)
  params(1)=A
  params(2)=B
  B_s=gauss_quad(E_s_in,0,params,w_arr,x_arr,[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.] ,[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.])
  B_s=B_s*((c**2.)/2.)
  write(*,*)"B_S", B_s

end function B_s
subroutine c_var(At,Z,w_arr11,x_arr11,w_arr10,x_arr10,files)
  implicit none
  integer,parameter  :: dp=kind(1.6)
  character(len=50) :: files
  real(kind=dp) :: E_c0,E_s0,e,pi,R,Chi,c,Eldm,B_c,B_s,tes
  real(kind=dp),dimension(10)::w_arr10,x_arr10
  real(kind=dp),dimension(11)::w_arr11,x_arr11
  integer :: At,Z,i
  pi=4.0*atan(1.0);e=1.6E-19
  R=0.53E-10*(At**(1.0/3))

  E_c0=9.*(Z*e)**2. /R
  E_s0=14.*(At**(2.0/3))
  Chi=(Z**2.)/At
  Chi=Chi/(45.)
  c=1.0
  print*, E_c0,E_s0,Chi
  
  open(1,file=files,status='new')
  do i=1,17
     print*,c,i
     Eldm=E_s0*(B_s(At,c,w_arr11,x_arr11)-1+2*Chi*(B_c(240,c,w_arr11,x_arr11,w_arr10,x_arr10)-1))
     write(1,*) c,Eldm
     print*,"done1"
     c=c+0.05
  end do
  print*, "done2"
  !close(1)
end subroutine c_var
