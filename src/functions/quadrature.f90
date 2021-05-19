module quadrature

use set_precision, only : prec
use set_constants, only : zero, one, two, three, four, half, fourth
use set_inputs, only : imax, i_low, i_high, ig_low, ig_high
use set_inputs, only : imax, i_low, i_high, ig_low, ig_high
use set_inputs, only : n_ghost
use soln_type, only : soln_t
use grid_type, only : grid_t

implicit none

private

public :: gauss_quad, gauss_pts

procedure( phi ), pointer :: fun

abstract interface

function phi(x,y)
  
  import :: prec
  real(prec) :: phi
  real(prec), intent(in) :: x, y
  
end function phi

end interface

contains

subroutine gauss_quad(P,xi,eta,w_xi,w_eta,fun,val)
  
  real(prec), dimension(4,2), intent(in) :: P
  real(prec), dimension(:), intent(in) :: xi, w_xi
  real(prec), dimension(:), intent(in) :: eta, w_eta
  real(prec), external :: fun
  real(prec), intent(out) :: val
  real(prec), dimension(size(xi),size(eta)) :: x, y, detJ
  integer :: i, j, i1, j1
  
  i1 = size(xi)
  j1 = size(eta)
  
  call bilinear_map(P,xi,eta,x,y,detJ)
  
  val = zero
  do i = 1,i1
    do j = 1,j1
      val = val + w_xi(i)*w_eta(j)*fun(x(i,j),y(i,j))*detJ(i,j)
    end do
  end do
  
end subroutine gauss_quad

subroutine bilinear_map(P,xi,eta,x,y,detJ)
  
  real(prec), dimension(4,2), intent(in) :: P
  real(prec), dimension(:), intent(in) :: xi, eta
  real(prec), dimension(size(xi),size(eta)), intent(out) :: x, y, detJ
  real(prec), dimension(4,4) :: A
  real(prec), dimension(4) :: alpha, beta
  real(prec), dimension(2,2) :: Jac
  integer :: i, j

A = fourth*reshape((/ one,-one,-one, one, &
                      one, one,-one,-one, &
                      one, one, one, one, &
                      one,-one, one,-one  /),shape(A))
alpha = matmul(A,P(:,1))
beta  = matmul(A,P(:,2))

do i = 1,size(xi)
  do j = 1,size(eta)
    x(i,j) = sum( alpha*(/ one, xi(i), eta(j), xi(i)*eta(j) /) )
    y(i,j) = sum(  beta*(/ one, xi(i), eta(j), xi(i)*eta(j) /) )
    Jac = transpose( matmul( reshape(                  &
          (/ -(1-eta(j)),1-eta(j),1+eta(j),-(1+eta(j)),&
             -(1-xi(i)),-(1+xi(i)),1+xi(i),1-xi(i) /), &
          (/ 2, 4 /) ), P ) )
    detJ(i,j) = abs( Jac(1,1)*Jac(2,2)-Jac(2,1)*Jac(1,2) )
  end do
end do

end subroutine bilinear_map

subroutine gauss_pts(xi,w,N)
  
  use quadrature_constants
  
  real(prec), dimension(:), intent(out) :: xi, w
  integer, intent(in) :: N
  
  select case(N)
  
  case(1)
    xi = a1
    w  = w1
  case(2)
    xi = a2
    w  = w2
  case(3)
    xi = a3
    w  = w3
  case(4)
    xi = a4
    w  = w4
  case(5)
    xi = a5
    w  = w5
  case default
  
  end select
  
end subroutine gauss_pts

end module quadrature
