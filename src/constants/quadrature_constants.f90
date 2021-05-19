module quadrature_constants

use set_precision, only : prec

implicit none

private

public :: w1, w2, w3, w4, w5
public :: a1, a2, a3, a4, a5

real(prec), dimension(1), parameter :: w1 = 2.0_prec
real(prec), dimension(2), parameter :: w2 = (/ 1.0_prec, 1.0_prec /)
real(prec), dimension(3), parameter :: w3 = (/ &
8.0_prec/9.0_prec, 5.0_prec/9.0_prec, 5.0_prec/9.0_prec /)
real(prec), dimension(4), parameter :: w4 = (/ &
(18.0_prec + sqrt(30.0_prec))/36.0_prec, &
(18.0_prec + sqrt(30.0_prec))/36.0_prec, &
(18.0_prec - sqrt(30.0_prec))/36.0_prec, &
(18.0_prec - sqrt(30.0_prec))/36.0_prec /)
real(prec), dimension(5), parameter :: w5 = (/ &
128.0_prec/225.0_prec, &
(322.0_prec + 13.0_prec*sqrt(70.0_prec))/900.0_prec, &
(322.0_prec + 13.0_prec*sqrt(70.0_prec))/900.0_prec, &
(322.0_prec - 13.0_prec*sqrt(70.0_prec))/900.0_prec, &
(322.0_prec - 13.0_prec*sqrt(70.0_prec))/900.0_prec /)


real(prec), dimension(1), parameter :: a1 = 0.0_prec
real(prec), dimension(2), parameter :: a2 = (/ &
1.0_prec/sqrt(3.0_prec), -1.0_prec/sqrt(3.0_prec) /)
real(prec), dimension(3), parameter :: a3 = (/ &
0.0_prec, sqrt(3.0_prec/5.0_prec), -sqrt(3.0_prec/5.0_prec) /)
real(prec), dimension(4), parameter :: a4 = (/ &
 sqrt( (3.0_prec/7.0_prec)-(2.0_prec/7.0_prec)*sqrt(6.0_prec/5.0_prec) ), &
-sqrt( (3.0_prec/7.0_prec)-(2.0_prec/7.0_prec)*sqrt(6.0_prec/5.0_prec) ), &
 sqrt( (3.0_prec/7.0_prec)+(2.0_prec/7.0_prec)*sqrt(6.0_prec/5.0_prec) ), &
-sqrt( (3.0_prec/7.0_prec)+(2.0_prec/7.0_prec)*sqrt(6.0_prec/5.0_prec) ) /)
real(prec), dimension(5), parameter :: a5 = (/ &
 0.0_prec, &
 (1.0_prec/3.0_prec)*sqrt( 5.0_prec-2.0_prec*sqrt(10.0_prec/7.0_prec) ), &
-(1.0_prec/3.0_prec)*sqrt( 5.0_prec-2.0_prec*sqrt(10.0_prec/7.0_prec) ), &
 (1.0_prec/3.0_prec)*sqrt( 5.0_prec+2.0_prec*sqrt(10.0_prec/7.0_prec) ), &
-(1.0_prec/3.0_prec)*sqrt( 5.0_prec+2.0_prec*sqrt(10.0_prec/7.0_prec) ) /)

end module quadrature_constants
