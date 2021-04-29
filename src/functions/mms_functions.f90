module mms_functions
contains
!=============================================================================80
!!!!MMS functions
pure elemental function rho_mms(length, x, y)

  use set_precision, only : prec
  use set_constants, only : pi, two
  use mms_constants, only : rho0, rhox, rhoy
  
  implicit none

  real(prec), intent(in) :: x
  real(prec), intent(in) :: y
  real(prec), intent(in) :: length
  real(prec) :: rho_mms

  rho_mms = rho0 + rhoy*cos((pi*y)/(two*length)) + rhox*sin((pi*x)/length)

end function rho_mms
!=============================================================================80
pure elemental function uvel_mms(length, x, y)

  use set_precision, only : prec
  use set_constants, only : pi, two, three, five
  use mms_constants, only : uvel0, uvelx, uvely
  
  implicit none

  real(prec), intent(in) :: x
  real(prec), intent(in) :: y
  real(prec), intent(in) :: length
  real(prec) :: uvel_mms
  
  uvel_mms = uvel0 + uvely*cos((three*pi*y)/(five*length)) +                   &
                 uvelx*sin((three*pi*x)/(two*length))

end function uvel_mms
!=============================================================================80
pure elemental function vvel_mms(length,x,y)

  use set_precision, only : prec
  use set_constants, only : pi, two, three
  use mms_constants, only : vvel0, vvelx, vvely
  
  implicit none

  real(prec), intent(in)  :: x
  real(prec), intent(in)  :: y
  real(prec), intent(in)  :: length
  real(prec) :: vvel_mms

  vvel_mms = vvel0 + vvelx*cos((pi*x)/(two*length)) +                          &
                 vvely*sin((two*pi*y)/(three*length))
end function vvel_mms
!=============================================================================80
pure elemental function press_mms(length,x,y)

  use set_precision, only : prec
  use set_constants, only : pi, two
  use mms_constants, only : press0, pressx, pressy
  
  implicit none

  real(prec), intent(in) :: x
  real(prec), intent(in) :: y
  real(prec), intent(in) :: length
  real(prec) :: press_mms

  press_mms = press0 + pressx*cos((two*pi*x)/length) + pressy*sin((pi*y)/length)

end function press_mms
!=============================================================================80
pure elemental function rmassconv(length,x,y)

  use set_precision, only : prec
  use set_constants, only : pi, two, three, five
  use mms_constants, only : rho0, rhox, rhoy, uvel0, uvelx, uvely,             &
                            vvel0, vvelx, vvely, press0, pressx, pressy
  
  implicit none

  real(prec), intent(in) :: x
  real(prec), intent(in) :: y
  real(prec), intent(in) :: length
  real(prec) :: rmassconv

  rmassconv = (three*pi*uvelx*cos((three*pi*x)/(two*length)) *                 &
    (rho0 + rhoy*cos((pi*y)/(two*length)) + rhox*sin((pi*x)/length))) /        &
    (two*length) + (two*pi*vvely*cos((two*pi*y)/(three*length)) *              &
    (rho0 + rhoy*cos((pi*y)/(two*length)) + rhox*sin((pi*x)/length))) /        &
    (three*length) + (pi*rhox*cos((pi*x)/length) *                             &
    (uvel0 + uvely*cos((three*pi*y)/(five*length)) + uvelx*sin((three*pi*x)/   &
    (two*length))))/length - (pi*rhoy*sin((pi*y)/(two*length)) *               &
    (vvel0 + vvelx*cos((pi*x)/(two*length)) + vvely*sin((two*pi*y) /           &
    (three*length))))/(two*length)

end function rmassconv
!=============================================================================80
pure elemental function xmtmconv(length,x,y)

  use set_precision, only : prec
  use set_constants, only : pi, two, three, five
  use mms_constants, only : rho0, rhox, rhoy, uvel0, uvelx, uvely,             &
                            vvel0, vvelx, vvely, press0, pressx, pressy
  
  implicit none

  real(prec), intent(in) :: x
  real(prec), intent(in) :: y
  real(prec), intent(in) :: length
  real(prec) :: xmtmconv

  xmtmconv = (three*pi*uvelx*cos((three*pi*x)/(two*length)) *                  &
    (rho0 + rhoy*cos((pi*y)/(two*length)) + rhox*sin((pi*x)/length)) *         &
    (uvel0 + uvely*cos((three*pi*y)/(five*length)) +                           &
    uvelx*sin((three*pi*x)/(two*length))))/length +                            &
    (two*pi*vvely*cos((two*pi*y) /                                             &
    (three*length))*(rho0 + rhoy*cos((pi*y)/(two*length)) +                    &
    rhox*sin((pi*x)/length))*(uvel0 + uvely*cos((three*pi*y) /                 &
    (five*length)) + uvelx*sin((three*pi*x)/(two*length))))/(three*length) +   &
    (pi*rhox*cos((pi*x)/length)*(uvel0 + uvely*cos((three*pi*y) /              &
    (five*length)) + uvelx*sin((three*pi*x)/(two*length)))**2)/length -        &
    (two*pi*pressx*sin((two*pi*x)/length))/length -                            &
    (pi*rhoy*(uvel0 + uvely*cos((three*pi*y)/(five*length)) +                  &
    uvelx*sin((three*pi*x)/(two*length)))*sin((pi*y)/(two*length))*            &
    (vvel0 + vvelx*cos((pi*x)/(two*length)) +                                  &
    vvely*sin((two*pi*y)/(three*length))))/(two*length) -                      &
    (three*pi*uvely*(rho0 + rhoy*cos((pi*y)/(two*length)) +                    &
    rhox*sin((pi*x)/length))*sin((three*pi*y)/(five*length))*(vvel0 + vvelx *  &
    cos((pi*x)/(two*length)) + vvely*sin((two*pi*y)/(three*length)))) /        &
    (five*length)

end function xmtmconv
!=============================================================================80
pure elemental function ymtmconv(length,x,y)

  use set_precision, only : prec
  use set_constants, only : pi, two, three, four, five
  use mms_constants, only : rho0, rhox, rhoy, uvel0, uvelx, uvely,             &
                            vvel0, vvelx, vvely, press0, pressx, pressy
  
  implicit none

  real(prec), intent(in) :: x
  real(prec), intent(in) :: y
  real(prec), intent(in) :: length
  real(prec) :: ymtmconv

  ymtmconv = (pi*pressy*cos((pi*y)/length))/length -                           &
    (pi*vvelx*sin((pi*x)/(two*length))*(rho0 + rhoy*cos((pi*y)/(two*length)) + &
    rhox*sin((pi*x)/length))*(uvel0 + uvely*cos((three*pi*y)/(five*length)) +  &
    uvelx*sin((three*pi*x)/(two*length))))/(two*length) +                      &
    (three*pi*uvelx*cos((three*pi*x)/(two*length)) *                           &
    (rho0 + rhoy*cos((pi*y)/(two*length)) + rhox*sin((pi*x)/length)) *         &
    (vvel0 + vvelx*cos((pi*x)/(two*length)) +                                  &
    vvely*sin((two*pi*y)/(three*length))))/(two*length) +                      &
    (four*pi*vvely*cos((two*pi*y) /                                            &
    (three*length))*(rho0 + rhoy*cos((pi*y)/(two*length)) +                    &
    rhox*sin((pi*x)/length))*(vvel0 + vvelx*cos((pi*x)/(two*length)) +         &
    vvely*sin((two*pi*y)/(three*length))))/(three*length) +                    &
    (pi*rhox*cos((pi*x)/length) *                                              &
    (uvel0 + uvely*cos((three*pi*y)/(five*length)) +                           &
    uvelx*sin((three*pi*x)/(two*length))) *                                    &
    (vvel0 + vvelx*cos((pi*x)/(two*length)) +                                  &
    vvely*sin((two*pi*y)/(three*length))))/length -                            &
    (pi*rhoy*sin((pi*y)/(two*length)) *                                        &
    (vvel0 + vvelx*cos((pi*x)/(two*length)) +                                  &
    vvely*sin((two*pi*y)/(three*length)))**2)/(two*length)

end function ymtmconv
!=============================================================================80
pure elemental function energyconv(gamma, length,x,y)

  use set_precision, only : prec
  use set_constants, only : pi, one, two, three, four, five, six
  use mms_constants, only : rho0, rhox, rhoy, uvel0, uvelx, uvely,             &
                            vvel0, vvelx, vvely, wvel0, wvelx, wvely,          &
                            press0, pressx, pressy
  
  implicit none

  real(prec), intent(in) :: gamma
  real(prec), intent(in) :: length
  real(prec), intent(in) :: x
  real(prec), intent(in) :: y
  real(prec) :: energyconv

  energyconv = (uvel0 + uvely*cos((three*pi*y)/(five*length)) +                &
    uvelx*sin((three*pi*x)/(two*length)))*((-two*pi*pressx*sin((two*pi*x) /    &
    length))/length + (rho0 + rhoy*cos((pi*y)/(two*length)) +                  &
    rhox*sin((pi*x)/length))*((-two*pi*pressx*sin((two*pi*x)/length))/         &
    ((-one + gamma)*length*(rho0 + rhoy*cos((pi*y)/(two*length)) +             &
    rhox*sin((pi*x)/length))) + ((three*pi*uvelx*cos((three*pi*x) /            &
    (two*length))*(uvel0 + uvely*cos((three*pi*y)/(five*length)) +             &
    uvelx*sin((three*pi*x)/(two*length))))/length - (pi*vvelx*sin((pi*x) /     &
    (two*length))*(vvel0 + vvelx*cos((pi*x)/(two*length)) +                    &
    vvely*sin((two*pi*y)/(three*length))))/length)/two - (pi*rhox*cos((pi*x) / &
    length)*(press0 + pressx*cos((two*pi*x)/length) +                          & 
    pressy*sin((pi*y)/length)))/((-one + gamma)*length*(rho0 + rhoy*cos((pi*y)/&
    (two*length)) + rhox*sin((pi*x)/length))**2)) +                            &
    (pi*rhox*cos((pi*x)/length)*((wvel0**2 + (uvel0 + uvely*cos((three*pi*y) / &
    (five*length)) + uvelx*sin((three*pi*x)/(two*length)))**2 +                &
    (vvel0 + vvelx*cos((pi*x)/(two*length)) + vvely*sin((two*pi*y) /           &
    (three*length)))**2)/two + (press0 + pressx*cos((two*pi*x)/length) +       &
    pressy*sin((pi*y)/length))/((-one + gamma) *                               &
    (rho0 + rhoy*cos((pi*y)/(two*length)) +                                    &
    rhox*sin((pi*x)/length)))))/length) +                                      &
    (three*pi*uvelx*cos((three*pi*x)/(two*length)) *                           &
    (press0 + pressx*cos((two*pi*x)/length) + pressy*sin((pi*y)/length) +      &
    (rho0 + rhoy*cos((pi*y)/(two*length)) + rhox*sin((pi*x)/length))*          &
    ((wvel0**2 + (uvel0 + uvely*cos((three*pi*y)/(five*length)) +              &
    uvelx*sin((three*pi*x)/(two*length)))**2 +                                 &
    (vvel0 + vvelx*cos((pi*x)/(two*length)) +                                  &
    vvely*sin((two*pi*y)/(three*length)))**2)/two +                            &
    (press0 + pressx*cos((two*pi*x)/length) +                                  &
    pressy*sin((pi*y)/length))/((-one + gamma) *                               &
    (rho0 + rhoy*cos((pi*y)/(two*length)) +                                    &
    rhox*sin((pi*x)/length))))))/(two*length) +                                &
    (two*pi*vvely*cos((two*pi*y)/(three*length)) *                             &
    (press0 + pressx*cos((two*pi*x)/length) +                                  &
    pressy*sin((pi*y)/length) + (rho0 + rhoy*cos((pi*y)/(two*length)) +        &
    rhox*sin((pi*x)/length))*((wvel0**2 +                                      &
    (uvel0 + uvely*cos((three*pi*y)/(five*length)) +                           &
    uvelx*sin((three*pi*x)/(two*length)))**2 +                                 &
    (vvel0 + vvelx*cos((pi*x)/(two*length)) +                                  &
    vvely*sin((two*pi*y)/(three*length)))**2)/two +                            &
    (press0 + pressx*cos((two*pi*x)/length) + pressy*sin((pi*y)/length)) /     &
    ((-one + gamma)*(rho0 + rhoy*cos((pi*y)/(two*length)) +                    &
    rhox*sin((pi*x)/length))))))/(three*length) + (vvel0 + vvelx*cos((pi*x) /  &
    (two*length)) + vvely*sin((two*pi*y)/(three*length))) *                    &
    ((pi*pressy*cos((pi*y)/length))/length - (pi*rhoy*sin((pi*y)/(two*length))*&
    ((wvel0**2 + (uvel0 + uvely*cos((three*pi*y)/(five*length)) +              &
    uvelx*sin((three*pi*x)/(two*length)))**2 + (vvel0 + vvelx *                &
    cos((pi*x)/(two*length)) + vvely*sin((two*pi*y)/(three*length)))**2)/two + &
    (press0 + pressx*cos((two*pi*x)/length) +                                  &
    pressy*sin((pi*y)/length))/((-one + gamma) *                               &
    (rho0 + rhoy*cos((pi*y)/(two*length)) +                                    &
    rhox*sin((pi*x)/length)))))/(two*length) +                                 &
    (rho0 + rhoy*cos((pi*y)/(two*length)) +                                    &
    rhox*sin((pi*x)/length))*((pi*pressy*cos((pi*y)/length)) /                 &
    ((-one + gamma)*length*(rho0 + rhoy*cos((pi*y)/(two*length)) +             &
    rhox*sin((pi*x)/length))) +                                                &
    ((-six*pi*uvely*(uvel0 + uvely*cos((three*pi*y) /                          &
    (five*length)) + uvelx*sin((three*pi*x)/(two*length))) *                   &
    sin((three*pi*y)/(five*length)))/(five*length) +                           &
    (four*pi*vvely*cos((two*pi*y) /                                            &
    (three*length))*(vvel0 + vvelx*cos((pi*x)/(two*length)) +                  &
    vvely*sin((two*pi*y)/(three*length))))/(three*length))/two +               &
    (pi*rhoy*sin((pi*y)/(two*length))*(press0 + pressx*cos((two*pi*x)/length) +&
    pressy*sin((pi*y)/length)))/(two*(-one + gamma)*length*                    &
    (rho0 + rhoy*cos((pi*y)/(two*length)) + rhox*sin((pi*x)/length))**2)))

end function energyconv

end module mms_functions
