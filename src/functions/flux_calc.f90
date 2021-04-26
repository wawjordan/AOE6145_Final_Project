module flux_calc
  
  use set_precision, only : prec
  use set_constants, only : zero, one, two, three, four, half, fourth
  use fluid_constants, only : gamma
  use set_inputs, only : neq, i_low, i_high, ig_low, ig_high, eps_roe
  use set_inputs, only : j_low, j_high, jg_low, jg_high
  use variable_conversion, only : cons2prim, speed_of_sound
  
  implicit none
  
  private
  
  public :: flux_fun, select_flux
  
  procedure( calc_flux ), pointer :: flux_fun
   
  abstract interface
    
  !================================ calc_flux ================================80
  !>
  !! Description:
  !!
  !! Inputs:      left_state  :
  !!              right_state : 
  !!
  !! Outputs:     F     :
  !<
  !===========================================================================80
  subroutine calc_flux(left_state, right_state, nx, ny, F)
    
    import :: prec, i_low, i_high, j_low, j_high, neq
    real(prec), dimension(:,:,:), intent(in) :: left_state, right_state
    real(prec), dimension(:,:), intent(in) :: nx, ny
    real(prec), dimension(i_low:i_high+1,j_low:j_high+1,neq), intent(out) :: F
    
  end subroutine calc_flux
    
  end interface
  
contains
  
  !================================ select_flux ==============================80
  !>
  !! Description:
  !<
  !===========================================================================80
  subroutine select_flux()
    
    use set_inputs, only : flux_scheme
    
    flux_fun => null()
    
    select case(flux_scheme)
    
    case(1)
      flux_fun => van_leer_flux
    case(2)
      flux_fun => roe_flux
    case default
    
    end select
  
  end subroutine select_flux
  
  !============================== van_leer_flux ==============================80
  !>
  !! Description:
  !!
  !! Inputs:      left  :
  !!              right : 
  !!
  !! Outputs:     F     :
  !<
  !===========================================================================80
  subroutine van_leer_flux(VL, VR, nx, ny, F)
    
    real(prec), dimension(:,:,:)               , intent(in)  :: VL, VR
    real(prec), dimension(:,:)                 , intent(in)  :: nx, ny
    real(prec), dimension(i_low:i_high+1,j_low:j_high+1,neq), intent(out) :: F
    real(prec), dimension(i_low:i_high+1,j_low:j_high+1)     :: aL, aR
    real(prec) :: rhoL, rhoR
    real(prec) ::   uL,   uR
    real(prec) ::   pL,   pR
    real(prec) ::   ML,   MR
    real(prec) ::  htL,  htR
    real(prec) :: M_plus, M_minus
    real(prec) :: p_plus, p_minus
    real(prec) :: c_plus, c_minus
    real(prec) :: d_plus, d_minus
    real(prec) :: alpha_plus, alpha_minus
    real(prec) :: beta_L, beta_R
    real(prec) :: Fc1, Fc2, Fc3, Fc4, Fp
    integer :: i, j
    
    call speed_of_sound(VL(:,:,4),VL(:,:,1),aL)
    call speed_of_sound(VR(:,:,4),VR(:,:,1),aR)
    do j = j_low-1,j_high
    do i = i_low-1,i_high
      rhoL = VL(i,j,1)
      rhoR = VR(i,j,1)
      uL   = VL(i,j,2)*nx(i,j) + VL(i,j,3)*ny(i,j)
      uR   = VR(i,j,2)*nx(i,j) + VR(i,j,3)*ny(i,j)
      pL   = VL(i,j,3)
      pR   = VR(i,j,3)
      
      ML = uL/aL(i,j)
      MR = uR/aR(i,j)
      M_plus  =  fourth*(ML+one)**2
      M_minus = -fourth*(MR-one)**2
      beta_L = -max(zero,one-int(abs(ML)))
      beta_R = -max(zero,one-int(abs(MR)))
      alpha_plus  = half*(one+sign(one,ML))
      alpha_minus = half*(one-sign(one,MR))
      c_plus  = alpha_plus*(one+beta_L)*ML - beta_L*M_plus
      c_minus = alpha_minus*(one+beta_R)*MR - beta_R*M_minus
      htL = aL(i,j)**2/(gamma-one) + half*uL**2
      htR = aR(i,j)**2/(gamma-one) + half*uR**2
      
      Fc1 = rhoL*aL(i,j)*c_plus + rhoR*aR(i,j)*c_minus
      Fc2 = rhoL*aL(i,j)*c_plus*VL(i,j,2) + rhoR*aR(i,j)*c_minus*VR(i,j,2)
      Fc3 = rhoL*aL(i,j)*c_plus*VL(i,j,3) + rhoR*aR(i,j)*c_minus*VR(i,j,3)
      Fc4 = rhoL*aL(i,j)*c_plus*htL + rhoR*aR(i,j)*c_minus*htR
      
      p_plus  = M_plus*(-ML + two)
      p_minus = M_minus*(-MR - two)
      d_plus  = alpha_plus*(one+beta_L) - beta_L*p_plus
      d_minus = alpha_minus*(one+beta_R) - beta_R*p_minus
      
      Fp = d_plus*pL + d_minus*pR
      
      F(i,j,1) = Fc1
      F(i,j,2) = Fc2 + Fp*nx(i,j)
      F(i,j,3) = Fc2 + Fp*ny(i,j)
      F(i,j,4) = Fc3
    end do
    end do
    
  end subroutine van_leer_flux
  
  !================================ roe_flux =================================80
  !>
  !! Description:
  !!
  !! Inputs:      left  :
  !!              right : 
  !!
  !! Outputs:     F     :
  !<
  !===========================================================================80
  subroutine roe_flux( VL, VR, nx, ny, F )
    
    real(prec), dimension(:,:,:)               , intent(in)  :: VL, VR
    real(prec), dimension(:,:)                 , intent(in)  :: nx, ny
    real(prec), dimension(i_low:i_high+1,j_low:j_high+1,neq), intent(out) :: F
    real(prec), dimension(i_low:i_high+1,j_low:j_high+1)     :: aL, aR
    real(prec), dimension(i_low:i_high+1,j_low:j_high+1)     :: R
    real(prec), dimension(4) :: FL, FR, rvec1, rvec2, rvec3, rvec4, lambda
    real(prec) ::  rhoL,  rhoR,  rho2
    real(prec) ::   uvL,   uvR,    u2
    real(prec) ::   vvL,   vvR,    v2
    real(prec) ::    pL,    pR,    a2
    real(prec) ::   htL,   htR,   ht2
    real(prec) :: uhatL, uhatR, uhat2
    real(prec) :: dw1, dw2, dw3, dw4
    integer :: i, j
    
    call speed_of_sound(VL(:,:,3),VL(:,:,1),aL)
    call speed_of_sound(VR(:,:,3),VR(:,:,1),aR)
    
    do j = j_low-1,j_high
    do i = i_low-1,i_high
      rhoL = VL(i,j,1)
      rhoR = VR(i,j,1)
      uvL  = VL(i,j,2)
      uvR  = VR(i,j,2)
      vvL  = VL(i,j,3)
      vvR  = VR(i,j,3)
      uhatL = uvL*nx(i,j) + vvL*ny(i,j)
      uhatR = uvR*nx(i,j) + vvR*ny(i,j)
      
      pL   = VL(i,j,4)
      pR   = VR(i,j,4)
      htL  = aL(i,j)**2/(gamma-one) + half*(uvL**2 + vvL**2)
      htR  = aR(i,j)**2/(gamma-one) + half*(uvR**2 + vvR**2)
      
      R(i,j) = sqrt(rhoR/rhoL)
      rho2 = R(i,j)*rhoL
      u2   = (R(i,j)*uvR + uvL)/(R(i,j) + one)
      v2   = (R(i,j)*vvR + vvL)/(R(i,j) + one)
      ht2  = (R(i,j)*htR + htL)/(R(i,j) + one)
      a2   = sqrt((gamma-one)*(ht2 - half*(u2**2 + v2**2)))
      uhat2 = u2*nx(i,j) + v2*ny(i,j)
      
      lambda = (/ uhat2, uhat2, uhat2 + a2, uhat2 - a2 /)
      
      rvec1 = (/ one, u2, v2, half*(u2**2 +v2**2) /)
      rvec2 = (/ zero, ny(i,j)*rho2, nx(i,j)*rho2, &
                     rho2*(ny(i,j)*u2-nx(i,j)*v2) /)
      rvec3 = half*(rho2/a2)*(/ one, u2+nx(i,j)*a2,&
                      v2+ny(i,j)*a2, ht2+uhat2*a2 /)
      rvec4 = half*(rho2/a2)*(/ one, u2-nx(i,j)*a2,&
                      v2-ny(i,j)*a2, ht2-uhat2*a2 /)
      
      lambda = abs(lambda)
      lambda = half*(one+sign(one,lambda-two*eps_roe*a2))*lambda &
           & + half*(one-sign(one,lambda-two*eps_roe*a2))*&
           & (lambda**2/(four*eps_roe*a2) + eps_roe*a2)
      
      dw1 = (rhoR - rhoL) - (pR - pL)/a2**2
      dw2 = ny(i,j)*(uvR - uvL) - nx(i,j)*(vvR - vvL)
      dw3 = nx(i,j)*(uvR - uvL) + ny(i,j)*(vvR - vvL) + (pR - pL)/(rho2*a2)
      dw4 = nx(i,j)*(uvR - uvL) + ny(i,j)*(vvR - vvL) - (pR - pL)/(rho2*a2)
     
      FL = (/ rhoL*uhatL,                  &
              rhoL*uvL*uhatL + pL*nx(i,j), &
              rhoL*vvL*uhatL + pR*ny(i,j), &
              rhoL*htL*uhatL /)
      FR = (/ rhoR*uhatR,                  &
              rhoR*uvR*uhatR + pR*nx(i,j), &
              rhoR*vvR*uhatR + pR*ny(i,j), &
              rhoR*htR*uhatR /)
      
      F(i,j,:) = half*( (FL+FR)      &
               - lambda(1)*dw1*rvec1 &
               - lambda(2)*dw2*rvec2 &
               - lambda(3)*dw3*rvec3 &
               - lambda(4)*dw4*rvec4 )
    end do
    end do
  end subroutine roe_flux
  
end module flux_calc
