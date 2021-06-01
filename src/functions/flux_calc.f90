module flux_calc
  
  use set_precision, only : prec
  use set_constants, only : zero, one, two, three, four, half, fourth
  use fluid_constants, only : gamma
  use set_inputs, only : neq, i_low, i_high, ig_low, ig_high, eps_roe
  use set_inputs, only : j_low, j_high, jg_low, jg_high, n_ghost
  use variable_conversion, only : cons2prim, speed_of_sound
  use other_subroutines, only : MUSCL_extrap
  use soln_type, only : soln_t
  use grid_type, only : grid_t
  
  implicit none
  
  private
  
  public :: flux_fun, select_flux, exact_flux, calc_flux_2D
  
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
    real(prec), dimension(:), intent(in) :: left_state, right_state
    real(prec), intent(in) :: nx, ny
    real(prec), dimension(:), intent(out) :: F
    
  end subroutine calc_flux
    
  end interface
  
contains

subroutine calc_flux_2D(soln,grid,Fnormal)

  type(soln_t), intent(in) :: soln
  type(grid_t), intent(in) :: grid
  real(prec), dimension(i_low:i_high+1,j_low:j_high+1,neq,2), &
                                                  intent(inout) :: Fnormal
  !real(prec), dimension(:,:,:,:), intent(inout) :: Fnormal
  real(prec), dimension(neq) :: left, right
  real(prec) :: nx, ny
  integer :: i,j,k
  integer, dimension(4) :: ind
  real(prec), dimension(4,neq) :: Vtmp, psiPtmp, psiMtmp
  Fnormal(:,:,:,:) = zero
  !write(*,*) lbound(Fnormal,1), ubound(Fnormal,1)
  !write(*,*) lbound(Fnormal,2), ubound(Fnormal,2)
  !write(*,*) lbound(Fnormal,3), ubound(Fnormal,3)
  !write(*,*) lbound(Fnormal,4), ubound(Fnormal,4)
  !stop
  !write(*,*) i_low,i_high,j_low,j_high 
  !call MUSCL_extrap(soln%V,soln%psi_plus,soln%psi_minus,left,right,ind1,ind2)
  
!  write(*,*) left(1,1,:,1)

  !write(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  do j = j_low,j_high+1
  !write(*,*) j
  do i = i_low,i_high+1
    nx = one
    ny = zero
    !nx = grid%n_xi(i,j,1)
    !ny = grid%n_xi(i,j,2)
    !ind = (/ ( k,k=i-2,i+1 ) /)
    !Vtmp = soln%V(ind,j,:)
    !psiPtmp = soln%psi_plus(ind,j,:,1)
    !psiMtmp = soln%psi_minus(ind,j,:,1)
    !call MUSCL_extrap(Vtmp, psiPtmp, psiMtmp, left, right)
    left =  soln%V(i,j,:)
    right = soln%V(i+1,j,:)
    call flux_fun(left,right,nx,ny,Fnormal(i,j,:,1))
    !write(*,*) Fnormal(i,j,:,1)
  end do
  end do
  !write(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  do j = j_low,j_high+1
  !write(*,*) j
  do i = i_low,i_high+1
    nx = zero
    ny = one
    !nx = grid%n_eta(i,j,1)
    !ny = grid%n_eta(i,j,2)
    !ind = (/ ( k,k=j-2,j+1 ) /)
    !Vtmp = soln%V(i,ind,:)
    !psiPtmp = soln%psi_plus(i,ind,:,2)
    !psiMtmp = soln%psi_minus(i,ind,:,2)
    !call MUSCL_extrap(Vtmp, psiPtmp, psiMtmp, left, right)
    left =  soln%V(i,j,:)
    right = soln%V(i,j+1,:)
    call flux_fun(left,right,nx,ny,Fnormal(i,j,:,2))
    !write(*,*) Fnormal(i,j,:,2)
  end do
  end do
  !stop
  
end subroutine calc_flux_2D


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
  
  subroutine exact_flux(V, nx, ny, F)
    
    real(prec), dimension(neq), intent(in) :: V
    real(prec), intent(in) :: nx, ny
    real(prec), dimension(neq), intent(out) :: F
    real(prec) :: rho, uvel, vvel, p, a, ht, un
    
    rho  = V(1)
    uvel = V(2)
    vvel = V(3)
    p    = V(4)
    call speed_of_sound(p,rho,a)
    un = uvel*nx + vvel*ny
    ht = a**2/(gamma-one) + half*un**2
    F(1) = rho*un
    F(2) = rho*uvel*un + nx*p
    F(3) = rho*vvel*un + ny*p
    F(4) = rho*ht*un
    
  end subroutine exact_flux
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
  subroutine van_leer_flux(left, right, nx, ny, F)
    
    real(prec), dimension(:), intent(in)  :: left, right
    real(prec), dimension(:), intent(out) :: F
    real(prec), intent(in) :: nx, ny
    real(prec) :: rhoL, uL, vL, pL, aL, htL, unL, ML, beta_L
    real(prec) :: rhoR, uR, vR, pR, aR, htR, unR, MR, beta_R
    real(prec) :: M_plus , alpha_plus , c_plus , p_plus , d_plus
    real(prec) :: M_minus, alpha_minus, c_minus, p_minus, d_minus
    real(prec) :: Fp
    
    rhoL = left(1)
    uL   = left(2)
    vL   = left(3)
    pL   = left(4)
    call speed_of_sound(pL,rhoL,aL) 
    rhoR = right(1)
    uR   = right(2)
    vR   = right(3)
    pR   = right(4)
    call speed_of_sound(pR,rhoR,aR)
    
    unL = uL*nx + vL*ny
    unR = uR*nx + vR*ny
    
    htL = aL**2/(gamma-one) + half*unL**2
    htR = aR**2/(gamma-one) + half*unR**2
    
    ML = unL/aL
    MR = unR/aR
    M_plus  =  fourth*(ML+one)**2
    M_minus = -fourth*(MR-one)**2
    
    beta_L = -max(zero,one-int(abs(ML)))
    beta_R = -max(zero,one-int(abs(MR)))
    alpha_plus  = half*(one+sign(one,ML))
    alpha_minus = half*(one-sign(one,MR))
    c_plus  = alpha_plus*(one+beta_L)*ML - beta_L*M_plus
    c_minus = alpha_minus*(one+beta_R)*MR - beta_R*M_minus
    
    p_plus  = M_plus*(-ML + two)
    p_minus = M_minus*(-MR - two)
    d_plus  = alpha_plus*(one+beta_L) - beta_L*p_plus
    d_minus = alpha_minus*(one+beta_R) - beta_R*p_minus
    
    Fp = d_plus*pL + d_minus*pR
    
    F(1) = rhoL*aL*c_plus     + rhoR*aR*c_minus
    F(2) = rhoL*aL*c_plus*uL  + rhoR*aR*c_minus*uR + Fp*nx
    F(3) = rhoL*aL*c_plus*vL  + rhoR*aR*c_minus*vR + Fp*ny
    F(4) = rhoL*aL*c_plus*htL + rhoR*aR*c_minus*htR
    
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
  subroutine roe_flux( left, right, nx, ny, F )
    
    real(prec), dimension(:), intent(in)  :: left, right
    real(prec), dimension(:), intent(out) :: F
    real(prec), intent(in) :: nx, ny
    real(prec), dimension(size(left)) :: rvec1, rvec2, rvec3, rvec4
    real(prec), dimension(size(left)) :: lambda, FL, FR
    real(prec) :: rhoL, uL, vL, pL, aL, htL, unL
    real(prec) :: rhoR, uR, vR, pR, aR, htR, unR
    real(prec) :: rho2, u2, v2, a2, ht2, un2
    real(prec) :: dw1, dw2, dw3, dw4, R
    
    rhoL = left(1)
    uL   = left(2)
    vL   = left(3)
    pL   = left(4)
    call speed_of_sound(pL,rhoL,aL)
    
    rhoR = right(1)
    uR   = right(2)
    vR   = right(3)
    pR   = right(4)
    call speed_of_sound(pR,rhoR,aR)
    
    unL = uL*nx + vL*ny
    unR = uR*nx + vR*ny
    
    htL = aL**2/(gamma-one) + half*unL**2
    htR = aR**2/(gamma-one) + half*unR**2
    
    R = sqrt(rhoR/rhoL)
    rho2 = R*rhoL
    u2   = (R*uR + uL)/(R + one)
    v2   = (R*vR + vL)/(R + one)
    ht2  = (R*htR + htL)/(R + one)
    a2   = sqrt((gamma-one)*(ht2 - half*(u2**2 + v2**2)))
    un2 = u2*nx + v2*ny
    !a2   = sqrt((gamma-one)*(ht2 - half*(un2**2)))
    
    rvec1 = (/ one, u2, v2, half*(u2**2 +v2**2) /)
    rvec2 = (/ zero, ny*rho2, nx*rho2, rho2*(ny*u2-nx*v2) /)
    rvec3 = half*(rho2/a2)*(/ one, u2+nx*a2, v2+ny*a2, ht2+un2*a2 /)
    rvec4 = half*(rho2/a2)*(/ one, u2-nx*a2, v2-ny*a2, ht2-un2*a2 /)
    lambda = (/ un2, un2, un2 + a2, un2 - a2 /)
     
    lambda = abs(lambda)
    lambda = half*(one+sign(one,lambda-two*eps_roe*a2))*lambda &
           & + half*(one-sign(one,lambda-two*eps_roe*a2))*&
           & (lambda**2/(four*eps_roe*a2) + eps_roe*a2)
    
    dw1 = (rhoR - rhoL) - (pR - pL)/a2**2
    dw2 = ny*(uR - uL) - nx*(vR - vL)
    dw3 = nx*(uR - uL) + ny*(vR - vL) + (pR - pL)/(rho2*a2)
    dw4 = nx*(uR - uL) + ny*(vR - vL) - (pR - pL)/(rho2*a2)
    
    FL = (/ rhoL*unL, rhoL*uL*unL + pL*nx, rhoL*vL*unL + pL*ny, rhoL*htL*unL /)
    FR = (/ rhoR*unR, rhoR*uR*unR + pR*nx, rhoR*vR*unR + pR*ny, rhoR*htR*unR /)
    
    F = half*( (FL+FR) - lambda(1)*dw1*rvec1 - lambda(2)*dw2*rvec2 &
                       - lambda(3)*dw3*rvec3 - lambda(4)*dw4*rvec4 )
    
  end subroutine roe_flux
  
end module flux_calc
