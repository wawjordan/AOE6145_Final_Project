module boundary_conds

  use set_precision, only : prec
  use set_constants, only : one, two
  use set_inputs,    only : imax, neq, eps, i_low, i_high, ig_low, ig_high
  use variable_conversion, only : prim2cons, cons2prim, speed_of_sound
  
  implicit none
  
  private
  
  public :: subsonic_dirichlet_BC, supersonic_inflow_BC, supersonic_outflow_BC
  public :: inviscid_wall_BC, wake_cut_BC
  public :: enforce_bndry
  
  procedure( enforce ), pointer :: enforce_BC
  
  abstract interface
!================================= enforce_bndry ============================80
!>
!! Description: 
!!
!! Inputs:      V :
!!
!! Outputs:     V : 
!<
!============================================================================80
subroutine enforce( bndry, soln )
  
  import :: bc_t, soln_t
  type(bc_t), intent(in) :: bndry
  type(soln_t), intent(inout) :: soln
  
end subroutine enforce
  
  end interface

  contains
  
!=========================== subsonic_dirichlet_BC  =========================80
!>
!! Description: 
!!
!! Inputs:      V_spec :
!!
!! Outputs:     V :
!<
!============================================================================80
subroutine subsonic_dirichlet_BC( bndry, soln )
  
  type(bc_t), intent(in) :: bndry
  type(soln_t), intent(inout) :: soln
  
  V = Vspec
  
end subroutine subsonic_dirichlet_BC

!============================ supersonic_inflow_BC ==========================80
!>
!! Description: 
!!
!! Inputs:      V :
!!
!! Outputs:     V : 
!<
!============================================================================80
subroutine supersonic_inflow_BC( bndry, soln )
  
  type(bc_t), intent(in) :: bndry
  type(soln_t), intent(inout) :: soln
  
  
end subroutine supersonic_inflow_BC

!=========================== supersonic_outflow_BC ==========================80
!>
!! Description: 
!!
!! Inputs:      U : 
!!              V :
!!
!! Outputs:     U :
!!              V : 
!<
!============================================================================80
subroutine supersonic_outflow_BC( bndry, soln )
  
  type(bc_t), intent(in) :: bndry
  type(soln_t), intent(inout) :: soln
  
  
end subroutine supersonic_outflow_BC

!================================= inviscid_wall_BC =========================80
!>
!! Description: 
!!
!! Inputs:      V :
!!
!! Outputs:     V : 
!<
!============================================================================80
subroutine inviscid_wall_BC( V_st, n_xi, n_eta, F )
  
  real(prec), dimension(:,:),  intent(in) :: V_st
  real(prec), dimension(2),    intent(in) :: n_xi, n_eta
  real(prec), dimension(:), intent(inout) :: F
  real(prec), dimension(lbound(V_st,1):ubound(V_st,1)) :: p
  real(prec) :: p_i
  
  p = V_st(:,4)
  
  p_i = p(1) - p_wall_extrap*half*(p(2)-p(1))
  
  F = (/ zero, (n_xi(1)+n_eta(1))*p_i, (n_xi(1)+n_eta(1))*p_i, zero /)
  
end subroutine inviscid_wall_BC

!================================= wake_cut_BC ==============================80
!>
!! Description: 
!!
!! Inputs:      V :
!!
!! Outputs:     V : 
!<
!============================================================================80
subroutine wake_cut_BC( )
  
end subroutine wake_cut_BC


end module boundary_conds
