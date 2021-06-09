module set_inputs

  use set_precision,   only : prec
  use set_constants,   only : zero, one, two, half, pi
  use fluid_constants, only : R_gas, gamma

  implicit none

  private

  public :: isMMS, isAxi, cart_grid, limiter_freeze, cons
  public :: flux_scheme, limiter_scheme
  public :: imax, jmax, neq, xmin, xmax, ymin, ymax, n_ghost
  public :: i_high, i_low, ig_high, ig_low
  public :: j_high, j_low, jg_high, jg_low
  public :: max_iter, soln_save, res_save, res_out , counter
  public :: CFL, eps, tol, eps_roe, beta_lim, epsM, kappaM
  public :: rho_inf, u_inf, p_inf, a_inf, T_inf, alpha, u0, v0, Lmms
  public :: grid_name, geometry_file
  public :: set_derived_inputs
   
  logical :: limiter_freeze = .false.
  logical :: cons           = .false.
  logical :: isMMS          = .true.
  logical :: cart_grid      = .false.
  logical :: isAxi          = .false.
  
  integer :: imax           = 129
  integer :: jmax           = 129
  integer :: i_low          = 0
  integer :: i_high         = 0
  integer :: ig_low         = 0
  integer :: ig_high        = 0
  integer :: j_low          = 0
  integer :: j_high         = 0
  integer :: jg_low         = 0
  integer :: jg_high        = 0
  integer :: neq            = 4
  integer :: n_ghost        = 2
  integer :: counter        = 1
  integer :: max_iter       = 150000
  integer :: soln_save      = 100
  integer :: res_save       = 10
  integer :: res_out        = 100
  integer :: flux_scheme    = 1
  integer :: limiter_scheme = 3

  real(prec) :: tol        = 1.0e-9_prec
  real(prec) :: eps        = 1.0e-3_prec
  real(prec) :: p_inf      = 1.0e5_prec
  real(prec) :: T_inf      = 600.0_prec
  real(prec) :: u_inf      = 200.0_prec
  real(prec) :: u0         = zero
  real(prec) :: v0         = zero
  real(prec) :: alpha      = zero
  real(prec) :: a_inf      = zero
  real(prec) :: rho_inf    = zero
  real(prec) :: xmin       = zero
  real(prec) :: xmax       = one
  real(prec) :: ymin       = zero
  real(prec) :: ymax       = one
  real(prec) :: CFL        = 0.1_prec
  real(prec) :: beta_lim   = two
  real(prec) :: Lmms       = one
  real(prec) :: eps_roe    = 0.1_prec
  real(prec) :: epsM       = one
  real(prec) :: kappaM     = zero

  character(200) :: grid_name = "../grids/curvilinear-grids/curv2d65.grd"
  !character(64) :: grid_name = "../grids/NACA64A006-grids/NACA64A006.extra-coarse.27x14.grd"
  !character(64) :: grid_name = "../grids/NACA64A006-grids/NACA64A006.fine.385x105.grd"
  character(200) :: geometry_file = "example.dat"
  contains

  !=========================== set_derived_inputs ============================80
  !>
  !! Description: Sets derived quantities and prints to STDOUT.
  !<
  !===========================================================================80
  subroutine set_derived_inputs
    
    a_inf   = sqrt(gamma*R_gas*T_inf)
    rho_inf = p_inf/(R_gas*T_inf)
    i_low = 1
    j_low = 1
    i_high = imax-1
    j_high = jmax-1
    ig_low  = i_low - n_ghost
    jg_low  = j_low - n_ghost
    ig_high = i_high + n_ghost
    jg_high = j_high + n_ghost
    
  end subroutine set_derived_inputs

end module set_inputs
