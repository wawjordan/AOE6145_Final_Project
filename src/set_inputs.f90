module set_inputs

  use set_precision,   only : prec
  use set_constants,   only : zero, one, half, pi
  use fluid_constants, only : R_gas, gamma

  implicit none

  private

  public :: imax, jmax, neq, xmin, xmax, ymin, ymax, n_ghost
  public :: i_high, i_low, ig_high, ig_low
  public :: j_high, j_low, jg_high, jg_low
  public :: CFL, k2, k4, eps, tol, eps_roe, beta_lim, epsM, kappaM
  public :: max_iter, max_newton_iter, newton_tol, counter
  public :: isMMS, isAxi, cart_grid, soln_save, res_save, res_out
  public :: p0, T0, a0, rho0, pb, p_ratio
  public :: set_derived_inputs, flux_scheme, limiter_scheme, cons
  public :: leftV, rightV, leftU, rightU, limiter_freeze, psi_plus, psi_minus
  public :: grid_name, geometry_file
   
  integer :: imax    = 5
  integer :: jmax    = 5
  integer :: i_low   = 0
  integer :: i_high  = 0
  integer :: ig_low  = 0
  integer :: ig_high = 0
  integer :: j_low   = 0
  integer :: j_high  = 0
  integer :: jg_low  = 0
  integer :: jg_high = 0
  integer :: neq  = 4
  integer :: max_iter = 150000
  integer :: n_ghost   = 2
  integer :: max_newton_iter = 1000

  real(prec) :: newton_tol = 1.0e-15_prec
  real(prec) :: tol = 1.0e-11_prec
  real(prec) :: eps        = 1.0e-3_prec
  real(prec) :: p0         = 300.0_prec
  real(prec) :: T0         = 600.0_prec
  real(prec) :: Astar      = 0.2_prec
  real(prec) :: a0         = zero
  real(prec) :: rho0       = zero
  real(prec) :: xmin       = zero
  real(prec) :: xmax       = one
  real(prec) :: ymin       = zero
  real(prec) :: ymax       = one
  real(prec) :: CFL        = 0.1_prec
  real(prec) :: p_ratio    = 0.5_prec
  real(prec) :: pb         = 150000_prec
  real(prec) :: k2         = 1.0_prec/2.0_prec
  real(prec) :: k4         = 1.0_prec/32.0_prec
  integer :: flux_scheme   = 1
  integer :: limiter_scheme = 2
  real(prec) :: beta_lim = 2
  integer :: soln_save     = 150000
  integer :: res_save      = 10
  integer :: res_out       = 1000
  real(prec) :: eps_roe    = 0.1_prec
  real(prec) :: epsM       = zero
  real(prec) :: kappaM     = -one
  logical :: limiter_freeze = .false.
  logical :: cons           = .true.
  logical :: isMMS = .true.
  logical :: cart_grid = .false.
  logical :: isAxi = .false.

  real(prec), dimension(:,:), allocatable :: leftV, rightV, leftU, rightU
  real(prec), dimension(:,:), allocatable :: psi_plus, psi_minus
  integer :: counter = 1
  character(64) :: grid_name = "../grids/curvilinear-grids/curv2d9.grd"
  !character(64) :: grid_name = "../grids/NACA64A006-grids/NACA64A006.extra-coarse.27x14.grd"
  !character(64) :: grid_name = "../grids/NACA64A006-grids/NACA64A006.fine.385x105.grd"
  character(64) :: geometry_file = "example.dat"
  contains

  !=========================== set_derived_inputs ============================80
  !>
  !! Description: Sets derived quantities and prints to STDOUT.
  !<
  !===========================================================================80
  subroutine set_derived_inputs
    
    a0   = sqrt(gamma*R_gas*T0)
    rho0 = 1000.0_prec*p0/(R_gas*T0)
    i_low = 1
    j_low = 1
    i_high = imax-1
    j_high = jmax-1
    ig_low  = i_low - n_ghost
    jg_low  = j_low - n_ghost
    ig_high = i_high + n_ghost
    jg_high = j_high + n_ghost
    
    
    pb = p_ratio*p0*1000_prec
    write(*,'(A8,F20.14,A13)') 'R     = ', R_gas, ' [J/(kmol*K)]'
    write(*,'(A8,F20.14)')     'gamma = ', gamma
    write(*,'(A8,F20.14,A6)')  'a_0   = ', a0, ' [m/s]'
    write(*,'(A8,F20.14,A9)')  'rho_0 = ', rho0, ' [kg/m^3]'
    write(*,'(A8,F20.14,A6)')  'P_0   = ', p0, ' [kPa]'
    write(*,'(A8,F20.14,A4)')  'T_0   = ', T0, ' [K]'
    write(*,'(A8,F20.14,A6)')  'A*    = ', Astar, ' [m^2]'

    !allocate(leftV(i_low-1:i_high,1:neq))
    !allocate(rightV(i_low-1:i_high,1:neq))
    !allocate(leftU(i_low-1:i_high,1:neq))
    !allocate(rightU(i_low-1:i_high,1:neq))
    !allocate(psi_plus(ig_low:ig_high,1:neq))
    !allocate(psi_minus(ig_low:ig_high,1:neq))
    
  end subroutine set_derived_inputs

end module set_inputs
