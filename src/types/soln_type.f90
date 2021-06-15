module soln_type

  use set_precision, only : prec
  use fluid_constants, only : gamma
  use set_constants, only : zero, one, half
  use set_inputs,    only : i_low, i_high, ig_low, ig_high
  use set_inputs,    only : j_low, j_high, jg_low, jg_high
  use set_inputs,    only : neq, max_iter, isMMS
  
  implicit none
  
  private
  
  public :: soln_t, allocate_soln, deallocate_soln, calc_mms

  type soln_t
    
    !sequence
    
    real(prec), allocatable, dimension(:,:,:) :: U ! conserved variables
    real(prec), allocatable, dimension(:,:,:) :: Fxi  ! normal fluxes
    real(prec), allocatable, dimension(:,:,:) :: Feta ! normal fluxes
    real(prec), allocatable, dimension(:,:,:) :: psi_p_xi  ! limiters
    real(prec), allocatable, dimension(:,:,:) :: psi_p_eta  ! limiters
    real(prec), allocatable, dimension(:,:,:) :: psi_m_xi  ! limiters
    real(prec), allocatable, dimension(:,:,:) :: psi_m_eta  ! limiters
    real(prec), allocatable, dimension(:,:,:) :: S ! source terms
    real(prec), allocatable, dimension(:,:,:) :: V ! primitive variables
    real(prec), allocatable, dimension(:,:,:) :: L ! eigenvalues
    real(prec), allocatable, dimension(:,:,:) :: R ! residuals
    real(prec), allocatable, dimension(:,:,:) :: Umms ! MMS conserved variables
    real(prec), allocatable, dimension(:,:,:) :: Vmms ! MMS primitive variables
    real(prec), allocatable, dimension(:,:,:) :: Smms ! MMS source terms
    real(prec), allocatable, dimension(:,:,:) :: DE ! discretization error
    real(prec), allocatable, dimension(:,:)   :: asnd
    real(prec), allocatable, dimension(:,:)   :: mach
    real(prec), allocatable, dimension(:,:)   :: temp
    real(prec), allocatable, dimension(:,:)   :: dt
    real(prec), allocatable, dimension(:,:)     :: DEnorm
    real(prec), allocatable, dimension(:)     :: rnorm
    real(prec), allocatable, dimension(:)     :: rold
    real(prec), allocatable, dimension(:)     :: rinit
    
  end type soln_t

  contains
  
  
  !============================= allocate_soln ===============================80
  !>
  !! Description: 
  !!
  !! Inputs:      soln : 
  !!
  !! Outputs:     soln : 
  !<
  !===========================================================================80
  subroutine allocate_soln( soln )

    type(soln_t), intent(inout) :: soln
    
    allocate( &
              soln%U( neq, ig_low:ig_high, jg_low:jg_high ), &
              soln%Fxi(  neq, i_low-1:i_high, j_low:j_high ), &
              soln%Feta( neq, i_low:i_high, j_low-1:j_high ), &
              soln%S( neq, ig_low:ig_high, jg_low:jg_high ), &
              soln%V( neq, ig_low:ig_high, jg_low:jg_high ), &
              soln%L( neq, ig_low:ig_high, jg_low:jg_high ), &
              soln%R(  neq,i_low:i_high,   j_low:j_high ), &
              soln%asnd( ig_low:ig_high, jg_low:jg_high ),   &
              soln%mach( ig_low:ig_high, jg_low:jg_high ),   &
              soln%temp( ig_low:ig_high, jg_low:jg_high ),   &
              soln%dt(   ig_low:ig_high, jg_low:jg_high ),   &
              soln%rnorm( neq ), &
              soln%rold( neq ),  &
              soln%rinit( neq ) )
    allocate( &
              soln%psi_p_xi(  neq, ig_low-1:ig_high, j_low:j_high ), &
              soln%psi_m_xi(  neq, ig_low-1:ig_high, j_low:j_high ), &
              soln%psi_p_eta( neq, i_low:i_high, jg_low-1:jg_high ), &
              soln%psi_m_eta( neq, i_low:i_high, jg_low-1:jg_high )  )
    
    if (isMMS) then
      allocate( soln%DE(   neq, ig_low:ig_high,  jg_low:jg_high ),  &
                soln%Vmms( neq, ig_low:ig_high,  jg_low:jg_high ),&
                soln%Umms( neq, ig_low:ig_high,  jg_low:jg_high ),&
                soln%Smms( neq, ig_low:ig_high,  jg_low:jg_high ),&
                soln%DEnorm( neq,3 ) )
      soln%DE     = zero
      soln%DEnorm = zero
    end if
    
    soln%U     = zero
    soln%Fxi   = zero
    soln%Feta  = zero
    soln%S     = zero
    soln%V     = zero
    soln%L     = one
    soln%R     = one
    soln%asnd  = one
    soln%mach  = one
    soln%temp  = one
    soln%dt    = one
    soln%rnorm = one
    soln%rold  = one
    soln%rinit = one
    
    soln%psi_p_xi  = one
    soln%psi_m_xi  = one
    soln%psi_p_eta = one
    soln%psi_m_eta = one

  end subroutine allocate_soln
  
  
  !=========================== deallocate_soln ===============================80
  !>
  !! Description: 
  !!
  !! Inputs:      soln : 
  !!
  !! Outputs:     soln : 
  !<
  !===========================================================================80
  subroutine deallocate_soln( soln )
  
    implicit none
    
    type(soln_t), intent(inout) :: soln
    
    deallocate( &
               soln%U,     &
               soln%Fxi,   &
               soln%Feta,  &
               soln%S,     &
               soln%V,     &
               soln%L,     &
               soln%R,     &
               soln%asnd,  &
               soln%mach,  &
               soln%temp,  &
               soln%dt,    &
               soln%rnorm, &
               soln%rold,  &
               soln%rinit   )
    deallocate( &
               soln%psi_p_xi,  soln%psi_m_xi, &
               soln%psi_p_eta, soln%psi_m_eta )
    
    if (isMMS) then
      deallocate( soln%DE, soln%Vmms, soln%Umms, soln%Smms, soln%DEnorm )
    end if
    
  end subroutine deallocate_soln
  
  subroutine calc_mms( grid, soln )
    
    use mms_functions, only : wrap_rho_mms, wrap_uvel_mms, &
                              wrap_vvel_mms, wrap_press_mms, &
                              wrap_rmassconv, wrap_xmtmconv, &
                              wrap_ymtmconv, wrap_energyconv
    use mms_functions, only : rho_mms, uvel_mms, vvel_mms, press_mms, &
                         rmassconv, xmtmconv, ymtmconv, energyconv
    use quadrature, only : gauss_quad, gauss_pts
    use grid_type, only : grid_t
    
    type(soln_t), intent(inout) :: soln
    type(grid_t), intent(inout) :: grid
    real(prec), dimension(:), allocatable :: xi, w
    real(prec), dimension(4,2) :: P
    real(prec) :: L
    integer :: N, i, j
    
    N = 5
    L = one
    
    allocate(xi(N),w(N))
    
    call gauss_pts(xi,w,N)
    
    soln%Smms(1,:,:) = rmassconv(L,grid%xc,grid%yc)
    soln%Smms(2,:,:) = xmtmconv(L,grid%xc,grid%yc)
    soln%Smms(3,:,:) = ymtmconv(L,grid%xc,grid%yc)
    soln%Smms(4,:,:) = energyconv(gamma,L,grid%xc,grid%yc)
    do j = jg_low, jg_high
      do i = ig_low, ig_high
        P = transpose( reshape( (/ &
          grid%x(i,j),     grid%y(i,j),      &
          grid%x(i+1,j),   grid%y(i+1,j),    &
          grid%x(i+1,j+1), grid%y(i+1,j+1),  &
          grid%x(i,j+1),   grid%y(i,j+1) /), (/2,4/) ) )
        
        !soln%Vmms(1,:,:) = rho_mms(L,grid%xc,grid%yc)
        !soln%Vmms(2,:,:) = uvel_mms(L,grid%xc,grid%yc)
        !soln%Vmms(3,:,:) = vvel_mms(L,grid%xc,grid%yc)
        !soln%Vmms(4,:,:) = press_mms(L,grid%xc,grid%yc)
        call gauss_quad(P,xi,xi,w,w,wrap_rho_mms,soln%Vmms(1,i,j))
        call gauss_quad(P,xi,xi,w,w,wrap_uvel_mms,soln%Vmms(2,i,j))
        call gauss_quad(P,xi,xi,w,w,wrap_vvel_mms,soln%Vmms(3,i,j))
        call gauss_quad(P,xi,xi,w,w,wrap_press_mms,soln%Vmms(4,i,j))
        soln%Vmms(:,i,j) = soln%Vmms(:,i,j)/grid%V(i,j)
        !call gauss_quad(P,xi,xi,w,w,wrap_rmassconv,soln%Smms(1,i,j))
        !call gauss_quad(P,xi,xi,w,w,wrap_xmtmconv,soln%Smms(2,i,j))
        !call gauss_quad(P,xi,xi,w,w,wrap_ymtmconv,soln%Smms(3,i,j))
        !call gauss_quad(P,xi,xi,w,w,wrap_energyconv,soln%Smms(4,i,j))
        !soln%Smms(:,i,j) = soln%Smms(:,i,j)/grid%V(i,j)
      end do
    end do
    
    deallocate(xi,w)
    
  end subroutine calc_mms
  
end module soln_type
