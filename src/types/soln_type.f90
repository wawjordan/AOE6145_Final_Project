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
    
    real(prec), allocatable, dimension(:,:,:) :: U ! conserved variables
    real(prec), allocatable, dimension(:,:,:,:) :: F ! normal fluxes
    real(prec), allocatable, dimension(:,:,:) :: S ! source terms
    real(prec), allocatable, dimension(:,:,:) :: V ! primitive variables
    real(prec), allocatable, dimension(:,:,:) :: L ! eigenvalues
    real(prec), allocatable, dimension(:,:,:) :: R ! residuals
    real(prec), allocatable, dimension(:,:,:,:) :: psi_plus  ! limiters
    real(prec), allocatable, dimension(:,:,:,:) :: psi_minus ! limiters
    real(prec), allocatable, dimension(:,:,:) :: Umms ! MMS conserved variables
    real(prec), allocatable, dimension(:,:,:) :: Vmms ! MMS primitive variables
    real(prec), allocatable, dimension(:,:,:) :: Smms ! MMS source terms
    real(prec), allocatable, dimension(:,:,:) :: DE ! discretization error
    real(prec), allocatable, dimension(:,:)   :: asnd
    real(prec), allocatable, dimension(:,:)   :: mach
    real(prec), allocatable, dimension(:,:)   :: temp
    real(prec), allocatable, dimension(:,:)   :: dt
    real(prec), allocatable, dimension(:)     :: DEnorm
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
              soln%U( ig_low:ig_high, jg_low:jg_high, neq ), &
              soln%F( i_low:i_high+1, j_low:j_high+1, neq, 2 ), &
              soln%S( ig_low:ig_high, jg_low:jg_high, neq ), &
              soln%V( ig_low:ig_high, jg_low:jg_high, neq ), &
              soln%L( ig_low:ig_high, jg_low:jg_high, neq ), &
              soln%R(  i_low:i_high,   j_low:j_high,  neq ), &
              soln%asnd( ig_low:ig_high, jg_low:jg_high ),   &
              soln%mach( ig_low:ig_high, jg_low:jg_high ),   &
              soln%temp( ig_low:ig_high, jg_low:jg_high ),   &
              soln%dt(   ig_low:ig_high, jg_low:jg_high ),   &
              soln%rnorm( neq ), &
              soln%rold( neq ),  &
              soln%rinit( neq ) )
    allocate( &
              soln%psi_plus(  ig_low:ig_high+1, jg_low:jg_high+1, neq,  2 ), &
              soln%psi_minus(  ig_low:ig_high+1, jg_low:jg_high+1, neq, 2 )  )
    
    if (isMMS) then
      allocate( soln%DE( i_low:i_high,  j_low:j_high, neq ),  &
                soln%Vmms( ig_low:ig_high,  jg_low:jg_high, neq ),&
                soln%Umms( ig_low:ig_high,  jg_low:jg_high, neq ),&
                soln%Smms( ig_low:ig_high,  jg_low:jg_high, neq ),&
                soln%DEnorm( neq ) )
      soln%DE     = zero
      soln%DEnorm = zero
    end if
    
    soln%U     = one
    soln%F     = one
    soln%S     = one
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
    
    soln%psi_plus  = one
    soln%psi_minus = one

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
               soln%F,     &
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
   ! deallocate( &
   !            soln%U,     &
   !            soln%F,     &
   !            soln%S,     &
   !            soln%V,     &
   !            soln%L,     &
   !            soln%R,     &
   !            soln%asnd,  &
   !            soln%mach,  &
   !            soln%temp,  &
   !            soln%dt,    &
   !            soln%rnorm, &
   !            soln%rold,  &
   !            soln%rinit   )
    deallocate( &
               soln%psi_plus,  &
               soln%psi_minus  )
    
    if (isMMS) then
      deallocate( soln%DE, soln%Vmms, soln%Umms, soln%Smms, soln%DEnorm )
    end if
    
  end subroutine deallocate_soln
  
  subroutine calc_mms( grid, soln )
    
    use mms_functions
    use quadrature
    use grid_type, only : grid_t
    
    type(soln_t), intent(inout) :: soln
    type(grid_t), intent(inout) :: grid
    real(prec), dimension(ig_low:ig_high,jg_low:jg_high) :: x,y
    real(prec) :: L
    integer :: N
    
    N = 5
    L = one
    
    call cv_averages(grid,N,wrap_rho_mms,soln%Vmms(:,:,1))
    call cv_averages(grid,N,wrap_uvel_mms,soln%Vmms(:,:,2))
    call cv_averages(grid,N,wrap_vvel_mms,soln%Vmms(:,:,3))
    call cv_averages(grid,N,wrap_press_mms,soln%Vmms(:,:,4))
    
    soln%Umms(:,:,1) = soln%Vmms(:,:,1)
    soln%Umms(:,:,2) = soln%Vmms(:,:,1)*soln%Vmms(:,:,2)
    soln%Umms(:,:,3) = soln%Vmms(:,:,1)*soln%Vmms(:,:,3)
    soln%Umms(:,:,4) = soln%Vmms(:,:,4)/(gamma-one) + &
                       half*soln%Vmms(:,:,1)*( &
                       soln%Vmms(:,:,2)*soln%Vmms(:,:,2) + &
                       soln%Vmms(:,:,3)*soln%Vmms(:,:,3) )
    
    x = grid%x(ig_low+1:ig_high+1,jg_low:jg_high) - &
        grid%x(ig_low:ig_high,jg_low:jg_high)
    y = grid%y(ig_low:ig_high,jg_low+1:jg_high+1) - &
        grid%y(ig_low:ig_high,jg_low:jg_high)
    
    call cv_averages(grid,1,wrap_rmassconv,soln%Smms(:,:,1))
    call cv_averages(grid,1,wrap_xmtmconv,soln%Smms(:,:,2))
    call cv_averages(grid,1,wrap_ymtmconv,soln%Smms(:,:,3))
    call cv_averages(grid,1,wrap_energyconv,soln%Smms(:,:,4))
!    soln%Smms(:,:,1) = rmassconv(L,x,y)
!    soln%Smms(:,:,2) = xmtmconv(L,x,y)
!    soln%Smms(:,:,3) = ymtmconv(L,x,y)
!    soln%Smms(:,:,4) = energyconv(gamma,L,x,y)

!    soln%Smms(:,:,1) = rmassconv(L,                      &
!                       grid%x(ig_low:ig_high,jg_low:jg_high),&
!                       grid%y(ig_low:ig_high,jg_low:jg_high))
!    soln%Smms(:,:,2) = xmtmconv(L,                       &
!                       grid%x(ig_low:ig_high,jg_low:jg_high),&
!                       grid%y(ig_low:ig_high,jg_low:jg_high))
!    soln%Smms(:,:,3) = ymtmconv(L,                       &
!                       grid%x(ig_low:ig_high,jg_low:jg_high),&
!                       grid%y(ig_low:ig_high,jg_low:jg_high))
!    soln%Smms(:,:,4) = energyconv(gamma, L,              &
!                       grid%x(ig_low:ig_high,jg_low:jg_high),&
!                       grid%y(ig_low:ig_high,jg_low:jg_high))
   
  end subroutine calc_mms
!  subroutine calc_mms( grid, soln )
!    
!    use fluid_constants, only : gamma
!    use mms_functions
!    use grid_type, only : grid_t
!    
!    type(soln_t), intent(inout) :: soln
!    type(grid_t), intent(inout) :: grid
!    
!    real(prec) :: L  = one
!    
!    soln%Vmms(:,:,1) = rho_mms(L,                        &
!                       grid%x(ig_low:ig_high,jg_low:jg_high),&
!                       grid%y(ig_low:ig_high,jg_low:jg_high))
!    soln%Vmms(:,:,2) = uvel_mms(L,                       &
!                       grid%x(ig_low:ig_high,jg_low:jg_high),&
!                       grid%y(ig_low:ig_high,jg_low:jg_high))
!    soln%Vmms(:,:,3) = vvel_mms(L,                       &
!                       grid%x(ig_low:ig_high,jg_low:jg_high),&
!                       grid%y(ig_low:ig_high,jg_low:jg_high))
!    soln%Vmms(:,:,4) = press_mms(L,                      &
!                       grid%x(ig_low:ig_high,jg_low:jg_high),&
!                       grid%y(ig_low:ig_high,jg_low:jg_high))
!    
!    soln%Smms(:,:,1) = rmassconv(L,                      &
!                       grid%x(ig_low:ig_high,jg_low:jg_high),&
!                       grid%y(ig_low:ig_high,jg_low:jg_high))
!    soln%Smms(:,:,2) = xmtmconv(L,                       &
!                       grid%x(ig_low:ig_high,jg_low:jg_high),&
!                       grid%y(ig_low:ig_high,jg_low:jg_high))
!    soln%Smms(:,:,3) = ymtmconv(L,                       &
!                       grid%x(ig_low:ig_high,jg_low:jg_high),&
!                       grid%y(ig_low:ig_high,jg_low:jg_high))
!    soln%Smms(:,:,4) = energyconv(gamma, L,              &
!                       grid%x(ig_low:ig_high,jg_low:jg_high),&
!                       grid%y(ig_low:ig_high,jg_low:jg_high))
!    
!  end subroutine calc_mms
end module soln_type
