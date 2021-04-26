module soln_type

  use set_precision, only : prec
  use set_constants, only : zero, one
  use set_inputs,    only : i_low, i_high, ig_low, ig_high
  use set_inputs,    only : j_low, j_high, jg_low, jg_high
  use set_inputs,    only : neq, max_iter, isMMS
  
  implicit none
  
  private
  
  public :: soln_t, allocate_soln, deallocate_soln

  type soln_t
    
    real(prec), allocatable, dimension(:,:,:) :: U ! conserved variables
    real(prec), allocatable, dimension(:,:,:) :: F ! normal fluxes
    real(prec), allocatable, dimension(:,:,:) :: S ! source terms
    real(prec), allocatable, dimension(:,:,:) :: V ! primitive variables
    real(prec), allocatable, dimension(:,:,:) :: L ! eigenvalues
    real(prec), allocatable, dimension(:,:,:) :: R ! residuals
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
              soln%F( ig_low:ig_high, jg_low:jg_high, neq ), &
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
    
    if (isMMS) then
      allocate( soln%DE( i_low:i_high,  j_low:j_high, neq ), &
                soln%DEnorm( neq ) )
      soln%DE     = zero
      soln%DEnorm = zero
    end if
    
    soln%U     = zero
    soln%F     = zero
    soln%S     = zero
    soln%V     = zero
    soln%L     = zero
    soln%R     = zero
    soln%asnd  = zero
    soln%mach  = zero
    soln%temp  = zero
    soln%dt    = zero
    soln%rnorm = zero
    soln%rold  = zero
    soln%rinit = zero

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
    
    if (isMMS) then
      deallocate( soln%DE, soln%DEnorm )
    end if
    
  end subroutine deallocate_soln

end module soln_type
