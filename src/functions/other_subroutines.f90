module other_subroutines
  
  use set_precision, only : prec
  use set_constants, only : zero, one, two, three, half, fourth
  use set_inputs, only : imax, neq, i_low, i_high, ig_low, ig_high
  use set_inputs, only : jmax, j_low, j_high, jg_low, jg_high
  use set_inputs, only : epsM, kappaM, isMMS, n_ghost
  use fluid_constants, only : gamma
  use variable_conversion, only : speed_of_sound
  use limiter_calc, only : limiter_fun, calc_consecutive_variations
  use soln_type, only : soln_t
  use grid_type, only : grid_t
  
  implicit none
  
  private
  
  public :: MUSCL_extrap, calc_de, calc_sources, Limit
  
  contains
  
  !============================= calculate_sources ===========================80
  !>
  !! Description: 
  !!
  !! Inputs:      V  : 
  !!              y  : 
  !!
  !! Outputs:     S  : 
  !<
  !===========================================================================80
  subroutine source_terms(V,y,isAxi,S)
    
    real(prec), dimension(:),   intent(in) :: V
    real(prec), dimension(:),   intent(out) :: S
    real(prec), intent(in) :: y
    logical, intent(in) :: isAxi
    real(prec) :: rho, uvel, vvel, p, a, ht
    
    rho  = V(1)
    uvel = V(2)
    vvel = V(3)
    p    = V(4)
    call speed_of_sound(p,rho,a)
    ht   = a**2/(gamma-one) + half*(uvel**2 + vvel**2)
    
    S = (/ rho*vvel, rho*uvel*vvel, rho*vvel**2, rho*vvel*ht /)
    S = merge(-one,zero,isAxi)/merge(y,1.0e-6_prec,(y>zero))*S
    
  end subroutine source_terms
  
  !================================ MUSCL_extrap =============================80
  !>
  !! Description: 
  !!
  !! Inputs:      V     : 
  !!
  !! Outputs:     left  : 
  !!              right :
  !<
  !===========================================================================80
  subroutine Limit(V,psi_plus,psi_minus)
    
    use set_inputs, only : limiter_freeze
    real(prec), dimension(:,:,:), intent(in)  :: V
    !real(prec), dimension(:,:,:,:), intent(inout)  :: psi_plus, psi_minus
    real(prec), dimension(neq,ig_low:ig_high, &
                          jg_low:jg_high,2),intent(inout) :: psi_plus, psi_minus
    !real(prec), dimension(lbound(psi_plus,1):ubound(psi_plus,1), &
    !                      lbound(psi_plus,2):ubound(psi_plus,2),neq) :: r_plus, r_minus
    real(prec), dimension(neq,ig_low:ig_high, &
                          jg_low:jg_high) :: r_plus, r_minus
    !write(*,*) lbound(psi_plus,1),ubound(psi_plus,1)
    !write(*,*) lbound(psi_plus,2),ubound(psi_plus,2)
    !stop
    r_plus = zero
    r_minus = zero
    if (limiter_freeze) then
      continue
    else
      call calc_consecutive_variations(V,r_plus,r_minus,1)
      call limiter_fun(r_plus,psi_plus(:,:,:,1))
      call limiter_fun(r_minus,psi_minus(:,:,:,1))
      
      call calc_consecutive_variations(V,r_plus,r_minus,2)
      call limiter_fun(r_plus,psi_plus(:,:,:,2))
      call limiter_fun(r_minus,psi_minus(:,:,:,2))
    end if
  
  end subroutine Limit

  subroutine MUSCL_extrap(stencil,psi_plus,psi_minus,left,right)

    real(prec), dimension(neq,4), intent(in)  :: stencil
    real(prec), dimension(neq,4), intent(in)  :: psi_plus, psi_minus
    real(prec), dimension(neq), intent(out) :: left, right
    integer :: i
    
    i = 2
    

    left  = stencil(:,i) + fourth*epsM*( &
         & (one-kappaM)*psi_plus(:,i-1)*(stencil(:,i)-stencil(:,i-1)) + &
         & (one+kappaM)*psi_minus(:,i)*(stencil(:,i+1)-stencil(:,i)) )
    right = stencil(:,i+1) - fourth*epsM*( &
         & (one+kappaM)*psi_minus(:,i+1)*(stencil(:,i+1)-stencil(:,i)) + &
         & (one-kappaM)*psi_plus(:,i)*(stencil(:,i+2)-stencil(:,i+1)) )
  end subroutine MUSCL_extrap
  
  !================================== calc_sources ==========================80
  !>
  !! Description: 
  !!
  !! Inputs:      soln : 
  !!              exact_soln : 
  !!              pnorm :
  !!
  !! Outputs:     DE     : 
  !!              DEnorm : 
  !<
  !===========================================================================80
  subroutine calc_sources( soln, grid )
    
    use set_inputs, only : isAxi
    
    type(soln_t), intent(inout) :: soln
    type(grid_t), intent(in) :: grid
    integer :: i, j
    
    do j = grid%j_low,grid%j_high
      do i = grid%i_low,grid%i_high
         call source_terms( soln%V(:,i,j), grid%y(i,j), isAxi, soln%S(:,i,j) )
      end do
    end do
    
  end subroutine calc_sources
     
  
  !================================== calc_de ==== ===========================80
  !>
  !! Description: 
  !!
  !! Inputs:      soln : 
  !!              exact_soln : 
  !!              pnorm :
  !!
  !! Outputs:     DE     : 
  !!              DEnorm : 
  !<
  !===========================================================================80
  subroutine calc_de( soln, DE, DEnorm, pnorm, cons )
    
    type(soln_t), intent(in) :: soln
    logical, intent(in) :: cons
    real(prec), dimension(:,:,:), intent(out) :: DE
    real(prec), dimension(3,neq), intent(out) :: DEnorm
    integer, intent(in) :: pnorm
    real(prec) :: Linv
    !Linv = one/real( (i_high-i_low)*(j_high-j_low) )
    !if (cons) then
    !  DE = soln%U(i_low:i_high,j_low:j_high,1:neq) &
    !   & - soln%Umms(i_low:i_high,j_low:j_high,1:neq)
    !else
    !  DE = soln%V(i_low:i_high,j_low:j_high,1:neq) &
    !   & - soln%Vmms(i_low:i_high,j_low:j_high,1:neq)
    !end if
    DE = zero
    Linv = one/real( (i_high-i_low)*(j_high-j_low) )
    if (cons) then
      DE = soln%U - soln%Umms
    else
      DE = soln%V - soln%Vmms
    end if
    
    if (pnorm == 0) then
      DEnorm(1,1) = maxval( abs( DE(1,:,:) ) )
      DEnorm(1,2) = maxval( abs( DE(2,:,:) ) )
      DEnorm(1,3) = maxval( abs( DE(3,:,:) ) )
      DEnorm(1,4) = maxval( abs( DE(4,:,:) ) )
    elseif (pnorm == 1) then
      DEnorm(2,1) = Linv*sum( abs( DE(1,:,:) ) )
      DEnorm(2,2) = Linv*sum( abs( DE(2,:,:) ) )
      DEnorm(2,3) = Linv*sum( abs( DE(3,:,:) ) )
      DEnorm(2,4) = Linv*sum( abs( DE(4,:,:) ) )
    elseif (pnorm == 2) then
      DEnorm(3,1) = sqrt( Linv*sum( DE(1,:,:)**2 ) )
      DEnorm(3,2) = sqrt( Linv*sum( DE(2,:,:)**2 ) )
      DEnorm(3,3) = sqrt( Linv*sum( DE(3,:,:)**2 ) )
      DEnorm(3,4) = sqrt( Linv*sum( DE(4,:,:)**2 ) )
    end if
    
  end subroutine calc_de
  
end module other_subroutines
