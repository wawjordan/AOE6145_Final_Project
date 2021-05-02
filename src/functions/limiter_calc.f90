module limiter_calc
  
  use set_precision, only : prec
  use set_constants, only : zero, one, two, three, four, half, fourth
  use set_inputs, only : neq, imax, i_low, i_high, ig_low, ig_high
  use set_inputs, only : beta_lim, n_ghost
  
  implicit none
  
  private
  
  public :: calc_consecutive_variations, limiter_fun, select_limiter
  
  procedure( calc_limiter ), pointer :: limiter_fun
   
  abstract interface
  !============================== calc_limiter ===============================80
  !>
  !! Description:
  !!
  !! Inputs:      r   : 
  !!
  !! Outputs:     psi : 
  !<
  !===========================================================================80  
  subroutine calc_limiter( r, psi )
    
    import :: prec
    real(prec), dimension(:,:), intent(in) :: r
    real(prec), dimension(:,:), intent(out) :: psi
    
  end subroutine calc_limiter
    
  end interface
  
contains
  
  !============================== select_limiter =============================80
  !>
  !! Description:
  !<
  !===========================================================================80  
  subroutine select_limiter()
    
    use set_inputs, only : limiter_scheme
    
    limiter_fun => null()
    
    select case(limiter_scheme)
    
    case(1)
      limiter_fun => van_leer_limiter
    case(2)
      limiter_fun => van_albada_limiter
    case(3)
      limiter_fun => minmod_limiter
    case(4)
      limiter_fun => beta_limiter
    case default
    
    end select
  
  end subroutine select_limiter
  
  subroutine calc_consecutive_variations(V,r_plus,r_minus)
    
    real(prec), dimension(:,:), intent(in)  :: V
    real(prec), dimension(:,:), intent(out) :: r_plus, r_minus
    real(prec), dimension(neq) :: den
    integer :: i, j, low, high

    write(*,*) lbound(V)
    write(*,*) ubound(V)
    write(*,*)
    write(*,*) lbound(r_plus)
    write(*,*) ubound(r_minus)
    low = lbound(V,1)+n_ghost
    high = ubound(V,1)-n_ghost
    do i = low,high
      den = V(i+1,:) - V(i,:)
      den = sign(one,den)*max(abs(den),1e-6_prec)
      r_plus(i,:)   = ( V(i+2,:) - V(i+1,:) )/den
      r_minus(i,:)  = ( V(i,:) - V(i-1,:) )/den
    end do
    r_plus(low-2,:) = r_plus(low-1,:)
    r_minus(low-2,:) = r_minus(low-1,:)
    
    r_plus(high+1,:) = r_plus(high,:)
    r_minus(high+1,:) = r_minus(high,:)
    do i = lbound(V,1),ubound(V,1)
      write(*,*) V(i,1), V(i,2), V(i,3), V(i,4)
    end do
    write(*,*) 'r_plus'
    do i = lbound(V,1),ubound(V,1)
      write(*,*) r_plus(i,1), r_plus(i,2), r_plus(i,3), r_plus(i,4)
    end do
    write(*,*) 'r_minus'
    do i = lbound(V,1),ubound(V,1)
      write(*,*) r_minus(i,1), r_minus(i,2), r_minus(i,3), r_minus(i,4)
    end do
!    low = lbound(V,1)+n_ghost
!    high = ubound(V,1)-n_ghost
!    do i = low,high
!      den = V(i+1,:) - V(i,:)
!      den = sign(one,den)*max(abs(den),1e-6_prec)
!      r_plus(i,:)   = ( V(i+2,:) - V(i+1,:) )/den
!      r_minus(i,:)  = ( V(i,:) - V(i-1,:) )/den
!    end do
!    r_plus(low-2,:) = r_plus(low-1,:)
!    r_minus(low-2,:) = r_minus(low-1,:)
!    
!    r_plus(high+1,:) = r_plus(high,:)
!    r_minus(high+1,:) = r_minus(high,:)
    
!    r_plus(i_low-2,:) = r_plus(i_low-1,:)
!    r_minus(i_low-2,:) = r_minus(i_low-1,:)
!    
!    r_plus(i_high+1,:) = r_plus(i_high,:)
!    r_minus(i_high+1,:) = r_minus(i_high,:)
    
  end subroutine calc_consecutive_variations

  !=========================== van_leer_limiter ==============================80
  !>
  !! Description:
  !!
  !! Inputs:      r   : 
  !!
  !! Outputs:     psi : 
  !<
  !===========================================================================80  
  subroutine van_leer_limiter( r, psi )
    
    real(prec), dimension(:,:), intent(in) :: r
    !real(prec), dimension(i_low-1:i_high,neq), intent(out)   :: psi
    real(prec), dimension(:,:), intent(out)   :: psi
    
    psi = (r + abs(r))/(one + r)
    psi = half*(one+sign(one,r))*psi
    
  end subroutine van_leer_limiter
  
  !======================== van_albada_limiter ===============================80
  !>
  !! Description:
  !!
  !! Inputs:      r   : 
  !!
  !! Outputs:     psi : 
  !<
  !===========================================================================80  
  subroutine van_albada_limiter( r, psi )
    
    real(prec), dimension(:,:), intent(in) :: r
    !real(prec), dimension(i_low-1:i_high,neq), intent(out)   :: psi
    real(prec), dimension(:,:), intent(out)   :: psi
    
    psi = (r**2 + r)/(one + r**2)
    psi = half*(one+sign(one,r))*psi
    
  end subroutine van_albada_limiter

  !============================== minmod_limiter =============================80
  !>
  !! Description:
  !!
  !! Inputs:      r   : 
  !!
  !! Outputs:     psi : 
  !<
  !===========================================================================80  
  subroutine minmod_limiter( r, psi )
    
    real(prec), dimension(:,:), intent(in) :: r
    !real(prec), dimension(i_low-1:i_high,neq), intent(out)   :: psi
    real(prec), dimension(:,:), intent(out)   :: psi
    
    psi = half*(one + sign(one,r))*min(r,one)
    
  end subroutine minmod_limiter
  
  !============================== beta_limiter ===============================80
  !>
  !! Description:
  !!
  !! Inputs:      r   : 
  !!
  !! Outputs:     psi : 
  !<
  !===========================================================================80  
  subroutine beta_limiter( r, psi )
    
    real(prec), dimension(:,:), intent(in) :: r
    !real(prec), dimension(i_low-1:i_high,neq), intent(out)   :: psi
    real(prec), dimension(:,:), intent(out)   :: psi
    
    psi = maxval((/ zero, min(beta_lim*r,one), min(r,beta_lim) /))
    psi = half*(one+sign(one,r))*psi
    
  end subroutine beta_limiter
  
end module limiter_calc
