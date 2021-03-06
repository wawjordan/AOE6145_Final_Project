module limiter_calc
  
  use set_precision, only : prec
  use set_constants, only : zero, one, two, three, four, half, fourth
  use set_inputs, only : i_low, i_high, ig_low, ig_high, imax
  use set_inputs, only : j_low, j_high, jg_low, jg_high, jmax
  use set_inputs, only : neq, beta_lim, n_ghost
  implicit none
  
  private
  
  public :: calculate_limiters, limiter_fun, select_limiter
  
  procedure( calc_limiter ), pointer :: limiter_fun => null()
   
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
  subroutine calc_limiter(r,psi)
    
    import :: prec
    real(prec), dimension(:,:,:), intent(in) :: r
    real(prec), dimension(:,:,:), intent(out) :: psi
    
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
    
    select case(limiter_scheme)
    
    case(0)
      limiter_fun => null_limiter
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
  
  !========================== calculate_limiter ==============================80
  !>
  !! Description:
  !!
  !! Inputs:      soln :
  !!
  !! Outputs:     soln :
  !<
  !===========================================================================80
  subroutine calculate_limiters(soln)
    
    use soln_type, only : soln_t
    type(soln_t), intent(inout) :: soln
    real(prec), dimension(neq,i_low-1:i_high,j_low:j_high) :: &
             psi_p_xi, psi_m_xi, r_plus_xi, r_minus_xi, den_xi
    real(prec), dimension(neq,i_low:i_high,j_low-1:j_high) :: &
             psi_p_eta, psi_m_eta, r_plus_eta, r_minus_eta, den_eta
    integer :: i
    
    den_xi = soln%V(:,i_low  :i_high+1,j_low:j_high)       &
           - soln%V(:,i_low-1:I_high  ,j_low:j_high)
    den_xi = sign(one,den_xi)*max(abs(den_xi),1e-6_prec)
    r_plus_xi  = ( soln%V(:,i_low+1:i_high+2,j_low:j_high) &
                 - soln%V(:,i_low  :i_high+1,j_low:j_high) )/den_xi
    
    r_minus_xi = ( soln%V(:,i_low-1:i_high  ,j_low:j_high) &
                 - soln%V(:,i_low-2:i_high-1,j_low:j_high) )/den_xi
    call limiter_fun(r_plus_xi,psi_p_xi)
    call limiter_fun(r_minus_xi,psi_m_xi)
    
    den_eta = soln%V(:,i_low:i_high,j_low  :j_high+1)      &
            - soln%V(:,i_low:i_high,j_low-1:j_high  )
    den_eta = sign(one,den_eta)*max(abs(den_eta),1e-6_prec)
    
    r_plus_eta  = ( soln%V(:,i_low:i_high,j_low+1:j_high+2) &
                  - soln%V(:,i_low:i_high,j_low  :j_high+1) )/den_eta
    
    r_minus_eta = ( soln%V(:,i_low:i_high,j_low-1:j_high  ) &
                  - soln%V(:,i_low:i_high,j_low-2:j_high-1) )/den_eta

    call limiter_fun(r_plus_eta,psi_p_eta)
    call limiter_fun(r_minus_eta,psi_m_eta)
    
    soln%psi_p_xi( :,i_low-1:i_high,j_low:j_high) = psi_p_xi
    soln%psi_m_xi( :,i_low-1:i_high,j_low:j_high) = psi_m_xi
    soln%psi_p_eta(:,i_low:i_high,j_low-1:j_high) = psi_p_eta
    soln%psi_m_eta(:,i_low:i_high,j_low-1:j_high) = psi_m_eta
    
    do i = 1,n_ghost
    soln%psi_p_xi( :,i_low-1-i,:) = &
         soln%psi_p_xi(:,i_low-i,:)
    soln%psi_m_xi( :,i_low-1-i,:) = &
         soln%psi_m_xi(:,i_low-i,:)
    soln%psi_p_xi( :,i_high+i,:) = &
         soln%psi_p_xi(:,i_high-1+i,:)
    soln%psi_m_xi( :,i_high+i,:) = &
         soln%psi_m_xi(:,i_high-1+i,:)
    
    soln%psi_p_eta( :,:,j_low-1-i) = &
         soln%psi_p_eta(:,:,j_low-i)
    soln%psi_m_eta( :,:,j_low-1-i) = &
         soln%psi_m_eta(:,:,j_low-i)
    soln%psi_p_eta( :,:,j_high+i) = &
         soln%psi_p_eta(:,:,j_high-1+i)
    soln%psi_m_eta( :,:,j_high+i) = &
         soln%psi_m_eta(:,:,j_high-1+i)
    end do
    !do i = 1,n_ghost
    !soln%psi_p_xi( :,i_low-1-i,:) = &
    !     two*soln%psi_p_xi(:,i_low-i,:) - soln%psi_p_xi(:,i_low+1-i,:)
    !soln%psi_m_xi( :,i_low-1-i,:) = &
    !     two*soln%psi_m_xi(:,i_low-i,:) - soln%psi_m_xi(:,i_low+1-i,:)
    !soln%psi_p_xi( :,i_high+i,:) = &
    !     two*soln%psi_p_xi(:,i_high-1+i,:) - soln%psi_p_xi(:,i_high-2+i,:)
    !soln%psi_m_xi( :,i_high+i,:) = &
    !     two*soln%psi_m_xi(:,i_high-1+i,:) - soln%psi_m_xi(:,i_high-2+i,:)
    !
    !soln%psi_p_eta( :,:,j_low-1-i) = &
    !     two*soln%psi_p_eta(:,:,j_low-i) - soln%psi_p_eta(:,:,j_low+1-i)
    !soln%psi_m_eta( :,:,j_low-1-i) = &
    !     two*soln%psi_m_eta(:,:,j_low-i) - soln%psi_m_eta(:,:,j_low+1-i)
    !soln%psi_p_eta( :,:,j_high+i) = &
    !     two*soln%psi_p_eta(:,:,j_high-1+i) - soln%psi_p_eta(:,:,j_high-2+i)
    !soln%psi_m_eta( :,:,j_high+i) = &
    !     two*soln%psi_m_eta(:,:,j_high-1+i) - soln%psi_m_eta(:,:,j_high-2+i)
    !end do

  end subroutine calculate_limiters
  
  !============================= null_limiter ================================80
  !>
  !! Description:
  !!
  !! Inputs:      r   :
  !!
  !! Outputs:     psi :
  !<
  !===========================================================================80
  subroutine null_limiter(r,psi)

    real(prec), dimension(:,:,:),  intent(in) :: r
    real(prec), dimension(:,:,:), intent(out) :: psi

    psi = one

  end subroutine null_limiter

  !=========================== van_leer_limiter ==============================80
  !>
  !! Description:
  !!
  !! Inputs:      r   :
  !!
  !! Outputs:     psi :
  !<
  !===========================================================================80
  subroutine van_leer_limiter(r,psi)

    real(prec), dimension(:,:,:),  intent(in) :: r
    real(prec), dimension(:,:,:), intent(out) :: psi
    
    psi = sign(one,one+r)*max(abs(one+r),1.0e-12_prec)
    psi = (r + abs(r))/psi
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
  subroutine van_albada_limiter(r,psi)

    real(prec), dimension(:,:,:),  intent(in) :: r
    real(prec), dimension(:,:,:), intent(out) :: psi

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
  subroutine minmod_limiter(r,psi)

    real(prec), dimension(:,:,:),  intent(in) :: r
    real(prec), dimension(:,:,:), intent(out) :: psi

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
  subroutine beta_limiter(r,psi)

    real(prec), dimension(:,:,:), intent(in) :: r
    real(prec), dimension(:,:,:), intent(out) :: psi

    psi = max(zero, min(beta_lim*r,one), min(r,beta_lim))
    psi = half*(one+sign(one,r))*psi

  end subroutine beta_limiter
  
end module limiter_calc
