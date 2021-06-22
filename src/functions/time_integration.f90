module time_integration
  
  use set_precision, only : prec
  use set_constants, only : zero, one, half, third, fourth
  use set_inputs,    only : neq, i_low, i_high, ig_low, ig_high
  use set_inputs,    only : j_low, j_high, jg_low, jg_high
  use grid_type,     only : grid_t
  use soln_type,     only : soln_t
  use variable_conversion, only : speed_of_sound
  
  implicit none
  
  contains
  
  subroutine SER(CFL,CFLmax,Rnorm,Rnorm2)
    real(prec), intent(inout) :: CFL
    real(prec), intent(in) :: CFLmax
    real(prec), dimension(4), intent(in) :: Rnorm, Rnorm2
    !if (all(Rnorm.LE.Rnorm2)) then
      CFL = min(minval(CFL*(Rnorm2/Rnorm)),CFLmax)
    !end if
    
  end subroutine SER
  
  subroutine RDM(CFL,CFLmax,beta,Rnorm,Rnorm2)
    real(prec), intent(inout) :: CFL
    real(prec), intent(in) :: CFLmax, beta
    real(prec), dimension(4), intent(in) :: Rnorm, Rnorm2
    real(prec) :: g
    !if (all(Rnorm.LE.Rnorm2)) then
      g = minval((Rnorm2-Rnorm)/Rnorm2)
    !else
    !  g = zero
    !end if
    CFL = min(CFL*beta**g,CFLmax)
    
  end subroutine RDM
  
  
  
  !=========================== calc_time_step ===============================80
  !>
  !! Description: 
  !!
  !! Inputs:      dx : 
  !!              V  : 
  !!
  !! Outputs:     lambda : 
  !!              dt     : 
  !<
  !==========================================================================80
  subroutine calc_time_step(grid,soln)
    use set_inputs, only : CFL, locTime
    type(grid_t), intent(in)    :: grid
    type(soln_t), intent(inout) :: soln
    real(prec), dimension(ig_low:ig_high,jg_low:jg_high) :: Lam_xi,Lam_eta
    
    call speed_of_sound(soln%V(4,:,:),soln%V(1,:,:),soln%asnd)
    
    Lam_xi  = abs(soln%V(2,:,:)*grid%n_xi_avg(:,:,1)   + &
                  soln%V(3,:,:)*grid%n_xi_avg(:,:,2))  + soln%asnd
    Lam_eta = abs(soln%V(2,:,:)*grid%n_eta_avg(:,:,1)  + &
                  soln%V(3,:,:)*grid%n_eta_avg(:,:,2)) + soln%asnd
    soln%dt = CFL*grid%V/( &
               Lam_xi*grid%A_xi(ig_low:ig_high,jg_low:jg_high) + &
               Lam_eta*grid%A_eta(ig_low:ig_high,jg_low:jg_high) )
    if (.not.locTime) then
      soln%dt = CFL*minval(soln%dt)
    end if
    
  end subroutine calc_time_step
!  subroutine calc_time_step( A_xi, A_eta, n_xi_avg, n_eta_avg, vol, V, dt )
!    
!    use set_inputs,          only : CFL
!    
!    real(prec), dimension(:,:,:), intent(in)  :: V
!    real(prec), dimension(:,:), intent(in)    :: A_xi, A_eta, vol
!    real(prec), dimension(:,:,:), intent(in)  :: n_xi_avg, n_eta_avg
!    real(prec), dimension(lbound(V,1):ubound(V,1),lbound(V,2):ubound(V,2)) &
!                                         :: asnd, lambda_xi, lambda_eta
!    real(prec), dimension(:,:),   intent(out) :: dt
!    integer :: low1, low2, high1, high2
!    
!    low1  = lbound(V,2)
!    high1 = ubound(V,2)
!    low2  = lbound(V,3)
!    high2 = ubound(V,3)
!
!
!    call speed_of_sound(V(4,:,:),V(1,:,:),asnd)
!    lambda_xi  = abs(V(2,:,:)*n_xi_avg(1,:,:) + &
!                     V(3,:,:)*n_xi_avg(2,:,:)) + asnd
!    lambda_eta = abs(V(2,:,:)*n_eta_avg(1,:,:) + &
!                     V(3,:,:)*n_eta_avg(2,:,:)) + asnd
!    dt = CFL*vol/(lambda_xi*A_xi(low1:high1,low2:high2) + &
!                  lambda_eta*A_eta(low1:high1,low2:high2))
!    dt = CFL*minval(dt)
!    
!  end subroutine calc_time_step
  
  !============================= explicit_RK ================================80
  !>
  !! Description: 
  !!
  !! Inputs:      grid : 
  !!              soln : 
  !!
  !! Outputs:     grid : 
  !!              soln : 
  !<
  !==========================================================================80

  subroutine explicit_RK( grid, soln, bnds )
    use flux_calc, only : calc_flux_2D
    use variable_conversion, only : update_states, prim2cons
    use bc_type, only : bc_t
    use set_inputs, only : num_BCs
    type(grid_t), intent(inout) :: grid
    type(soln_t), intent(inout) :: soln
    type(bc_t), dimension(:), intent(inout) :: bnds
    real(prec), dimension(4) :: k
    integer :: i, j, n
    
    k = (/ fourth, third, half, one /)
    !k = one
    do j = 1,4
    do n = 1,num_BCs
      call bnds(n)%enforce(soln)
    end do
    call prim2cons(soln%U,soln%V)
    call calc_flux_2D(grid,soln)
    call calc_residual(grid,soln)
    
    do i = 1,4
    soln%U(i,i_low:i_high,j_low:j_high) = &
         soln%U(i,i_low:i_high,j_low:j_high) - &
       k(j)*( soln%R(i,i_low:i_high,j_low:j_high)*  &
        soln%dt(i_low:i_high,j_low:j_high) )/  &
         grid%V(i_low:i_high,j_low:j_high)
    end do
        call update_states(soln)
    end do
    
  end subroutine explicit_RK
  
  !============================= calc_residual ==============================80
  !>
  !! Description: 
  !!
  !! Inputs:      grid : 
  !!              soln : 
  !!
  !! Outputs:     grid : 
  !!              soln : 
  !<
  !==========================================================================80
  subroutine calc_residual(grid,soln)
    
    type(grid_t), intent(inout) :: grid
    type(soln_t), intent(inout) :: soln
    integer :: i,j
    do j = j_low,j_high
    do i = i_low,i_high
      soln%R(:,i,j) = grid%A_xi(i+1,j)*soln%Fxi(:,i,j) &
                    - grid%A_xi(i,j)*soln%Fxi(:,i-1,j)  &
                    + grid%A_eta(i,j+1)*soln%Feta(:,i,j) &
                    - grid%A_eta(i,j)*soln%Feta(:,i,j-1) &
                    - grid%V(i,j)*soln%S(:,i,j)
    end do
    end do
  end subroutine calc_residual
  
  !============================= residual_norms ==============================80
  !>
  !! Description: 
  !!
  !! Inputs:      R : 
  !!
  !! Outputs:     Rnorm : 
  !<
  !===========================================================================80
  subroutine residual_norms( R, Rnorm, rinit )
    
    real(prec), dimension(neq,i_low:i_high,j_low:j_high), intent(in) :: R
    real(prec), dimension(neq), intent(in)  :: rinit
    real(prec), dimension(neq), intent(out) :: Rnorm
    real(prec) :: Linv
    integer :: i
    
    Linv = one/real(size(R(1,:,:)),prec)
    
    do i = 1,neq
      Rnorm(i) = sqrt(Linv*sum(R(i,:,:)**2))
    end do
    Rnorm = Rnorm/rinit
     
  end subroutine residual_norms

end module time_integration
