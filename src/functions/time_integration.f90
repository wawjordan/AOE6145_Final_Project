module time_integration
  
  use set_precision, only : prec
  use set_constants, only : one, half, third, fourth
  use set_inputs,    only : neq, i_low, i_high, ig_low, ig_high
  use set_inputs,    only : j_low, j_high, jg_low, jg_high
  use grid_type
  use soln_type
  use variable_conversion
  
  implicit none
  
  contains
  
  !============================= calc_time_step ==============================80
  !>
  !! Description: 
  !!
  !! Inputs:      dx : 
  !!              V  : 
  !!
  !! Outputs:     lambda : 
  !!              dt     : 
  !<
  !===========================================================================80
  subroutine calc_time_step( A_xi, A_eta, n_xi_avg, n_eta_avg, vol, V, dt )
    
    use set_inputs,          only : CFL
    
    real(prec), dimension(:,:,:), intent(in)  :: V
    real(prec), dimension(:,:), intent(in)    :: A_xi, A_eta, vol
    real(prec), dimension(:,:,:), intent(in)  :: n_xi_avg, n_eta_avg
    real(prec), dimension(lbound(V,1):ubound(V,1),lbound(V,2):ubound(V,2)) &
                                         :: asnd, lambda_xi, lambda_eta
    real(prec), dimension(:,:),   intent(out) :: dt
    integer :: low1, low2, high1, high2
    
    low1  = lbound(V,1)
    high1 = ubound(V,1)
    low2  = lbound(V,2)
    high2 = ubound(V,2)


    call speed_of_sound(V(:,:,4),V(:,:,1),asnd)
    lambda_xi  = abs(V(:,:,2)*n_xi_avg(:,:,1) + &
                     V(:,:,3)*n_xi_avg(:,:,2)) + asnd
    lambda_eta = abs(V(:,:,2)*n_eta_avg(:,:,1) + &
                     V(:,:,3)*n_eta_avg(:,:,2)) + asnd
    dt = CFL*vol/(lambda_xi*A_xi(low1:high1,low2:high2) + &
                  lambda_eta*A_eta(low1:high1,low2:high2))
    dt(:,:) = minval(dt)
    
  end subroutine calc_time_step
  
  !============================= explicit_euler ==============================80
  !>
  !! Description: 
  !!
  !! Inputs:      grid : 
  !!              S    : 
  !!              F    :
  !!              dt   : 
  !!              U    :
  !!
  !! Outputs:     lambda : 
  !!              dt     : 
  !!              R      : 
  !!              U      : 
  !<
  !===========================================================================80

  subroutine explicit_RK( grid, S, dt, F, U, R, N )
    
    type(grid_t), intent(in) :: grid
    real(prec), dimension(ig_low:ig_high,jg_low:jg_high,neq),&
                                               intent(inout) :: U
    real(prec), dimension(i_low:i_high,j_low:j_high,neq),&
                                               intent(out) :: R
    real(prec), dimension(i_low:i_high+1,j_low:j_high+1,neq,2),&
                                               intent(in) :: F
    real(prec), dimension(i_low:i_high,j_low:j_high,neq),&
                                               intent(in) :: S
    real(prec), dimension(i_low:i_high,j_low:j_high), intent(in) :: dt
    integer, intent(in) :: N
    real(prec), dimension(4) :: k
    integer :: i, j
    
    k = (/ fourth, third, half, one /)
    
    call calc_residual(grid%A_xi,grid%A_eta,grid%V,S,F,R)
    do j = 1,N
    do i = 1,neq
    U(i_low:i_high,j_low:j_high,i) = U(i_low:i_high,j_low:j_high,i) &
    - k(j)*R(:,:,i)*dt/grid%V(i_low:i_high,j_low:j_high)
    end do
    end do
    
  end subroutine explicit_RK
  
  subroutine calc_residual(A_xi,A_eta,V,S,F,R)
    
    real(prec), dimension(ig_low:ig_high+1,jg_low:jg_high),&
                                                  intent(in) :: A_xi
    real(prec), dimension(ig_low:ig_high,jg_low:jg_high+1),&
                                                  intent(in) :: A_eta
    real(prec), dimension(i_low:i_high,j_low:j_high),&
                                                  intent(in) :: V
    real(prec), dimension(i_low:i_high,j_low:j_high,neq),&
                                                  intent(in) :: S
    real(prec), dimension(i_low:i_high+1,j_low:j_high+1,neq,2),&
                                                  intent(in) :: F
    real(prec), dimension(i_low:i_high,j_low:j_high,neq),&
                                                 intent(out) :: R
    
    integer :: i,j
    
    do j = j_low,j_high
    do i = i_low,i_high
      R(i,j,:) = A_xi(i+1,j)*F(i+1,j,:,1)  + A_xi(i,j)*F(i,j,:,1)  &
               + A_eta(i,j+1)*F(i,j+1,:,2) + A_eta(i,j)*F(i,j,:,2) &
               - V(i,j)*S(i,j,:)
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
  subroutine residual_norms( R, Rnorm, pnorm, rinit )
    
    real(prec), dimension(i_low:i_high,1:neq), intent(in) :: R
    real(prec), dimension(1,1:neq), intent(in)  :: rinit
    real(prec),  intent(inout) :: Rnorm(1,1:neq)
    integer :: pnorm
    
    if (pnorm == 0) then
      Rnorm(1,1:neq) = maxval(abs(R),1)
    elseif (pnorm == 1) then
      Rnorm(1,1) = (one/real(size(R,1)))*sum(abs(R(:,1)),1)
      Rnorm(1,2) = (one/real(size(R,1)))*sum(abs(R(:,2)),1)
      Rnorm(1,3) = (one/real(size(R,1)))*sum(abs(R(:,3)),1)
    elseif (pnorm == 2) then 
      Rnorm(1,1) = sqrt((one/real(size(R,1)))*sum(R(:,1)**2,1))
      Rnorm(1,2) = sqrt((one/real(size(R,1)))*sum(R(:,2)**2,1))
      Rnorm(1,3) = sqrt((one/real(size(R,1)))*sum(R(:,3)**2,1))
    else
      Rnorm(1,1) = maxval(abs(R(:,1)))
      Rnorm(1,2) = maxval(abs(R(:,2)))
      Rnorm(1,3) = maxval(abs(R(:,3)))
    end if
     Rnorm(1,1) = Rnorm(1,1)/rinit(1,1)
     Rnorm(1,2) = Rnorm(1,2)/rinit(1,2)
     Rnorm(1,3) = Rnorm(1,3)/rinit(1,3)
     
  end subroutine residual_norms

end module time_integration
