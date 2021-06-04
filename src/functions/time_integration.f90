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
    
    low1  = lbound(V,2)
    high1 = ubound(V,2)
    low2  = lbound(V,3)
    high2 = ubound(V,3)


    call speed_of_sound(V(4,:,:),V(1,:,:),asnd)
    lambda_xi  = abs(V(2,:,:)*n_xi_avg(1,:,:) + &
                     V(3,:,:)*n_xi_avg(2,:,:)) + asnd
    lambda_eta = abs(V(2,:,:)*n_eta_avg(1,:,:) + &
                     V(3,:,:)*n_eta_avg(2,:,:)) + asnd
    dt = CFL*vol/(lambda_xi*A_xi(low1:high1,low2:high2) + &
                  lambda_eta*A_eta(low1:high1,low2:high2))
    dt = CFL*minval(dt)
    
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

  subroutine explicit_RK( grid, soln)
    
    type(grid_t), intent(inout) :: grid
    type(soln_t), intent(inout) :: soln
    integer :: i, j
    
    call calc_residual(grid,soln)
    !call calc_residual(grid%A_xi(i_low:i_high+1,j_low:j_high),&
    !                  grid%A_eta(i_low:i_high,j_low:j_high+1),&
    !                  grid%V(i_low:i_high,j_low:j_high),&
    !                  soln%S(i_low:i_high,j_low:j_high,:),&
    !                  soln%F(i_low:i_high+1,j_low:j_high+1,:,:),&
    !                  soln%R(i_low:i_high,j_low:j_high,:))
    do i = 1,4
    soln%U(i,i_low:i_high,j_low:j_high) = &
         soln%U(i,i_low:i_high,j_low:j_high) - &
       ( soln%R(i,i_low:i_high,j_low:j_high)*  &
        soln%dt(i_low:i_high,j_low:j_high) )/  &
         grid%V(i_low:i_high,j_low:j_high)
    end do
    
  end subroutine explicit_RK
  !subroutine explicit_RK( grid, soln )
  !  
  !  type(grid_t), intent(inout) :: grid
  !  type(soln_t), intent(inout) :: soln
  !  !real(prec), dimension(ig_low:ig_high,jg_low:jg_high,neq),&
  !  !                                           intent(inout) :: U
  !  !real(prec), dimension(i_low:i_high,j_low:j_high,neq),&
  !  !                                           intent(inout) :: R
  !  !real(prec), dimension(i_low:i_high+1,j_low:j_high+1,neq,2),&
  !  !                                           intent(in) :: F
  !  !real(prec), dimension(i_low:i_high,j_low:j_high,neq),&
  !  !                                           intent(in) :: S
  !  !real(prec), dimension(i_low:i_high,j_low:j_high), intent(in) :: dt
  !  !integer, intent(in) :: N
  !  real(prec), dimension(4) :: k
  !  integer :: i, j
  !  
  !  k = (/ fourth, third, half, one /)
  !  do j = 1,4
  !  call calc_residual(grid%A_xi(i_low:i_high+1,j_low:j_high),&
  !                    grid%A_eta(i_low:i_high,j_low:j_high+1),&
  !                    grid%V(i_low:i_high,j_low:j_high),&
  !                    soln%S(i_low:i_high,j_low:j_high,:),&
  !                    soln%F(i_low:i_high+1,j_low:j_high+1,:,:),&
  !                    soln%R(i_low:i_high,j_low:j_high,:))
  !  do i = 1,4
  !  soln%U(i_low:i_high,j_low:j_high,i) = &
  !       soln%U(i_low:i_high,j_low:j_high,i) - &
  !     k(j)*( soln%R(i_low:i_high,j_low:j_high,i)*  &
  !      soln%dt(i_low:i_high,j_low:j_high) )/  &
  !       grid%V(i_low:i_high,j_low:j_high)
  !  end do
  !  end do
  !  
  !end subroutine explicit_RK
  
  subroutine calc_residual(grid,soln)
    
    type(grid_t), intent(inout) :: grid
    type(soln_t), intent(inout) :: soln
    integer :: i,j
    do j = j_low,j_high
    do i = i_low,i_high
      soln%R(i,j,:) = grid%A_xi(i+1,j)*soln%Fxi(:,i+1,j) &
                    - grid%A_xi(i,j)*soln%Fxi(:,i,j)  &
                    + grid%A_eta(i,j+1)*soln%Feta(:,i,j+1) &
                    - grid%A_eta(i,j)*soln%Feta(:,i,j) &
                    - grid%V(i,j)*soln%S(:,i,j)
    end do
    end do
  end subroutine calc_residual
  !***subroutine calc_residual(A_xi,A_eta,V,S,F,R)
  !  
  !  real(prec), dimension(:,:), intent(in) :: A_xi
  !  real(prec), dimension(:,:), intent(in) :: A_eta
  !  real(prec), dimension(:,:), intent(in) :: V
  !  real(prec), dimension(:,:,:), intent(in) :: S
  !  real(prec), dimension(:,:,:,:), intent(in) :: F
  !  real(prec), dimension(:,:,:), intent(inout) :: R
  !  integer :: i,j,lowi,highi,lowj,highj
  !  lowi  = lbound(R,1)
  !  highi = ubound(R,1)
  !  lowj  = lbound(R,2)
  !  highj = ubound(R,2)
  !  do j = lowj,highj
  !  do i = lowi,highi
  !    !R(i,j,:) = A_xi(i+1,j)*F(i+1,j,:,1) - A_xi(i,j)*F(i,j,:,1)  &
  !    !         + A_eta(i,j+1)*F(i,j+1,:,2) - A_eta(i,j)*F(i,j,:,2) &
  !    !         - V(i,j)*S(i,j,:)
  !    R(i,j,:) = A_xi(i+1,j)*F(i+1,j,:,1) - A_xi(i,j)*F(i,j,:,1)  &
  !             + A_eta(i,j+1)*F(i,j+1,:,2) - A_eta(i,j)*F(i,j,:,2) &
  !             - V(i,j)*S(i,j,:)
  !  end do
  !  end do
  !end subroutine calc_residual
  !subroutine calc_residual(A_xi,A_eta,V,S,F,R)
  !  
  !  real(prec), dimension(ig_low:ig_high+1,jg_low:jg_high),&
  !                                                intent(in) :: A_xi
  !  real(prec), dimension(ig_low:ig_high,jg_low:jg_high+1),&
  !                                                intent(in) :: A_eta
  !  real(prec), dimension(ig_low:ig_high,jg_low:jg_high),&
  !                                                intent(in) :: V
  !  real(prec), dimension(i_low:i_high,j_low:j_high,neq),&
  !                                                intent(in) :: S
  !  real(prec), dimension(i_low:i_high+1,j_low:j_high+1,neq,2),&
  !                                                intent(in) :: F
  !  real(prec), dimension(i_low:i_high,j_low:j_high,neq),&
  !                                               intent(inout) :: R
  !  
  !  integer :: i,j
  !  
  !  do j = j_low,j_high
  !  do i = i_low,i_high
  !    R(i,j,:) = A_xi(i+1,j)*F(i+1,j,:,1) - A_xi(i,j)*F(i,j,:,1)  &
  !             + A_eta(i,j+1)*F(i,j+1,:,2) - A_eta(i,j)*F(i,j,:,2) &
  !             - V(i,j)*S(i,j,:)
  !  end do
  !  end do
  !end subroutine calc_residual
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
    
    real(prec), dimension(neq,i_low:i_high,j_low:j_high), intent(in) :: R
    real(prec), dimension(neq,3), intent(in)  :: rinit
    real(prec), dimension(neq,3), intent(out) :: Rnorm
    integer, intent(in) :: pnorm
    real(prec) :: Linv
    integer :: i
    
    Linv = one/real(size(R(1,:,:)),prec)
    
    do i = 1,neq
      Rnorm(i,1) = maxval(abs(R(i,:,:)))
      Rnorm(i,2) = Linv*sum(abs(R(i,:,:)))
      Rnorm(i,3) = sqrt(Linv*sum(R(i,:,:)**2))
    end do
    Rnorm = Rnorm/rinit
     
  end subroutine residual_norms

end module time_integration
