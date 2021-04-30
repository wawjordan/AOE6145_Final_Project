module time_integration
  
  use set_precision, only : prec
  use set_constants, only : one, half
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
  subroutine calc_time_step( dx, V, asnd, lambda, dt )
    
    use set_inputs,          only : CFL
    
    real(prec), dimension(ig_low:ig_high,neq), intent(in)  :: V
    real(prec), intent(in) :: dx
    real(prec), dimension(ig_low:ig_high),   intent(out) :: lambda
    real(prec), dimension(i_low:i_high),   intent(out) :: dt
    
    real(prec), dimension(ig_low:ig_high), intent(in)    :: asnd
    
    lambda(:) = abs(V(:,2)) + asnd
    dt(:) = CFL*dx/lambda(i_low:i_high)
    !dt(:) = minval(CFL*dx/lambda(i_low:i_high))
    
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

  subroutine explicit_euler( grid, src, dt, F, U, R )
    
    type(grid_t), intent(in) :: grid
    
    real(prec), dimension(ig_low:ig_high,jg_low:jg_high,neq),&
                                               intent(inout) :: U
    real(prec), dimension(i_low:i_high,j_low:j_high,neq),&
                                               intent(out) :: R
    real(prec), dimension(i_low:i_high+1,j_low:j_high+1,neq),&
                                               intent(in) :: F
    real(prec), dimension(i_low:i_high,j_low:j_high,neq),&
                                               intent(in) :: src
    real(prec), dimension(i_low:i_high,j_low:j_high), intent(in) :: dt
    
    integer :: i
    
  end subroutine explicit_euler
  
  subroutine calc_residual(A_xi,A_eta,V,n_xi,n_eta,S,F,U,R)
    
    real(prec), dimension(ig_low:ig_high+1,jg_low:jg_high),&
                                                  intent(in) :: A_xi
    real(prec), dimension(ig_low:ig_high,jg_low:jg_high+1),&
                                                  intent(in) :: A_eta
    real(prec), dimension(i_low:i_high,j_low:j_high),&
                                                  intent(in) :: V
    real(prec), dimension(ig_low:ig_high+1,jg_low:jg_high,2),&
                                                  intent(in) :: n_xi
    real(prec), dimension(ig_low:ig_high,jg_low:jg_high+1,2),&
                                                  intent(in) :: n_eta
    real(prec), dimension(i_low:i_high,j_low:j_high,neq),&
                                                  intent(in) :: S
    real(prec), dimension(i_low:i_high+1,j_low:j_high+1,neq),&
                                                  intent(in) :: F
    real(prec), dimension(ig_low:ig_high,jg_low:jg_high,neq),&
                                               intent(inout) :: U
    real(prec), dimension(i_low:i_high,j_low:j_high,neq),&
                                                 intent(out) :: R
    
    integer :: i,j
    
    do j = j_low,j_high
    do i = i_low,i_high
      R(i,j,:) = A_xi(i+1,j)*F(i+1,j,:)  + A_xi(i,j)*F(i,j,:)  &
               + A_eta(i,j+1)*F(i,j+1,:) + A_eta(i,j)*F(i,j,:) &
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
