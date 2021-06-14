module other_subroutines
  
  use set_precision, only : prec
  use set_constants, only : zero, one, two, three, half, fourth
  use set_inputs, only : imax, neq, i_low, i_high, ig_low, ig_high
  use set_inputs, only : jmax, j_low, j_high, jg_low, jg_high
  use set_inputs, only : epsM, kappaM, isMMS, n_ghost
  use fluid_constants, only : gamma
  use variable_conversion, only : speed_of_sound
  use limiter_calc, only : limiter_fun
  use soln_type, only : soln_t
  use grid_type, only : grid_t
  
  implicit none
  
  private
  
  public :: MUSCL_extrap, surface_MUSCL_extrap, calc_de, calc_sources, &
            airfoil_forces
  
  contains
  
  !============================= calculate_sources ==========================80
  !>
  !! Description: 
  !!
  !! Inputs:      V  : 
  !!              y  : 
  !!
  !! Outputs:     S  : 
  !<
  !==========================================================================80
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
  
  !================================ MUSCL_extrap ============================80
  !>
  !! Description: 
  !!
  !! Inputs:      V     : 
  !!
  !! Outputs:     left  : 
  !!              right :
  !<
  !==========================================================================80
  subroutine MUSCL_extrap(soln,Lxi,Rxi,Leta,Reta)
    type(soln_t), intent(inout) :: soln
    real(prec), dimension(neq,i_low-1:i_high,j_low:j_high) :: Lxi, Rxi
    real(prec), dimension(neq,i_low:i_high,j_low-1:j_high) :: Leta, Reta
    integer :: i
    
    Lxi  = soln%V(:,i_low-1:i_high  ,j_low:j_high) + fourth*epsM*(      &
           (one-kappaM)*soln%psi_p_xi(:,i_low-2:i_high-1,j_low:j_high)* &
                             ( soln%V(:,i_low-1:i_high  ,j_low:j_high)  &
                             - soln%V(:,i_low-2:i_high-1,j_low:j_high) )&
         + (one+kappaM)*soln%psi_m_xi(:,i_low-1:i_high  ,j_low:j_high)* &
                             ( soln%V(:,i_low  :i_high+1,j_low:j_high)  &
                             - soln%V(:,i_low-1:i_high  ,j_low:j_high) ))
    Rxi  = soln%V(:,i_low  :i_high+1,j_low:j_high) + fourth*epsM*(      &
           (one-kappaM)*soln%psi_p_xi(:,i_low  :i_high+1,j_low:j_high)* &
                             ( soln%V(:,i_low+1:i_high+2,j_low:j_high)  &
                             - soln%V(:,i_low  :i_high+1,j_low:j_high) )&
         + (one+kappaM)*soln%psi_m_xi(:,i_low-1:i_high  ,j_low:j_high)* &
                             ( soln%V(:,i_low  :i_high+1,j_low:j_high)  &
                             - soln%V(:,i_low-1:i_high  ,j_low:j_high) ))

    Leta = soln%V(:,i_low:i_high,j_low-1:j_high  ) + fourth*epsM*(       &
           (one-kappaM)*soln%psi_p_eta(:,i_low:i_high,j_low-2:j_high-1)* &
                              ( soln%V(:,i_low:i_high,j_low-1:j_high  )  &
                              - soln%V(:,i_low:i_high,j_low-2:j_high-1) )&
         + (one+kappaM)*soln%psi_m_eta(:,i_low:i_high,j_low-1:j_high  )* &
                              ( soln%V(:,i_low:i_high,j_low  :j_high+1)  &
                              - soln%V(:,i_low:i_high,j_low-1:j_high  ) ))
    Reta = soln%V(:,i_low:i_high,j_low  :j_high+1) + fourth*epsM*(      &
           (one-kappaM)*soln%psi_p_eta(:,i_low:i_high,j_low  :j_high+1)* &
                              ( soln%V(:,i_low:i_high,j_low+1:j_high+2)  &
                              - soln%V(:,i_low:i_high,j_low  :j_high+1) )&
         + (one+kappaM)*soln%psi_m_eta(:,i_low:i_high,j_low-1:j_high  )* &
                              ( soln%V(:,i_low:i_high,j_low  :j_high+1)  &
                              - soln%V(:,i_low:i_high,j_low-1:j_high  ) ))
    
  end subroutine MUSCL_extrap
  
  subroutine surface_MUSCL_extrap(soln,Left,Right,indices)
    type(soln_t), intent(in) :: soln
    real(prec), dimension(:,:) :: Left, right
    integer, dimension(4), intent(in) :: indices
    integer :: low, high, ind, dir
    low  = indices(1)
    high = indices(2)
    ind  = indices(3)
    dir  = indices(4)
    if (dir==2) then ! eta-face (j-constant)
      Left  = soln%V(:,low:high,ind-1) + fourth*epsM*(    &
           (one-kappaM)*soln%psi_p_eta(:,low:high,ind-2)* &
                              ( soln%V(:,low:high,ind-1)  &
                              - soln%V(:,low:high,ind-2) )&
         + (one+kappaM)*soln%psi_m_eta(:,low:high,ind-1)* &
                              ( soln%V(:,low:high,ind  )  &
                              - soln%V(:,low:high,ind-1) ))
      Right = soln%V(:,low:high,ind  ) + fourth*epsM*(    &
           (one-kappaM)*soln%psi_p_eta(:,low:high,ind  )* &
                              ( soln%V(:,low:high,ind+1)  &
                              - soln%V(:,low:high,ind  ) )&
         + (one+kappaM)*soln%psi_m_eta(:,low:high,ind-1)* &
                              ( soln%V(:,low:high,ind  )  &
                              - soln%V(:,low:high,ind-1) ))
    elseif (dir==1) then ! xi-face (i-constant)
      Left  = soln%V(:,ind-1,low:high) + fourth*epsM*(    &
           (one-kappaM)*soln%psi_p_xi(:,ind-2,low:high)*  &
                             ( soln%V(:,ind-1,low:high)   &
                             - soln%V(:,ind-2,low:high))  &
         + (one+kappaM)*soln%psi_m_xi(:,ind-1,low:high)*  &
                             ( soln%V(:,ind,  low:high)   &
                             - soln%V(:,ind-1,low:high) ) )
      Right = soln%V(:,ind  ,low:high) + fourth*epsM*(    &
           (one-kappaM)*soln%psi_p_xi(:,ind  ,low:high)*  &
                             ( soln%V(:,ind+1,low:high)   &
                             - soln%V(:,ind  ,low:high))  &
         + (one+kappaM)*soln%psi_m_xi(:,ind-1,low:high)*  &
                             ( soln%V(:,ind  ,low:high)   &
                             - soln%V(:,ind-1,low:high) ) )
    end if
    
  end subroutine surface_MUSCL_extrap
  
  subroutine airfoil_forces(soln,grid,indices,Cp,Cl,Cd)
    use set_constants, only : pi 
    use set_inputs, only : p_inf, rho_inf, u_inf, alpha
    type(soln_t), intent(in) :: soln
    type(grid_t), intent(in) :: grid
    integer, dimension(4), intent(in) :: indices
    real(prec), dimension(:,:), allocatable, intent(out) :: Cp
    real(prec), dimension(:,:), allocatable :: Left, Right
    real(prec), intent(out) :: Cl, Cd
    real(prec) :: chord, Fx, Fy, Lprime, Dprime, h
    integer :: i, low, high, ind, dir
    
    
    low  = indices(1)
    high = indices(2)
    ind  = indices(3)
    dir  = indices(4)
    chord = maxval(abs(grid%x(low:high,ind) - grid%x(low,ind) ) )
    allocate( Left(neq,low:high), Right(neq,low:high), Cp(3,low:high) )
    call surface_MUSCL_extrap(soln,Left,Right,indices)
    Fx = zero
    Fy = zero
    do i = low,high
      h = (grid%x(i+1,ind) - grid%x(i,ind))
      Cp(1,i) = half*(grid%x(i+1,ind) + grid%x(i,ind))
      Cp(2,i) = Right(4,i)
      Cp(3,i) = two*(Right(4,i) - p_inf)/(rho_inf*u_inf**2)
      Fx = Fx + h*(grid%n_eta(i,ind,1)*grid%A_eta(i,ind)*Right(4,i))
      Fy = Fy + h*(grid%n_eta(i,ind,2)*grid%A_eta(i,ind)*Right(4,i))
    end do
    Lprime = Fy*cos((pi/180.0_prec)*alpha) - Fx*sin((pi/180.0_prec)*alpha)
    Dprime = Fy*sin((pi/180.0_prec)*alpha) + Fx*cos((pi/180.0_prec)*alpha)
    Cl = two*Lprime/(chord*rho_inf*u_inf**2)
    Cd = two*Dprime/(chord*rho_inf*u_inf**2)
    deallocate(left, right)
  end subroutine airfoil_forces
  
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
  !==========================================================================80
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
     
  
  !================================== calc_de ===============================80
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
  !==========================================================================80
  subroutine calc_de( soln, DE, DEnorm, cons )
    
    type(soln_t), intent(in) :: soln
    logical, intent(in) :: cons
    real(prec), dimension(:,:,:), intent(out) :: DE
    real(prec), dimension(neq,3), intent(out) :: DEnorm
    real(prec) :: Linv
    integer :: i
    !if (cons) then
    !  DE = soln%U(i_low:i_high,j_low:j_high,1:neq) &
    !   & - soln%Umms(i_low:i_high,j_low:j_high,1:neq)
    !else
    !  DE = soln%V(i_low:i_high,j_low:j_high,1:neq) &
    !   & - soln%Vmms(i_low:i_high,j_low:j_high,1:neq)
    !end if
    DE = zero
    Linv = one/real( (i_high-i_low)*(j_high-j_low),prec )
    if (cons) then
      DE = soln%U - soln%Umms
    else
      DE = soln%V - soln%Vmms
    end if
    do i = 1,neq
      DEnorm(i,1) = Linv*sum( abs( DE(i,i_low:i_high,j_low:j_high) ) )
      DEnorm(i,2) = sqrt( Linv*sum( DE(i,i_low:i_high,j_low:j_high)**2 ) )
      DEnorm(i,3) = maxval( abs( DE(i,i_low:i_high,j_low:j_high) ) )
    end do
    
  end subroutine calc_de
  
end module other_subroutines
