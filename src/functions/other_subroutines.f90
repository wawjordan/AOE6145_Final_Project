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
  
  public :: output_file_headers, output_soln, MUSCL_extrap, calc_de, output_res
  public :: calc_sources, output_exact_soln, Limit
  
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
!  subroutine MUSCL_extrap( V, psi_plus, psi_minus, left, right, i, j )
    
!    use set_inputs, only : limiter_freeze
!    real(prec), dimension(:,:,:), intent(in)  :: V
!    real(prec), dimension(:,:,:,:), intent(inout)  :: psi_plus, psi_minus
!    real(prec), dimension(:,:,:,:), intent(out) :: left, right
!    real(prec), dimension(lbound(V,1):ubound(V,1), &
!                          lbound(V,2):ubound(V,2),neq) :: r_plus, r_minus
!    integer, dimension(:), intent(in) :: i,j
    !integer :: i,j, low, high
    
    !low = lbound(V,1)+n_ghost
    !high = ubound(V,1)-n_ghost
    
!    r_plus = zero
!    r_minus = zero
!    if (limiter_freeze) then
!      continue
!    else
!      call calc_consecutive_variations(V,r_plus,r_minus,1)
!      call limiter_fun(r_plus,psi_plus(:,:,:,1))
!      call limiter_fun(r_minus,psi_minus(:,:,:,1))
!      
!      call calc_consecutive_variations(V,r_plus,r_minus,2)
!      call limiter_fun(r_plus,psi_plus(:,:,:,2))
!      call limiter_fun(r_minus,psi_minus(:,:,:,2))
!    end if
    
!    left(:,:,:,1)  = V(i,j,:) + fourth*epsM*( &
!         & (one-kappaM)*psi_plus(i-1,j,:,1)*(V(i,j,:)-V(i-1,j,:)) + &
!         & (one+kappaM)*psi_minus(i,j,:,1)*(V(i+1,j,:)-V(i,j,:)) )
!    right(:,:,:,1) = V(i+1,j,:) - fourth*epsM*( &
!         & (one+kappaM)*psi_minus(i+1,j,:,1)*(V(i+1,j,:)-V(i,j,:)) + &
!         & (one-kappaM)*psi_plus(i,j,:,1)*(V(i+2,j,:)-V(i+1,j,:)) )
!    
!    left(:,:,:,2)  = V(i,j,:) + fourth*epsM*( &
!         & (one-kappaM)*psi_plus(i,j-1,:,2)*(V(i,j,:)-V(i,j-1,:)) + &
!         & (one+kappaM)*psi_minus(i,j,:,2)*(V(i,j+1,:)-V(i,j,:)) )
!    right(:,:,:,2) = V(i,j+1,:) - fourth*epsM*( &
!         & (one+kappaM)*psi_minus(i,j+1,:,2)*(V(i,j+1,:)-V(i,j,:)) + &
!         & (one-kappaM)*psi_plus(i,j,:,2)*(V(i,j+2,:)-V(i,j+1,:)) )
!left(:,:,:,1) = V(i-1,j,:) 
!right(:,:,:,1) = V(i,j,:) 
!left(:,:,:,2) = V(i,j-1,:) 
!right(:,:,:,2) = V(i,j,:)


!    do i = low-1,high
!      j = i-low+2
!      left(j,:) = V(i,:) + fourth*epsM*( &
!         & (one-kappaM)*psi_plus(i-1,:)*(V(i,:)-V(i-1,:)) + &
!         & (one+kappaM)*psi_minus(i,:)*(V(i+1,:)-V(i,:)) )
!      right(j,:) = V(i+1,:) - fourth*epsM*( &
!         & (one+kappaM)*psi_minus(i+1,:)*(V(i+1,:)-V(i,:)) + &
!         & (one-kappaM)*psi_plus(i,:)*(V(i+2,:)-V(i+1,:)) )
!    end do

!  end subroutine MUSCL_extrap
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

  !========================== output_file_headers ============================80
  !>
  !! Description: 
  !!
  !! Inputs:      grid     : 
  !!              soln     : 
  !!              num_iter : 
  !<
  !===========================================================================80
  subroutine output_file_headers
    
    use set_inputs, only : flux_scheme, limiter_scheme, beta_lim, cfl
    use set_inputs, only : grid_name
    
    character(len=1024) :: dirname_hist
    character(len=1024) :: dirname_field
    character(len=1024) :: filename
    character(len=64) ::   ncells_str
    character(len=64) ::   MMS_str
    character(len=64) ::   CFL_str
    character(len=64) ::   flux_str
    character(len=64) ::   order_str
    character(len=64) ::   limiter_str
    character(len=64) ::   kappa_str
    integer :: resunit, fldunit
    
    resunit = 30
    fldunit = 40
    
    ! Set up output directories
!    write (ncells_str  , "(A5,I0.3,A6,I0.3)") "imax=",imax,"_jmax=",jmax
!    write (CFL_str  , "(A4,F0.4)") "CFL=",cfl
!    
!    select case(flux_scheme)
!    case(1)
!      write (flux_str,"(A3)") "V-L"
!    case(2)
!      write (flux_str,"(A3)") "ROE"
!    case default
!    end select
!    if (nint(epsM).eq.0) then
!      write (order_str,"(A9)") "1st-order"
!    else
!      write (order_str,"(A9)") "2nd-order"
!    end if
!    
!    select case(limiter_scheme)
!    case(1)
!      write (limiter_str,"(A3)") "_VL"
!    case(2)
!      write (limiter_str,"(A3)") "_VA"
!    case(3)
!      write (limiter_str,"(A3)") "_MM"
!    case(4)
!      write (limiter_str,"(A2,I0.2,A1)") "_B",int(10*beta_lim),"_"
!    case default
!    end select
!    
!    write (kappa_str, "(A2,I3.2,A1)") "_K"  , int(10*kappaM),"_"
!    write (dirname_hist, *) adjustl(trim(MMS_str)),"/hist/"
!    write (dirname_field, *) adjustl(trim(MMS_str)),"/field/"
    !write (filename,*)  trim(ncells_str)//  &
    !&                   trim(flux_str)//    &
    !&                   trim(order_str)//   &
    !&                   trim(limiter_str)// &
    !&                   trim(kappa_str)//   &
    !&                   trim(CFL_str)
    
    !call execute_command_line ('mkdir -p ../results/' // &
    !& adjustl(trim(dirname_hist)))
    
    !call execute_command_line ('mkdir -p ../results/' // &
    !& adjustl(trim(dirname_field)))
    ! Set up output files (history and solution)
    !open(30,file='history.dat',status='unknown')
    !open(resunit,file= '../results/'//trim(adjustl(dirname_hist))//  &
    !&               trim(adjustl(filename))//'_history.dat',status='unknown')
    open(resunit,file= '../results/example_res.dat')
    if(isMMS) then
      write(resunit,*) 'TITLE = "MMS ('// &
      & trim(adjustl(grid_name))//')"'
      write(resunit,*) 'variables="Iteration"   &
      &   "R<sub>1</sub>""R<sub>2</sub>"   &
      &   "R<sub>3</sub>""R<sub>4</sub>"   &
      &   "DE<sub>1</sub>""DE<sub>2</sub>" &
      &   "DE<sub>3</sub>""DE<sub>4</sub>"'
    else
      write(resunit,*) 'TITLE = "'// &
      & trim(adjustl(grid_name))//'"'
      write(resunit,*) 'variables="Iteration"   &
      &   "R<sub>1</sub>""R<sub>2</sub>"   &
      &   "R<sub>3</sub>""R<sub>4</sub>"'
    end if
    
    open(fldunit,file= '../results/example_out.dat')
    write(fldunit,*) 'TITLE = "blahblahblah"'
    if(isMMS) then
      write(fldunit,*) 'variables="x""y""<greek>r</greek>""u""v""p""M"&
      & "DE<sub>1</sub>""DE<sub>2</sub>""DE<sub>3</sub>""DE<sub>4</sub>"'
    else
      write(fldunit,*) 'variables="x""y""<greek>r</greek>""u""v""p""M"'
    endif
    
  end subroutine output_file_headers
  
  !======================= output_exact_soln ================================80
  !>
  !! Description: 
  !!
  !! Inputs:      grid     : 
  !!              soln     : 
  !!              num_iter : 
  !<
  !==========================================================================80
  subroutine output_exact_soln(grid,soln)
    
    type( grid_t ), intent(in) :: grid
    type( soln_t ), intent(in) :: soln
    
    integer :: i, j, fldunit
    fldunit = 20
    
    open(fldunit,file= '../results/exact_soln.dat')
    write(fldunit,*) 'TITLE = "Manufactured Solution"'
    write(fldunit,*) 'variables="x""y"&
    & "<greek>r</greek><sub>exact</sub>"&
    & "u<sub>exact</sub>""v<sub>exact</sub>""p<sub>exact</sub>"&
    & "S<sub>1</sub>""S<sub>2</sub>""S<sub>3</sub>""S<sub>4</sub>"'
    write(fldunit,*) 'ZONE'
    write(fldunit,*) 'T= "1"'
    write(fldunit,*) 'I=',ig_high-ig_low+2,' J=',jg_high-jg_low+2
    write(fldunit,*) 'DT = ( DOUBLE, DOUBLE, DOUBLE, DOUBLE, &
                           & DOUBLE, DOUBLE, DOUBLE, DOUBLE, &
                           & DOUBLE, DOUBLE )'
    write(fldunit,*) 'DATAPACKING = BLOCK'
    write(fldunit,*) 'VARLOCATION = ([3-10]=CELLCENTERED)'
    
    write(fldunit,*) ((grid%x(i,j),i=ig_low,ig_high+1),    & !  1
                     NEW_LINE('a'), j=jg_low,jg_high+1)    
    write(fldunit,*) ((grid%y(i,j),i=ig_low,ig_high+1),    & !  2
                     NEW_LINE('a'), j=jg_low,jg_high+1)    
    write(fldunit,*) ((soln%Vmms(1,i,j),i=ig_low,ig_high), & !  3
                     NEW_LINE('a'), j=jg_low,jg_high)      
    write(fldunit,*) ((soln%Vmms(2,i,j),i=ig_low,ig_high), & !  4
                     NEW_LINE('a'), j=jg_low,jg_high)      
    write(fldunit,*) ((soln%Vmms(3,i,j),i=ig_low,ig_high), & !  5
                     NEW_LINE('a'), j=jg_low,jg_high)      
    write(fldunit,*) ((soln%Vmms(4,i,j),i=ig_low,ig_high), & !  6
                     NEW_LINE('a'), j=jg_low,jg_high)      
    write(fldunit,*) ((soln%Smms(1,i,j),i=ig_low,ig_high), & !  7
                     NEW_LINE('a'), j=jg_low,jg_high)      
    write(fldunit,*) ((soln%Smms(2,i,j),i=ig_low,ig_high), & !  8
                     NEW_LINE('a'), j=jg_low,jg_high)      
    write(fldunit,*) ((soln%Smms(3,i,j),i=ig_low,ig_high), & !  9
                     NEW_LINE('a'), j=jg_low,jg_high)      
    write(fldunit,*) ((soln%Smms(4,i,j),i=ig_low,ig_high), & ! 10
                     NEW_LINE('a'), j=jg_low,jg_high)      
    
  end subroutine output_exact_soln
  !subroutine output_exact_soln(grid,soln)
  !  
  !  type( grid_t ), intent(in) :: grid
  !  type( soln_t ), intent(in) :: soln
  !  
  !  integer :: i, j, fldunit
  !  fldunit = 20
  !  
  !  open(fldunit,file= '../results/exact_soln.dat')
  !  write(fldunit,*) 'TITLE = "Manufactured Solution"'
  !  write(fldunit,*) 'variables="x""y"&
  !  & "<greek>r</greek><sub>exact</sub>"&
  !  & "u<sub>exact</sub>""v<sub>exact</sub>""p<sub>exact</sub>"&
  !  & "S<sub>1</sub>""S<sub>2</sub>""S<sub>3</sub>""S<sub>4</sub>"'
  !  write(fldunit,*) 'ZONE'
  !  write(fldunit,*) 'T= "1"'
  !  write(fldunit,*) 'I=',i_high-i_low+2,' J=',j_high-j_low+2
  !  write(fldunit,*) 'DT = ( DOUBLE, DOUBLE, DOUBLE, DOUBLE, &
  !                         & DOUBLE, DOUBLE, DOUBLE, DOUBLE, &
  !                         & DOUBLE, DOUBLE )'
  !  write(fldunit,*) 'DATAPACKING = BLOCK'
  !  write(fldunit,*) 'VARLOCATION = ([3-10]=CELLCENTERED)'
  !  
  !  write(fldunit,*) ((grid%x(i,j),i=i_low,i_high+1),    & !  1
  !                   NEW_LINE('a'), j=j_low,j_high+1)    
  !  write(fldunit,*) ((grid%y(i,j),i=i_low,i_high+1),    & !  2
  !                   NEW_LINE('a'), j=j_low,j_high+1)    
  !  write(fldunit,*) ((soln%Vmms(i,j,1),i=i_low,i_high), & !  3
  !                   NEW_LINE('a'), j=j_low,j_high)      
  !  write(fldunit,*) ((soln%Vmms(i,j,2),i=i_low,i_high), & !  4
  !                   NEW_LINE('a'), j=j_low,j_high)      
  !  write(fldunit,*) ((soln%Vmms(i,j,3),i=i_low,i_high), & !  5
  !                   NEW_LINE('a'), j=j_low,j_high)      
  !  write(fldunit,*) ((soln%Vmms(i,j,4),i=i_low,i_high), & !  6
  !                   NEW_LINE('a'), j=j_low,j_high)      
  !  write(fldunit,*) ((soln%Smms(i,j,1),i=i_low,i_high), & !  7
  !                   NEW_LINE('a'), j=j_low,j_high)      
  !  write(fldunit,*) ((soln%Smms(i,j,2),i=i_low,i_high), & !  8
  !                   NEW_LINE('a'), j=j_low,j_high)      
  !  write(fldunit,*) ((soln%Smms(i,j,3),i=i_low,i_high), & !  9
  !                   NEW_LINE('a'), j=j_low,j_high)      
  !  write(fldunit,*) ((soln%Smms(i,j,4),i=i_low,i_high), & ! 10
  !                   NEW_LINE('a'), j=j_low,j_high)      
  !  
  !end subroutine output_exact_soln
  
  !============================= output_soln ================================80
  !>
  !! Description: 
  !!
  !! Inputs:      grid     : 
  !!              soln     : 
  !!              num_iter : 
  !<
  !==========================================================================80
  subroutine output_soln(grid,soln,num_iter)
    
    use set_inputs, only : counter
    
    type( grid_t ), intent(in) :: grid
    type( soln_t ), intent(in) :: soln
    integer,        intent(in) :: num_iter
    
    integer :: i, j, fldunit
    fldunit = 40
    
    open(fldunit,status='unknown')
    write(fldunit,*) 'ZONE'
    write(fldunit,*) 'T= "',counter,'"'
    write(fldunit,*) 'I=',ig_high-ig_low+2,' J=',jg_high-jg_low+2
    if(isMMS) then
      write(fldunit,*) 'DT = ( DOUBLE, DOUBLE, DOUBLE, DOUBLE, &
                             & DOUBLE, DOUBLE, DOUBLE, DOUBLE, &
                             & DOUBLE, DOUBLE, DOUBLE)'
      write(fldunit,*) 'DATAPACKING = BLOCK'
      write(fldunit,*) 'VARLOCATION = ([3-11]=CELLCENTERED)'
      
      write(fldunit,*) ((grid%x(i,j),i=ig_low,ig_high+1),    & !  1
                       NEW_LINE('a'), j=jg_low,jg_high+1)    
      write(fldunit,*) ((grid%y(i,j),i=ig_low,ig_high+1),    & !  2
                       NEW_LINE('a'), j=jg_low,jg_high+1)    
      write(fldunit,*) ((soln%V(1,i,j),i=ig_low,ig_high),    & !  3
                       NEW_LINE('a'), j=jg_low,jg_high)      
      write(fldunit,*) ((soln%V(2,i,j),i=ig_low,ig_high),    & !  4
                       NEW_LINE('a'), j=jg_low,jg_high)      
      write(fldunit,*) ((soln%V(3,i,j),i=ig_low,ig_high),    & !  5
                       NEW_LINE('a'), j=jg_low,jg_high)      
      write(fldunit,*) ((soln%V(4,i,j),i=ig_low,ig_high),    & !  6
                       NEW_LINE('a'), j=jg_low,jg_high)      
      write(fldunit,*) ((soln%mach(i,j),i=ig_low,ig_high),   & !  7
                       NEW_LINE('a'), j=jg_low,jg_high)      
      write(fldunit,*) ((soln%DE(1,i,j),i=ig_low,ig_high),   & !  8
                       NEW_LINE('a'), j=jg_low,jg_high)      
      write(fldunit,*) ((soln%DE(2,i,j),i=ig_low,ig_high),   & !  9
                       NEW_LINE('a'), j=jg_low,jg_high)      
      write(fldunit,*) ((soln%DE(3,i,j),i=ig_low,ig_high),   & ! 10
                       NEW_LINE('a'), j=jg_low,jg_high)      
      write(fldunit,*) ((soln%DE(4,i,j),i=ig_low,ig_high),   & ! 11
                       NEW_LINE('a'), j=jg_low,jg_high)      
      
    else
      
      write(fldunit,*) 'DT = (DOUBLE, DOUBLE, DOUBLE, DOUBLE,&
                           &  DOUBLE, DOUBLE, DOUBLE)'
      write(fldunit,*) 'DATAPACKING = BLOCK'
      write(fldunit,*) 'VARLOCATION = ([3-7]=CELLCENTERED'
      
      write(fldunit,*) ((grid%x(i,j),i=ig_low,ig_high+1),    & !  1
                       NEW_LINE('a'), j=jg_low,jg_high+1)    
      write(fldunit,*) ((grid%y(i,j),i=ig_low,ig_high+1),    & !  2
                       NEW_LINE('a'), j=jg_low,jg_high+1)    
      write(fldunit,*) ((soln%V(1,i,j),i=ig_low,ig_high),    & !  3
                       NEW_LINE('a'), j=jg_low,jg_high)      
      write(fldunit,*) ((soln%V(2,i,j),i=ig_low,ig_high),    & !  4
                       NEW_LINE('a'), j=jg_low,jg_high)      
      write(fldunit,*) ((soln%V(3,i,j),i=ig_low,ig_high),    & !  5
                       NEW_LINE('a'), j=jg_low,jg_high)      
      write(fldunit,*) ((soln%V(4,i,j),i=ig_low,ig_high),    & !  6
                       NEW_LINE('a'), j=jg_low,jg_high)      
      write(fldunit,*) ((soln%mach(i,j),i=ig_low,ig_high),   & !  7
                       NEW_LINE('a'), j=jg_low,jg_high)      
      
    endif
    
  end subroutine output_soln
 ! subroutine output_soln(grid,soln,num_iter)
 !   
 !   use set_inputs, only : counter
 !   
 !   type( grid_t ), intent(in) :: grid
 !   type( soln_t ), intent(in) :: soln
 !   integer,        intent(in) :: num_iter
 !   
 !   integer :: i, j, fldunit
 !   fldunit = 40
 !   
 !   open(fldunit,status='unknown')
 !   write(fldunit,*) 'ZONE'
 !   write(fldunit,*) 'T= "',counter,'"'
 !   write(fldunit,*) 'I=',i_high-i_low+2,' J=',j_high-j_low+2
 !   if(isMMS) then
 !     write(fldunit,*) 'DT = ( DOUBLE, DOUBLE, DOUBLE, DOUBLE, &
 !                            & DOUBLE, DOUBLE, DOUBLE, DOUBLE, &
 !                            & DOUBLE, DOUBLE, DOUBLE)'
 !     write(fldunit,*) 'DATAPACKING = BLOCK'
 !     write(fldunit,*) 'VARLOCATION = ([3-11]=CELLCENTERED)'
 !     
 !     write(fldunit,*) ((grid%x(i,j),i=i_low,i_high+1),    & !  1
 !                      NEW_LINE('a'), j=j_low,j_high+1)    
 !     write(fldunit,*) ((grid%y(i,j),i=i_low,i_high+1),    & !  2
 !                      NEW_LINE('a'), j=j_low,j_high+1)    
 !     write(fldunit,*) ((soln%V(i,j,1),i=i_low,i_high),    & !  3
 !                      NEW_LINE('a'), j=j_low,j_high)      
 !     write(fldunit,*) ((soln%V(i,j,2),i=i_low,i_high),    & !  4
 !                      NEW_LINE('a'), j=j_low,j_high)      
 !     write(fldunit,*) ((soln%V(i,j,3),i=i_low,i_high),    & !  5
 !                      NEW_LINE('a'), j=j_low,j_high)      
 !     write(fldunit,*) ((soln%V(i,j,4),i=i_low,i_high),    & !  6
 !                      NEW_LINE('a'), j=j_low,j_high)      
 !     write(fldunit,*) ((soln%mach(i,j),i=i_low,i_high),   & !  7
 !                      NEW_LINE('a'), j=j_low,j_high)      
 !     write(fldunit,*) ((soln%DE(i,j,1),i=i_low,i_high),   & !  8
 !                      NEW_LINE('a'), j=j_low,j_high)      
 !     write(fldunit,*) ((soln%DE(i,j,2),i=i_low,i_high),   & !  9
 !                      NEW_LINE('a'), j=j_low,j_high)      
 !     write(fldunit,*) ((soln%DE(i,j,3),i=i_low,i_high),   & ! 10
 !                      NEW_LINE('a'), j=j_low,j_high)      
 !     write(fldunit,*) ((soln%DE(i,j,4),i=i_low,i_high),   & ! 11
 !                      NEW_LINE('a'), j=j_low,j_high)      
 !     
 !   else
 !     
 !     write(fldunit,*) 'DT = (DOUBLE, DOUBLE, DOUBLE, DOUBLE,&
 !                          &  DOUBLE, DOUBLE, DOUBLE)'
 !     write(fldunit,*) 'DATAPACKING = BLOCK'
 !     write(fldunit,*) 'VARLOCATION = ([3-7]=CELLCENTERED'
 !     
 !     write(fldunit,*) ((grid%x(i,j),i=i_low,i_high+1),    & !  1
 !                      NEW_LINE('a'), j=j_low,j_high+1)    
 !     write(fldunit,*) ((grid%y(i,j),i=i_low,i_high+1),    & !  2
 !                      NEW_LINE('a'), j=j_low,j_high+1)    
 !     write(fldunit,*) ((soln%V(i,j,1),i=i_low,i_high),    & !  3
 !                      NEW_LINE('a'), j=j_low,j_high)      
 !     write(fldunit,*) ((soln%V(i,j,2),i=i_low,i_high),    & !  4
 !                      NEW_LINE('a'), j=j_low,j_high)      
 !     write(fldunit,*) ((soln%V(i,j,3),i=i_low,i_high),    & !  5
 !                      NEW_LINE('a'), j=j_low,j_high)      
 !     write(fldunit,*) ((soln%V(i,j,4),i=i_low,i_high),    & !  6
 !                      NEW_LINE('a'), j=j_low,j_high)      
 !     write(fldunit,*) ((soln%mach(i,j),i=i_low,i_high),   & !  7
 !                      NEW_LINE('a'), j=j_low,j_high)      
 !     
 !   endif
 !   
 ! end subroutine output_soln
  
  !============================= output_res ==================================80
  !>
  !! Description: 
  !!
  !! Inputs:      grid     : 
  !!              soln     : 
  !!              num_iter : 
  !<
  !===========================================================================80
  subroutine output_res(soln,num_iter)
    
    type( soln_t), intent(in) :: soln
    integer, intent(in) :: num_iter
    integer :: i, resunit
    
    resunit = 30
    
    if(isMMS) then
      write(resunit,*) num_iter,(soln%rnorm(i),i=1,neq),(soln%DEnorm(i),i=1,neq)
    else
      write(resunit,*) num_iter,(soln%rnorm(i),i=1,neq)
    end if
    
  end subroutine output_res
  
end module other_subroutines
