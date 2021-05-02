module other_subroutines
  
  use set_precision, only : prec
  use set_constants, only : zero, one, two, three, half, fourth
  use set_inputs, only : imax, neq, i_low, i_high, ig_low, ig_high
  use set_inputs, only : jmax, j_low, j_high, jg_low, jg_high
  use set_inputs, only : epsM, kappaM, isMMS, n_ghost
  use fluid_constants, only : gamma
  use variable_conversion
  use limiter_calc, only : limiter_fun, calc_consecutive_variations
  use soln_type, only : soln_t
  use grid_type, only : grid_t
  
  implicit none
  
  private
  
  public :: output_file_headers, output_soln, MUSCL_extrap
  
  contains
  
  !============================= calculate_sources ===========================80
  !>
  !! Description: 
  !!
  !! Inputs:      V  : 
  !!              dA : 
  !!
  !! Outputs:     S  : 
  !<
  !===========================================================================80
  subroutine calculate_sources(P,dA,S)
    
    real(prec), dimension(i_low:i_high),   intent(in) :: dA
    real(prec), dimension(ig_low:ig_high),   intent(in) :: P
    real(prec), dimension(i_low:i_high),   intent(out) :: S
    
    S(i_low:i_high) = P(i_low:i_high)*dA(i_low:i_high)
    
  end subroutine calculate_sources
  
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
  subroutine MUSCL_extrap( V, left, right )
    
    use set_inputs, only : limiter_freeze
    real(prec), dimension(:,:), intent(in)  :: V
    real(prec), dimension(:,:), intent(out) :: left, right
    real(prec), dimension(lbound(V,1):ubound(V,1),neq) :: r_plus, r_minus
    real(prec), dimension(lbound(V,1):ubound(V,1),neq):: psi_plus, psi_minus
    integer :: i,j, low, high
    
    low = lbound(V,1)+n_ghost
    high = ubound(V,1)-n_ghost
    
    if (limiter_freeze) then
      continue
    else
      call calc_consecutive_variations(V,r_plus,r_minus)
      do i = lbound(V,1),ubound(V,1)
        write(*,*) r_plus(i,1), r_minus(i,2)
        write(*,*) V(i,1), V(i,2)
      end do
      call limiter_fun(r_plus,psi_plus)
      !stop
      call limiter_fun(r_minus,psi_minus)
    end if
    
    do i = low-1,high
      j = i - low +2
      left(j,:) = V(i,:) + fourth*epsM*( &
         & (one-kappaM)*psi_plus(i-1,:)*(V(i,:)-V(i-1,:)) + &
         & (one+kappaM)*psi_minus(i,:)*(V(i+1,:)-V(i,:)) )
      right(j,:) = V(i+1,:) - fourth*epsM*( &
         & (one+kappaM)*psi_minus(i+1,:)*(V(i+1,:)-V(i,:)) + &
         & (one-kappaM)*psi_plus(i,:)*(V(i+2,:)-V(i+1,:)) )
    end do

  end subroutine MUSCL_extrap
  
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
  subroutine calc_de( soln, exact_soln, DE, DEnorm, pnorm, cons )
    
    type(soln_t), intent(inout) :: soln
    type(soln_t), intent(in) :: exact_soln
    logical, intent(in) :: cons
    real(prec), dimension(:,:,:), intent(out) :: DE
    real(prec), dimension(1,1:neq), intent(out) :: DEnorm
    integer, intent(in) :: pnorm
    real(prec) :: Linvi, Linvj
    Linvi = one/real(i_high-i_low)
    Linvj = one/real(j_high-j_low)
    if (cons) then
      DE = soln%U(i_low:i_high,j_low:j_high,1:neq) &
       & - soln%Umms(i_low:i_high,j_low:j_high,1:neq)
    else
      DE = soln%V(i_low:i_high,j_low:j_high,1:neq) &
       & - soln%Vmms(i_low:i_high,j_low:j_high,1:neq)
    end if
    
    if (pnorm == 0) then
      DEnorm(1,1) = Linvj*maxval(sum(abs(DE(:,:,1)),2))
      DEnorm(1,2) = Linvj*maxval(sum(abs(DE(:,:,2)),2))
      DEnorm(1,3) = Linvj*maxval(sum(abs(DE(:,:,3)),2))
      DEnorm(1,4) = Linvj*maxval(sum(abs(DE(:,:,4)),2))
    elseif (pnorm == 1) then
      DEnorm(1,1) = Linvi*maxval(sum(abs(DE(:,:,1)),1))
      DEnorm(1,2) = Linvi*maxval(sum(abs(DE(:,:,2)),1))
      DEnorm(1,3) = Linvi*maxval(sum(abs(DE(:,:,3)),1))
      DEnorm(1,4) = Linvi*maxval(sum(abs(DE(:,:,4)),1))
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
    write (ncells_str  , "(A5,I0.3,A6,I0.3)") "imax=",imax,"_jmax=",jmax
    write (CFL_str  , "(A4,F0.4)") "CFL=",cfl
    
    select case(flux_scheme)
    case(1)
      write (flux_str,"(A3)") "V-L"
    case(2)
      write (flux_str,"(A3)") "ROE"
    case default
    end select
    if (nint(epsM).eq.0) then
      write (order_str,"(A9)") "1st-order"
    else
      write (order_str,"(A9)") "2nd-order"
    end if
    
    select case(limiter_scheme)
    case(1)
      write (limiter_str,"(A3)") "_VL"
    case(2)
      write (limiter_str,"(A3)") "_VA"
    case(3)
      write (limiter_str,"(A3)") "_MM"
    case(4)
      write (limiter_str,"(A2,I0.2,A1)") "_B",int(10*beta_lim),"_"
    case default
    end select
    
    write (kappa_str, "(A2,I3.2,A1)") "_K"  , int(10*kappaM),"_"
    write (dirname_hist, *) adjustl(trim(MMS_str)),"/hist/"
    write (dirname_field, *) adjustl(trim(MMS_str)),"/field/"
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
!    if(isMMS) then
!      write(resunit,*) 'TITLE = "MMS ('// &
!      & trim(adjustl(grid_name))//')"'
!      write(resunit,*) 'variables="Iteration"   &
!      &   "R<sub>1</sub>""R<sub>2</sub>"   &
!      &   "R<sub>3</sub>""R<sub>4</sub>"   &
!      &   "DE<sub>1</sub>""DE<sub>2</sub>" &
!      &   "DE<sub>3</sub>""DE<sub>4</sub>"'
!    else
!      write(resunit,*) 'TITLE = "'// &
!      & trim(adjustl(grid_name))//'"'
!      write(resunit,*) 'variables="Iteration"   &
!      &   "R<sub>1</sub>""R<sub>2</sub>"   &
!      &   "R<sub>3</sub>""R<sub>4</sub>"'
!    end if
    
    open(fldunit,file= '../results/example_out.dat')
!    open(fldunit,file= '../results/'//trim(adjustl(dirname_field))//  &
!    &               trim(adjustl(filename))//'_field.dat',status='unknown')
    write(fldunit,*) 'TITLE = "blahblahblah"'
    if(isMMS) then
      write(fldunit,*) 'variables="x""y""<greek>r</greek>""u""v""p""M"&
      & "<greek>r</greek><sub>exact</sub>"&
      & "u<sub>exact</sub>""v<sub>exact</sub>""p<sub>exact</sub>"&
      & "DE<sub>1</sub>""DE<sub>2</sub>""DE<sub>3</sub>""DE<sub>4</sub>"'
    else
      write(fldunit,*) 'variables="x""y""<greek>r</greek>""u""v""p""M"'
    endif
    
  end subroutine output_file_headers
  
  !============================= output_soln =================================80
  !>
  !! Description: 
  !!
  !! Inputs:      grid     : 
  !!              soln     : 
  !!              num_iter : 
  !<
  !===========================================================================80
  subroutine output_soln(grid,soln,num_iter)
    
    use set_inputs, only : counter
    
    type( grid_t ), intent(in) :: grid
    type( soln_t ), intent(in) :: soln
    integer,        intent(in) :: num_iter
    
    integer :: i, j, fldunit
    fldunit = 40
    
    open(fldunit,status='unknown')
    ! Repeat the following each time you want to write out the solution
    !write(40,*) 'zone T="',num_iter,'" '
    write(fldunit,*) 'ZONE'
    write(fldunit,*) 'T= "',counter,'"'
    write(fldunit,*) 'I=',i_high-i_low+2,' J=',j_high-j_low+2
    if(isMMS) then
      write(fldunit,*) 'DT = (DOUBLE, DOUBLE, DOUBLE, DOUBLE,&
                           &  DOUBLE, DOUBLE, DOUBLE, DOUBLE, &
                           &  DOUBLE, DOUBLE, DOUBLE, DOUBLE, &
                           &  DOUBLE, DOUBLE, DOUBLE)'
      write(fldunit,*) 'DATAPACKING = BLOCK'
      write(fldunit,*) 'VARLOCATION = ([3-15]=CELLCENTERED)'
      
      write(fldunit,*) ((grid%x(i,j),i=i_low,i_high+1), NEW_LINE('a'), j=j_low,j_high+1)
      write(fldunit,*) ((grid%y(i,j),i=i_low,i_high+1), NEW_LINE('a'), j=j_low,j_high+1)
      write(fldunit,*) ((soln%V(i,j,1),i=i_low,i_high), NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%V(i,j,2),i=i_low,i_high), NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%V(i,j,3),i=i_low,i_high), NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%V(i,j,4),i=i_low,i_high), NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%mach(i,j),i=i_low,i_high), NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%Vmms(i,j,1),i=i_low,i_high), NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%Vmms(i,j,2),i=i_low,i_high), NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%Vmms(i,j,3),i=i_low,i_high), NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%Vmms(i,j,4),i=i_low,i_high), NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%DE(i,j,1),i=i_low,i_high), NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%DE(i,j,2),i=i_low,i_high), NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%DE(i,j,3),i=i_low,i_high), NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%DE(i,j,4),i=i_low,i_high), NEW_LINE('a'), j=j_low,j_high)
      
    else
      
      write(fldunit,*) 'DT = (DOUBLE, DOUBLE, DOUBLE, DOUBLE,&
                           &  DOUBLE, DOUBLE, DOUBLE)'
      write(fldunit,*) 'DATAPACKING = BLOCK'
      write(fldunit,*) 'VARLOCATION = ([3-7]=CELLCENTERED'
      
      write(fldunit,*) ((grid%x(i,j),i=i_low,i_high+1), j=j_low,j_high+1)
      write(fldunit,*) ((grid%y(i,j),i=i_low,i_high+1), j=j_low,j_high+1)
      write(fldunit,*) ((soln%V(i,j,1),i=i_low,i_high), j=j_low,j_high)
      write(fldunit,*) ((soln%V(i,j,2),i=i_low,i_high), j=j_low,j_high)
      write(fldunit,*) ((soln%V(i,j,3),i=i_low,i_high), j=j_low,j_high)
      write(fldunit,*) ((soln%V(i,j,4),i=i_low,i_high), j=j_low,j_high)
      write(fldunit,*) ((soln%mach(i,j),i=i_low,i_high), j=j_low,j_high)
   
    endif
    
  end subroutine output_soln
  
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
