module file_handling
  
  use set_precision, only : prec
  use set_constants, only : zero, one, two, half
  use set_inputs, only : neq, n_ghost, isMMS
  use set_inputs, only : imax, i_low, i_high, ig_low, ig_high
  use set_inputs, only : jmax, j_low, j_high, jg_low, jg_high
  use grid_type,  only : grid_t, allocate_grid
  use soln_type,  only : soln_t
  
  implicit none
  
  private
  
  public :: grid_in, grid_out, output_soln, output_exact_soln, output_res, &
            output_file_headers
  
contains
  !=============================== grid_in  =================================80
  !>
  !! Description:
  !!
  !! Inputs:      filename :
  !!
  !! Outputs:     grid :
  !<
  !==========================================================================80
  subroutine grid_in(filename,grid)

    type(grid_t), intent(out) :: grid
    character(len=*), intent(in) :: filename
    integer :: fstat
    integer :: i, j, k, kmax, nzones
    real(prec) :: zztemp

    open(unit=12,file=filename,status='old',action='read',iostat=fstat)
    if (fstat > 0) then
      write(*,*) 'ERROR: error reading', trim(filename)
      stop
    elseif ( fstat == 0 ) then
      read(12,*) nzones
      read(12,*) imax, jmax, kmax
      i_low   = 1
      j_low   = 1
      i_high  = imax-1
      j_high  = jmax-1
      ig_low  = i_low  - n_ghost
      jg_low  = j_low  - n_ghost
      ig_high = i_high + n_ghost
      jg_high = j_high + n_ghost

      grid%i_low   = i_low
      grid%j_low   = j_low
      grid%i_high  = i_high
      grid%j_high  = j_high
      grid%ig_low  = ig_low
      grid%jg_low  = jg_low
      grid%ig_high = ig_high
      grid%jg_high = jg_high

      call allocate_grid(grid)
      if (kmax>2) then
        !! Read in x-coordinate
        do k = 1, kmax
          do j = j_low, j_high+1
            do i = i_low, i_high+1
              read(12,*) grid%x(i,j)
            enddo
          enddo
      enddo
      !! Read in y-coordinate
      do k = 1, kmax
        do j = j_low, j_high+1
          do i = i_low, i_high+1
            read(12,*) grid%y(i,j)
          enddo
        enddo
      enddo

      else
        read(12,*) &
         (((grid%x(i,j),i=i_low,i_high+1),j=j_low,j_high+1),k=1,kmax), &
         (((grid%y(i,j),i=i_low,i_high+1),j=j_low,j_high+1),k=1,kmax), &
         (((zztemp,i=1,imax),j=1,jmax),k=1,kmax)
      end if

    end if

    close(12)
    imax = i_high-i_low+2
    jmax = j_high-j_low+2

  end subroutine grid_in
  
subroutine grid_out(filename,grid)
  
  type(grid_t), intent(in) :: grid
  character(len=*), intent(in) :: filename
  integer :: i, j, funit
   
  funit = 50

  open(unit=funit,file=filename,status='unknown')
  write(funit,*) 'TITLE = "2D Curvilinear Mesh Data"'
  write(funit,*)'variables="x(m)","y(m)","n<sub><greek>x</greek>,x</sub>",&
    & "n<sub><greek>x</greek>,y</sub>", "n<sub><greek>n</greek>,x</sub>",&
    & "n<sub><greek>n</greek>,y</sub>", "A<sub><greek>x</greek></sub>",&
    & "A<sub><greek>n</greek></sub>","V"'

  write(funit,*) 'ZONE'
  write(funit,*) 'T = "Meshy-Mesh"'
  write(funit,*) 'I=',ig_high-ig_low+2,' J=',jg_high-jg_low+2
  write(funit,*) 'DT = (DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE,&
                   & DOUBLE, DOUBLE, DOUBLE, DOUBLE)'
  write(funit,*) 'DATAPACKING = BLOCK'
  write(funit,*) 'VARLOCATION = ([3-9]=CELLCENTERED)'
  
  
  write(funit,*) ((grid%x(i,j),i=ig_low,ig_high+1),      &
                 NEW_LINE('a'),  j=jg_low,jg_high+1)
  write(funit,*) ((grid%y(i,j),i=ig_low,ig_high+1),      &
                 NEW_LINE('a'), j=jg_low,jg_high+1)
  write(funit,*) ((grid%n_xi(i,j,1),i=ig_low,ig_high), &
                 NEW_LINE('a'),  j=jg_low,jg_high)
  write(funit,*) ((grid%n_xi(i,j,2),i=ig_low,ig_high), &
                 NEW_LINE('a'),  j=jg_low,jg_high)
  write(funit,*) ((grid%n_eta(i,j,1),i=ig_low,ig_high),&
                 NEW_LINE('a'),  j=jg_low,jg_high)
  write(funit,*) ((grid%n_eta(i,j,2),i=ig_low,ig_high),&
                 NEW_LINE('a'),  j=jg_low,jg_high)
  write(funit,*) ((grid%A_xi(i,j),i=ig_low,ig_high),     &
                 NEW_LINE('a'),  j=jg_low,jg_high)
  write(funit,*) ((grid%A_eta(i,j),i=ig_low,ig_high),    &
                 NEW_LINE('a'),  j=jg_low,jg_high)
  write(funit,*) ((grid%V(i,j),i=ig_low,ig_high),        &
                 NEW_LINE('a'),  j=jg_low,jg_high)
  
  close(funit)

  end subroutine grid_out
  
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
    integer :: i, j, resunit

    resunit = 30
    i = 1
    open(resunit,status='unknown')
    if(isMMS) then
      write(resunit,*) num_iter,(soln%rnorm(i),j=1,neq),(soln%DEnorm(j,i),j=1,neq)
    else
      write(resunit,*) num_iter,(soln%rnorm(i),j=1,neq)
    end if

  end subroutine output_res
  
end module file_handling
