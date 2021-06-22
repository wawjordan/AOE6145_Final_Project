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
  
  public :: grid_in, grid_out, output_file_headers, output_airfoil_forces, &
           output_soln, output_exact_soln, output_flux, output_res, &
           output_pressure_loss
  
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
  subroutine grid_in(grid,directory,filename)

    type(grid_t), intent(out) :: grid
    character(len=*), intent(in) :: directory
    character(len=*), intent(in) :: filename
    integer :: fstat
    integer :: i, j, k, kmax, nzones
    real(prec) :: zztemp
    
    open(unit=12,file=adjustl(trim(directory))//&
                      adjustl(trim(filename)), &
         status='old',action='read',iostat=fstat)
    if (fstat > 0) then
      write(*,*) 'ERROR: error reading', trim(directory//filename)
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
  
subroutine grid_out(grid,dirname,base_filename)
  use set_inputs, only : cart_grid, grid_name
  type(grid_t), intent(in) :: grid
  character(*), intent(in) :: dirname, base_filename
  integer :: i, j, funit
  
  funit = 50
  open(funit,file= '../results/'//trim(adjustl(dirname))//&
                 & trim(adjustl(base_filename))//'_geometry.dat')
  
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
  subroutine output_file_headers(dirname, base_filename)

    use set_inputs, only : flux_scheme, limiter_scheme, beta_lim, cfl
    use set_inputs, only : grid_name, alpha, m_inf, T_inf, rho_inf, p_inf
    use set_inputs, only : epsM, kappaM, inlet, C_grid, flux_out

    character(*), intent(in) :: dirname, base_filename
    integer :: flxunit, resunit, fldunit, foilunit, inletunit
    
    call execute_command_line ('mkdir -p ../results/'//adjustl(trim(dirname)))
    
    flxunit   = 10
    resunit   = 30
    fldunit   = 40
    foilunit  = 37
    inletunit = 39

    if (flux_out) then
      open(flxunit,file= '../results/'//trim(adjustl(dirname))//&
                     & trim(adjustl(base_filename))//'_flux.dat')
      !open(flxunit,file= '../results/flux_out.dat')
      write(flxunit,*) 'TITLE = "bloody flux"'
      write(flxunit,*) 'variables="x""y" &
        & "F<sub><greek>x</greek>1</sub>"&
        & "F<sub><greek>x</greek>2</sub>"&
        & "F<sub><greek>x</greek>3</sub>"&
        & "F<sub><greek>x</greek>4</sub>"&
        & "F<sub><greek>e</greek>1</sub>"&
        & "F<sub><greek>e</greek>2</sub>"&
        & "F<sub><greek>e</greek>3</sub>"&
        & "F<sub><greek>e</greek>4</sub>"&
        & "<greek>r</greek>""u""v""p"'
        call write_metadata(flxunit)
    end if
    open(resunit,file= '../results/'//trim(adjustl(dirname))//&
                     & trim(adjustl(base_filename))//'_res.dat')
    !open(resunit,file= '../results/example_res.dat')
    if(isMMS) then
      write(resunit,*) 'TITLE = "MMS ('// &
      & trim(adjustl(grid_name))//')"'
      write(resunit,*) 'variables="Iteration"   &
      &   "R<sub>1</sub>""R<sub>2</sub>"   &
      &   "R<sub>3</sub>""R<sub>4</sub>"   &
      &   "DE<sub>1</sub><sup>1<sup>"      &
      &   "DE<sub>2</sub><sup>1<sup>"      &
      &   "DE<sub>3</sub><sup>1<sup>"      &
      &   "DE<sub>4</sub><sup>1<sup>"      &
      &   "DE<sub>1</sub><sup>2<sup>"      &
      &   "DE<sub>2</sub><sup>2<sup>"      &
      &   "DE<sub>3</sub><sup>2<sup>"      &
      &   "DE<sub>4</sub><sup>2<sup>"      &
      &   "DE<sub>1</sub><sup>3<sup>"      &
      &   "DE<sub>2</sub><sup>3<sup>"      &
      &   "DE<sub>3</sub><sup>3<sup>"      &
      &   "DE<sub>4</sub><sup>3<sup>""CFL"'
    elseif(C_grid) then
      write(resunit,*) 'TITLE = "'// &
      & trim(adjustl(grid_name))//'"'
      write(resunit,*) 'variables="Iteration"   &
      &   "R<sub>1</sub>""R<sub>2</sub>"   &
      &   "R<sub>3</sub>""R<sub>4</sub>"   &
      &   "C<sub>L</sub>""C<sub>D</sub>""CFL"'
    elseif(inlet) then
      write(resunit,*) 'TITLE = "'// &
      & trim(adjustl(grid_name))//'"'
      write(resunit,*) 'variables="Iteration"   &
      &   "R<sub>1</sub>""R<sub>2</sub>"   &
      &   "R<sub>3</sub>""R<sub>4</sub>"   &
      &   "Total Pressure Loss""CFL"'
    else
      write(resunit,*) 'TITLE = "'// &
      & trim(adjustl(grid_name))//'"'
      write(resunit,*) 'variables="Iteration"   &
      &   "R<sub>1</sub>""R<sub>2</sub>"   &
      &   "R<sub>3</sub>""R<sub>4</sub>""CFL"'
    end if
    call write_metadata(resunit)

    open(fldunit,file= '../results/'//trim(adjustl(dirname))//&
                     & trim(adjustl(base_filename))//'_field.dat')
    !open(fldunit,file= '../results/example_out.dat')
    write(fldunit,*) 'TITLE = "blahblahblah"'
    if(isMMS) then
      write(fldunit,*) 'variables="x""y""<greek>r</greek>""u""v""p""M"&
      & "DE<sub>1</sub>""DE<sub>2</sub>""DE<sub>3</sub>""DE<sub>4</sub>"'
    else
      write(fldunit,*) 'variables="x""y""<greek>r</greek>""u""v""p""M"'
    end if
    call write_metadata(fldunit)
    
    if (inlet) then
      open(inletunit,file= '../results/'//trim(adjustl(dirname))//&
                       & trim(adjustl(base_filename))//'_losses.dat')
      write(inletunit,*) 'TITLE = "Happy Inlet"'
      write(inletunit,*) 'variables="p<sub>0</sub>""y"'
      call write_metadata(inletunit)
    end if
    
    if (C_grid) then
      open(foilunit,file= '../results/'//trim(adjustl(dirname))//&
                       & trim(adjustl(base_filename))//'_forces.dat')
      write(foilunit,*) 'TITLE = "Happy Airfoil"'
      write(foilunit,*) 'variables="x""p""C<sub>p</sub>"'
      call write_metadata(foilunit)
    end if
    
  end subroutine output_file_headers
  
  subroutine write_metadata(fileunit)
    use set_inputs, only : flux_scheme, limiter_scheme, beta_lim, CFL
    use set_inputs, only : grid_name, alpha, m_inf, u_inf, rho_inf, p_inf, T_inf
    use set_inputs, only : epsM, kappaM, eps_roe, inlet, C_grid, flux_out
    integer, intent(in) :: fileunit
    write(fileunit,"(A24,I0.1,A1)") 'DATASETAUXDATA Flux = "',flux_scheme,'"'
    if (flux_scheme==2) then
      write(fileunit,"(A26,F6.4,A1)") 'DATASETAUXDATA RoeFix = "',eps_roe,'"'
    end if
    if (epsM > zero) then
      write(fileunit,"(A24,I0.1,A1)") 'DATASETAUXDATA Order = "',2,'"'
      write(fileunit,"(A24,F4.1,A1)") 'DATASETAUXDATA Kappa = "',kappaM,'"'
      write(fileunit,"(A26,I0.1,A1)") 'DATASETAUXDATA Limiter = "', &
                                                  & limiter_scheme,'"'
      if (limiter_scheme==4) then
        write(fileunit,"(A23,F3.1,A1)") 'DATASETAUXDATA beta = "',beta_lim,'"'
      end if
    else
      write(fileunit,"(A24,I0.1,A1)") 'DATASETAUXDATA Order = "',1,'"'
    end if
    write(fileunit,"(A25,F7.4,A1)") 'DATASETAUXDATA rhoInf = "',rho_inf,'"'
    write(fileunit,"(A23,F9.4,A1)") 'DATASETAUXDATA vInf = "',u_inf,'"'
    write(fileunit,"(A24,F7.4,A1)") 'DATASETAUXDATA alpha = "',alpha,'"'
    write(fileunit,"(A23,F12.4,A1)") 'DATASETAUXDATA pInf = "',p_inf,'"'
    write(fileunit,"(A23,F9.4,A1)") 'DATASETAUXDATA tInf = "',T_inf,'"'
    write(fileunit,"(A23,F6.4,A1)") 'DATASETAUXDATA mInf = "',M_inf,'"'
    write(fileunit,"(A22,F6.4,A1)") 'DATASETAUXDATA CFL = "',CFL,'"'
  end subroutine write_metadata
   
  subroutine output_pressure_loss(grid,soln,num_iter)
    use other_subroutines, only : inlet_pressure_loss
    type(grid_t), intent(in) :: grid
    type(soln_t), intent(in) :: soln
    integer,      intent(in) :: num_iter
    real(prec), allocatable, dimension(:,:) :: p0
    real(prec) :: loss, p0_inf
    integer :: i, j, inletunit, low, high
    low = j_low
    high = j_high
    allocate( p0(2,low:high) )
    call inlet_pressure_loss(soln,grid,(/low,high,i_high,1/),loss,p0_inf,p0)
    inletunit = 39
    
    write(inletunit,*) 'ZONE'
    write(inletunit,*) 'T= "',num_iter,'"'
    write(inletunit,*) 'I=',high-low+1
    write(inletunit,*) 'AUXDATA Loss= "',loss,'"'
    write(inletunit,*) 'AUXDATA p0inf= "',p0_inf,'"'
    write(inletunit,*) 'DT = ( DOUBLE, DOUBLE )'
    write(inletunit,*) 'DATAPACKING = BLOCK'
    write(inletunit,*) (p0(2,i),i=low,high)
    write(inletunit,*) (p0(1,i),i=low,high)
    deallocate(p0)
  end subroutine output_pressure_loss
  
  subroutine output_airfoil_forces(grid,soln,num_iter)
    use set_inputs, only : neq, imax, index1, index2
    use other_subroutines, only : airfoil_forces
    type(grid_t), intent(in) :: grid
    type(soln_t), intent(in) :: soln
    integer,      intent(in) :: num_iter
    real(prec), allocatable, dimension(:,:) :: Cp
    real(prec) :: Cl, Cd
    integer :: i, j, foilunit, low, high
    low = index2+1
    high = imax-index2-1
    allocate( Cp(3,low:high) )
    call airfoil_forces(soln,grid,(/low,high,j_low,2/),Cp,Cl,Cd)
    foilunit = 37
    
    write(foilunit,*) 'ZONE'
    write(foilunit,*) 'T= "',num_iter,'"'
    write(foilunit,*) 'I=',high-low+1
    write(foilunit,*) 'AUXDATA Cl= "',Cl,'"'
    write(foilunit,*) 'AUXDATA Cd= "',Cd,'"'
    write(foilunit,*) 'DT = ( DOUBLE, DOUBLE, DOUBLE )'
    write(foilunit,*) 'DATAPACKING = BLOCK'
    write(foilunit,*) (Cp(1,i),i=low,high)
    write(foilunit,*) (Cp(2,i),i=low,high)
    write(foilunit,*) (Cp(3,i),i=low,high)
    deallocate(Cp)
  end subroutine output_airfoil_forces
  !======================= output_exact_soln ================================80
  !>
  !! Description:
  !!
  !! Inputs:      grid     :
  !!              soln     :
  !!              num_iter :
  !<
  !==========================================================================80
  subroutine output_exact_soln(grid,soln,dirname,base_filename)

    type( grid_t ), intent(in) :: grid
    type( soln_t ), intent(in) :: soln
    character(*), intent(in) :: dirname, base_filename
    integer :: i, j, fldunit
    fldunit = 20

    open(fldunit,file= '../results/'//trim(adjustl(dirname))//&
                     & trim(adjustl(base_filename))//'_exact.dat')
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
    write(fldunit,*) 'I=',i_high-i_low+2,' J=',j_high-j_low+2
    if(isMMS) then
      write(fldunit,*) 'DT = ( DOUBLE, DOUBLE, DOUBLE, DOUBLE, &
                             & DOUBLE, DOUBLE, DOUBLE, DOUBLE, &
                             & DOUBLE, DOUBLE, DOUBLE)'
      write(fldunit,*) 'DATAPACKING = BLOCK'
      write(fldunit,*) 'VARLOCATION = ([3-11]=CELLCENTERED)'

      write(fldunit,*) ((grid%x(i,j),i=i_low,i_high+1),    & !  1
                       NEW_LINE('a'), j=j_low,j_high+1)
      write(fldunit,*) ((grid%y(i,j),i=i_low,i_high+1),    & !  2
                       NEW_LINE('a'), j=j_low,j_high+1)
      write(fldunit,*) ((soln%V(1,i,j),i=i_low,i_high),    & !  3
                       NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%V(2,i,j),i=i_low,i_high),    & !  4
                       NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%V(3,i,j),i=i_low,i_high),    & !  5
                       NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%V(4,i,j),i=i_low,i_high),    & !  6
                       NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%mach(i,j),i=i_low,i_high),   & !  7
                       NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%DE(1,i,j),i=i_low,i_high),   & !  8
                       NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%DE(2,i,j),i=i_low,i_high),   & !  9
                       NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%DE(3,i,j),i=i_low,i_high),   & ! 10
                       NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%DE(4,i,j),i=i_low,i_high),   & ! 11
                       NEW_LINE('a'), j=j_low,j_high)

    else
    
    write(fldunit,*) 'DT = (DOUBLE, DOUBLE, DOUBLE, DOUBLE,&
                           &  DOUBLE, DOUBLE, DOUBLE)'
      write(fldunit,*) 'DATAPACKING = BLOCK'
      write(fldunit,*) 'VARLOCATION = ([3-7]=CELLCENTERED)'

      write(fldunit,*) ((grid%x(i,j),i=i_low,i_high+1),    & !  1
                       NEW_LINE('a'), j=j_low,j_high+1)
      write(fldunit,*) ((grid%y(i,j),i=i_low,i_high+1),    & !  2
                       NEW_LINE('a'), j=j_low,j_high+1)
      write(fldunit,*) ((soln%V(1,i,j),i=i_low,i_high),    & !  3
                       NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%V(2,i,j),i=i_low,i_high),    & !  4
                       NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%V(3,i,j),i=i_low,i_high),    & !  5
                       NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%V(4,i,j),i=i_low,i_high),    & !  6
                       NEW_LINE('a'), j=j_low,j_high)
      write(fldunit,*) ((soln%mach(i,j),i=i_low,i_high),   & !  7
                       NEW_LINE('a'), j=j_low,j_high)

    endif

  end subroutine output_soln
!  subroutine output_soln(grid,soln,num_iter)
!
!    use set_inputs, only : counter
!
!    type( grid_t ), intent(in) :: grid
!    type( soln_t ), intent(in) :: soln
!    integer,        intent(in) :: num_iter
!
!    integer :: i, j, fldunit
!    fldunit = 40
!
!    open(fldunit,status='unknown')
!    write(fldunit,*) 'ZONE'
!    write(fldunit,*) 'T= "',counter,'"'
!    write(fldunit,*) 'I=',ig_high-ig_low+2,' J=',jg_high-jg_low+2
!    if(isMMS) then
!      write(fldunit,*) 'DT = ( DOUBLE, DOUBLE, DOUBLE, DOUBLE, &
!                             & DOUBLE, DOUBLE, DOUBLE, DOUBLE, &
!                             & DOUBLE, DOUBLE, DOUBLE)'
!      write(fldunit,*) 'DATAPACKING = BLOCK'
!      write(fldunit,*) 'VARLOCATION = ([3-11]=CELLCENTERED)'
!
!      write(fldunit,*) ((grid%x(i,j),i=ig_low,ig_high+1),    & !  1
!                       NEW_LINE('a'), j=jg_low,jg_high+1)
!      write(fldunit,*) ((grid%y(i,j),i=ig_low,ig_high+1),    & !  2
!                       NEW_LINE('a'), j=jg_low,jg_high+1)
!      write(fldunit,*) ((soln%V(1,i,j),i=ig_low,ig_high),    & !  3
!                       NEW_LINE('a'), j=jg_low,jg_high)
!      write(fldunit,*) ((soln%V(2,i,j),i=ig_low,ig_high),    & !  4
!                       NEW_LINE('a'), j=jg_low,jg_high)
!      write(fldunit,*) ((soln%V(3,i,j),i=ig_low,ig_high),    & !  5
!                       NEW_LINE('a'), j=jg_low,jg_high)
!      write(fldunit,*) ((soln%V(4,i,j),i=ig_low,ig_high),    & !  6
!                       NEW_LINE('a'), j=jg_low,jg_high)
!      write(fldunit,*) ((soln%mach(i,j),i=ig_low,ig_high),   & !  7
!                       NEW_LINE('a'), j=jg_low,jg_high)
!      write(fldunit,*) ((soln%DE(1,i,j),i=ig_low,ig_high),   & !  8
!                       NEW_LINE('a'), j=jg_low,jg_high)
!      write(fldunit,*) ((soln%DE(2,i,j),i=ig_low,ig_high),   & !  9
!                       NEW_LINE('a'), j=jg_low,jg_high)
!      write(fldunit,*) ((soln%DE(3,i,j),i=ig_low,ig_high),   & ! 10
!                       NEW_LINE('a'), j=jg_low,jg_high)
!      write(fldunit,*) ((soln%DE(4,i,j),i=ig_low,ig_high),   & ! 11
!                       NEW_LINE('a'), j=jg_low,jg_high)
!
!    else
!    
!    write(fldunit,*) 'DT = (DOUBLE, DOUBLE, DOUBLE, DOUBLE,&
!                           &  DOUBLE, DOUBLE, DOUBLE)'
!      write(fldunit,*) 'DATAPACKING = BLOCK'
!      write(fldunit,*) 'VARLOCATION = ([3-7]=CELLCENTERED)'
!
!      write(fldunit,*) ((grid%x(i,j),i=ig_low,ig_high+1),    & !  1
!                       NEW_LINE('a'), j=jg_low,jg_high+1)
!      write(fldunit,*) ((grid%y(i,j),i=ig_low,ig_high+1),    & !  2
!                       NEW_LINE('a'), j=jg_low,jg_high+1)
!      write(fldunit,*) ((soln%V(1,i,j),i=ig_low,ig_high),    & !  3
!                       NEW_LINE('a'), j=jg_low,jg_high)
!      write(fldunit,*) ((soln%V(2,i,j),i=ig_low,ig_high),    & !  4
!                       NEW_LINE('a'), j=jg_low,jg_high)
!      write(fldunit,*) ((soln%V(3,i,j),i=ig_low,ig_high),    & !  5
!                       NEW_LINE('a'), j=jg_low,jg_high)
!      write(fldunit,*) ((soln%V(4,i,j),i=ig_low,ig_high),    & !  6
!                       NEW_LINE('a'), j=jg_low,jg_high)
!      write(fldunit,*) ((soln%mach(i,j),i=ig_low,ig_high),   & !  7
!                       NEW_LINE('a'), j=jg_low,jg_high)
!
!    endif
!
!  end subroutine output_soln
  !============================= output_flux ================================80
  !>
  !! Description:
  !!
  !! Inputs:      grid     :
  !!              soln     :
  !!              num_iter :
  !<
  !==========================================================================80
  subroutine output_flux(grid,soln,num_iter)

    type( grid_t ), intent(in) :: grid
    type( soln_t ), intent(in) :: soln
    integer,        intent(in) :: num_iter
    
    integer :: i, j, flxunit
    flxunit = 10
    
    open(flxunit,status='unknown')
    write(flxunit,*) 'ZONE'
    write(flxunit,*) 'T= "',num_iter,'"'
    write(flxunit,*) 'I=',i_high-i_low+2,' J=',j_high-j_low+2
    write(flxunit,*) 'DT = ( DOUBLE, DOUBLE, DOUBLE, DOUBLE, &
                           & DOUBLE, DOUBLE, DOUBLE, DOUBLE, &
                           & DOUBLE, DOUBLE, DOUBLE, DOUBLE, &
                           & DOUBLE, DOUBLE )'
    write(flxunit,*) 'DATAPACKING = BLOCK'
    write(flxunit,*) 'VARLOCATION = ([3-14]=CELLCENTERED)'
   


 
    write(flxunit,*) ((grid%x(i,j),i=i_low,i_high+1),      & !  1
                     NEW_LINE('a'), j=j_low,j_high+1)      
    write(flxunit,*) ((grid%y(i,j),i=i_low,i_high+1),      & !  2
                     NEW_LINE('a'), j=j_low,j_high+1)      
    write(flxunit,*) ((soln%Fxi(1,i-1,j),i=i_low,i_high),& !  3
                     NEW_LINE('a'), j=j_low,j_high)        
    write(flxunit,*) ((soln%Fxi(2,i-1,j),i=i_low,i_high),& !  4
                     NEW_LINE('a'), j=j_low,j_high)        
    write(flxunit,*) ((soln%Fxi(3,i-1,j),i=i_low,i_high),& !  5
                     NEW_LINE('a'), j=j_low,j_high)        
    write(flxunit,*) ((soln%Fxi(4,i-1,j),i=i_low,i_high),& !  6
                     NEW_LINE('a'), j=j_low,j_high)        
    write(flxunit,*) ((soln%Feta(1,i,j-1),i=i_low,i_high), & !  7
                     NEW_LINE('a'), j=j_low,j_high)      
    write(flxunit,*) ((soln%Feta(2,i,j-1),i=i_low,i_high), & !  8
                     NEW_LINE('a'), j=j_low,j_high)      
    write(flxunit,*) ((soln%Feta(3,i,j-1),i=i_low,i_high), & !  9
                     NEW_LINE('a'), j=j_low,j_high)      
    write(flxunit,*) ((soln%Feta(4,i,j-1),i=i_low,i_high), & ! 10
                     NEW_LINE('a'), j=j_low,j_high)      
    write(flxunit,*) ((soln%V(1,i,j),i=i_low,i_high),      & ! 11
                     NEW_LINE('a'), j=j_low,j_high)        
    write(flxunit,*) ((soln%V(2,i,j),i=i_low,i_high),      & ! 12
                     NEW_LINE('a'), j=j_low,j_high)        
    write(flxunit,*) ((soln%V(3,i,j),i=i_low,i_high),      & ! 13
                     NEW_LINE('a'), j=j_low,j_high)        
    write(flxunit,*) ((soln%V(4,i,j),i=i_low,i_high),      & ! 14
                     NEW_LINE('a'), j=j_low,j_high)        
    
  end subroutine output_flux
  
  !============================= output_res ==================================80
  !>
  !! Description:
  !!
  !! Inputs:      grid     :
  !!              soln     :
  !!              num_iter :
  !<
  !===========================================================================80
  subroutine output_res(grid,soln,num_iter)
    
    use set_inputs, only : neq, imax, index2, CFL, C_grid, inlet
    use other_subroutines, only : airfoil_forces, inlet_pressure_loss
    use set_inputs, only : CFL
    type( soln_t), intent(in) :: soln
    type( grid_t), intent(in) :: grid
    integer, intent(in) :: num_iter
    integer :: i, j, resunit, low, high
    real(prec), allocatable, dimension(:,:) :: dummy
    real(prec) :: Cl, Cd, loss, p0_inf
    
    resunit = 30
    open(resunit,status='unknown')
    if(isMMS) then
      write(resunit,*) num_iter,(soln%rnorm(j),j=1,neq),&
                      ((soln%DEnorm(j,i),j=1,neq),i=1,3),CFL
    elseif (C_grid) then
      low = index2+1
      high = imax-index2-1
      allocate( dummy(3,low:high) )
      call airfoil_forces(soln,grid,(/low,high,j_low,2/),dummy,Cl,Cd)
      deallocate(dummy)
      write(resunit,*) num_iter,(soln%rnorm(j),j=1,neq), Cl, Cd, CFL
    elseif (inlet) then
      low = j_low
      high = j_high
      allocate( dummy(2,low:high) )
      call inlet_pressure_loss(soln,grid,(/low,high,i_high,1/),&
                                            loss,p0_inf,dummy)
      write(resunit,*) num_iter,(soln%rnorm(j),j=1,neq), loss, CFL
      deallocate(dummy)
    else
      write(resunit,*) num_iter,(soln%rnorm(j),j=1,neq), CFL
    end if

  end subroutine output_res
  
end module file_handling
