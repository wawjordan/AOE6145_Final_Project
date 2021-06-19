module namelist
  
  use fluid_constants, only : R_gas, gamma
  use set_inputs, only : grid_dir, grid_name, cart_grid, C_grid 
  use set_inputs, only : index1, index2, imax, jmax, n_ghost
  use set_inputs, only : xmin, xmax, ymin, ymax, isAxi, Lmms
  use set_inputs, only : isMMS, u_inf, alpha, p_inf, T_inf, M_inf, u0, v0
  use set_inputs, only : CFL, max_iter, eps, tol, locTime, flux_out
  use set_inputs, only : flux_scheme, limiter_scheme, eps_roe, beta_lim
  use set_inputs, only : geometry_file, soln_save, res_save,res_out, cons
  use set_inputs, only : epsM, kappaM, limiter_freeze
  use set_inputs, only : num_BCs, bounds
  
  implicit none
  private
  
  public :: read_namelist
  
contains
  
  subroutine read_namelist(file_path)
    
    !use set_precision, only : prec
    !use set_constants, only : zero, one
    
    character(len=*), intent(in) :: file_path
    integer :: fstat
    integer :: funit
    logical :: fexist
    logical :: fopen = .false.
    !character(len=20) :: file_path = 'q1d.nml'
    namelist /grid/ grid_dir, grid_name, cart_grid, C_grid, &
                    index1, index2, imax, jmax, n_ghost
    namelist /geometry/ xmin, xmax, ymin, ymax, isAxi, Lmms
    namelist /constants/ R_gas, gamma, num_BCs
    namelist /initial/ isMMS, p_inf, u0, v0, u_inf, alpha, T_inf, M_inf
    namelist /numerical/ CFL, max_iter, eps, tol, locTime
    namelist /boundary/ bounds
    namelist /flux/ flux_scheme, limiter_scheme, eps_roe, beta_lim
    namelist /output/ geometry_file, soln_save, res_save, res_out, cons, &
                      flux_out
    namelist /reconstruction/ epsM, kappaM, limiter_freeze
    
    inquire( file=file_path,exist=fexist )
    if ( .not. fexist ) then
      write(*,*) 'ERROR: ',trim(file_path),' does not exist.'
      return
    end if
    
    open(file=file_path,status='old',action='read',iostat=fstat,newunit=funit, &
         access='sequential')
    fopen = .true.
    ! grid
    rewind(funit)
    read(funit,nml=grid,iostat=fstat)
    if (fstat > 0) then
      write(*,*) 'ERROR: error in namelist "grid".'
      stop
    end if
    
    ! geometry
    rewind(funit)
    read(funit,nml=geometry,iostat=fstat)
    if (fstat>0) then
      write(*,*) 'ERROR: error in namelist "geometry".'
      stop
    end if
    
    ! constants 
    rewind(funit)
    read(funit,nml=constants,iostat=fstat)
    if (fstat>0) then
      write(*,*) 'ERROR: error in namelist "constants".'
      stop
    end if
    allocate(bounds(6,num_BCs))
    
    ! initial
    rewind(funit)
    read(funit,nml=initial,iostat=fstat)
    if (fstat>0) then
      write(*,*) 'ERROR: error in namelist "initial".'
      stop
    end if
    
    
    ! numerical
    rewind(funit)
    read(funit,nml=numerical,iostat=fstat)
    if (fstat>0) then
      write(*,*) 'ERROR: error in namelist "numerical".'
      stop
    end if
    
    ! boundary
    rewind(funit)
    read(funit,nml=boundary,iostat=fstat)
    if (fstat > 0) then
      write(*,*) 'ERROR: error in namelist "boundary".'
      stop
    end if
    
    ! flux
    rewind(funit)
    read(funit,nml=flux,iostat=fstat)
    if (fstat>0) then
      write(*,*) 'ERROR: error in namelist "flux".'
      stop
    end if
    
    ! output
    rewind(funit)
    read(funit,nml=output,iostat=fstat)
    if (fstat>0) then
      write(*,*) 'ERROR: error in namelist "output".'
      stop
    end if
    
    ! reconstruction
    rewind(funit)
    read(funit,nml=reconstruction,iostat=fstat)
    if (fstat>0) then
      write(*,*) 'ERROR: error in namelist "reconstruction".'
      stop
    end if
    
    close(funit)
    
  end subroutine read_namelist
      


end module namelist
