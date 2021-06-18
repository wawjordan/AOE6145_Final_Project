module geometry

  use set_inputs,      only : geometry_file, grid_dir, grid_name
  use soln_type,       only : soln_t, allocate_soln, deallocate_soln 
  use grid_type,       only : grid_t, allocate_grid, deallocate_grid, &
                              ghost_shape, ghost_shape_C, cell_geometry

  implicit none
  
  private
  
  public :: setup_geometry, teardown_geometry
  
  contains
  
  !============================= setup geometry  =============================80
  !>
  !! Description: 
  !!
  !! Inputs:      grid : 
  !!              soln : 
  !!
  !! Outputs:     grid :
  !!              soln :
  !<
  !===========================================================================80
  subroutine setup_geometry( grid, soln )
    
    use set_inputs, only : cart_grid, C_grid, index1, index2, save_grid
    use file_handling, only : grid_in, grid_out
    
    type( soln_t ), intent(inout) :: soln
    type( grid_t ), intent(inout) :: grid
    
    if (cart_grid) then
      call cartesian_grid(grid)
    else
      call grid_in(grid,grid_dir,grid_name)
    end if
    if (C_grid) then
      call ghost_shape_C(grid,index1,index2)
    else
      call ghost_shape(grid)
    end if
    call cell_geometry(grid)
    !call grid_out(grid,geometry_file)
    if (save_grid) then
      call grid_out(grid)
    end if
    call allocate_soln( soln )
    
  end subroutine setup_geometry
  
  !========================== teardown geometry  =============================80
  !>
  !! Description: 
  !!
  !! Inputs:      grid : 
  !!              soln : 
  !!
  !! Outputs:     grid :
  !!              soln :
  !<
  !===========================================================================80
  subroutine teardown_geometry( grid, soln )
    
    type( soln_t ) :: soln
    type( grid_t ) :: grid
    
    call deallocate_grid( grid )
    call deallocate_soln( soln )
    
  end subroutine teardown_geometry
  
  subroutine cartesian_grid( grid )
    
    use set_precision, only : prec
    use set_constants, only : one
    use set_inputs, only : imax, jmax, xmin, xmax, ymin, ymax, n_ghost
    use set_inputs, only : i_low, i_high, j_low, j_high
    use set_inputs, only : ig_low, ig_high, jg_low, jg_high
    
    type( grid_t ), intent(inout) :: grid
    integer :: i, j
    
    call allocate_grid(grid)
    grid%i_low   = i_low
    grid%j_low   = j_low
    grid%i_high  = i_high
    grid%j_high  = j_high
    grid%ig_low  = ig_low
    grid%jg_low  = jg_low
    grid%ig_high = ig_high
    grid%jg_high = jg_high
    
    do j = j_low,j_high+1
      do i = i_low,i_high+1
        grid%x(i,j) = xmin + real(i-1,prec)/real(imax-1,prec)*(xmax-xmin)
        grid%y(i,j) = ymin + real(j-1,prec)/real(jmax-1,prec)*(ymax-ymin)
      end do
    end do
    
  end subroutine cartesian_grid
  
end module geometry
