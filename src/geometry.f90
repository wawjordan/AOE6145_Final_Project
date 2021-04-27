module geometry

  use set_inputs,      only : geometry_file, grid_name
  use soln_type,       only : soln_t, allocate_soln, deallocate_soln 
  use grid_type,       only : grid_t, allocate_grid, deallocate_grid, &
                              ghost_shape, cell_geometry

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
    
    use file_handling, only : grid_in, grid_out
    
    type( soln_t ), intent(inout) :: soln
    type( grid_t ), intent(inout) :: grid
    
    call grid_in(grid_name,grid)
    call ghost_shape(grid)
    call cell_geometry(grid)
    call grid_out(geometry_file,grid)
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

end module geometry
