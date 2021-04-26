module geometry

  use set_precision,   only : prec
  use set_constants,   only : zero, one, two, half, pi
  use fluid_constants, only : R_gas, gamma
  use set_inputs,      only : area, darea, eps, geometry_file, grid_name
  use soln_type 
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
    
    !integer :: i
    
    !call allocate_soln( soln )
    
    call grid_in(grid_name,grid)
    call ghost_shape(grid)
    call cell_geometry(grid)
    call grid_out(geometry_file,grid)
    
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
    !call deallocate_soln( soln )
    
  end subroutine teardown_geometry

end module geometry
