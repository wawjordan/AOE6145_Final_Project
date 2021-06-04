module point_type
  
  use set_precision, only : prec
  
  implicit none
  
  private
  
  public :: point_t, scalar_assignment, array_to_point, point_to_array, &
            point_add, add_scalar, point_subtract, subtract_scalar, &
            multiply_by_scalar, divide_by_scalar, dist, &
            assignment(=), operator(+), operator(-), operator(*), operator(/)
  
  type :: point_t
    real(prec) :: x
    real(prec) :: y
  end type point_t

  interface assignment (=)
    module procedure scalar_assignment
    module procedure array_to_point
    module procedure point_to_array
  end interface
  
  interface operator (+)
    module procedure point_add
    module procedure add_scalar
  end interface
  
  interface operator (-)
    module procedure point_subtract
    module procedure subtract_scalar
  end interface
  
  interface operator (*)
    module procedure multiply_by_scalar
  end interface
  
  interface operator (/)
    module procedure divide_by_scalar
  end interface
  
contains
  pure elemental subroutine scalar_assignment(point1,val)
    type(point_t), intent(out) :: point1
    real(prec), intent(in) :: val
    point1%x = val
    point1%y = val
  end subroutine scalar_assignment
  
  subroutine array_to_point( point1, array )
    type(point_t), intent(inout) :: point1
    real(prec), dimension(2), intent(in) :: array
    point1%x = array(1)
    point1%y = array(2)
  end subroutine array_to_point
  
  subroutine point_to_array( array, point1 )
    real(prec), dimension(2), intent(out) :: array
    type(point_t), intent(in) :: point1
    array(1) = point1%x
    array(2) = point1%y
  end subroutine point_to_array
  
  pure elemental function point_add(point1,point2)
    type(point_t) :: point_add
    type(point_t), intent(in) :: point1, point2
    point_add%x = point1%x + point2%x
    point_add%y = point1%y + point2%y
  end function point_add
  
  function add_scalar(point1,val)
    type(point_t) :: add_scalar
    type(point_t), intent(in) :: point1
    real(prec), intent(in)  :: val
    add_scalar%x = point1%x + val
    add_scalar%y = point1%y + val
  end function add_scalar
  
  pure elemental function point_subtract(point1,point2)
    type(point_t) :: point_subtract
    type(point_t), intent(in) :: point1, point2
    point_subtract%x = point1%x - point2%x
    point_subtract%y = point1%y - point2%y
  end function point_subtract
  
  function subtract_scalar(point1,val)
    type(point_t) :: subtract_scalar
    type(point_t), intent(in) :: point1
    real(prec), intent(in)  :: val
    subtract_scalar%x = point1%x - val
    subtract_scalar%y = point1%y - val
  end function subtract_scalar
  
  function multiply_by_scalar(point1,val)
    type(point_t) :: multiply_by_scalar
    type(point_t), intent(in) :: point1
    real(prec), intent(in)  :: val
    multiply_by_scalar%x = point1%x*val
    multiply_by_scalar%y = point1%y*val
  end function multiply_by_scalar
  
  function divide_by_scalar(point1,val)
    type(point_t) :: divide_by_scalar
    type(point_t), intent(in) :: point1
    real(prec), intent(in)  :: val
    divide_by_scalar%x = point1%x/val
    divide_by_scalar%y = point1%y/val
  end function divide_by_scalar
  
  function dist( point1, point2 )
    real(prec) :: dist
    type(point_t), intent(in) :: point1, point2
    real(prec), dimension(2) :: temp
    
    call point_to_array( temp, (point1-point2) )
    
    dist = sqrt(sum(temp**2))
  end function dist
    
end module point_type
