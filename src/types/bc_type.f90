module bc_type
  
  use set_precision, only : prec
  use set_constants, only : zero, one
  use soln_type, only : soln_t
  use grid_type, only : grid_t
  !use set_inputs,    only : i_low, i_high, ig_low, ig_high
  !use set_inputs,    only : j_low, j_high, jg_low, jg_high
  !use set_inputs,    only : neq, isMMS, isAxi
  
  implicit none
  
  private
  
  public :: bc_t, subsonic_mms_bc
!!  
  type bc_t
    
    integer, allocatable, dimension(:) :: i1, i2
    real(prec), allocatable, dimension(:,:,:) :: Vspec
    
  contains
    
    procedure :: initialize => init_bc
    procedure :: enforce => enforce_primitive_bc
    
  end type bc_t
!!  
  type, extends(bc_t) :: subsonic_mms_bc
    
    real(prec) :: length
    
  contains
    
    procedure :: init_mms_bc
    
  end type subsonic_mms_bc
!!  
  contains
    
  subroutine init_bc(this,i_low,i_high,j_low,j_high,Vspec)
    
    use set_inputs, only : neq
    
    class(bc_t) :: this
    integer :: i_low, i_high, j_low, j_high
    real(prec), dimension(:,:,:) :: Vspec
    integer :: s1, s2, k
    s1 = i_high-i_low
    s2 = j_high-j_low
    
    allocate(this%i1(s1), this%i2(s2), this%Vspec(s1,s2,neq))
    
    do k = i_low,i_high
      this%i1 = k
    end do
    
    do k = j_low,j_high
      this%i2 = k
    end do
    
    this%Vspec = Vspec
    
  end subroutine init_bc
  
  subroutine init_mms_bc(this,grid,i_low,i_high,j_low,j_high,length)
    
    use set_inputs, only : neq
    use grid_type, only : grid_t
    use mms_functions, only : rho_mms, uvel_mms, vvel_mms, press_mms
    
    class(subsonic_mms_bc) :: this
    type(grid_t), intent(in) :: grid
    integer, intent(in) :: i_low, i_high, j_low, j_high
    real(prec), optional :: length
    real(prec), dimension(:,:), allocatable :: x,y
    integer :: s1, s2, k
    s1 = i_high-i_low
    s2 = j_high-j_low
    
    allocate( x(i_low:i_high,j_low:j_high), y(i_low:i_high,j_low:j_high) )
    allocate( this%i1(s1), this%i2(s2), this%Vspec(s1,s2,neq) )
    
    do k = i_low,i_high
      this%i1 = k
    end do
    
    do k = j_low,j_high
      this%i2 = k
    end do
    
    if (present(length)) then
      this%length = length
    else
      this%length = one
    end if
    x = grid%x(this%i1,this%i2)
    y = grid%y(this%i1,this%i2)
    
    this%Vspec(:,:,1) = rho_mms(this%length,x,y)
    this%Vspec(:,:,2) = uvel_mms(this%length,x,y)
    this%Vspec(:,:,3) = vvel_mms(this%length,x,y)
    this%Vspec(:,:,4) = press_mms(this%length,x,y)
    
  end subroutine init_mms_bc
  
  subroutine enforce_primitive_bc( this, soln )
    
    use soln_type, only : soln_t
    
    class(bc_t) :: this
    type(soln_t), intent(inout) :: soln
    
    soln%V(this%i1,this%i2,:) = this%Vspec
    
  end subroutine enforce_primitive_bc
  
end module bc_type
