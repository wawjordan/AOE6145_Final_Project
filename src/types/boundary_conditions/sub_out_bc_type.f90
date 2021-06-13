module sub_out_bc_type

  use set_precision, only : prec
  use set_constants, only : pi
  use soln_type,     only : soln_t
  use grid_type,     only : grid_t
  use set_inputs,    only : neq, rho_inf, u_inf, p_inf, alpha

  implicit none

  private

  type, public :: sub_out_bc_t
    integer, dimension(2) :: i1, j1, dir
    integer :: bc_id
    real(prec), allocatable, dimension(:,:) :: rho, uvel, vvel, press
  contains
    procedure, public :: set_bc  => set_bc_sub
    procedure, public :: set_val => set_val_sub
    procedure, public :: enforce => enforce_sub
  end type dir_bc_t

  private :: set_bc_sub, set_val_sub, enforce_sub
contains

  subroutine set_bc_sub(this,grid,ID,i_low,i_high,j_low,j_high)
    class(dir_bc_t), intent(inout) :: this
    type(grid_t),    intent(in)    :: grid
    integer,         intent(in)    :: ID
    integer,         intent(in)    :: i_low, i_high, j_low, j_high
    this%i1 = (/ i_low, i_high /)
    this%j1 = (/ j_low, j_high /)
    this%bc_id = ID
    allocate( this%rho(  i_low:i_high,j_low:j_high),&
              this%uvel( i_low:i_high,j_low:j_high),&
              this%vvel( i_low:i_high,j_low:j_high),&
              this%press(i_low:i_high,j_low:j_high) )
  end subroutine set_bc_sub

  subroutine set_val_sub(this,soln)
    class(dir_bc_t), intent(inout) :: this
    type(soln_t),    intent(in)    :: soln
      this%rho(this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = 
      this%uvel(this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = &
                                        u_inf*cos((pi/180.0_prec)*alpha)
      this%vvel(this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = &
                                        u_inf*sin((pi/180.0_prec)*alpha)
      this%press(this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = p_inf
  end subroutine set_val_sub

  subroutine enforce_sub(this,soln)
    class(dir_bc_t), intent(in)    :: this
    type(soln_t),    intent(inout) :: soln

    soln%V(1,this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = this%rho
    soln%V(2,this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = this%uvel
    soln%V(3,this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = this%vvel
    soln%V(4,this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = this%press
  end subroutine enforce_sub

end module dir_bc_type
