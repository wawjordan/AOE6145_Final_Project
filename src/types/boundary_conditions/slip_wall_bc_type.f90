module slip_wall_bc_type

  use set_precision, only : prec
  use set_constants, only : pi, zero
  use soln_type,     only : soln_t
  use grid_type,     only : grid_t
  use set_inputs,    only : neq, epsM

  implicit none

  private

  type, public :: slip_wall_bc_t
    integer, dimension(2) :: i1, j1
    integer :: bc_id, dir, offset
    real(prec), allocatable, dimension(:,:,:) :: F
    real(prec), allocatable, dimension(:,:)   :: nx, ny
  contains
    procedure, public :: set_bc  => set_bc_xi
    procedure, public :: set_val => set_val_xi
    procedure, public :: enforce => enforce_xi
  end type slip_wall_bc_t
  
  type, extends(slip_wall_bc_t), public :: xi_wall_t
  end type xi_wall_t
  
  type, extends(slip_wall_bc_t), public :: eta_wall_t
  contains
    procedure, public :: set_bc  => set_bc_eta
    procedure, public :: set_val => set_val_eta
    procedure, public :: enforce => enforce_eta
  end type eta_wall_t

  private :: set_bc_xi,  set_val_xi,  enforce_xi, &
             set_bc_eta, set_val_eta, enforce_eta
contains

  subroutine set_bc_xi(this,grid,ID,i_low,i_high,j_low,j_high)
    class(slip_wall_bc_t), intent(inout) :: this
    type(grid_t),          intent(in)    :: grid
    integer,               intent(in)    :: ID
    integer,               intent(in)    :: i_low, i_high, j_low, j_high
    this%i1 = (/ i_low, i_high /)
    this%j1 = (/ j_low, j_high /)
    this%bc_id = ID
    allocate( this%F( neq, i_low:i_high, j_low:j_high ), &
                  this%nx( i_low:i_high, j_low:j_high ), &
                  this%ny( i_low:i_high, j_low:j_high ) )
    if (i_low<=grid%i_low) then
      this%dir = 1
      this%offset = 1
      this%nx = grid%n_xi(i_low:i_high,j_low:j_high,1)
      this%ny = grid%n_xi(i_low:i_high,j_low:j_high,2)
    else
      this%dir = -1
      this%offset = 0
      this%nx = grid%n_xi(i_low+1:i_high+1,j_low:j_high,1)
      this%ny = grid%n_xi(i_low+1:i_high+1,j_low:j_high,2)
    end if
  end subroutine set_bc_xi
  
  subroutine set_bc_eta(this,grid,ID,i_low,i_high,j_low,j_high)
    class(eta_wall_t), intent(inout) :: this
    type(grid_t),      intent(in)    :: grid
    integer,           intent(in)    :: ID
    integer,           intent(in)    :: i_low, i_high, j_low, j_high
    this%i1 = (/ i_low, i_high /)
    this%j1 = (/ j_low, j_high /)
    this%bc_id = ID
    allocate( this%F( neq, i_low:i_high, j_low:j_high ), &
                  this%nx( i_low:i_high, j_low:j_high ), &
                  this%ny( i_low:i_high, j_low:j_high ) )
    if (j_low<=grid%j_low) then
      this%dir = 1
      this%offset = 1
      this%nx = grid%n_eta(i_low:i_high,j_low:j_high,1)
      this%ny = grid%n_eta(i_low:i_high,j_low:j_high,2)
    else
      this%dir = -1
      this%offset = 0
      this%nx = grid%n_eta(i_low:i_high,j_low+1:j_high+1,1)
      this%ny = grid%n_eta(i_low:i_high,j_low+1:j_high+1,2)
    end if
  end subroutine set_bc_eta
  
  subroutine set_val_xi(this,soln)
    class(slip_wall_bc_t), intent(inout) :: this
    type(soln_t),          intent(in)    :: soln
    integer :: i, j
    real(prec) :: press
    do j = this%j1(1), this%j1(2)
      do i = this%i1(1), this%i1(2)
        press = soln%V(4,i,j) + epsM*(soln%V(4,i,j) - soln%V(4,i+this%dir,j))
        this%F(:,i,j) = &
              (/ zero, this%nx(i,j)*press, this%ny(i,j)*press, zero/)
      end do
    end do
  end subroutine set_val_xi
  
  subroutine set_val_eta(this,soln)
    class(eta_wall_t), intent(inout) :: this
    type(soln_t),      intent(in)    :: soln
    integer :: i, j
    real(prec) :: press
    do j = this%j1(1), this%j1(2)
      do i = this%i1(1), this%i1(2)
        press = soln%V(4,i,j) + epsM*(soln%V(4,i,j) - soln%V(4,i,j+this%dir))
        this%F(:,i,j) = &
              (/ zero, this%nx(i,j)*press, this%ny(i,j)*press, zero/)
      end do
    end do
  end subroutine set_val_eta
  
  subroutine enforce_xi(this,soln)
    class(slip_wall_bc_t), intent(in)    :: this
    type(soln_t),          intent(inout) :: soln
    soln%Fxi(:,this%i1(1)-1:this%i1(2)-1,this%j1(1):this%j1(2)) = this%F
  end subroutine enforce_xi
  
  subroutine enforce_eta(this,soln)
    class(eta_wall_t), intent(in)    :: this
    type(soln_t),      intent(inout) :: soln
    integer :: i
    i = this%offset
    soln%Feta(:,this%i1(1):this%i1(2),this%j1(1)-i:this%j1(2)-i) = this%F
    !do i = this%i1(1), this%i1(2)
    !write(*,*) this%F(2,i,:), soln%Feta(2,i,this%j1(2)-1)
    !end do
    !stop
  end subroutine enforce_eta

end module slip_wall_bc_type
