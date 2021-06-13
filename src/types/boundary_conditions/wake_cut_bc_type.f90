module wake_cut_bc_type

  use set_precision, only : prec
  use set_constants, only : pi, zero
  use soln_type,     only : soln_t
  use grid_type,     only : grid_t

  implicit none

  private

  type, public :: wake_cut_bc_t
    integer :: bc_id, M, N
    integer, allocatable, dimension(:,:) :: i1, i2, j1, j2
  contains
    procedure, public :: set_bc  => set_bc_xi
    procedure, public :: enforce => enforce_xi
  end type wake_cut_bc_t
  
  type, extends(wake_cut_bc_t), public :: xi_cut_t
  end type xi_cut_t
  
  type, extends(wake_cut_bc_t), public :: eta_cut_t
  contains
    procedure, public :: set_bc  => set_bc_eta
    procedure, public :: enforce => enforce_eta
  end type eta_cut_t

  private :: set_bc_xi,  enforce_xi, &
             set_bc_eta, enforce_eta
contains

  subroutine set_bc_xi(this,grid,ID,low,high,n_ghost)
    use set_inputs, only : imax
    class(wake_cut_bc_t), intent(inout) :: this
    type(grid_t),         intent(in)    :: grid
    integer,              intent(in)    :: ID
    integer,              intent(in)    :: low, high, n_ghost
    integer :: i, j, d1
    this%bc_id = ID
    this%N = n_ghost
    d1 = high-low + 1
    this%M = d1
    allocate( this%i1(d1,n_ghost), this%j1(d1,n_ghost), &
              this%i2(d1,n_ghost), this%j2(d1,n_ghost)  )
    do j = 1, n_ghost
      do i = 1, d1
        this%i1(i,j) = low - 1 + i
        this%i2(i,j) = imax - low + 1 - i
        this%j1(i,j) = 1 - j
        this%j2(i,j) = j
      end do
    end do
  end subroutine set_bc_xi
  
  subroutine set_bc_eta(this,grid,ID,low,high,n_ghost)
    use set_inputs, only : jmax
    class(eta_cut_t), intent(inout) :: this
    type(grid_t),     intent(in)    :: grid
    integer,          intent(in)    :: ID
    integer,          intent(in)    :: low, high, n_ghost
    integer :: i, j, d1
    this%bc_id = ID
    this%N = n_ghost
    d1 = high-low + 1
    this%M = d1
    allocate( this%i1(n_ghost,d1), this%j1(n_ghost,d1), &
              this%i2(n_ghost,d1), this%j2(n_ghost,d1)  )
    do j = 1, d1
      do i = 1, n_ghost
        this%i1(i,j) = 1 - i
        this%i2(i,j) = i
        this%j1(i,j) = low - 1 + j
        this%j2(i,j) = jmax - low + 1 - j
      end do
    end do
  end subroutine set_bc_eta
  
  subroutine enforce_xi(this,soln)
    class(wake_cut_bc_t), intent(in)    :: this
    type(soln_t),         intent(inout) :: soln
    integer :: i, j
    
    do j = 1, this%N
      do i = 1, this%M
        soln%V(:,this%i1(i,j),this%j1(i,j)) = &
              soln%V(:,this%i2(i,j),this%j2(i,j))
        soln%V(:,this%i2(i,j),this%j1(i,j)) = &
              soln%V(:,this%i1(i,j),this%j2(i,j))
      end do
    end do
  end subroutine enforce_xi
  
  subroutine enforce_eta(this,soln)
    class(eta_cut_t), intent(in)    :: this
    type(soln_t),         intent(inout) :: soln
    integer :: i, j
    
    do j = 1, this%M
      do i = 1, this%N
        soln%V(:,this%i1(i,j),this%j1(i,j)) = &
              soln%V(:,this%i2(i,j),this%j2(i,j))
        soln%V(:,this%i2(i,j),this%j1(i,j)) = &
              soln%V(:,this%i1(i,j),this%j2(i,j))
      end do
    end do
  end subroutine enforce_eta
  

end module wake_cut_bc_type
