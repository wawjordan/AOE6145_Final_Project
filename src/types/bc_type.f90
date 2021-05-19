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
    
    integer, allocatable, dimension(:,:) :: ind
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
    
  subroutine init_bc(this,ind,Vspec)
    
    use set_inputs, only : neq
    
    class(bc_t) :: this
    integer, dimension(:,:) :: ind
    real(prec), dimension(:,:,:) :: Vspec
    integer :: s1, s2
    s1 = size(Vspec,1)
    s2 = size(Vspec,2)
    
    allocate(this%ind(s1,s2), this%Vspec(s1,s2,neq))
    
    this%ind = ind
    this%Vspec = Vspec
    
  end subroutine init_bc
  
  subroutine init_mms_bc(this,ind,grid,length)
    
    use set_inputs, only : neq
    use grid_type, only : grid_t
    use mms_functions, only : rho_mms, uvel_mms, vvel_mms, press_mms
    
    class(subsonic_mms_bc) :: this
    integer, dimension(:,:) :: ind
    type(grid_t) :: grid
    real(prec), optional :: length
    real(prec), dimension(size(ind,1),1) :: x,y
    integer :: i, j, s1, s2
    s1 = size(ind,1)
    
    allocate(this%ind(s1,2), this%Vspec(s1,1,neq))
    
    this%ind = ind
    
    if (present(length)) then
      this%length = length
    else
      this%length = one
    end if
    !write(*,*) lbound(ind,1), ubound(ind,1)
    !write(*,*)
    do i = 1,s1
      x(i,1) = grid%x(ind(i,1),ind(i,2))
      y(i,1) = grid%y(ind(i,1),ind(i,2))
      write(*,*) ind(i,1),ind(i,2)
      !write(*,*) 'i=',ind(i,1),' j=',ind(i,2),' x=',x(i,1),' y=',y(i,1)
    end do
    !x = grid%x(ind(:,1),ind(:,2))
    !y = grid%y(ind(:,1),ind(:,2))
    !write(*,*) s1
    !do i = 1,s1
    !  write(*,*) 'x= ',x(i,1),'  y= ',y(i,1)
    !end do
    
    this%Vspec(:,:,1) = rho_mms(this%length,x,y)
    this%Vspec(:,:,2) = uvel_mms(this%length,x,y)
    this%Vspec(:,:,3) = vvel_mms(this%length,x,y)
    this%Vspec(:,:,4) = press_mms(this%length,x,y)
    
  end subroutine init_mms_bc
  
  subroutine enforce_primitive_bc( this, soln )
    
    use soln_type, only : soln_t
    
    class(bc_t) :: this
    type(soln_t), intent(inout) :: soln
    
    soln%V(this%ind(:,1),this%ind(:,2),:) = this%Vspec
    
  end subroutine enforce_primitive_bc
  
end module bc_type
