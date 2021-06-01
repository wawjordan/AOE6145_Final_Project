module bc_type
  
  use set_precision, only : prec
  use set_constants, only : zero, one
  use soln_type, only : soln_t
  use grid_type, only : grid_t
  !use set_inputs,    only : i_low, i_high, ig_low, ig_high
  !use set_inputs,    only : j_low, j_high, jg_low, jg_high
  use set_inputs,    only : neq
  
  implicit none
  
  private
  
  public :: dirichlet_mms_bc_t
!!
  
  type bc_t
    integer, allocatable, dimension(:) :: i1, i2
  contains
    procedure :: init_bc
  end type bc_t
!!
  type, extends(bc_t) :: cell_center_bc_t
    real(prec), allocatable, dimension(:,:,:) :: val
  contains
    procedure :: init => init_cell_center_bc
    procedure :: enforce_prim => enforce_prim_bc
    procedure :: enforce_cons => enforce_cons_bc
  end type cell_center_bc_t
!!  
  type, extends(bc_t) :: flux_bc_t
    real(prec), allocatable, dimension(:,:,:,:) :: F
  contains
    procedure :: init => init_flux_bc
  end type flux_bc_t
!!
  type, extends(cell_center_bc_t) :: dirichlet_mms_bc_t
    real(prec) :: length
  contains
    procedure :: init_dirichlet_mms_bc
  end type dirichlet_mms_bc_t
!!  
!  type, extends(cell_center_bc_t) :: ss_outflow_bc_t
!    real(prec) :: length
!  contains
!    procedure :: init__bc
!  end type ss_outflow_bc_t
!!
  contains
  
  subroutine init_bc(this,i_low,i_high,j_low,j_high)
    
    class(bc_t) :: this
    integer, intent(in) :: i_low, i_high, j_low, j_high
    integer :: s1, s2, k
    s1 = i_high-i_low
    s2 = j_high-j_low
    
    allocate( this%i1(s1), this%i2(s2) )
    
    do k = i_low,i_high
      this%i1 = k
    end do
    
    do k = j_low,j_high
      this%i2 = k
    end do
  
  end subroutine init_bc
  
  subroutine init_cell_center_bc(this,i_low,i_high,j_low,j_high,val)
    
    class(cell_center_bc_t) :: this
    integer, intent(in) :: i_low, i_high, j_low, j_high
    real(prec), dimension(:,:,:), optional :: val
    call init_bc(this,i_low,i_high,j_low,j_high)
    allocate( this%val(i_low:i_high,j_low:j_high,neq) )
    
    if (present(val)) then
      this%val = val
    else
      this%val = zero
    end if
    
  end subroutine init_cell_center_bc
  
  subroutine init_dirichlet_mms_bc(this,grid,i_low,i_high,j_low,j_high,length)
    
    use mms_functions, only : rho_mms, uvel_mms, vvel_mms, press_mms
    class(dirichlet_mms_bc_t) :: this
    type(grid_t), intent(in) :: grid
    integer, intent(in) :: i_low, i_high, j_low, j_high
    real(prec), optional :: length
    
    call init_cell_center_bc(this,i_low,i_high,j_low,j_high)
    
    if (present(length)) then
      this%length = length
    else
      this%length = one
    end if
    
    this%val(:,:,1) = rho_mms(this%length,grid%xc(this%i1,this%i2),&
                                     grid%yc(this%i1,this%i2) )
    this%val(:,:,2) = uvel_mms(this%length,grid%xc(this%i1,this%i2),&
                                     grid%yc(this%i1,this%i2) )
    this%val(:,:,3) = vvel_mms(this%length,grid%xc(this%i1,this%i2),&
                                     grid%yc(this%i1,this%i2) )
    this%val(:,:,4) = press_mms(this%length,grid%xc(this%i1,this%i2),&
                                     grid%yc(this%i1,this%i2) )
    
  end subroutine init_dirichlet_mms_bc
  
  subroutine init_flux_bc(this,i_low,i_high,j_low,j_high,F)
    
    class(flux_bc_t) :: this
    integer, intent(in) :: i_low, i_high, j_low, j_high
    real(prec), dimension(:,:,:,:), intent(in) :: F
    call init_bc(this,i_low,i_high,j_low,j_high)
    allocate( this%F(i_low:i_high,j_low:j_high,neq,2) )
    
    this%F = F
    
  end subroutine init_flux_bc
  
  
  subroutine enforce_prim_bc( this, soln )
    
    use soln_type, only : soln_t
    
    class(cell_center_bc_t) :: this
    type(soln_t), intent(inout) :: soln
    
    soln%V(this%i1,this%i2,:) = this%val
    
  end subroutine enforce_prim_bc
  
  subroutine enforce_cons_bc( this, soln )
    
    use soln_type, only : soln_t
    
    class(cell_center_bc_t) :: this
    type(soln_t), intent(inout) :: soln
    
    soln%U(this%i1,this%i2,:) = this%val
    
  end subroutine enforce_cons_bc
  
  subroutine enforce_flux_bc( this, soln )
    
    use soln_type, only : soln_t
    
    class(flux_bc_t) :: this
    type(soln_t), intent(inout) :: soln
    
    soln%F(this%i1,this%i2,:,:) = this%F
    
  end subroutine enforce_flux_bc
  
end module bc_type
