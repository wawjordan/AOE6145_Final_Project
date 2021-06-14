module bc_type
  
  use set_precision, only : prec
  use fluid_constants, only : R_gas, gamma
  use set_constants, only : pi, zero, one, two, half
  use soln_type, only : soln_t
  use grid_type, only : grid_t
  use set_inputs,    only : neq, pb, rho_inf, u_inf, p_inf, alpha
  
  implicit none
  
  private
  
  type, public :: bc_t
    integer, dimension(2) :: i1, j1, dir
    integer :: bc_id
    real(prec), allocatable, dimension(:,:,:) :: Fxi, Feta
    real(prec), allocatable, dimension(:,:) :: nx, ny, &
                              rho, uvel, vvel, press, mach
  contains
    procedure, public :: set_bc => set_bc_sub
    procedure, public :: set_val => set_val_sub
    procedure, public :: enforce => enforce_sub
  end type bc_t
  
  private :: set_bc_sub, set_val_sub, enforce_sub
  
contains
  
  subroutine set_bc_sub(this,grid,ID,dir,i_low,i_high,j_low,j_high)
    implicit none
    class(bc_t) :: this
    type(grid_t), intent(in) :: grid
    integer, dimension(2), intent(in) :: dir
    integer, intent(in) :: ID
    integer, intent(in) :: i_low, i_high, j_low, j_high
    integer :: i,j
    this%i1 = (/ i_low, i_high /)
    this%j1 = (/ j_low, j_high /)
    this%bc_id = ID
    this%dir = dir
    allocate( this%nx(i_low:i_high,j_low:j_high), &
              this%ny(i_low:i_high,j_low:j_high) )
    if (dir(1)==0) then ! i-index constant
    !write(*,*) 'i-index constant'
      if (i_low<grid%i_low) then ! "left"-side
      !write(*,*) 'left'
        do i = i_low, i_high
          this%nx(i,:) = grid%n_xi(i_high,j_low:j_high,1)
          this%ny(i,:) = grid%n_xi(i_high,j_low:j_high,2)
        end do
      else                       ! "right"-side
      !write(*,*) 'right'
        do i = i_low, i_high
          this%nx(i,:) = grid%n_xi(i_low,j_low:j_high,1)
          this%ny(i,:) = grid%n_xi(i_low,j_low:j_high,2)
        end do
      end if
    else               ! j-index constant
    !write(*,*) 'j-index constant'
      if (j_low<grid%j_low) then ! "bottom"
      !write(*,*) 'bottom'
        do j = j_low, j_high
          this%nx(:,j) = grid%n_eta(i_low:i_high,j_high,1)
          this%ny(:,j) = grid%n_eta(i_low:i_high,j_high,2)
        end do
      else                       ! "top"
      !write(*,*) 'top'
        do j = j_low, j_high
          this%nx(:,j) = grid%n_eta(i_low:i_high,j_low,1)
          this%ny(:,j) = grid%n_eta(i_low:i_high,j_low,2)
        end do
        !write(*,*) 'i, nx,ny'
        !do i = i_low, i_high
        !  write(*,*) i, this%nx(i,j_low), this%ny(i,j_low)
        !end do
      end if
    end if
    
    select case(ID)
    case(1:5)
      allocate( this%rho(  i_low:i_high,j_low:j_high),&
                this%uvel( i_low:i_high,j_low:j_high),&
                this%vvel( i_low:i_high,j_low:j_high),&
                this%press(i_low:i_high,j_low:j_high) )
    case(6)
      allocate( this%Fxi (neq,i_low:i_high,j_low:j_high),&
                this%Feta(neq,i_low:i_high,j_low:j_high) )
    case(7)
      allocate( this%mach(i_low:i_high,j_low:j_high) )
    case default
    end select
    
    
  end subroutine set_bc_sub
  
  subroutine set_val_sub(this,soln)
    use set_inputs, only : epsM
    implicit none
    class(bc_t), intent(inout) :: this
    type(soln_t), intent(in) :: soln
    integer :: i, j, d1, d2
    select case(this%bc_id)
    case(1)  ! MMS dirichlet
      this%rho(this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = &
           soln%Vmms(1,this%i1(1):this%i1(2),this%j1(1):this%j1(2))
      this%uvel(this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = &
           soln%Vmms(2,this%i1(1):this%i1(2),this%j1(1):this%j1(2))
      this%vvel(this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = &
           soln%Vmms(3,this%i1(1):this%i1(2),this%j1(1):this%j1(2))
      this%press(this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = &
           soln%Vmms(4,this%i1(1):this%i1(2),this%j1(1):this%j1(2))
    case(2)  ! Far field dirichlet
      this%rho(this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = rho_inf
      this%uvel(this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = &
                                        u_inf*cos((pi/180.0_prec)*alpha)
      this%vvel(this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = &
                                        u_inf*sin((pi/180.0_prec)*alpha)
      this%press(this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = p_inf
!    case(3)  ! subsonic inflow
!      d1 = this%dir(1)
!      d2 = this%dir(2)
!      do j = this%j1(1), this%j1(2)
!        do i = this%i1(1), this%i1(2)
!          this%rho(i,j)   = soln%V(1,i+d1,j+d2) + epsM*( &
!                            soln%V(1,i+d1,j+d2) - soln%V(1,i+2*d1,j+2*d2) )
!          this%uvel(i,j)  = soln%V(2,i+d1,j+d2) + epsM*( &
!                            soln%V(2,i+d1,j+d2) - soln%V(2,i+2*d1,j+2*d2) )
!          this%vvel(i,j)  = soln%V(3,i+d1,j+d2) + epsM*( &
!                            soln%V(3,i+d1,j+d2) - soln%V(3,i+2*d1,j+2*d2) )
!          this%press(i,j) = pb
!        end do
!      end do
!    case(4)  ! supersonic inflow
!      d1 = this%dir(1)
!      d2 = this%dir(2)
!      do j = this%j1(1), this%j1(2)
!        do i = this%i1(1), this%i1(2)
!          this%rho(i,j)   = soln%V(1,i+d1,j+d2) + epsM*( &
!                            soln%V(1,i+d1,j+d2) - soln%V(1,i+2*d1,j+2*d2) )
!          this%uvel(i,j)  = soln%V(2,i+d1,j+d2) + epsM*( &
!                            soln%V(2,i+d1,j+d2) - soln%V(2,i+2*d1,j+2*d2) )
!          this%vvel(i,j)  = soln%V(3,i+d1,j+d2) + epsM*( &
!                            soln%V(3,i+d1,j+d2) - soln%V(3,i+2*d1,j+2*d2) )
!          this%press(i,j) = pb
!        end do
!      end do
    case(3)  ! subsonic outflow
      d1 = this%dir(1)
      d2 = this%dir(2)
      do j = this%j1(1), this%j1(2)
        do i = this%i1(1), this%i1(2)
          this%rho(i,j)   = soln%V(1,i+d1,j+d2) + epsM*( &
                            soln%V(1,i+d1,j+d2) - soln%V(1,i+2*d1,j+2*d2) )
          this%uvel(i,j)  = soln%V(2,i+d1,j+d2) + epsM*( &
                            soln%V(2,i+d1,j+d2) - soln%V(2,i+2*d1,j+2*d2) )
          this%vvel(i,j)  = soln%V(3,i+d1,j+d2) + epsM*( &
                            soln%V(3,i+d1,j+d2) - soln%V(3,i+2*d1,j+2*d2) )
          this%press(i,j) = p_inf
        end do
      end do
    case(4) ! supersonic outflow
      d1 = this%dir(1)
      d2 = this%dir(2)
      do j = this%j1(1), this%j1(2)
        do i = this%i1(1), this%i1(2)
          this%rho(i,j)   = soln%V(1,i+d1,j+d2) + epsM*( &
                            soln%V(1,i+d1,j+d2) - soln%V(1,i+2*d1,j+2*d2) )
          this%uvel(i,j)  = soln%V(2,i+d1,j+d2) + epsM*( &
                            soln%V(2,i+d1,j+d2) - soln%V(2,i+2*d1,j+2*d2) )
          this%vvel(i,j)  = soln%V(3,i+d1,j+d2) + epsM*( &
                            soln%V(3,i+d1,j+d2) - soln%V(3,i+2*d1,j+2*d2) )
          this%press(i,j) = soln%V(4,i+d1,j+d2) + epsM*( &
                            soln%V(4,i+d1,j+d2) - soln%V(4,i+2*d1,j+2*d2) )
        end do
      end do
    case(5)  ! slip wall w/ ghost cells
      d1 = this%dir(1)
      d2 = this%dir(2)
      
      do j = this%j1(1), this%j1(2)
        do i = this%i1(1), this%i1(2)
          call reflect_vec( soln%V(2,i+d1,j+d2), soln%V(3,i+d1,j+d2), &
                              this%nx(i,j), this%ny(i,j), &
                              this%uvel(i,j), this%vvel(i,j) )
          !write(*,*) soln%V(2,i+d1,j+d2), soln%V(3,i+d1,j+d2) ,this%uvel(i,j),this%vvel(i,j)
          this%press(i,j) = soln%V(4,i+d1,j+d2) + epsM*( &
                            soln%V(4,i+d1,j+d2) - soln%V(4,i+2*d1,j+2*d2) )
          this%rho(i,j)   = this%press(i,j)/( R_gas*soln%temp(i+d1,j+d2) )
        end do
      end do
    case default
    end select
  end subroutine set_val_sub 

  subroutine enforce_sub(this,soln)
    class(bc_t), intent(in)   :: this
    type(soln_t), intent(inout):: soln
    
    soln%V(1,this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = this%rho
    soln%V(2,this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = this%uvel
    soln%V(3,this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = this%vvel
    soln%V(4,this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = this%press
  end subroutine enforce_sub
  
  elemental subroutine reflect_vec(u1,v1,nx,ny,u0,v0)
    real(prec), intent(in) :: u1,v1,nx,ny
    real(prec), intent(inout) :: u0, v0
    real(prec) :: nxny
    nxny = nx/ny
    u0 = ( (one-nxny**2)*u1 - two*nxny*v1 )/(one+nxny**2)
    v0 = -nxny*(u1+u0) - v1
  end subroutine
  
end module bc_type
