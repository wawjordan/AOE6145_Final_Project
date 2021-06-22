module bc_type
  
  use set_precision, only : prec
  use fluid_constants, only : R_gas, gamma
  use set_constants, only : pi, zero, one, two, half
  use soln_type, only : soln_t
  use grid_type, only : grid_t
  use variable_conversion, only : speed_of_sound
  use set_inputs,    only : neq, pb, rho_inf, u_inf, p_inf, alpha
  
  implicit none
  
  private
  public :: reflect_vec
  
  type, public :: bc_t
    integer :: bc_id, i_offset, j_offset
    !integer, dimension(2) :: i1, j1
    integer, dimension(3) :: loop1, loop2
    integer, dimension(2) :: i0, j0
    integer, allocatable, dimension(:,:) :: i1, j1, i2, j2
!    real(prec), allocatable, dimension(:,:,:) :: Fxi, Feta
    real(prec), allocatable, dimension(:,:) :: nx, ny
!                              rho, uvel, vvel, press, mach
  contains
    procedure, public :: set_bc  => set_bc_sub
!    procedure, public :: set_val => set_val_sub
    procedure, public :: enforce => enforce_sub
  end type bc_t
  
  private :: set_bc_sub, enforce_sub !, set_val_sub
  
contains
  
  subroutine set_bc_sub(this,grid,ID,i_low,i_high,j_low,j_high,ij)
    use set_inputs, only : imax, jmax
    implicit none
    class(bc_t) :: this
    type(grid_t), intent(in) :: grid
    integer, intent(in) :: ID
    integer, intent(in) :: i_low, i_high, j_low, j_high
    integer, intent(in) :: ij
    integer :: i,j
    !this%i1 = (/ i_low, i_high /)
    !this%j1 = (/ j_low, j_high /)
    this%bc_id = ID
    select case(ID)
    case(1:5)
      allocate( this%nx(i_low:i_high,j_low:j_high), &
                this%ny(i_low:i_high,j_low:j_high)  )
      this%i0 = (/ i_low, i_high /)
      this%j0 = (/ j_low, j_high /)
      
      select case(ij)
      case(1) ! "left" (xi / i == constant)
        do i = i_low, i_high
          this%nx(i,:) = grid%n_xi(i_high,j_low:j_high,1)
          this%ny(i,:) = grid%n_xi(i_high,j_low:j_high,2)
        end do
        this%i_offset = -1
        this%j_offset = 0
        this%loop1 = (/j_low,j_high,1/)
        this%loop2 = (/i_high,i_low,-1/)
      case(2) ! "right" (xi / i == constant)
        do i = i_low, i_high
          this%nx(i,:) = grid%n_xi(i_low,j_low:j_high,1)
          this%ny(i,:) = grid%n_xi(i_low,j_low:j_high,2)
        end do
        this%i_offset = 1
        this%j_offset = 0
        this%loop1 = (/j_low,j_high,1/)
        this%loop2 = (/i_low,i_high,1/)
      case(3) ! "bottom" (eta / j == constant)
        do j = j_low, j_high
         this%nx(:,j) = grid%n_eta(i_low:i_high,j_high+1,1)
         this%ny(:,j) = grid%n_eta(i_low:i_high,j_high+1,2)
         !this%nx(:,j) = grid%n_eta(i_low:i_high,j_high,1)
         !this%ny(:,j) = grid%n_eta(i_low:i_high,j_high,2)
        end do
        this%i_offset = 0
        this%j_offset = -1
        this%loop1 = (/j_high,j_low,-1/)
        this%loop2 = (/i_low,i_high,1/)
      case(4) ! "top" (eta / j == constant)
        do j = j_low, j_high
          this%nx(:,j) = grid%n_eta(i_low:i_high,j_low,1)
          this%ny(:,j) = grid%n_eta(i_low:i_high,j_low,2)
        end do
        this%i_offset = 0
        this%j_offset = 1
        this%loop1 = (/j_low,j_high,1/)
        this%loop2 = (/i_low,i_high,1/)
      case default
      end select
    case(6)
      select case(ij)
      case(3:4)
        this%i_offset = i_high - i_low + 1
        this%j_offset = j_high - j_low + 1
        this%loop1 = (/1,this%j_offset,1/)
        this%loop2 = (/1,this%i_offset,1/)
        
        allocate( this%i1(this%i_offset,this%j_offset), &
                  this%j1(this%i_offset,this%j_offset), &
                  this%i2(this%i_offset,this%j_offset), &
                  this%j2(this%i_offset,this%j_offset)  )
        do j = this%loop1(1),this%loop1(2),this%loop1(3)
          do i = this%loop2(1),this%loop2(2),this%loop2(3)
            this%i1(i,j) = j_low - 1 + i
            this%i2(i,j) = imax - j_low + 1 - i
            this%j1(i,j) = 1 - j
            this%j2(i,j) = j
            !write(*,*) this%i1(i,j),this%i2(i,j),this%j1(i,j),this%j2(i,j)
          end do
        end do
      case(1:2)
        this%i_offset = i_high - i_low + 1
        this%j_offset = j_high - j_low + 1
        this%loop1 = (/1,this%j_offset,1/)
        this%loop2 = (/1,this%i_offset,1/)
       
        allocate( this%i1(this%i_offset,this%j_offset), &
                  this%j1(this%i_offset,this%j_offset), &
                  this%i2(this%i_offset,this%j_offset), &
                  this%j2(this%i_offset,this%j_offset)  )
        do j = this%loop1(1),this%loop1(2),this%loop1(3)
          do i = this%loop2(1),this%loop2(2),this%loop2(3)
            this%i1(i,j) = 1 - i
            this%i2(i,j) = i
            this%j1(i,j) = i_low - 1 + j
            this%j2(i,j) = jmax - i_low + 1 - j
            !write(*,*) this%i1(i,j),this%i2(i,j),this%j1(i,j),this%j2(i,j)
          end do
        end do
      case default
      end select
    end select
!    select case(ID)
!    case(1:5)
!      allocate( this%rho(  i_low:i_high,j_low:j_high),&
!                this%uvel( i_low:i_high,j_low:j_high),&
!                this%vvel( i_low:i_high,j_low:j_high),&
!                this%press(i_low:i_high,j_low:j_high) )
!    case(6)
!      allocate( this%Fxi (neq,i_low:i_high,j_low:j_high),&
!                this%Feta(neq,i_low:i_high,j_low:j_high) )
!    case(7)
!      allocate( this%mach(i_low:i_high,j_low:j_high) )
!    case default
!    end select
    
    
  end subroutine set_bc_sub
  
  subroutine enforce_sub(this,soln)
    use set_inputs, only : epsM
    implicit none
    class(bc_t), intent(inout) :: this
    type(soln_t), intent(inout) :: soln
    real(prec) :: asnd,mach,vmag,vx,vy
    integer :: i, j, d1, d2, i2, j2
    select case(this%bc_id)
    case(1)  ! MMS dirichlet
      soln%V(1,this%i0(1):this%i0(2),this%j0(1):this%j0(2)) = &
           soln%Vmms(1,this%i0(1):this%i0(2),this%j0(1):this%j0(2))
      soln%V(2,this%i0(1):this%i0(2),this%j0(1):this%j0(2)) = &
           soln%Vmms(2,this%i0(1):this%i0(2),this%j0(1):this%j0(2))
      soln%V(3,this%i0(1):this%i0(2),this%j0(1):this%j0(2)) = &
           soln%Vmms(3,this%i0(1):this%i0(2),this%j0(1):this%j0(2))
      soln%V(4,this%i0(1):this%i0(2),this%j0(1):this%j0(2)) = &
           soln%Vmms(4,this%i0(1):this%i0(2),this%j0(1):this%j0(2))
      !this%rho(this%i0(1):this%i0(2),this%j0(1):this%j0(2)) = &
      !     soln%Vmms(1,this%i0(1):this%i0(2),this%j0(1):this%j0(2))
      !this%uvel(this%i0(1):this%i0(2),this%j0(1):this%j0(2)) = &
      !     soln%Vmms(2,this%i0(1):this%i0(2),this%j0(1):this%j0(2))
      !this%vvel(this%i0(1):this%i0(2),this%j0(1):this%j0(2)) = &
      !     soln%Vmms(3,this%i0(1):this%i0(2),this%j0(1):this%j0(2))
      !this%press(this%i0(1):this%i0(2),this%j0(1):this%j0(2)) = &
      !     soln%Vmms(4,this%i0(1):this%i0(2),this%j0(1):this%j0(2))
    case(2)  ! Far field dirichlet
      soln%V(1,this%i0(1):this%i0(2),this%j0(1):this%j0(2)) = rho_inf
      soln%V(2,this%i0(1):this%i0(2),this%j0(1):this%j0(2)) = &
                                        u_inf*cos((pi/180.0_prec)*alpha)
      soln%V(3,this%i0(1):this%i0(2),this%j0(1):this%j0(2)) = &
                                        u_inf*sin((pi/180.0_prec)*alpha)
      soln%V(4,this%i0(1):this%i0(2),this%j0(1):this%j0(2)) = p_inf
      !this%rho(this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = rho_inf
      !this%uvel(this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = &
      !                                  u_inf*cos((pi/180.0_prec)*alpha)
      !this%vvel(this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = &
      !                                  u_inf*sin((pi/180.0_prec)*alpha)
      !this%press(this%i1(1):this%i1(2),this%j1(1):this%j1(2)) = p_inf
    case(3)  ! subsonic outflow
      d1 = this%i_offset
      d2 = this%j_offset
      do j = this%loop1(1),this%loop1(2),this%loop1(3)
        do i = this%loop2(1),this%loop2(2),this%loop2(3)
          soln%V(1,i,j)   = soln%V(1,i-d1,j-d2) + epsM*( &
                            soln%V(1,i-d1,j-d2) - soln%V(1,i-2*d1,j-2*d2) )
          soln%V(2,i,j)  = soln%V(2,i-d1,j-d2) + epsM*( &
                            soln%V(2,i-d1,j-d2) - soln%V(2,i-2*d1,j-2*d2) )
          soln%V(3,i,j)  = soln%V(3,i-d1,j-d2) + epsM*( &
                            soln%V(3,i-d1,j-d2) - soln%V(3,i-2*d1,j-2*d2) )
          soln%V(4,i,j) = p_inf
          !this%rho(i,j)   = soln%V(1,i-d1,j-d2) + epsM*( &
          !                  soln%V(1,i-d1,j-d2) - soln%V(1,i-2*d1,j-2*d2) )
          !this%uvel(i,j)  = soln%V(2,i-d1,j-d2) + epsM*( &
          !                  soln%V(2,i-d1,j-d2) - soln%V(2,i-2*d1,j-2*d2) )
          !this%vvel(i,j)  = soln%V(3,i-d1,j-d2) + epsM*( &
          !                  soln%V(3,i-d1,j-d2) - soln%V(3,i-2*d1,j-2*d2) )
          !this%press(i,j) = p_inf
        end do
      end do
    case(4) ! supersonic outflow
      d1 = this%i_offset
      d2 = this%j_offset
      do j = this%loop1(1),this%loop1(2),this%loop1(3)
        do i = this%loop2(1),this%loop2(2),this%loop2(3)
          soln%V(1,i,j)   = soln%V(1,i-d1,j-d2) + epsM*( &
                            soln%V(1,i-d1,j-d2) - soln%V(1,i-2*d1,j-2*d2) )
          soln%V(2,i,j)  = soln%V(2,i-d1,j-d2) + epsM*( &
                            soln%V(2,i-d1,j-d2) - soln%V(2,i-2*d1,j-2*d2) )
          soln%V(3,i,j)  = soln%V(3,i-d1,j-d2) + epsM*( &
                            soln%V(3,i-d1,j-d2) - soln%V(3,i-2*d1,j-2*d2) )
          soln%V(4,i,j) = soln%V(4,i-d1,j-d2) + epsM*( &
                            soln%V(4,i-d1,j-d2) - soln%V(4,i-2*d1,j-2*d2) )
          !this%rho(i,j)   = soln%V(1,i-d1,j-d2) + epsM*( &
          !                  soln%V(1,i-d1,j-d2) - soln%V(1,i-2*d1,j-2*d2) )
          !this%uvel(i,j)  = soln%V(2,i-d1,j-d2) + epsM*( &
          !                  soln%V(2,i-d1,j-d2) - soln%V(2,i-2*d1,j-2*d2) )
          !this%vvel(i,j)  = soln%V(3,i-d1,j-d2) + epsM*( &
          !                  soln%V(3,i-d1,j-d2) - soln%V(3,i-2*d1,j-2*d2) )
          !this%press(i,j) = soln%V(4,i-d1,j-d2) + epsM*( &
          !                  soln%V(4,i-d1,j-d2) - soln%V(4,i-2*d1,j-2*d2) )
        end do
      end do
    case(5)  ! slip wall w/ ghost cells
      d1 = this%i_offset
      d2 = this%j_offset 
     !write(*,*) d1, d2
     !write(*,*) this%i1(1), this%i1(2) 
     !write(*,*) this%j1(1), this%j1(2)
     ! do j = this%j1(1), this%j1(2)
     !   do i = this%i1(2), this%i1(1),-1
     !     write(*,*) this%nx(i,j), this%ny(i,j)
     !   end do
     ! end do
     !stop
      i2 = 0
      j2 = 0
      do j = this%loop1(1),this%loop1(2),this%loop1(3)
        j2 = j2+1
        do i = this%loop2(1),this%loop2(2),this%loop2(3)
          i2 = i2+1
          call reflect_vec( soln%V(2,i-d1*(2*i2-1),j-d2*(2*j2-1)), &
                            soln%V(3,i-d1*(2*i2-1),j-d2*(2*j2-1)), &
                              this%nx(i,j), this%ny(i,j), &
                              soln%V(2,i,j), soln%V(3,i,j) )
          soln%V(4,i,j) = soln%V(4,i-d1,j-d2) + epsM*( &
                            soln%V(4,i-d1,j-d2) - soln%V(4,i-2*d1,j-2*d2) )
          soln%V(1,i,j) = soln%V(4,i,j)/&
              ( R_gas*soln%temp(i-d1,j-d2) )
          !write(*,*) i,j,'-->',i-d1*(2*i2-1),j-d2*(2*j2-1)
          !call reflect_vec( soln%V(2,i-d1,j-d2), soln%V(3,i-d1,j-d2), &
          !                    this%nx(i,j), this%ny(i,j), &
          !                    soln%V(2,i,j), soln%V(3,i,j) )
          !soln%V(4,i,j) = soln%V(4,i-d1,j-d2) + epsM*( &
          !                  soln%V(4,i-d1,j-d2) - soln%V(4,i-2*d1,j-2*d2) )
          !soln%V(1,i,j) = soln%V(4,i,j)/( R_gas*soln%temp(i-d1,j-d2) )
          !call reflect_vec( soln%V(2,i-d1,j-d2), soln%V(3,i-d1,j-d2), &
          !                    this%nx(i,j), this%ny(i,j), &
          !                    this%uvel(i,j), this%vvel(i,j) )
          !this%press(i,j) = soln%V(4,i-d1,j-d2) + zero*epsM*( &
          !                  soln%V(4,i-d1,j-d2) - soln%V(4,i-2*d1,j-2*d2) )
          !this%rho(i,j)   = this%press(i,j)/( R_gas*soln%temp(i-d1,j-d2) )
        end do
      end do
    case(6)  ! wake cut
      !write(*,*)
      !write(*,*) this%loop1(1),this%loop1(2) 
      !write(*,*) this%loop2(1),this%loop2(2)
      !stop
      do j = this%loop1(1),this%loop1(2),this%loop1(3)
        do i = this%loop2(1),this%loop2(2),this%loop2(3)
          soln%V(:,this%i1(i,j),this%j1(i,j)) = &
                soln%V(:,this%i2(i,j),this%j2(i,j))
          soln%V(:,this%i2(i,j),this%j1(i,j)) = &
                soln%V(:,this%i1(i,j),this%j2(i,j))
        end do
      end do  
    case default
    end select
  end subroutine enforce_sub 

!  subroutine enforce_sub(this,soln)
!    class(bc_t), intent(in)   :: this
!    type(soln_t), intent(inout):: soln
!    
!    soln%V(1,this%i0(1):this%i0(2),this%j0(1):this%j0(2)) = this%rho
!    soln%V(2,this%i0(1):this%i0(2),this%j0(1):this%j0(2)) = this%uvel
!    soln%V(3,this%i0(1):this%i0(2),this%j0(1):this%j0(2)) = this%vvel
!    soln%V(4,this%i0(1):this%i0(2),this%j0(1):this%j0(2)) = this%press
!  end subroutine enforce_sub
  
  elemental subroutine reflect_vec(u1,v1,nx,ny,u0,v0)
    real(prec), intent(in) :: u1,v1,nx,ny
    real(prec), intent(inout) :: u0, v0
    real(prec) :: nxny
    nxny = nx/ny
    u0 = ( (one-nxny**2)*u1 - two*nxny*v1 )/(one+nxny**2)
    v0 = -nxny*(u1+u0) - v1
  end subroutine
  
end module bc_type
