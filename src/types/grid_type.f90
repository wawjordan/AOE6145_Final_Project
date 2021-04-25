module grid_type
  
  use set_precision, only : prec
  use set_constants, only : zero, one, two, half
  
  implicit none
  
  type grid_t

    integer :: imax, jmax, n_ghost
    integer :: i_low, i_high, ig_low, ig_high
    integer :: j_low, j_high, jg_low, jg_high
    
    
    real(prec), allocatable, dimension(:,:) :: x, y
    real(prec), allocatable, dimension(:,:) :: A_xi, A_eta
    real(prec), allocatable, dimension(:,:) :: n_xi_x, n_xi_y
    real(prec), allocatable, dimension(:,:) :: n_eta_x, n_eta_y
    real(prec), allocatable, dimension(:,:) :: V
    
  end type grid_t

  contains
  
  !============================= allocate_grid ===============================80
  !>
  !! Description: 
  !!
  !! Inputs:      grid :  
  !!
  !! Outputs:     grid : 
  !<
  !===========================================================================80
  subroutine allocate_grid( grid )
    
    use set_inputs   , only : grid_name, n_ghost
    
    type( grid_t ), intent( inout ) :: grid
    integer :: i, j, k, imax, jmax, kmax, nzones
    integer :: i_low, i_high, ig_low, ig_high
    integer :: j_low, j_high, jg_low, jg_high
    real(prec) :: zztemp
    
    open(unit=12,file=grid_name)
    read(12,*) nzones
    read(12,*) imax, jmax, kmax
    i_low   = 1
    j_low   = 1
    i_high  = imax-1
    j_high  = jmax-1
    ig_low  = i_low  - n_ghost
    jg_low  = j_low  - n_ghost
    ig_high = i_high + n_ghost
    jg_high = j_high + n_ghost
    
    grid%i_low   = i_low
    grid%j_low   = j_low
    grid%i_high  = i_high
    grid%j_high  = j_high
    grid%ig_low  = ig_low
    grid%jg_low  = jg_low
    grid%ig_high = ig_high
    grid%jg_high = jg_high
    
    allocate( grid%x(       ig_low:ig_high+1,jg_low:jg_high+1), &
              grid%y(       ig_low:ig_high+1,jg_low:jg_high+1), &
              grid%A_xi(    ig_low:ig_high+1,jg_low:jg_high), &
              grid%A_eta(   ig_low:ig_high,jg_low:jg_high+1), &
              grid%n_xi_x(  ig_low:ig_high+1,jg_low:jg_high), &
              grid%n_xi_y(  ig_low:ig_high+1,jg_low:jg_high), &
              grid%n_eta_x( ig_low:ig_high,jg_low:jg_high+1), &
              grid%n_eta_y( ig_low:ig_high,jg_low:jg_high+1), &
              grid%V(     ig_low:ig_high  ,jg_low:jg_high)   )
    grid%x = zero
    grid%y = zero
    read(12,*) (((grid%x(i,j),i=i_low,i_high+1),j=j_low,j_high+1),k=1,kmax),  &
               (((grid%y(i,j),i=i_low,i_high+1),j=j_low,j_high+1),k=1,kmax),  &
               (((zztemp,i=1,imax),j=1,jmax),k=1,kmax)
    
    !======ALTERNATE FORM (for reading 2D Grid from 3D File======
    !open(unit=12,file=grid_name)
    !read(12,*) nzones
    !read(12,*) imax, jmax, kmax
    !! Read in x-coordinate
    !do k = 1, kmax
    !do j = j_low, j_high+1
    !do i = i_low, i_high+1
    !  read(12,*) grid%x(i,j)
    !enddo
    !enddo
    !enddo
    !! Read in y-coordinate
    !do k = 1, kmax
    !do j = j_low, j_high+1
    !do i = i_low, i_high+1
    !  read(12,*) grid%y(i,j)
    !enddo
    !enddo
    !enddo
    
    close(unit=12)
    
    ! extrapolate coordinates to form ghost cells
    do j = j_low-1, jg_low, -1
    do i = i_low, i_high+1
      grid%x(i,j) = two*grid%x(i,j+1) - grid%x(i,j+2)
      grid%y(i,j) = two*grid%y(i,j+1) - grid%y(i,j+2)
    end do
    end do
    
    do j = j_high+2, jg_high+1
    do i = i_low, i_high+1
      grid%x(i,j) = two*grid%x(i,j-1) - grid%x(i,j-2)
      grid%y(i,j) = two*grid%y(i,j-1) - grid%y(i,j-2)
    end do
    end do
   !! 
    do j = j_low, j_high+1
    do i = i_low-1, ig_low, -1
      grid%x(i,j) = two*grid%x(i+1,j) - grid%x(i+2,j)
      grid%y(i,j) = two*grid%y(i+1,j) - grid%y(i+2,j)
    end do
    end do
    
    do j = j_low, j_high+1
    do i = i_high+2, ig_high+1
      grid%x(i,j) = two*grid%x(i-1,j) - grid%x(i-2,j)
      grid%y(i,j) = two*grid%y(i-1,j) - grid%y(i-2,j)
    end do
    end do
   !!!! 
    do j = j_low-1, jg_low,-1
    do i = i_low-1, ig_low,-1
      grid%x(i,j) = grid%x(i+1,j+1)
      grid%y(i,j) = grid%y(i+1,j+1)
    end do
    end do
    
    do j = j_high+2, jg_high+1
    do i = i_low-1, ig_low,-1
      grid%x(i,j) = grid%x(i+1,j-1)
      grid%y(i,j) = grid%y(i+1,j-1)
    end do
    end do
   !! 
    do j = j_low-1, jg_low,-1
    do i = i_high+2, ig_high+1
      grid%x(i,j) = grid%x(i-1,j+1)
      grid%y(i,j) = grid%y(i-1,j+1)
    end do
    end do
    
    do j = j_high+2, jg_high+1
    do i = i_high+2, ig_high+1
      grid%x(i,j) = grid%x(i-1,j-1)
      grid%y(i,j) = grid%y(i-1,j-1)
    end do
    end do
    
    ! calculate cell face areas, normals, and volumes
    !do j = jg_low, jg_high
    !do i = ig_low, ig_high
    !  grid%A_xi(i,j) =  sqrt( (grid%x(i,j+1)-grid%x(i,j))**2 &
    !                        + (grid%y(i,j+1)-grid%y(i,j))**2 )
    !  grid%A_eta(i,j) = sqrt( (grid%x(i+1,j)-grid%x(i,j))**2 &
    !                        + (grid%y(i+1,j)-grid%y(i,j))**2 )
    !  grid%n_xi_x(i,j)   =  ( (grid%y(i,j+1)-grid%y(i,j)) )!/grid%A_xi(i,j) 
    !  grid%n_xi_y(i,j)   = -( (grid%x(i,j+1)-grid%x(i,j)) )!/grid%A_xi(i,j)
    !  grid%n_eta_x(i,j)  =  ( (grid%y(i+1,j)-grid%y(i,j)) )!/grid%A_eta(i,j)
    !  grid%n_eta_y(i,j)  = -( (grid%x(i+1,j)-grid%x(i,j)) )!/grid%A_eta(i,j)
    !  call cell_volume( (/ grid%x(i,j), grid%y(i,j) /), &
    !                (/ grid%x(i+1,j), grid%y(i+1,j) /), &
    !            (/ grid%x(i+1,j+1), grid%y(i+1,j+1) /), &
    !                (/ grid%x(i,j+1), grid%y(i,j+1) /), &
    !                                  grid%V(i,j)       )
    !end do
    !end do
    
    !j = jg_high+1
    !do i = ig_low, ig_high
    !  grid%A_eta(i,j) = sqrt( (grid%x(i+1,j)-grid%x(i,j))**2 &
    !                        + (grid%y(i+1,j)-grid%y(i,j))**2 )
    !  grid%n_eta_x(i,j)  =  ( (grid%y(i+1,j)-grid%y(i,j)) )!/grid%A_eta(i,j)
    !  grid%n_eta_y(i,j)  = -( (grid%x(i+1,j)-grid%x(i,j)) )!/grid%A_eta(i,j)
    !end do
    !
    !i = ig_high+1
    !do j = jg_low, jg_high
    !  grid%A_xi(i,j) = sqrt( (grid%x(i,j+1)-grid%x(i,j))**2 &
    !                       + (grid%y(i,j+1)-grid%y(i,j))**2 )
    !  grid%n_xi_x(i,j) =   ( (grid%y(i,j+1)-grid%y(i,j)) )/grid%A_xi(i,j)
    !  grid%n_xi_y(i,j) =  -( (grid%x(i,j+1)-grid%x(i,j)) )/grid%A_xi(i,j)
    !end do
 
!==========================SAMPLE SOLUTION OUTPUT==========================
!======== This writes out a Tecplot output file in "point" format =========

  ! Note: write this once at the beginning
  open(40,file='2dEuler-curvilinear.dat',status='unknown')
  write(40,*) 'TITLE = "2D Curvilinear Euler Equations Field Data"'
  write(40,*)'variables="x(m)""y(m)"'

  ! Write this every time you want to output the solution
  !    => assumes n has been passed as the iteration #
  write(40,*) 'zone T="',1,'" '
  write(40,*) 'I=',imax+2*n_ghost,' J=',jmax+2*n_ghost
  write(40,*) 'DATAPACKING=POINT'  ! Tecplot "point" format
  do j = jg_low, jg_high+1
  do i = ig_low, ig_high+1
    ! Note: V is assumed to the vector of primitive variables V
    !       = [rho, u, v, p]^T interpolated to the grid nodes
    write(40,*)grid%x(i,j),grid%y(i,j)
  enddo
  enddo
   
    
  end subroutine allocate_grid
  
  subroutine cell_volume(A,B,C,D,V)
    
    real(prec), dimension(2), intent(in)  :: A, B, C, D
    real(prec), intent(out) :: V
    real(prec), dimension(2) :: AC, BD
    
    AC = C - A
    BD = D - B
    V = half*abs( AC(1)*BD(2) - AC(2)*BD(1) )
    
  end subroutine cell_volume
  !============================= deallocate_grid =============================80
  !>
  !! Description: 
  !!
  !! Inputs:      grid :
  !!
  !! Outputs:     grid : 
  !<
  !===========================================================================80
  subroutine deallocate_grid( grid )
    
    type( grid_t ), intent( inout ) :: grid
    
    deallocate( grid%x, grid%y, grid%A_xi, grid%A_eta, grid%V, &
                grid%n_xi_x, grid%n_xi_y, grid%n_eta_x, grid%n_eta_y )
    
  end subroutine deallocate_grid
  
end module grid_type

