module grid_type
  
  use set_precision, only : prec
  use set_constants, only : zero, one, two, half
  use set_inputs,    only : i_low, i_high, ig_low, ig_high
  use set_inputs,    only : j_low, j_high, jg_low, jg_high
  
  implicit none
  
  private
  
  public :: grid_t, allocate_grid, deallocate_grid
  public :: ghost_shape, ghost_shape_C, cell_geometry
  
  type grid_t
    
    sequence
    
    integer :: imax, jmax, n_ghost
    integer :: i_low, i_high, ig_low, ig_high
    integer :: j_low, j_high, jg_low, jg_high
    
    
    real(prec), allocatable, dimension(:,:) :: x, y
    real(prec), allocatable, dimension(:,:) :: xc, yc
    real(prec), allocatable, dimension(:,:) :: A_xi, A_eta
    real(prec), allocatable, dimension(:,:,:) :: n_xi, n_eta
    real(prec), allocatable, dimension(:,:,:) :: n_xi_avg
    real(prec), allocatable, dimension(:,:,:) :: n_eta_avg
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
    
    use set_inputs, only : grid_name, n_ghost
    
    type( grid_t ), intent( inout ) :: grid
    
    allocate( grid%x(       ig_low:ig_high+1,jg_low:jg_high+1), &
              grid%y(       ig_low:ig_high+1,jg_low:jg_high+1), &
              grid%xc(      ig_low:ig_high,jg_low:jg_high), &
              grid%yc(      ig_low:ig_high,jg_low:jg_high), &
              grid%A_xi(    ig_low:ig_high+1,jg_low:jg_high), &
              grid%A_eta(   ig_low:ig_high,jg_low:jg_high+1), &
              grid%n_xi(    ig_low:ig_high+1,jg_low:jg_high,2),   &
              grid%n_eta(   ig_low:ig_high,jg_low:jg_high+1,2),   &
              grid%n_xi_avg(  ig_low:ig_high,jg_low:jg_high,2), &
              grid%n_eta_avg( ig_low:ig_high,jg_low:jg_high,2), &
              grid%V(       ig_low:ig_high,jg_low:jg_high)   )
   grid%x = zero
   grid%y = zero
    
  end subroutine allocate_grid
  
  subroutine ghost_shape_C(grid,low,high)
    use set_inputs, only : n_ghost, imax
    type(grid_t), intent(inout) :: grid
    integer, intent(inout) :: low, high
    real(prec) :: x1, x2, x3, x4, y1, y2, y3, y4, m1, m2
    integer :: i, j, d1
    integer, allocatable, dimension(:,:) :: i1, i2, j1, j2
    

    ! extrapolate coordinates to form ghost cells
    ! "bottom"
    do j = j_low-1, jg_low, -1
    do i = i_low, i_high+1
      grid%x(i,j) = two*grid%x(i,j+1) - grid%x(i,j+2)
      grid%y(i,j) = two*grid%y(i,j+1) - grid%y(i,j+2)
    end do
    end do
    !do j = j_low-1, jg_low, -1
    !do i = i_low, i_high+1
    !  grid%x(i,j) = half*(3.0_prec*grid%x(i,j+1) - grid%x(i,j+2))
    !  grid%y(i,j) = half*(3.0_prec*grid%y(i,j+1) - grid%y(i,j+2))
    !end do
    !end do


    ! "top"
    do j = j_high+2, jg_high+1
    do i = i_low, i_high+1
      grid%x(i,j) = two*grid%x(i,j-1) - grid%x(i,j-2)
      grid%y(i,j) = two*grid%y(i,j-1) - grid%y(i,j-2)
    end do
    end do
    
    ! "left"
    do j = j_low, j_high+1
    do i = i_low-1, ig_low, -1
      grid%x(i,j) = two*grid%x(i+1,j) - grid%x(i+2,j)
      grid%y(i,j) = two*grid%y(i+1,j) - grid%y(i+2,j)
    end do
    end do
    
    ! "right"
    do j = j_low, j_high+1
    do i = i_high+2, ig_high+1
      grid%x(i,j) = two*grid%x(i-1,j) - grid%x(i-2,j)
      grid%y(i,j) = two*grid%y(i-1,j) - grid%y(i-2,j)
    end do
    end do
   !!!!
    
    
    ! "top left"
    do j = j_high+2, jg_high+1
    do i = i_low-1, ig_low,-1
      x1 = grid%x(i+2,j)
      x2 = grid%x(i+1,j)
      x3 = grid%x(i,j-2)
      x4 = grid%x(i,j-1)
      
      y1 = grid%y(i+2,j)
      y2 = grid%y(i+1,j)
      y3 = grid%y(i,j-2)
      y4 = grid%y(i,j-1)
      
      call extrap_corner( (/x1,y1/), (/x2,y2/), (/x3,y3/), (/x4,y4/), &
                          grid%x(i,j), grid%y(i,j) )
    end do
    end do
     
    
    ! "top right"
    do j = j_high+2, jg_high+1
    do i = i_high+2, ig_high+1
      x1 = grid%x(i-2,j)
      x2 = grid%x(i-1,j)
      x3 = grid%x(i,j-2)
      x4 = grid%x(i,j-1)
      
      y1 = grid%y(i-2,j)
      y2 = grid%y(i-1,j)
      y3 = grid%y(i,j-2)
      y4 = grid%y(i,j-1)
      
      call extrap_corner( (/x1,y1/), (/x2,y2/), (/x3,y3/), (/x4,y4/), &
                          grid%x(i,j), grid%y(i,j) )
      
    end do
    end do
    
    
    high = high
    low = low-n_ghost
    d1 = high-low +1
    allocate( i1(d1,n_ghost), j1(d1,n_ghost), &
              i2(d1,n_ghost), j2(d1,n_ghost)  )
    do j = 1, n_ghost
      do i = 1, d1
        i1(i,j) = low - 1 + i
        i2(i,j) = imax - low + 2 - i
        j1(i,j) = 1 - j
        j2(i,j) = j + 1
      end do
    end do
    
    do j = 1, n_ghost
    do i = 1, d1
      grid%x(i1(i,j),j1(i,j)) = grid%x(i2(i,j),j2(i,j))
      grid%y(i1(i,j),j1(i,j)) = grid%y(i2(i,j),j2(i,j))

      grid%x(i2(i,j),j1(i,j)) = grid%x(i1(i,j),j2(i,j))
      grid%y(i2(i,j),j1(i,j)) = grid%y(i1(i,j),j2(i,j))
    end do
    end do
    
    deallocate(i1,i2,j1,j2)
    
    ! "bottom left"
    do j = j_low-1, jg_low,-1
    do i = i_low-1, ig_low,-1
      x1 = grid%x(i+2,j)
      x2 = grid%x(i+1,j)
      x3 = grid%x(i,j+2)
      x4 = grid%x(i,j+1)
      
      y1 = grid%y(i+2,j)
      y2 = grid%y(i+1,j)
      y3 = grid%y(i,j+2)
      y4 = grid%y(i,j+1)
    !  call extrap_corner( (/x1,y1/), (/x2,y2/), (/x3,y3/), (/x4,y4/), &
    !                      grid%x(i,j), grid%y(i,j) )
      grid%x(i,j) = two*x4 - x3
      grid%y(i,j) = two*y4 - y3
    end do
    end do

    
    ! "bottom right"
    do j = j_low-1, jg_low,-1
    do i = i_high+2, ig_high+1
      x1 = grid%x(i-2,j)
      x2 = grid%x(i-1,j)
      x3 = grid%x(i,j+2)
      x4 = grid%x(i,j+1)
      
      y1 = grid%y(i-2,j)
      y2 = grid%y(i-1,j)
      y3 = grid%y(i,j+2)
      y4 = grid%y(i,j+1)
      
     ! call extrap_corner( (/x1,y1/), (/x2,y2/), (/x3,y3/), (/x4,y4/), &
     !                     grid%x(i,j), grid%y(i,j) )
      grid%x(i,j) = two*x4 - x3
      grid%y(i,j) = two*y4 - y3
    end do
    end do
    
    
    
    
  end subroutine ghost_shape_C
   
  subroutine ghost_shape(grid)
    
    type(grid_t), intent(inout) :: grid
    real(prec) :: x1, x2, x3, x4, y1, y2, y3, y4, m1, m2
    integer :: i, j
    
    ! extrapolate coordinates to form ghost cells
    ! "bottom"
    do j = j_low-1, jg_low, -1
    do i = i_low, i_high+1
      grid%x(i,j) = two*grid%x(i,j+1) - grid%x(i,j+2)
      grid%y(i,j) = two*grid%y(i,j+1) - grid%y(i,j+2)
    end do
    end do
    
    ! "top"
    do j = j_high+2, jg_high+1
    do i = i_low, i_high+1
      grid%x(i,j) = two*grid%x(i,j-1) - grid%x(i,j-2)
      grid%y(i,j) = two*grid%y(i,j-1) - grid%y(i,j-2)
    end do
    end do
    
    ! "left"
    do j = j_low, j_high+1
    do i = i_low-1, ig_low, -1
      grid%x(i,j) = two*grid%x(i+1,j) - grid%x(i+2,j)
      grid%y(i,j) = two*grid%y(i+1,j) - grid%y(i+2,j)
    end do
    end do
    
    ! "right"
    do j = j_low, j_high+1
    do i = i_high+2, ig_high+1
      grid%x(i,j) = two*grid%x(i-1,j) - grid%x(i-2,j)
      grid%y(i,j) = two*grid%y(i-1,j) - grid%y(i-2,j)
    end do
    end do
   !!!!
    
    ! "bottom left"
    do j = j_low-1, jg_low,-1
    do i = i_low-1, ig_low,-1
      x1 = grid%x(i+2,j)
      x2 = grid%x(i+1,j)
      x3 = grid%x(i,j+2)
      x4 = grid%x(i,j+1)
      
      y1 = grid%y(i+2,j)
      y2 = grid%y(i+1,j)
      y3 = grid%y(i,j+2)
      y4 = grid%y(i,j+1)
      
      call extrap_corner( (/x1,y1/), (/x2,y2/), (/x3,y3/), (/x4,y4/), &
                          grid%x(i,j), grid%y(i,j) )
    end do
    end do
    
    ! "top left"
    do j = j_high+2, jg_high+1
    do i = i_low-1, ig_low,-1
      x1 = grid%x(i+2,j)
      x2 = grid%x(i+1,j)
      x3 = grid%x(i,j-2)
      x4 = grid%x(i,j-1)
      
      y1 = grid%y(i+2,j)
      y2 = grid%y(i+1,j)
      y3 = grid%y(i,j-2)
      y4 = grid%y(i,j-1)
      
      call extrap_corner( (/x1,y1/), (/x2,y2/), (/x3,y3/), (/x4,y4/), &
                          grid%x(i,j), grid%y(i,j) )
    end do
    end do
     
    ! "bottom right"
    do j = j_low-1, jg_low,-1
    do i = i_high+2, ig_high+1
      x1 = grid%x(i-2,j)
      x2 = grid%x(i-1,j)
      x3 = grid%x(i,j+2)
      x4 = grid%x(i,j+1)
      
      y1 = grid%y(i-2,j)
      y2 = grid%y(i-1,j)
      y3 = grid%y(i,j+2)
      y4 = grid%y(i,j+1)
      
      call extrap_corner( (/x1,y1/), (/x2,y2/), (/x3,y3/), (/x4,y4/), &
                          grid%x(i,j), grid%y(i,j) )
    end do
    end do
    
    ! "top right"
    do j = j_high+2, jg_high+1
    do i = i_high+2, ig_high+1
      x1 = grid%x(i-2,j)
      x2 = grid%x(i-1,j)
      x3 = grid%x(i,j-2)
      x4 = grid%x(i,j-1)
      
      y1 = grid%y(i-2,j)
      y2 = grid%y(i-1,j)
      y3 = grid%y(i,j-2)
      y4 = grid%y(i,j-1)
      
      call extrap_corner( (/x1,y1/), (/x2,y2/), (/x3,y3/), (/x4,y4/), &
                          grid%x(i,j), grid%y(i,j) )
      
    end do
    end do
    

  end subroutine ghost_shape
  
  subroutine extrap_corner(P1,P2,P3,P4,x,y)
    
    real(prec), dimension(2), intent(in) :: P1,P2,P3,P4
    real(prec), intent(out) :: x,y
    real(prec) :: D, x1,x2,x3,x4, y1,y2,y3,y4
    
    x1 = P1(1)
    x2 = P2(1)
    x3 = P3(1)
    x4 = P4(1)
    
    y1 = P1(2)
    y2 = P2(2)
    y3 = P3(2)
    y4 = P4(2)
    
    D = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)
      
    x = ((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/D
    y = ((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/D
    
  end subroutine extrap_corner

  subroutine cell_geometry(grid)
    
    type(grid_t), intent(inout) :: grid 
    integer :: i, j
    
    ! calculate cell face areas, normals, and volumes
    do j = jg_low, jg_high
    do i = ig_low, ig_high
      grid%A_xi(i,j) =  sqrt( (grid%x(i,j+1)-grid%x(i,j))**2 &
                            + (grid%y(i,j+1)-grid%y(i,j))**2 )
      grid%A_eta(i,j) = sqrt( (grid%x(i+1,j)-grid%x(i,j))**2 &
                            + (grid%y(i+1,j)-grid%y(i,j))**2 )
      grid%n_xi(i,j,1)   = ( (grid%y(i,j+1)-grid%y(i,j)) )/grid%A_xi(i,j) 
      grid%n_xi(i,j,2)   = -( (grid%x(i,j+1)-grid%x(i,j)) )/grid%A_xi(i,j)
      grid%n_eta(i,j,1)  = -( (grid%y(i+1,j)-grid%y(i,j)) )/grid%A_eta(i,j)
      grid%n_eta(i,j,2)  = ( (grid%x(i+1,j)-grid%x(i,j)) )/grid%A_eta(i,j)
     call volume_calc( (/ grid%x(i,j), grid%x(i+1,j), grid%x(i+1,j+1), &
                                         grid%x(i,j+1),grid%x(i,j) /), &
                       (/ grid%y(i,j), grid%y(i+1,j), grid%y(i+1,j+1), &
                                         grid%y(i,j+1),grid%y(i,j) /), &
                       grid%xc(i,j), grid%yc(i,j), grid%V(i,j) )
    end do
    end do
    
    j = jg_high+1
    do i = ig_low, ig_high
      grid%A_eta(i,j) = sqrt( (grid%x(i+1,j)-grid%x(i,j))**2 &
                            + (grid%y(i+1,j)-grid%y(i,j))**2 )
      grid%n_eta(i,j,1)  =  -( (grid%y(i+1,j)-grid%y(i,j)) )/grid%A_eta(i,j)
      grid%n_eta(i,j,2)  =  ( (grid%x(i+1,j)-grid%x(i,j)) )/grid%A_eta(i,j)
    end do
    
    i = ig_high+1
    do j = jg_low, jg_high
      grid%A_xi(i,j) = sqrt( (grid%x(i,j+1)-grid%x(i,j))**2 &
                           + (grid%y(i,j+1)-grid%y(i,j))**2 )
      grid%n_xi(i,j,1) =   ( (grid%y(i,j+1)-grid%y(i,j)) )/grid%A_xi(i,j)
      grid%n_xi(i,j,2) =   -( (grid%x(i,j+1)-grid%x(i,j)) )/grid%A_xi(i,j)
    end do
    
    do j = jg_low, jg_high
      do i = ig_low, ig_high
        grid%n_xi_avg(i,j,1) = half*(grid%n_xi(i+1,j,1)+grid%n_xi(i,j,1))
        grid%n_xi_avg(i,j,2) = half*(grid%n_xi(i+1,j,2)+grid%n_xi(i,j,2))
        grid%n_eta_avg(i,j,1) = half*(grid%n_eta(i,j+1,1)+grid%n_eta(i,j,1))
        grid%n_eta_avg(i,j,2) = half*(grid%n_eta(i,j+1,2)+grid%n_eta(i,j,2))
      end do
    end do
  
  end subroutine cell_geometry
  
!  subroutine cell_volume(A,B,C,D,V)
!    
!    real(prec), dimension(2), intent(in)  :: A, B, C, D
!    real(prec), intent(out) :: V
!    real(prec), dimension(2) :: AC, BD
!    
!    AC = C - A
!    BD = D - B
!    V = half*abs( AC(1)*BD(2) - AC(2)*BD(1) )
!    
!  end subroutine cell_volume
  
  subroutine volume_calc(x,y,cx,cy,V)
    
    real(prec), dimension(:), intent(in) :: x,y
    real(prec), intent(out) :: cx, cy, V
    integer :: i, N
    
    cx = zero
    cy = zero
    V  = zero
    N = size(x)
    
    do i = 1,N-1
      V = V + x(i)*y(i+1) - x(i+1)*y(i)
    end do
    V = half*V
    
    do i = 1,N-1
      cx = cx + (x(i)+x(i+1))*(x(i)*y(i+1)-x(i+1)*y(i))
      cy = cy + (y(i)+y(i+1))*(x(i)*y(i+1)-x(i+1)*y(i))
    end do
    cx = cx/(6.0_prec*V)
    cy = cy/(6.0_prec*V)
    
   V = abs(V)
   
  end subroutine volume_calc
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
    
    deallocate( grid%x, grid%y, grid%xc, grid%yc, grid%A_xi, grid%A_eta, &
                grid%n_xi, grid%n_eta, grid%n_xi_avg, grid%n_eta_avg, grid%V )
    
  end subroutine deallocate_grid
  
end module grid_type

