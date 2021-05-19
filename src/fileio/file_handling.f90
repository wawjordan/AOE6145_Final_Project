module file_handling
  
  use set_precision, only : prec
  use set_constants, only : zero, one, two, half
  !use fluid_constants, only :
  use set_inputs, only : n_ghost
  use set_inputs, only : imax, i_low, i_high, ig_low, ig_high
  use set_inputs, only : jmax, j_low, j_high, jg_low, jg_high
  use grid_type,  only : grid_t, allocate_grid
  
  implicit none
  
  private
  
  public :: grid_in, grid_out
  
contains
  
subroutine grid_in(filename,grid)
  
  type(grid_t), intent(out) :: grid
  character(len=*), intent(in) :: filename
  integer :: fstat
  integer :: i, j, k, kmax, nzones
  real(prec) :: zztemp
  
  open(unit=12,file=filename,status='old',action='read',iostat=fstat)
  if ( fstat == 0 ) then  
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
    
    call allocate_grid(grid)
    
    if (kmax>2) then
      !! Read in x-coordinate
      do k = 1, kmax
        do j = j_low, j_high+1
          do i = i_low, i_high+1
            read(12,*) grid%x(i,j)
          enddo
        enddo
    enddo
    !! Read in y-coordinate
    do k = 1, kmax
      do j = j_low, j_high+1
        do i = i_low, i_high+1
          read(12,*) grid%y(i,j)
        enddo
      enddo
    enddo

    else
      read(12,*) &
       (((grid%x(i,j),i=i_low,i_high+1),j=j_low,j_high+1),k=1,kmax), &
       (((grid%y(i,j),i=i_low,i_high+1),j=j_low,j_high+1),k=1,kmax), &
       (((zztemp,i=1,imax),j=1,jmax),k=1,kmax)
    end if
    
    elseif (fstat > 0) then
      write(*,*) 'ERROR: error reading', trim(filename)
      stop
    end if
  
  close(12)
  imax = i_high-i_low+2
  jmax = j_high-j_low+2
  
end subroutine grid_in

subroutine grid_out(filename,grid)
  
  type(grid_t), intent(in) :: grid
  character(len=*), intent(in) :: filename
  integer :: i, j, funit
   
  funit = 50

  open(unit=funit,file=filename,status='unknown')
  write(funit,*) 'TITLE = "2D Curvilinear Mesh Data"'
  write(funit,*)'variables="x(m)","y(m)","n<sub><greek>x</greek>,x</sub>",&
    & "n<sub><greek>x</greek>,y</sub>", "n<sub><greek>n</greek>,x</sub>",&
    & "n<sub><greek>n</greek>,y</sub>", "A<sub><greek>x</greek></sub>",&
    & "A<sub><greek>n</greek></sub>","V"'

  write(funit,*) 'ZONE'
  write(funit,*) 'T = "Meshy-Mesh"'
  write(funit,*) 'I=',ig_high-ig_low+2,' J=',jg_high-jg_low+2
  write(funit,*) 'DT = (DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE,&
                   & DOUBLE, DOUBLE, DOUBLE, DOUBLE)'
  write(funit,*) 'DATAPACKING = BLOCK'
  write(funit,*) 'VARLOCATION = ([3-9]=CELLCENTERED)'
  
  
  write(funit,*) ((grid%x(i,j),i=ig_low,ig_high+1),      &
                 NEW_LINE('a'),  j=jg_low,jg_high+1)
  write(funit,*) ((grid%y(i,j),i=ig_low,ig_high+1),      &
                 NEW_LINE('a'), j=jg_low,jg_high+1)
  write(funit,*) ((grid%n_xi(i,j,1),i=ig_low,ig_high), &
                 NEW_LINE('a'),  j=jg_low,jg_high)
  write(funit,*) ((grid%n_xi(i,j,2),i=ig_low,ig_high), &
                 NEW_LINE('a'),  j=jg_low,jg_high)
  write(funit,*) ((grid%n_eta(i,j,1),i=ig_low,ig_high),&
                 NEW_LINE('a'),  j=jg_low,jg_high)
  write(funit,*) ((grid%n_eta(i,j,2),i=ig_low,ig_high),&
                 NEW_LINE('a'),  j=jg_low,jg_high)
  write(funit,*) ((grid%A_xi(i,j),i=ig_low,ig_high),     &
                 NEW_LINE('a'),  j=jg_low,jg_high)
  write(funit,*) ((grid%A_eta(i,j),i=ig_low,ig_high),    &
                 NEW_LINE('a'),  j=jg_low,jg_high)
  write(funit,*) ((grid%V(i,j),i=ig_low,ig_high),        &
                 NEW_LINE('a'),  j=jg_low,jg_high)
  
  close(funit)

  end subroutine grid_out
!subroutine grid_out(filename,grid)
!  
!  type(grid_t), intent(in) :: grid
!  character(len=*), intent(in) :: filename
!  integer :: i, j, funit
!   
!  funit = 50
!
!  open(unit=funit,file=filename,status='unknown')
!  write(funit,*) 'TITLE = "2D Curvilinear Mesh Data"'
!  write(funit,*)'variables="x(m)","y(m)","n<sub><greek>x</greek>,x</sub>",&
!    & "n<sub><greek>x</greek>,y</sub>", "n<sub><greek>n</greek>,x</sub>",&
!    & "n<sub><greek>n</greek>,y</sub>", "A<sub><greek>x</greek></sub>",&
!    & "A<sub><greek>n</greek></sub>","V"'
!
!  write(funit,*) 'ZONE'
!  write(funit,*) 'T = "Meshy-Mesh"'
!  write(funit,*) 'I=',i_high-i_low+2,' J=',j_high-j_low+2
!  write(funit,*) 'DT = (DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE,&
!                   & DOUBLE, DOUBLE, DOUBLE, DOUBLE)'
!  write(funit,*) 'DATAPACKING = BLOCK'
!  write(funit,*) 'VARLOCATION = ([7-9]=CELLCENTERED)'
!  
!  
!  write(funit,*) ((grid%x(i,j),i=i_low,i_high+1),      &
!                 NEW_LINE('a'),  j=j_low,j_high+1)
!  write(funit,*) ((grid%y(i,j),i=i_low,i_high+1),      &
!                 NEW_LINE('a'), j=j_low,j_high+1)
!  write(funit,*) ((grid%n_xi(i,j,1),i=i_low,i_high+1), &
!                 NEW_LINE('a'),  j=j_low,j_high+1)
!  write(funit,*) ((grid%n_xi(i,j,2),i=i_low,i_high+1), &
!                 NEW_LINE('a'),  j=j_low,j_high+1)
!  write(funit,*) ((grid%n_eta(i,j,1),i=i_low,i_high+1),&
!                 NEW_LINE('a'),  j=j_low,j_high+1)
!  write(funit,*) ((grid%n_eta(i,j,2),i=i_low,i_high+1),&
!                 NEW_LINE('a'),  j=j_low,j_high+1)
!!  write(funit,*) ((grid%n_xi_x(i,j),i=i_low,i_high), &
!!                 NEW_LINE('a'),  j=j_low,j_high)
!!  write(funit,*) ((grid%n_xi_y(i,j),i=i_low,i_high), &
!!                 NEW_LINE('a'),  j=j_low,j_high)
!!  write(funit,*) ((grid%n_eta_x(i,j),i=i_low,i_high),&
!!                 NEW_LINE('a'),  j=j_low,j_high)
!!  write(funit,*) ((grid%n_eta_y(i,j),i=i_low,i_high),&
!!                 NEW_LINE('a'),  j=j_low,j_high)
!  write(funit,*) ((grid%A_xi(i,j),i=i_low,i_high),     &
!                 NEW_LINE('a'),  j=j_low,j_high)
!  write(funit,*) ((grid%A_eta(i,j),i=i_low,i_high),    &
!                 NEW_LINE('a'),  j=j_low,j_high)
!  write(funit,*) ((grid%V(i,j),i=i_low,i_high),        &
!                 NEW_LINE('a'),  j=j_low,j_high)
!  
!  close(funit)
!
!  end subroutine grid_out

end module file_handling
