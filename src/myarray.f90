module myarray
!-------------------------------------------------------------------------------
! Module my_array:
! A module creating an easy to work with array object.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!type my_array
! A data type containing the number in an array and the array itself.
! Used for doing multiple initial conditions simultaneously.
!
!subroutine my_array_init:
!  -  Creates an instance of my_array
! input:
!   -  N = The number of elements
! output:
!   -  my_a_out = The instance of my_array with N blocks allocated.
!
!subroutine my_array_scaffold:
!  -  Creates an instance of my_array of the same size
! input:
!   -  my_a_in = The already created my_array
! output:
!   -  my_a_out = The instance of my_array that is the same size as the input blocks allocated.
!
!subroutine my_array_clone:
!  -  Copies the content of my_array into a second my_array
! input:
!   -  my_a_in = The already created my_array
! output:
!   -  my_a_out = The instance of my_array that is that has the same content (and size)
!
!subroutine my_array_to_str
!  -  Converts the array to a string
! input:
!   my_a_in = The already created and filled my_array
!   delimeter (optional) = The 1 character delimiter. Default is " ".
! output:
!   output_str = The contents of the array as a str.
!subroutine free_my_array:
!  -  Garbage collects my_array
! input:
!   -  my_a_in = The already created my_array
! output:
!   -  None
!
implicit none

type my_array !Contains the details of a single instance of similar data.
    sequence
    integer :: N ! The number elements
	real*8, dimension(:), allocatable :: a ! The array itself
end type my_array

contains

function new_my_array(N)
  integer, intent(in) :: N
  type(my_array) :: new_my_array
  new_my_array%N = N
  allocate(new_my_array%a(N))
end function new_my_array

subroutine clone_my_array(my_a_in,my_a_out)
  type(my_array), intent(in) :: my_a_in 
  type(my_array), intent(out) :: my_a_out 

  my_a_out = new_my_array(my_a_in%N)
  my_a_out%a(:) = my_a_in%a(:)
end subroutine clone_my_array

subroutine single_value_my_array(my_a_inout,val)
  type(my_array), intent(inout) :: my_a_inout
  real*8, intent(in) :: val 

  my_a_inout%a(:) = val
end subroutine single_value_my_array

subroutine lin_range_my_array(my_a_inout,min_val,max_val,stat)
  type(my_array), intent(inout) :: my_a_inout
  real*8, intent(in) :: min_val, max_val 
  integer, intent(out) :: stat
  real*8 :: step
  integer :: i
  
  if (my_a_inout%N < 2) then
    stat = 1
  else
    step = (max_val - min_val)/float(my_a_inout%N-1)
    do i = 1, my_a_inout%N, 1
      my_a_inout%a(i) = min_val + float(i-1)*step
    end do
  end if
end subroutine lin_range_my_array

function table_column_into_my_array(a2d,Nrow,Ncol,colN)
  real*8, dimension(Nrow,Ncol), intent(in) :: a2d
  integer, intent(in) :: Nrow, Ncol, colN
  type(my_array) :: table_column_into_my_array
  integer i

  table_column_into_my_array = new_my_array(Nrow)
  do i = 1, Nrow, 1
    table_column_into_my_array%a(i) = a2d(i,colN)
  end do
end function table_column_into_my_array

subroutine report_my_arrays(label_array,N,my_arrays,delimiter)
  type(my_array), intent(in) :: label_array
  integer, intent(in) :: N
  type(my_array), intent(in), dimension(N) :: my_arrays
  character(16) :: temp_str, myformat
  character(17*(N+1)) :: line_str
  character(1), intent(in), optional :: delimiter
  character(1) :: d
  integer :: i, j

  if ( present(delimiter) ) then
    d = delimiter
  else
    d = " "
  end if

  myformat = '(E13.6E3)'
  do i = 1, label_array%N, 1
    write(line_str,myformat) label_array%a(i)
    do j = 1, N, 1
      write(temp_str,myformat) my_arrays(j)%a(i)
      write(line_str,*) adjustl(trim(line_str)),d,adjustl(trim(temp_str))
    end do
    write(*,*) adjustl(trim(line_str))
  end do 
end subroutine report_my_arrays

subroutine stringify_my_array(my_a_in,output_str,delimeter)
  type(my_array), intent(in) :: my_a_in 
  character(my_a_in%N*17), intent(out) :: output_str
  character(1), intent(in), optional :: delimeter
  character(1) :: d
  character(16) :: current_element, myformat
  integer :: i
  
  if ( present(delimeter) ) then
    d = delimeter
  else
    d = " "
  end if

  myformat = '(E13.6E3)'
  write(output_str,myformat) my_a_in%a(1)
  if (my_a_in%N > 1) then
    do i = 2, my_a_in%N, 1
      write(current_element,myformat) my_a_in%a(i)
      write(output_str,*) trim(output_str),d,trim(current_element)
      !write(*,*) output_str
    end do
  endif
end subroutine stringify_my_array
   
subroutine free_my_array(my_a_in)
  type(my_array), intent(inout) :: my_a_in 
  integer :: ierr
  deallocate(my_a_in%a, stat = ierr)
  if (ierr /= 0) then
    write(*,*) "my_array could not be freed. Stat =", ierr
  end if
end subroutine free_my_array

end module myarray
