! RMS envelope integration code.
! Integration is conducted by the leap frog alogorithm
! as written by Brandon Zerbe.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Main Program
!	- Evolves the provided statistics according to a range of longitudinal 
!     chirps to the desired time.  Returns the statistics at that time for
!     each chirp.
!-------------------------------------------------------------------------------
program main

  use constants
  use myarray
  use spheroidal_envelope
  use f90argparser

  !Declare the variables
  implicit none
  type(my_array), dimension(4) :: s
  type(my_array), dimension(2) :: em
  type(my_array) :: muz0
  type(my_constants) :: consts
  integer :: N, Nchirps
  real*8 ::  N_real
  real*8 :: ti, tf, frac, q, m, gam
  real*8, dimension(2) :: s0, em0
  real*8 :: mux0, muzmin, muzmax
  integer :: i, stat
  type(opt_container):: opt_c
  logical :: debug

  !Define what to expect as input and provide the user help.

  call opt_container_init(opt_c,17)
  opt_c%description = 'Simultaneously runs the evolution of the spheroidal envelope equation for a range of longitudinal chirps.'
  call add_option( opt_c, &
                   'The initial transverse standard deviation in m.', &
                   has_arg = .true. )
  call add_option( opt_c, &
                   'The initial longitudinal standard deviation in m.', &
                   has_arg = .true. )
  call add_option( opt_c, &
                   'The initial time derivative of the transverse standard deviation in m/s. Use n for negative.', &
                   has_arg = .true. )
  call add_option( opt_c, &
                   'The minimum initial time derivative of the longitudinal standard deviation in m/s.  Use n for negative.', &
                   has_arg = .true. )
  call add_option( opt_c, &
                   'The maximum initial time derivative of the longitudinal standard deviation in m/s.  Use n for negative.', &
                   has_arg = .true. )
  call add_option( opt_c, &
                   'The total number of chirps to scan.', &
                   has_arg = .true. )
  call add_option( opt_c, &
                   'The conserved transverse emittance in m.', &
                   has_arg = .true. )
  call add_option( opt_c, &
                   'The conserved longitudinal emittance in m.', &
                   has_arg = .true. )
  call add_option( opt_c, 'The final time.', has_arg = .true. )
  call add_option( opt_c, 'The number of particles in the uniform distribution.', has_arg = .true. )
  call add_option( opt_c, 'The fraction of total envelope per time step.', long='frac',has_arg = .true. ,val='0.0001')
  call add_option( opt_c, 'The charge, in Coulombs, of each particle.', long = 'charge', short='q', &
                   has_arg = .true., val = '1.602176634e-19')
  call add_option( opt_c, 'The mass, in kg, of each particle.', long = 'mass', short='m', &
                   has_arg = .true., val = '9.10938356e-31')
  call add_option( opt_c, 'The Lorentz factor of the beam.', long = 'gamma', short='g', &
                   has_arg = .true., val = '1.')
  call add_option( opt_c, 'Specifies the initial time of the simulation.', &
                   long='t_initial', short='t', has_arg = .true., val = '0' )
  call add_option( opt_c, 'Flag that turns on debugging.', &
                   long='debug', short="d", has_arg = .false. )

  !Fill the options container from the command line
  call process_cli( opt_c )
  
  !Store the command line and/or default values to the appropriate variables.
  
  read(opt_c%opts(1)%val,*) s0(1)
  read(opt_c%opts(2)%val,*) s0(2)
  if (opt_c%opts(3)%val(1:1) == 'n') then
    read(opt_c%opts(3)%val(2:),*) mux0
    mux0 = -mux0
  else
    read(opt_c%opts(3)%val,*) mux0
  end if
  if (opt_c%opts(4)%val(1:1) == 'n') then
    read(opt_c%opts(4)%val(2:),*) muzmin
    muzmin = -muzmin
  else
    read(opt_c%opts(4)%val,*) muzmin
  end if
  if (opt_c%opts(5)%val(1:1) == 'n') then
    read(opt_c%opts(5)%val(2:),*) muzmax
    mux0 = -muzmax
  else
    read(opt_c%opts(5)%val,*) muzmax
  end if
  read(opt_c%opts(6)%val,*) Nchirps
  read(opt_c%opts(7)%val,*) em0(1)
  read(opt_c%opts(8)%val,*) em0(2)
  read(opt_c%opts(9)%val,*) tf
  read(opt_c%opts(10)%val,*) N_real
  read(opt_c%opts(11)%val,*) frac
  read(opt_c%opts(12)%val,*) q
  read(opt_c%opts(13)%val,*) m
  read(opt_c%opts(14)%val,*) gam
  read(opt_c%opts(15)%val,*) ti
  if (opt_c%opts(16)%val == "") then
    debug = .false.
  else
    debug = .true.
  end if

  if (debug) then
    write(*,*) "Through options."
  end if

  N = int(N_real)

  !Initialize arrays for scanning
  do i = 1, 2, 1
    s(i) = new_my_array(Nchirps)
    s(i+2) = new_my_array(Nchirps)
    em(i) = new_my_array(Nchirps)
  end do

  if (debug) then
    write(*,*) "Array initialized."
  end if

  !Put in initial values
  do i = 1, 2, 1
    s(i)%a(:) = s0(i)
    em(i)%a(:) = em0(i)
  end do
  s(3)%a(:) = mux0
  call lin_range_my_array(s(4),muzmin,muzmax,stat)
  call clone_my_array(s(4),muz0)

  if (debug) then
    write(*,*) "Arrays filled."
  end if

  !Run the integration
  consts = new_my_constants(gam,q,m)
  call run_adaptive_envelope(N,s,em,ti,tf,frac,consts,debug)
  call report_my_arrays(muz0,4,s)

  if (debug) then
    write(*,*) "Before loop."
  end if

  do i = 1, 2, 1
    call free_my_array(s(i))
    call free_my_array(s(i+2))
    call free_my_array(em(i))
  end do

end program main
