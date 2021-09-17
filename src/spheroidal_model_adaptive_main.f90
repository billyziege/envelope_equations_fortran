! RMS envelope integration code.
! Integration is conducted by the leap frog alogorithm
! as written by Brandon Zerbe.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Main Program
!	- Evolves the provided statistics according to the RMS envelope equation
!     with assumed cylindrically symmetry
!-------------------------------------------------------------------------------
program main

  use constants
  use myarray
  use spheroidal_envelope
  use f90argparser

  !Declare the variables
  implicit none
  type(my_array), dimension(4) :: s !s(1)=sx, s(2)=sz,s(3)=dsx/dt,s(4)=dsz/dt
  type(my_array), dimension(2) :: em
  type(my_array) :: ts, report_t
  type(my_constants) :: consts
  integer :: every, N, stat, i, j
  real*8 :: every_real, N_real
  real*8 :: ti, tf , q, m, gam, frac
  type(opt_container):: opt_c
  logical :: debug

  !Define what to expect as input and provide the user help.

  call opt_container_init(opt_c,15)
  opt_c%description = 'Runs the evolution of the spheroidal envelope equation reporting the standard deviaitons at every steps.'
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
                   'The initial time derivative of the longitudinal standard deviation in m/s.  Use n for negative.', &
                   has_arg = .true. )
  call add_option( opt_c, &
                   'The conserved transverse emittance in m.', &
                   has_arg = .true. )
  call add_option( opt_c, &
                   'The conserved longitudinal emittance in m.', &
                   has_arg = .true. )
  call add_option( opt_c, 'The final time.', has_arg = .true. )
  call add_option( opt_c, 'The number of particles in the uniform distribution.', has_arg = .true. )
  call add_option( opt_c, 'The fraction used for the adaptive time step.', has_arg = .true., long='frac', val='0.001' )
  call add_option( opt_c, 'The charge, in Coulombs, of each particle.', long = 'charge', short='q', &
                   has_arg = .true., val = '1.602176634e-19')
  call add_option( opt_c, 'The mass, in kg, of each particle.', long = 'mass', short='m', &
                   has_arg = .true., val = '9.10938356e-31')
  call add_option( opt_c, 'The Lorentz factor of the beam.', long = 'gamma', short='g', &
                   has_arg = .true., val = '1.')
  call add_option( opt_c, 'Specifies the regular interval at which the reporting mechanism triggers.  Default is 100 times.', &
                   long='every', has_arg = .true., val = '100' )
  call add_option( opt_c, 'Specifies the initial time of the simulation.', &
                   long='t_initial', short='t', has_arg = .true., val = '0' )
  call add_option( opt_c, 'Flag that turns on debugging.', &
                   long='debug', short="d", has_arg = .false. )

  !Fill the options container from the command line
  call process_cli( opt_c )
  
  !Store the command line and/or default values to the appropriate variables.
  
  do j = 1, 2, 1
    s(j) = new_my_array(1)
    s(j+2) = new_my_array(1)
    em(j) = new_my_array(1)
  end do
  report_t = new_my_array(1)
  read(opt_c%opts(1)%val,*) s(1)%a(1)
  read(opt_c%opts(2)%val,*) s(2)%a(1)
  if (opt_c%opts(3)%val(1:1) == 'n') then
    read(opt_c%opts(3)%val(2:),*) s(3)%a(1)
    s(3)%a(1) = -1*s(3)%a(1)
  else
    read(opt_c%opts(3)%val,*) s(3)%a(1)
  end if
  if (opt_c%opts(4)%val(1:1) == 'n') then
    read(opt_c%opts(4)%val(2:),*) s(4)%a(1)
    s(4)%a(1) = -1*s(4)%a(1)
  else
    read(opt_c%opts(4)%val,*) s(4)%a(1)
  end if
  read(opt_c%opts(5)%val,*) em(1)%a(1)
  read(opt_c%opts(6)%val,*) em(2)%a(1)
  read(opt_c%opts(7)%val,*) tf
  read(opt_c%opts(8)%val,*) N_real
  read(opt_c%opts(9)%val,*) frac
  read(opt_c%opts(10)%val,*) q
  read(opt_c%opts(11)%val,*) m
  read(opt_c%opts(12)%val,*) gam
  read(opt_c%opts(13)%val,*) every_real
  read(opt_c%opts(14)%val,*) ti
  if (opt_c%opts(15)%val == "") then
    debug = .false.
  else
    debug = .true.
  end if

  if (debug) then
    write(*,*) "Through options"
  end if

  N = int(N_real)
  every = int(every_real)
  ts = new_my_array(every)
  call lin_range_my_array(ts,ti,tf,stat)
  

  !Run the integration
  consts = new_my_constants(gam,q,m)
  report_t%a(1) = ti
  call report_my_arrays(report_t,4,s)
  do i = 2, ts%N, 1
    call run_adaptive_envelope(N,s,em,ts%a(i-1),ts%a(i),frac,consts,debug)
    report_t%a(1) = ts%a(i)
    !call report_my_arrays(report_t,4,s)
  end do


  do j = 1, 2, 1
    call free_my_array(s(j))
    call free_my_array(s(j+2))
    call free_my_array(em(j))
  end do

end program main
