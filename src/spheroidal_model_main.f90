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
  type(my_array) :: sx, sz, mux, muz
  type(my_array) :: emx, emz
  type(my_constants) :: consts
  integer :: Nt, every, N
  real*8 :: Nt_real, every_real, N_real
  real*8 :: t, dt, q, m, gam
  type(opt_container):: opt_c
  logical :: debug

  !Define what to expect as input and provide the user help.

  call opt_container_init(opt_c,16)
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
  call add_option( opt_c, 'The time step (dt) in seconds.', has_arg = .true. )
  call add_option( opt_c, 'The number of time steps to run.', has_arg = .true. )
  call add_option( opt_c, 'The number of particles in the uniform distribution.', has_arg = .true. )
  call add_option( opt_c, 'The charge, in Coulombs, of each particle.', long = 'charge', short='q', &
                   has_arg = .true., val = '1.602176634e-19')
  call add_option( opt_c, 'The mass, in kg, of each particle.', long = 'mass', short='m', &
                   has_arg = .true., val = '9.10938356e-31')
  call add_option( opt_c, 'The Lorentz factor of the beam.', long = 'gamma', short='g', &
                   has_arg = .true., val = '1.')
  call add_option( opt_c, 'Specifies at what interval of steps the reporting mechanism triggers.  Default is 10 every steps.', &
                   long='every', has_arg = .true., val = '10' )
  call add_option( opt_c, 'Specifies the initial time of the simulation.', &
                   long='t_initial', short='t', has_arg = .true., val = '0' )
  call add_option( opt_c, 'Flag that turns on debugging.', &
                   long='debug', short="d", has_arg = .false. )

  !Fill the options container from the command line
  call process_cli( opt_c )
  
  !Store the command line and/or default values to the appropriate variables.
  
  call init_my_array(1,sx)
  call init_my_array(1,sz)
  call init_my_array(1,mux)
  call init_my_array(1,muz)
  call init_my_array(1,emx)
  call init_my_array(1,emz)
  read(opt_c%opts(1)%val,*) sx%a(1)
  read(opt_c%opts(2)%val,*) sz%a(1)
  if (opt_c%opts(3)%val(1:1) == 'n') then
    read(opt_c%opts(3)%val(2:),*) mux%a(1)
    mux%a(1) = -1*mux%a(1)
  else
    read(opt_c%opts(3)%val,*) mux%a(1)
  end if
  if (opt_c%opts(4)%val(1:1) == 'n') then
    read(opt_c%opts(4)%val(2:),*) muz%a(1)
    muz%a(1) = -1*muz%a(1)
  else
    read(opt_c%opts(4)%val,*) muz%a(1)
  end if
  read(opt_c%opts(5)%val,*) emx%a(1)
  read(opt_c%opts(6)%val,*) emz%a(1)
  read(opt_c%opts(7)%val,*) dt
  read(opt_c%opts(7)%val,*) Nt_real
  read(opt_c%opts(9)%val,*) N_real
  read(opt_c%opts(10)%val,*) q
  read(opt_c%opts(11)%val,*) m
  read(opt_c%opts(12)%val,*) gam
  read(opt_c%opts(13)%val,*) every_real
  read(opt_c%opts(14)%val,*) t
  if (opt_c%opts(15)%val == "") then
    debug = .false.
  else
    debug = .true.
  end if

  if (debug) then
    write(*,*) "Through options"
  end if

  N = int(N_real)
  Nt = int(Nt_real)
  every = int(every_real)

  !Run the integration
  call init_my_constants(consts,gam,q,m)
  call run_envelope(N,sx,sz,mux,muz,emx,emz,every,dt,Nt,t,consts,debug)

  call free_my_array(sx)
  call free_my_array(sz)
  call free_my_array(mux)
  call free_my_array(muz)
  call free_my_array(emx)
  call free_my_array(emz)

end program main
