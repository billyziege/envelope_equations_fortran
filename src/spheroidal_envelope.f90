module spheroidal_envelope
!-------------------------------------------------------------------------------
! Module spheroidal envelope:
! A module implementing the spheroidal envelope 
! equations integration and application.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!subroutine calc_dmudt:
!  - Calculates the second order derivative of spatial standard deviation.
! input:
!	- sx = transverse standard deviation (s_x)
!	- sz = longitudinal standard deviation (s_z)
!	- Vsqx = s_x^2 s_(v_x)^2 - s_(x,v_x)^2 (where s_(v_x) is the standard deviation of the x-velocity)
!            and s_(x,v_x) is the covariance.
!	- Vsqz = s_z^2 s_(v_z)^2 - s_(z,v_z)^2
!	- kellip = Elliptical constant including squared initial plasma frequency.
! output:
!	- dmuxdt = The second derivative of s_x
!   - dmuzdt = The second derivative of s_z
!
!Need to notate run_evelope
use constants
use myarray
implicit none

contains

subroutine geometric_parameters_one(alpha,G1,G2)
  real*8, intent(in) :: alpha  
  real*8, intent(out) :: G1, G2

  G1 = 0.5
  G2 = 0.5
  if (alpha > 1) then
    G1 =  -log(alpha + sqrt(alpha**2 -1))/(alpha**2-1)**1.5 +alpha/(alpha**2-1)
    G2 =  2*alpha**2/(alpha**2-1)*(alpha*log(alpha + sqrt(alpha**2-1))/sqrt(alpha**2-1) -1)
  else if (alpha < 1) then
    if (alpha > 1e-10) then
      G1 =  asin(sqrt(1-alpha**2))/(1-alpha**2)**1.5 -alpha/(1-alpha**2)
      G2 =  2*alpha**2/(alpha**2-1)*(alpha*asin(sqrt(1-alpha**2))/sqrt(1-alpha**2) -1)
    else
      G1 =  acos(-1.0)/2
      G2 =  0
    end if
  end if
end subroutine

subroutine geometric_parameters(alpha,G)
  type(my_array), intent(in) :: alpha  
  type(my_array), intent(inout), dimension(2) :: G

  !G1(:) = 1/sqrt(2.)
  !G2(:) = 1/(2.*sqrt(2.))
  !where (alpha<1)
  !  where (alpha > 1e-10)
  !    G1 = alpha/(alpha**2-1)*(alpha-acos(alpha)/sqrt(1.-alpha**2))
  !    G2 =  1/(alpha**2-1)*(alpha*asin(sqrt(1.-alpha**2))/sqrt(1.-alpha**2) -1)
  !  else where
  !    G1 =  acos(-1.0)
  !    G2 =  0
  !  end where
  !else where (alpha>1) 
  !  G1 =  alpha/(alpha**2-1)*(-acosh(alpha)/(alpha**2-1)**.5 +alpha)
  !  G2 = 1/(alpha**2-1)*(alpha*acosh(alpha)/sqrt(alpha**2-1.)-1.)
  !end where
    
  !axself(:) = wsq/(2*sx%a(:)*sz%a(:))*G1(:)
  !azself(:) = wsq/(sx%a(:)**2)*G2(:)

  G(1)%a(:) = 0.5
  G(2)%a(:) = 0.5
  where (alpha%a > 1) 
    G(1)%a =  -log(alpha%a + sqrt(alpha%a**2 -1))/(alpha%a**2-1)**1.5 +alpha%a/(alpha%a**2-1)
    G(2)%a =  2*alpha%a**2/(alpha%a**2-1)*(alpha%a*log(alpha%a + sqrt(alpha%a**2-1))/sqrt(alpha%a**2-1) -1)
  else where (alpha%a < 1)
    where (alpha%a > 1e-10) 
      G(1)%a =  asin(sqrt(1-alpha%a**2))/(1-alpha%a**2)**1.5 -alpha%a/(1-alpha%a**2)
      G(2)%a =  2*alpha%a**2/(alpha%a**2-1)*(alpha%a*asin(sqrt(1-alpha%a**2))/sqrt(1-alpha%a**2) -1)
    else where
      G(1)%a =  acos(-1.0)/2
      G(2)%a =  0
    end where
  end where
    
end subroutine

subroutine determine_time_step(s,dmudt,frac,dt,threshold_in)
  type(my_array), intent(in), dimension(4) :: s
  type(my_array), intent(in), dimension(2) :: dmudt
  real*8, intent(in) :: frac
  type(my_array), intent(inout) :: dt
  real*8, intent(in), optional :: threshold_in
  real*8 :: threshold
  real*8, dimension(2) :: dt_temp, inside
  integer :: i, j
  if ( present(threshold_in) ) then
    threshold = threshold_in
  else
    threshold = 1e-18
  end if
  do i = 1, s(1)%N, 1
    do j= 1, 2, 1
      if (s(2+j)%a(i) < 0) then
        inside(j) = s(2+j)%a(i)**2 - 2*frac*s(j)%a(i)*dmudt(j)%a(i)
        if (inside(j) > 0) then
          dt_temp(j) = -s(j+2)%a(i)/dmudt(j)%a(i) - sqrt(inside(j))/dmudt(j)%a(i) 
        else
          dt_temp(j) = -s(j+2)%a(i)/dmudt(j)%a(i)
        end if
      else
        inside(j) = s(j+2)%a(i)**2 + 2*frac*s(j)%a(i)*dmudt(j)%a(i)
        dt_temp(j) = -s(j+2)%a(i)/dmudt(j)%a(1) + sqrt(inside(j))/dmudt(j)%a(1) 
      end if
    end do
    dt%a(i) = min(dt_temp(1),dt_temp(2))
    if (dt%a(i) < threshold) then
      dt%a(i) = threshold
    end if
  end do
  !write(*,*) dt%a(1)
end subroutine

subroutine check_time_step(t,tf,dt)
  type(my_array), intent(in) :: t
  real*8, intent(in) :: tf
  type(my_array), intent(inout) :: dt
  where (dt%a >= tf - t%a)
    dt%a = tf - t%a
  end where
end subroutine

!subroutine calc_dmudt_one_at_a_time(sx,sz,alpha,Vsqx,Vsqz,wsq,dmuxdt,dmuzdt)
!  real*8, intent(in) :: sx, sz, alpha, Vsqx, Vsqz, wsq
!  real*8, intent(inout) :: dmuxdt, dmuzdt
!  real*8 :: G1, G2
!  real*8 :: axself, azself, axeff, azeff
  
!  call geometric_parameters(s,alpha,G1,G2)
!  axself(:) = wsq/sx%a(:)**2*G1(:)
!  azself(:) = wsq/sz%a(:)**2*G2(:)
!!
!  if (Vsqx == 0) then
!    axeff = 0
!  else 
!    axeff = Vsqx/sx**3
!  end if
!  if (Vsqz%a == 0) then
!    azeff = 0
!  else
!    azeff = Vsqz/sz**3
!  end if
!  dmuxdt = axself + axeff
!  dmuzdt = azself + azeff
!end subroutine

subroutine calc_dmudt(s,alpha,Vsq, wsq,dmudt)
  type(my_array), intent(in), dimension(4) :: s
  type(my_array), intent(in), dimension(2) :: Vsq
  type(my_array), intent(in) :: alpha  
  real*8, intent(in) :: wsq
  type(my_array), intent(inout), dimension(2) :: dmudt
  type(my_array), dimension(2) :: G
  integer :: i
  
  G(1) = new_my_array(alpha%N)
  G(2) = new_my_array(alpha%N)
  call geometric_parameters(alpha,G)
  do i = 1, 2, 1 
    dmudt(i)%a(:) = wsq/s(i)%a(:)**2*G(i)%a(:) + Vsq(i)%a(:)/s(i)%a(:)**3
  end do
  call free_my_array(G(1))
  call free_my_array(G(2))
end subroutine

subroutine run_adaptive_envelope(N,s,em,ti,tf,frac,constants,debug)
  integer, intent(in) :: N !This is the number of particles, not the size of the arrays.
  type(my_array), intent(inout), dimension(4) :: s
  type(my_array), intent(in), dimension(2) :: em
  real*8, intent(in) :: ti, tf, frac
  type(my_constants), intent(in) :: constants
  logical, intent(in), optional :: debug
  type(my_array), dimension(2) :: dmudt, Vsq
  real*8 :: wsq
  type(my_array) :: t, dt, alpha
  character(s(1)%N*64 + 16) :: report_str
  character(s(1)%N*16) :: current_str
  integer :: i

  if (debug) then
    write(*,*) "Inside run_envelope"
  end if
  wsq = constants%q**2*N/(20.*constants%m*sqrt(5.)*acos(-1.0)/3.)/(2.*constants%gam*constants%epsilon_0)!plasma frequency without volume
  if (debug) then
    write(*,*) wsq
    write(*,*) constants%q
    write(*,*) constants%m
  end if
  
  do i = 1, 2, 1 
    Vsq(i) = new_my_array(em(i)%N)
    Vsq(i)%a(:) = constants%c**2*em(i)%a(:)**2
    dmudt(i) = new_my_array(s(1)%N)
  end do
  if (debug) then
    write(*,*) "Initial setup."
  end if
  !Assigning additional parameters 
  t = new_my_array(s(1)%N)
  dt = new_my_array(s(1)%N)
  alpha = new_my_array(s(1)%N)
  t%a(:) = ti
  alpha%a(:) = s(2)%a(:)/s(1)%a(:)
  call calc_dmudt(s, alpha, Vsq, wsq, dmudt)
  call determine_time_step(s,dmudt,frac,dt)
  call check_time_step(t,tf,dt)
  

  if (debug) then
    write(*,*) "Before loop."
  end if
  do while ( any(t%a /= tf) )
    if (debug) then
      write(*,*) "Inside loop."
    end if
    !write(*,*) dt%a(1),t%a(1)
    t%a(:) = t%a(:) + dt%a(:)
    do i = 1, 2, 1 
      if (i == 100) then
        write(*,*) "Checking dmuzdt"
        write(*,*) dt%a(1)
        write(*,*) s(3)%a(1)
        write(*,*) s(4)%a(1)
        write(*,*) dmudt(1)%a(1)
        write(*,*) dmudt(2)%a(1)
      end if
      s(i)%a(:) = s(i)%a(:) + s(i+2)%a(:)*dt%a(:) + 0.5*dmudt(i)%a(:)*dt%a(:)**2
    end do
    if (debug) then
      write(*,*) "Before where."
    end if
    do i = 1, 2, 1
      where (s(i)%a < 0) 
        s(i)%a = s(i)%a - s(i+2)%a*dt%a - 0.5*dmudt(i)%a*dt%a**2
        s(i+2)%a = -s(i+2)%a
      else where
        s(i+2)%a = s(i+2)%a + dmudt(i)%a*dt%a
      end where
    end do
    alpha%a(:) = s(2)%a(:)/s(1)%a(:)
    call calc_dmudt(s, alpha, Vsq, wsq, dmudt)
    call determine_time_step(s,dmudt,frac,dt)
    call check_time_step(t,tf,dt)
  end do
  if (debug) then
    write(*,*) "Loop complete."
  end if
  if (debug) then
    write(*,*) "Run complete."
  end if
  call free_my_array(t)
  call free_my_array(dt)
end subroutine

subroutine run_envelope(N,s,em,every,dt,Nt,t,constants,debug)
  type(my_array), intent(inout), dimension(4) :: s
  type(my_array), intent(in), dimension(2) :: em
  integer, intent(in) :: every, Nt, N
  real*8, intent(in) :: dt
  real*8, intent(inout) :: t
  type(my_array) :: report_t
  type(my_constants), intent(in) :: constants
  logical, intent(in), optional :: debug
  type(my_array) , dimension(2):: dmudt, Vsq
  type(my_array) :: alpha
  real*8 :: wsq, hdt
  character(s(1)%N*64 + 16) :: report_str
  character(s(1)%N*16) :: current_str
  integer :: i, j

  if (debug) then
    write(*,*) "Inside run_envelope"
  end if
  wsq = constants%q**2*N/(20.*constants%m*sqrt(5.)*acos(-1.0)/3.)/(2.*constants%gam*constants%epsilon_0)!plasma frequency without volume
  if (debug) then
    write(*,*) wsq
    write(*,*) constants%q
    write(*,*) constants%m
  end if
  
  do j = 1, 2, 1
    Vsq(i) = new_my_array(em(i)%N)
    Vsq(i)%a(:) = constants%c**2*em(i)%a(:)**2
    dmudt(i) = new_my_array(s(1)%N)
  end do
  report_t = new_my_array(1)
  if (debug) then
    write(*,*) "Initial setup."
  end if

  !Assing additional parameters 
  hdt = 0.5*dt
  alpha = new_my_array(s(1)%N)
  alpha%a(:) = s(2)%a(:)/s(1)%a(:)
  call calc_dmudt(s, alpha, Vsq, wsq, dmudt)
  if ( every > 0 ) then
    report_t%a(1) = t
    call report_my_arrays(report_t,4,s)
  end if
  if (debug) then
    write(*,*) "Before loop."
  end if

  do i = 1, Nt, 1
    if (debug) then
      write(*,*) "Inside loop."
    end if
    t = t + dt
    do j = 1, 2, 1
      s(j)%a(:) = s(j)%a(:) + (s(j+2)%a(:) + dmudt(j)%a(:)*hdt)*dt
      where (s(j)%a < 0)
        s(j)%a = s(j)%a - (s(j+2)%a + dmudt(j)%a*hdt)*dt
        s(j+2)%a = -s(j+2)%a
      else where
        s(j+2)%a = s(j+2)%a + dmudt(j)%a*dt
      end where
    end do
    if ( every > 0 .AND. mod(i+1,every) == 0) then
      if (debug) then
        write(*,*) "Reporting."
      end if
      report_t%a(1) = t
      call report_my_arrays(report_t,4,s)
    end if
    if (debug) then
      write(*,*) "After report."
    end if
    alpha%a(:) = s(2)%a(:)/s(1)%a(:)
    call calc_dmudt(s, alpha, Vsq, wsq, dmudt)
    if (debug) then
      write(*,*) i
    end if
  end do
  if (debug) then
    write(*,*) "Loop complete."
  end if
  if (every == 0 .OR. mod(i+1,every) /= 0) then
    report_t%a(1) = t
    call report_my_arrays(report_t,4,s)
  end if
  if (debug) then
    write(*,*) "Run complete."
  end if
end subroutine

end module spheroidal_envelope
