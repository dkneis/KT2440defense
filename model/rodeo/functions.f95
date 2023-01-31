module functions
implicit none
contains

! Monod model with an offset parameter
real(8) function monod(x_current, x_half, x_zero)
  real(8), intent(in):: x_current, x_half, x_zero
  real(8), parameter:: ZERO=0.0_8, ONE=1.0_8, TWO=2.0_8
  real(8), parameter:: SMALL_FRACTION=0.001_8
  real(8):: pwr, x
  if (x_zero .ge. x_half) then
    call rexit("condition 'x_zero < x_half' not met")
  end if
  ! Handle offset
  x = x_current - x_zero
  if (x .lt. ZERO) then
    monod = ZERO
  else
    ! If x << h, we introduce a power term to avoid negative concentrations
    ! in the ODE solving process. The value of the power term continuously
    ! increases to a maximum value of 2 as x -> 0.
    if (x .lt. SMALL_FRACTION * x_half) then
      pwr = -ONE/(SMALL_FRACTION * x_half) * x + TWO
    else
      pwr = ONE
    end if
    monod = x**pwr / (x**pwr + x_half**pwr)
  end if
end function

! MIC model with parabola interpolation close to the MIC
! The parameters L and R define the range around the MIC
! where the parabola is used (as multipliers of MIC).
! Note: The parabola is defined by four boundary counditions,
!   two values and two derivatives. It appears that some of
!   the conditions are fulfilled magically, even though the
!   parabola has only a single free parameter.
real(8) function inh(x, mic)
  real(8), intent(in):: x, mic
  real(8), parameter:: small=0.01_8, &
    ZERO=0.0_8, ONE=1.0_8, TWO=2.0_8
  real(8):: a, L, H
  L = ONE - small
  H = ONE + small
  a = -ONE / (TWO * (L*mic**TWO - H*mic**TWO))
  if (x .le. L*mic) then
    inh = ONE - x/mic
  else if (x .ge. H*mic) then
    inh = ZERO
  else
    inh = a*(x-H*mic)**TWO
  end if
end function

!Sigmoidal functions
real(8) function jumpUp (x, lower, upper)
  real(8), intent(in):: x, lower, upper
  real(8):: mid
  real(8), parameter:: ZERO=0.0_8, ONE=1.0_8, TWO=2.0_8, HALF=0.5_8
  mid = (lower + upper) / TWO
  if ((x .gt. lower) .and. (x .le. mid)) then
    jumpUp = (x - lower)**TWO / (mid - lower)**TWO * HALF
  else if ((x .lt. upper) .and. (x .ge. mid)) then
    jumpUp = ONE - (upper - x)**TWO / (upper - mid)**TWO * HALF
  else if (x <= lower) then
    jumpUp = ZERO
  else
    jumpUp = ONE
  end if
end function

real(8) function jumpDown (x, lower, upper)
  real(8), intent(in):: x, lower, upper
  real(8), parameter:: ONE=1.0_8
  jumpDown = -ONE * jumpUp(x, lower, upper) + ONE
end function

! Simplified versions of the sigmoidal functions with fixed width of interval
real(8) function on (x, t)
  real(8), intent(in):: x, t
  real(8), parameter:: LOWER=0.9_8, UPPER=1.1_8
  on = jumpUp(x, t*LOWER, t*UPPER)
end function

real(8) function off (x, t)
  real(8), intent(in):: x, t
  real(8), parameter:: LOWER=0.9_8, UPPER=1.1_8
  off = jumpDown(x, t*LOWER, t*UPPER)
end function

end module
