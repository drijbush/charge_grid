Module charge_ranger_module

  ! Sets up abstact data type to do with finding range of a gaussian. Allows
  ! multiple implementations sharing much of the same code base depending on
  ! how it is decided to cut off the tails of the gaussian.

  Use numbers_module, Only : wp

  Implicit None

  Type, Abstract, Public :: charge_ranger
     Contains
       Procedure                                  , Public  :: find_range
       ! The function determining where to cut the gaussian off. Note
       ! this should be implemented for alpha = 1 in Exp(-alpha*r^2)
       Procedure( charge_integ ), NoPass, Deferred, Private :: integ_val 
  End type charge_ranger

  Private

  Interface
     Pure Function charge_integ( v ) Result( r )
       Import :: wp
       Implicit None
       Real( wp ) :: r
       Real( wp ), Intent( In ) :: v
     End Function charge_integ
  End Interface

Contains

  Pure Function find_range( q_ranger, alpha_sq, tol ) Result( range )

    ! Find the range of the gaussian Exp ( - alpha_sq * r * r ). This is done
    ! by solving the equation F(Exp ( - alpha_sq * r * r ))-tol = 0. Typically the function
    ! is to integrate the gaussian out to a certain distance, so showing how much
    ! charge would be "lost" if the tails is ignored

    ! Simple bisection is used. The integral of a gaussian is essentially (the precise details
    ! depend on how exactly the bounds are specified) an error function. As we want liitle charge left
    ! we are very close the limiting value of erf, 1. Here the gradien is essentially constant and
    ! zero, so we gain little by using non-linear equation solvers with gradient information, so KISS
    ! and use bisection - it converges very quickly for any sensible value of tol.

    Real( wp ) :: range

    Real( wp ), Parameter :: scale = 1.5_wp

    Class( charge_ranger ), Intent( In ) :: q_ranger
    Real( wp )            , Intent( In ) :: alpha_sq
    Real( wp )            , Intent( In ) :: tol

    Real( wp ) :: hi, lo
    Real( wp ) :: mid
    Real( wp ) :: val

    ! Bound solution
    lo = 0.1_wp
    Do While( q_ranger%integ_val( lo ) - tol < 0.0_wp )
       lo = lo * scale
    End Do
    hi = lo
    lo = lo / scale

    ! Solve by bisection
    val = -1.0_wp
    Do While( hi - lo > 1e-5_wp .Or. Abs( val ) > 1e-15_wp )
       mid = 0.5_wp * ( hi + lo )
       ! For the current implementation (depening on the dynamic type of q_ranger)
       ! get the current value of the function
       val =  q_ranger%integ_val( mid ) - tol
       If( val > 0.0_wp ) Then
          hi = mid
       Else
          lo = mid
       End If
    End Do
    mid = 0.5_wp * ( hi + lo )

    ! Adjust for actual value of alpha
    range = mid / Sqrt( alpha_sq )

  End Function find_range

End Module charge_ranger_module
