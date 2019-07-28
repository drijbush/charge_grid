Module charge_ranger_radial_module

  ! Implementation of a charge ranger assuming spherical integration about 
  ! the centre of a normalised gaussian

  Use numbers_module      , Only : wp
  Use charge_ranger_module, Only : charge_ranger

  Implicit None

  Type, Extends( charge_ranger ), Public :: charge_ranger_radial
   Contains
     Procedure, NoPass, Public :: integ_val => radial_integ_val
  End type charge_ranger_radial

  Private

Contains

  Pure Function radial_integ_val( v ) Result( r )

    ! This returns 4 * pi * N * integral(0,v) dr r^2 g(-r^2)
    ! Note alpha is 1

    Real( wp ) :: r

    Real( wp ), Intent( In ) :: v
    Real( wp ), Parameter :: sqrt_pi = Sqrt( 3.1415926535897932384626_wp )

    r = Erf( v ) - 2.0_wp * v * Exp( - v * v ) / sqrt_pi

  End Function radial_integ_val

End Module charge_ranger_radial_module
