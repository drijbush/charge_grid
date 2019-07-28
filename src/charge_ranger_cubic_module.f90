Module charge_ranger_cubic_module

  ! Implementation of a charge ranger assuming integration over a cube 
  ! centred on a normalised gaussian

  Use numbers_module      , Only : wp
  Use charge_ranger_module, Only : charge_ranger

  Implicit None

  Type, Extends( charge_ranger ), Public :: charge_ranger_cubic
   Contains
     Procedure, NoPass, Public :: integ_val => cubic_integ_val
  End type charge_ranger_cubic

  Private

Contains

  Pure Function cubic_integ_val( v ) Result( r )

    ! Return the value of the integration over the cube centered on
    ! the gaussian with unit alpha

    Real( wp ) :: r

    Real( wp ), Intent( In ) :: v

    r = Erf( v ) ** 3

  End Function cubic_integ_val

End Module charge_ranger_cubic_module
