
Module vec_operators_module

  Use numbers_module, Only : wp
  
  Implicit None

  Public :: Operator( .vec. )
  Public :: Operator( .dot. )
  
  Private

  Interface Operator( .vec. )
     Procedure vec_prod
  End Interface Operator( .vec. )  
  
  Interface Operator( .dot. )
     Procedure dot_prod
  End Interface Operator( .dot. )
  
Contains
  
  Pure Function vec_prod( a, b ) Result( r )
    Implicit None
    Real( wp ), Dimension( 1:3 ) :: r
    Real( wp ), Dimension( 1:3 ), Intent( In ) :: a
    Real( wp ), Dimension( 1:3 ), Intent( In ) :: b
    r( 1 ) =   ( a( 2 ) * b( 3 ) - b( 2 ) * a( 3 ) )
    r( 2 ) = - ( a( 1 ) * b( 3 ) - b( 1 ) * a( 3 ) )
    r( 3 ) =   ( a( 1 ) * b( 2 ) - b( 1 ) * a( 2 ) )
  End Function vec_prod

  Pure Function dot_prod( a, b ) Result( r )
    Implicit None
    Real( wp ) :: r
    Real( wp ), Dimension( 1: ), Intent( In ) :: a
    Real( wp ), Dimension( 1: ), Intent( In ) :: b
    r = Dot_Product( a, b )
  End Function dot_prod

End Module vec_operators_module
