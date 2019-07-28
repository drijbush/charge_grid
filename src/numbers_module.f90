Module numbers_module

  Implicit None

  Integer, Parameter, Public :: wp = Selected_real_kind( 12, 70 )
  Integer, Parameter, Public :: li = Selected_int_kind( 18 )

  Real( wp ), Parameter, Public :: pi = 3.14159265358979323846264338327950288419716939937510_wp

  Private

End Module numbers_module
