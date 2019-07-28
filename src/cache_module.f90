Module cache_module

  Implicit None

  ! Assumes a 256 KByte Cache
  Integer, Parameter, Private :: cache_l2_size_default =  2**18 / 2 ** 3 

  ! Number of 64 bit reals can fit in cache
  Integer, Protected, Public  :: cache_l2_size = cache_l2_size_default

  Public :: cache_set_size
  Public :: cache_determine_size

  Private

Contains

  Subroutine cache_set_size( n_reals )

    Integer, Intent( In ) :: n_reals

    cache_l2_size = n_reals

  End Subroutine cache_set_size

  Subroutine cache_determine_size

    ! Need to work out how to call sysconf or similar
    ! to find the size of (the largest?) non-shared cache
    ! At Moment (2019) usualy L2
    Stop "Not implemented"

  End Subroutine cache_determine_size

End Module cache_module
