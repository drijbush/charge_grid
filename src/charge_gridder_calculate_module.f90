Module charge_gridder_calculate_module

  ! Abstract type module to cover base code for implementations
  ! of gridding charge

  Use charge_ranger_module , Only : charge_ranger
  Use charge_gridder_module, Only : charge_gridder

  Implicit None

  Type, Extends(  charge_gridder ), Public :: charge_gridder_calculate
   Contains
     Procedure, NoPass, Private :: get_cache_block      => cache_block_calculate
     Procedure,         Private :: grid_charge_in_block => grid_block_calculate
     Procedure,         Private :: set_params           => set_params_calculate
  End type charge_gridder_calculate

  Private

Contains

  Pure Subroutine cache_block_calculate( n_range, n_charge, nc_block )

    Use cache_module, Only : cache_l2_size

    Implicit None

    Integer, Dimension( 1:3 )        , Intent( In    ) :: n_range
    Integer,                           Intent( In    ) :: n_charge
    Integer, Dimension( 1:3  )       , Intent(   Out ) :: nc_block

    Integer :: i

    nc_block = 0
    i = 0
    Do 
       nc_block( i + 1 ) = nc_block( i + 1 ) + 1
       If( Product( nc_block ) + 4 * n_charge + 0 * Product( n_range ) >= cache_l2_size ) Exit
       i = Mod( i + 1, 3 )
    End Do
    nc_block( i + 1 ) = nc_block( i + 1 ) - 1

  End Subroutine cache_block_calculate
     
  Subroutine grid_block_calculate( q_gridder, n_grid, grid_lb, origin, start, range, finish, &
       gvecs_dir, alpha_sq, r_charge, n_centre, q, grid )

    Use numbers_module, Only : wp, pi

    Implicit None

    Class( charge_gridder_calculate ),                                    Intent( InOut ) :: q_gridder
    Integer   , Dimension( 1:3    ),                                      Intent( In    ) :: n_grid    
    Integer   , Dimension( 1:3    ),                                      Intent( In    ) :: grid_lb   
    Integer   , Dimension( 1:3    ),                                      Intent( In    ) :: origin    
    Integer   , Dimension( 1:3    ),                                      Intent( In    ) :: start     
    Integer   , Dimension( 1:3    ),                                      Intent( In    ) :: range
    Integer   , Dimension( 1:3    ),                                      Intent( In    ) :: finish     
    Real( wp ), Dimension( 1:, 1: ),                                      Intent( In    ) :: gvecs_dir 
    Real( wp ),                                                           Intent( In    ) :: alpha_sq     
    Real( wp ), Dimension( 1:3    ),                                      Intent( In    ) :: r_charge  
    Integer   , Dimension( 1:3    ),                                      Intent( In    ) :: n_centre
    Real( wp ),                                                           Intent( In    ) :: q         
    Real( wp ), Dimension( grid_lb( 1 ):, grid_lb( 2 ):, grid_lb( 3 ): ), Intent( InOut ) :: grid      

    Real( wp ), Dimension( 1:3 ) :: r3, r32, r321, r

    Real( wp ) :: qnorm

    Integer, Dimension( 1:3 ) :: ilo, ihi

    Integer :: i1, i2, i3

    !!!! NEED TO WORK THROUGH ORIGIN PROPERLY !!!!

    qnorm = q * alpha_sq * Sqrt( alpha_sq ) / ( pi ** 1.5_wp )

    ilo = Max( start , n_centre - range )
    ihi = Min( finish, n_centre + range )

    Do i3 = ilo( 3 ), ihi( 3 )
       r3 = i3 * gvecs_dir( :, 3 )
       Do i2 = ilo( 2 ), ihi( 2 )
          r32 = r3 + i2 * gvecs_dir( :, 2 )
          Do i1 = ilo( 1 ), ihi( 1 )
             r321 = r32 + i1 * gvecs_dir( :, 1 )
             r = r_charge - r321
             grid( i1, i2, i3 ) = grid( i1, i2, i3 ) + qnorm * Exp( - alpha_sq * Dot_product( r, r ) )
          End Do
       End Do
    End Do

  End Subroutine grid_block_calculate

  Subroutine set_params_calculate( q_gridder, range, gvecs_dir, gvecs_inv, alpha_sq )

    Use numbers_module, Only : wp

    Implicit None

    Class( charge_gridder_calculate ), Intent( InOut ) :: q_gridder 
    Integer   , Dimension( 1:3    )  , Intent( In    ) :: range
    Real( wp ), Dimension( 1:, 1: )  , Intent( In    ) :: gvecs_dir 
    Real( wp ), Dimension( 1:, 1: )  , Intent( In    ) :: gvecs_inv 
    Real( wp ),                        Intent( In    ) :: alpha_sq     

    ! No extra Params to set for this emthod

    Return

    ! Shut the compiler up about unused arguments
    Write( *, * ) range, gvecs_dir, gvecs_inv, alpha_sq

  End Subroutine set_params_calculate

End Module charge_gridder_calculate_module
