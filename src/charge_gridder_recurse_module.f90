Module charge_gridder_recurse_module

  ! Type to cover gridding of charge by recursing from one point, i.e. writing
  ! Exp( -alpha*(r-dr).(r-dr) ) = Exp( - alpha r.r ) * Exp( + 2.0 * alpha * r.dr ) * Exp( -alpha dr.dr )
  !
  ! Also need for next point
  ! Exp( + 2.0 * alpha * ( r - dr ).dr ) =  Exp( + 2.0 * alpha * r.dr ) * Exp( - 2.0_wp * alpha * dr.dr )

  Use numbers_module       , Only : wp
  Use charge_ranger_module , Only : charge_ranger
  Use charge_gridder_module, Only : charge_gridder

  Implicit None

  ! Unrolling factors for loops in charge gridding - only innermost (i.e. 1) currently supported
!!$    Integer, Dimension( 1:3 ), Parameter :: n_unroll = [ 2, 1, 1 ]
!!$    Integer, Dimension( 1:3 ), Parameter :: n_unroll = [ 4, 1, 1 ]
    Integer, Dimension( 1:3 ), Parameter :: n_unroll = [ 4, 2, 1 ]

  Type, Extends(  charge_gridder ), Public :: charge_gridder_recurse
     ! Array to hold the terms of the form Exp(-alpha*n*dr.dr) that we need
     Real( wp ), Dimension( 0:2 * ( Maxval( n_unroll ) ** 2 ), 1:3, 1:3 ) :: Exp_dd 
   Contains
     Procedure, NoPass, Private :: get_cache_block      => cache_block_recurse
     Procedure,         Private :: grid_charge_in_block => grid_block_recurse
     Procedure,         Private :: set_params           => set_params_recurse
  End type charge_gridder_recurse

  Private

Contains

  Pure Subroutine cache_block_recurse( n_range, n_charge, nc_block )

    Use cache_module, Only : cache_l2_size

    Implicit None

    Integer, Dimension( 1:3 )      , Intent( In    ) :: n_range
    Integer,                         Intent( In    ) :: n_charge
    Integer, Dimension( 1:3  )     , Intent(   Out ) :: nc_block

    Integer :: i

    nc_block = 0
    i = 0
    Do 
       nc_block( i + 1 ) = nc_block( i + 1 ) + 1
       If( Product( nc_block ) + 4 * n_charge + 0 * Product( n_range ) >= cache_l2_size ) Exit
       i = Mod( i + 1, 3 )
    End Do
    nc_block( i + 1 ) = nc_block( i + 1 ) - 1

  End Subroutine cache_block_recurse
     
  Subroutine grid_block_recurse( q_gridder, n_grid, grid_lb, origin, start, range, finish, &
       gvecs_dir, alpha_sq, r_charge, n_centre, q, grid )

    Use numbers_module, Only : wp, pi

    Implicit None

    Class( charge_gridder_recurse ),                                      Intent( InOut ) :: q_gridder
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

    Real( wp ), Dimension( 0:n_unroll( 1 ) - 1, 0:n_unroll( 2 ) - 1 ) :: Exp_rr_1, Exp_rd_1

    Real( wp ) :: qnorm
    Real( wp ) :: Exp_1dd_1, Exp_2dd_1
    Real( wp ) :: Exp_update_rr_1, Exp_update_rd_1
    Real( wp ) :: Exp_2dd_12, Exp_2dd_23, Exp_2dd_13
    Real( wp ) :: Exp_1dd_2, Exp_2dd_2
    Real( wp ) :: Exp0_rd_1_start, Exp0_rd_12_start
    Real( wp ) :: Exp0_rd_2_start
    Real( wp ) :: Exp0_rr_2, Exp0_rd_2
    Real( wp ) :: Exp0_rr_3, Exp0_rd_3
    Real( wp ) :: Exp_1dd_3, Exp_2dd_3
    Real( wp ) :: e12

    Integer, Dimension( 1:3 ) :: ilo, ihi, n_left, n_top

    Integer :: i1, i2, i3
    Integer :: i1_unroll, i2_unroll

    !!!! NEED TO WORK THROUGH ORIGIN PROPERLY !!!!

    Exp_1dd_1 = q_gridder%Exp_dd( 1, 1, 1 )
    Exp_2dd_1 = q_gridder%Exp_dd( 2, 1, 1 )
    Exp_update_rr_1 = q_gridder%Exp_dd(     n_unroll( 1 ) ** 2, 1, 1 )
    Exp_update_rd_1 = q_gridder%Exp_dd( 2 * n_unroll( 1 ) ** 2, 1, 1 )

    Exp_1dd_2 = q_gridder%Exp_dd( 1, 2, 2 )
    Exp_2dd_2 = q_gridder%Exp_dd( 2, 2, 2 )

    Exp_1dd_3 = q_gridder%Exp_dd( 1, 3, 3 )
    Exp_2dd_3 = q_gridder%Exp_dd( 2, 3, 3 )

    Exp_2dd_12 = q_gridder%Exp_dd( 2, 1, 2 )
    Exp_2dd_23 = q_gridder%Exp_dd( 2, 2, 3 )
    Exp_2dd_13 = q_gridder%Exp_dd( 2, 1, 3 )

    qnorm = q * alpha_sq * Sqrt( alpha_sq ) / ( pi ** 1.5_wp )

    ilo = Max( start , n_centre - range )
    ihi = Min( finish, n_centre + range )
    n_left = Mod( ihi - ilo + 1, n_unroll )
    n_left = Merge( n_unroll, n_left, n_left == 0 )
    n_top = ihi - n_left 

    r3 = ilo( 3 ) * gvecs_dir( :, 3 )
    r32 = r3 + ilo( 2 ) * gvecs_dir( :, 2 )
    r321 = r32 + ilo( 1 ) * gvecs_dir( :, 1 )
    r = r_charge - r321
    Exp0_rr_3       = qnorm * Exp( -          alpha_sq * Dot_product( r, r ) )
    Exp0_rd_3       =         Exp(   2.0_wp * alpha_sq * Dot_product( r, gvecs_dir( :, 3 ) ) )
    Exp0_rd_2_start =         Exp(   2.0_wp * alpha_sq * Dot_product( r, gvecs_dir( :, 2 ) ) )
    Exp0_rd_1_start =         Exp(   2.0_wp * alpha_sq * Dot_product( r, gvecs_dir( :, 1 ) ) )
    Do i3 = ilo( 3 ), ihi( 3 )
       Exp0_rr_2        = Exp0_rr_3
       Exp0_rd_2        = Exp0_rd_2_start
       Exp0_rd_12_start = Exp0_rd_1_start
!!$       Do i2 = ilo( 2 ), ihi( 2 )
       Do i2 = ilo( 2 ), n_top( 2 ), n_unroll( 2 )
          Exp_rr_1( 0, 0 ) = Exp0_rr_2
          Exp_rd_1( 0, 0 ) = Exp0_rd_12_start
          ! Calculate The subsequent unrolling loops by first stepping one x at a time
          Do i1_unroll = 1, n_unroll( 1 ) - 1
             Exp_rr_1( i1_unroll, 0 ) = Exp_rr_1( i1_unroll - 1, 0 ) * Exp_rd_1( i1_unroll - 1, 0 ) * Exp_1dd_1
             Exp_rd_1( i1_unroll, 0 ) = Exp_rd_1( i1_unroll - 1, 0 ) * Exp_2dd_1
          End Do
          Do i2_unroll = 1, n_unroll( 2 ) - 1
             e12 = Exp0_rd_2
             Do i1_unroll = 0, n_unroll( 1 ) - 1
                Exp_rr_1( i1_unroll, i2_unroll ) = &
                     Exp_rr_1( i1_unroll, i2_unroll - 1 ) * e12 * Exp_1dd_2
                Exp_rd_1( i1_unroll, i2_unroll ) = Exp_rd_1( i1_unroll, i2_unroll - 1 ) * Exp_2dd_2
                e12 = e12 * Exp_2dd_12
             End Do
             ! Now step each element 1 y at a time filling in the rest of the elements
             ! Update the 2nd vector quatites while the power above is being evaluated
             Exp0_rr_2 = Exp0_rr_2 * Exp0_rd_2 * Exp_1dd_2
             Exp0_rd_2 = Exp0_rd_2 * Exp_2dd_2
             Exp0_rd_12_start = Exp0_rd_12_start * Exp_2dd_12
          End Do
          ! Update the 2nd vector quatites while the power above is being evaluated
          Exp0_rr_2 = Exp0_rr_2 * Exp0_rd_2 * Exp_1dd_2
          Exp0_rd_2 = Exp0_rd_2 * Exp_2dd_2
          Exp0_rd_12_start = Exp0_rd_12_start * Exp_2dd_12
          ! Now update the rd steps to move the appropriate number of unrolling steps
          ! at a time
          Exp_rd_1 = Exp_rd_1 ** n_unroll( 1 )
          Do i1 = ilo( 1 ), n_top( 1 ), n_unroll( 1 )
             grid( i1:i1 + n_unroll( 1 ) - 1, i2:i2 + n_unroll( 2 ) - 1, i3 ) = &
                  grid( i1:i1 + n_unroll( 1 ) - 1, i2:i2 + n_unroll( 2 ) - 1, i3 ) + Exp_rr_1
             Exp_rr_1 = Exp_rr_1 * Exp_rd_1 * Exp_update_rr_1
             Exp_rd_1 = Exp_rd_1 * Exp_update_rd_1
          End Do
          grid( n_top( 1 ) + 1:ihi( 1 ), i2:i2 + n_unroll( 2 ) - 1, i3 ) = &
               grid( n_top( 1 ) + 1:ihi( 1 ), i2:i2 + n_unroll( 2 ) - 1, i3 ) + &
               Exp_rr_1( 0:n_left( 1 ) - 1, : )
       End Do
       Do i2 = n_top( 2 ) + 1, ihi( 2 ), n_unroll( 2 )
          Exp_rr_1( 0, 0 ) = Exp0_rr_2
          Exp_rd_1( 0, 0 ) = Exp0_rd_12_start
          ! Calculate The subsequent unrolling loops by first stepping one x at a time
          Do i1_unroll = 1, n_unroll( 1 ) - 1
             Exp_rr_1( i1_unroll, 0 ) = Exp_rr_1( i1_unroll - 1, 0 ) * Exp_rd_1( i1_unroll - 1, 0 ) * Exp_1dd_1
             Exp_rd_1( i1_unroll, 0 ) = Exp_rd_1( i1_unroll - 1, 0 ) * Exp_2dd_1
          End Do
          Do i2_unroll = 1, n_left( 2 ) - 1
             e12 = Exp0_rd_2
             Do i1_unroll = 0, n_unroll( 1 ) - 1
                Exp_rr_1( i1_unroll, i2_unroll ) = &
                     Exp_rr_1( i1_unroll, i2_unroll - 1 ) * e12 * Exp_1dd_2
                Exp_rd_1( i1_unroll, i2_unroll ) = Exp_rd_1( i1_unroll, i2_unroll - 1 ) * Exp_2dd_2
                e12 = e12 * Exp_2dd_12
             End Do
             ! Update the 2nd vector quatites while the power above is being evaluated
             Exp0_rr_2 = Exp0_rr_2 * Exp0_rd_2 * Exp_1dd_2
             Exp0_rd_2 = Exp0_rd_2 * Exp_2dd_2
             Exp0_rd_12_start = Exp0_rd_12_start * Exp_2dd_12
          End Do
          ! Update the 2nd vector quatites while the power above is being evaluated
          Exp0_rr_2 = Exp0_rr_2 * Exp0_rd_2 * Exp_1dd_2
          Exp0_rd_2 = Exp0_rd_2 * Exp_2dd_2
          Exp0_rd_12_start = Exp0_rd_12_start * Exp_2dd_12
          ! Now update the rd steps to move the appropriate number of unrolling steps
          ! at a time
          Exp_rd_1 = Exp_rd_1 ** n_unroll( 1 )
          Do i1 = ilo( 1 ), n_top( 1 ), n_unroll( 1 )
             grid( i1:i1 + n_unroll( 1 ) - 1, n_top( 2 ) + 1:ihi( 2 ), i3 ) = &
                  grid( i1:i1 + n_unroll( 1 ) - 1, n_top( 2 ) + 1:ihi( 2 ), i3 ) + &
                  Exp_rr_1( :, 0:n_left( 2 ) - 1 )
             Exp_rr_1( :, 0:n_left( 2 ) - 1 ) = Exp_rr_1( :, 0:n_left( 2 ) - 1 ) * &
                  Exp_rd_1( :, 0:n_left( 2 ) - 1 ) * Exp_update_rr_1
             Exp_rd_1( :, 0:n_left( 2 ) - 1 ) = Exp_rd_1( :, 0:n_left( 2 ) - 1 ) * Exp_update_rd_1
          End Do
       End Do
       grid( n_top( 1 ) + 1:ihi( 1 ), n_top( 2 ) + 1:ihi( 2 ), i3 ) = &
            grid( n_top( 1 ) + 1:ihi( 1 ), n_top( 2 ) + 1:ihi( 2 ), i3 ) + &
            Exp_rr_1( 0:n_left( 1 ) - 1, 0:n_left( 2 ) - 1 )
       Exp0_rr_3 = Exp0_rr_3 * Exp0_rd_3 * Exp_1dd_3
       Exp0_rd_3 = Exp0_rd_3 * Exp_2dd_3
       Exp0_rd_2_start = Exp0_rd_2_start * Exp_2dd_23
       Exp0_rd_1_start = Exp0_rd_1_start * Exp_2dd_13
    End Do

  End Subroutine grid_block_recurse

  Subroutine set_params_recurse( q_gridder, range, gvecs_dir, gvecs_inv, alpha_sq )

    Use numbers_module, Only : wp

    Implicit None

    Class( charge_gridder_recurse ), Intent( InOut ) :: q_gridder 
    Integer   , Dimension( 1:3    ), Intent( In    ) :: range
    Real( wp ), Dimension( 1:, 1: ), Intent( In    ) :: gvecs_dir 
    Real( wp ), Dimension( 1:, 1: ), Intent( In    ) :: gvecs_inv 
    Real( wp ),                      Intent( In    ) :: alpha_sq     

    Integer :: i, idir1, idir2

    Do idir1 = 1, 3
       Do idir2 = idir1, 3
          Do i = Lbound( q_gridder%Exp_dd, Dim = 1 ), Ubound( q_gridder%Exp_dd, Dim = 1 )
             q_gridder%Exp_dd( i, idir2, idir1 ) = &
                  Exp( - i * alpha_sq * Dot_product( gvecs_dir( :, idir2 ), gvecs_dir( :, idir1 ) ) )
             q_gridder%Exp_dd( i, idir1, idir2 ) = q_gridder%Exp_dd( i, idir2, idir1 )
          End Do
       End Do
    End Do

    Return

    ! Shut the compiler up about unused arguments
    Write( *, * ) range, gvecs_dir, gvecs_inv, alpha_sq

  End Subroutine set_params_recurse

End Module charge_gridder_recurse_module
