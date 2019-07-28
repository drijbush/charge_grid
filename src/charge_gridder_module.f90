Module charge_gridder_module

  ! Abstract type module to cover base code for implementations
  ! of gridding charge

  Use numbers_module      , Only : wp, li
  Use charge_ranger_module, Only : charge_ranger

  Implicit None

  Type, Abstract, Public :: charge_gridder
     Class( charge_ranger ), Allocatable, Private :: q_ranger
     Integer( li )                      , Private :: count_points = 0_li
   Contains
     Procedure,                                  Public  :: set_charge_ranger
     Procedure,                                  Public  :: grid_charge
     Procedure,                                  Public  :: reset_count
     Procedure,                                  Public  :: get_count
     Procedure( cache_block ), Deferred, NoPass, Private :: get_cache_block
     Procedure( grid_block  ), Deferred,         Private :: grid_charge_in_block
     Procedure( set_parms   ), Deferred,         Private :: set_params
  End type charge_gridder

  Private

  Interface

     Pure Subroutine cache_block( n_range, n_charge, nc_block )
       Implicit None
       Integer, Dimension( 1:3 ), Intent( In    ) :: n_range
       Integer,                   Intent( In    ) :: n_charge
       Integer, Dimension( 1:3 ), Intent(   Out ) :: nc_block
     End Subroutine cache_block
     
     Subroutine grid_block( q_gridder, n_grid, grid_lb, origin, start, range, finish, &
          gvecs_dir, alpha_sq, r_charge, n_centre, q, grid )
       Import :: wp, charge_gridder
       Implicit None
       Class( charge_gridder ),                                              Intent( InOut ) :: q_gridder
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

     End Subroutine grid_block

     Subroutine set_parms( q_gridder, range, gvecs_dir, gvecs_inv, alpha_sq )
       Import :: wp, charge_gridder
       Implicit None
       Class( charge_gridder ),                                              Intent( InOut ) :: q_gridder
       Integer   , Dimension( 1:3    ),                                      Intent( In    ) :: range
       Real( wp ), Dimension( 1:, 1: ),                                      Intent( In    ) :: gvecs_dir 
       Real( wp ), Dimension( 1:, 1: ),                                      Intent( In    ) :: gvecs_inv 
       Real( wp ),                                                           Intent( In    ) :: alpha_sq     
     End Subroutine set_parms

  End Interface

Contains

  Pure Subroutine set_charge_ranger( q_gridder, q_ranger )

    ! Set the method to determine the range of the gaussian

    Class( charge_gridder ), Intent( InOut ) :: q_gridder
    Class( charge_ranger  ), Intent( In    ) :: q_ranger  

    If( Allocated( q_gridder%q_ranger ) ) Then
       Deallocate( q_gridder%q_ranger )
    End If

    Allocate( q_gridder%q_ranger, Source = q_ranger )

  End Subroutine set_charge_ranger

  Subroutine grid_charge( q_gridder, n_grid, grid_lb, origin, start, finish, gvecs_dir, gvecs_inv, &
       q_tol, alpha, r_charge, q, grid )

    ! Driver routine for gridding a set of charge distributions

    Use charge_ranger_cubic_module, Only : charge_ranger_cubic

    Class( charge_gridder ),                                              Intent( InOut ) :: q_gridder
    Integer   , Dimension( 1:3    ),                                      Intent( In    ) :: n_grid    !! The total number of points in the global grid
    Integer   , Dimension( 1:3    ),                                      Intent( In    ) :: grid_lb   !! lower bounds grid
    Integer   , Dimension( 1:3    ),                                      Intent( In    ) :: origin    !! Where (0,0,0) is on the grid
    Integer   , Dimension( 1:3    ),                                      Intent( In    ) :: start     !! start point for calculation
    Integer   , Dimension( 1:3    ),                                      Intent( In    ) :: finish    !! finish point for calculation
    Real( wp ), Dimension( 1:, 1: ),                                      Intent( In    ) :: gvecs_dir !! The direct grid vectors
    Real( wp ), Dimension( 1:, 1: ),                                      Intent( In    ) :: gvecs_inv !! The invers grid vectors
    Real( wp ),                                                           Intent( In    ) :: q_tol     !! Tolerance for ranging - how much charge ignore
    Real( wp ),                                                           Intent( In    ) :: alpha     !! Alpha ** 2 in gaussian
    Real( wp ), Dimension( 1:, 1: ),                                      Intent( In    ) :: r_charge  !! The position of the charges
    Real( wp ), Dimension( 1:      ),                                     Intent( In    ) :: q         !! The value of the charge
    Real( wp ), Dimension( grid_lb( 1 ):, grid_lb( 2 ):, grid_lb( 3 ): ), Intent(   Out ) :: grid      !! The charge gird

    Real( wp ) :: r_gauss
    Real( wp ) :: alpha_sq

    Integer, Dimension( :, : ), Allocatable :: n_centre

    Integer, Dimension( 1:3 ) :: nc_block
    Integer, Dimension( 1:3 ) :: n_range

    Integer :: i_block_3, i_block_2, i_block_1
    Integer :: i2, i3, i1
    Integer :: n_charge
    Integer :: lo3, hi3
    Integer :: lo2, hi2
    Integer :: lo1, hi1
    Integer :: i

    ! Save a few bazillion multiplies
    alpha_sq = alpha * alpha
    
    If( .Not. Allocated( q_gridder%q_ranger ) ) Then
       Allocate( charge_ranger_cubic :: q_gridder%q_ranger )
    End If

    r_gauss = q_gridder%q_ranger%find_range( alpha_sq, 1.0_wp - q_tol )
    ! Find number of points in each direction
    Do i = 1, 3
       n_range( i ) = Nint( r_gauss / Sqrt( Dot_product( gvecs_dir( :, i ), gvecs_dir( :, i ) ) ) )
    End Do
    
    ! Number of charges
    n_charge = Size( r_charge, Dim = 2 )

    ! Find centre point for each charge distribution
    n_centre = Nint( Matmul( gvecs_inv, r_charge ) )

    ! Get the cache blocking parameters
    Call q_gridder%get_cache_block( n_range, n_charge, nc_block )

    ! Set any parameter peculiar to this gridder
    Call q_gridder%set_params( n_range, gvecs_dir, gvecs_inv, alpha_sq )

    !$omp do collapse( 3 ) schedule( dynamic )
    Do i_block_3 = start( 3 ), finish( 3 ), nc_block( 3 )
       Do i_block_2 = start( 2 ), finish( 2 ), nc_block( 2 )
          Do i_block_1 = start( 1 ), finish( 1 ), nc_block( 1 )

             ! This is the OMP task. Each block is over a different part of th grid
             ! thus avoiding a reduction
             lo3 = i_block_3
             lo2 = i_block_2
             lo1 = i_block_1
             hi3 = Min( lo3 + nc_block( 3 ) - 1, finish( 3 ) )
             hi2 = Min( lo2 + nc_block( 2 ) - 1, finish( 2 ) )
             hi1 = Min( lo1 + nc_block( 1 ) - 1, finish( 1 ) )
             Do i3 = lo3, hi3
                Do i2 = lo2, hi2
                   Do i1 = lo1, hi1
                      grid( i1, i2, i3 ) = 0.0_wp
                   End Do
                End Do
             End Do
             Do i = 1, n_charge
                ! Check this charge centre has some points in the current block
                If( n_centre( 1, i ) + n_range( 1 ) >= lo1 .And. n_centre( 1, i ) - n_range( 1 ) <= hi1 ) Then
                   If( n_centre( 2, i ) + n_range( 2 ) >= lo2 .And. n_centre( 2, i ) - n_range( 1 ) <= hi2 ) Then
                      If( n_centre( 3, i ) + n_range( 3 ) >= lo3 .And. n_centre( 3, i ) - n_range( 1 ) <= hi3 ) Then
                         q_gridder%count_points = q_gridder%count_points + &
                              Product( Int( [ hi1, hi2, hi3 ] - [ lo1, lo2, lo3 ] + 1, li ) )
                         Call q_gridder%grid_charge_in_block( n_grid, &
                              grid_lb, origin, [ lo1, lo2, lo3 ], n_range, [ hi1, hi2, hi3 ], &
                              gvecs_dir, alpha_sq, &
                              r_charge( :, i ), n_centre( :, i ), q( i ), &
                              grid )
                      End If
                   End If
                End If
             End Do

          End Do
       End Do
    End Do
    !$omp end do

  End Subroutine grid_charge

  Subroutine reset_count( q_gridder )

    Class( charge_gridder ), Intent( InOut ) :: q_gridder

    q_gridder%count_points = 0_li

  End Subroutine reset_count

  Function get_count( q_gridder ) Result( count )

    Integer( li ) :: count

    Class( charge_gridder ), Intent( In ) :: q_gridder

    ! Do it this way to make sure get the right result for OMP
    ! runs
    count = 0
    !$omp critical
    count = count + q_gridder%count_points
    !$omp end critical
    !$omp barrier

  End Function get_count

End Module charge_gridder_module
