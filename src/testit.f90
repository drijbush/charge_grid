
Program testit

  Use numbers_module, Only : wp, pi
  Use charge_gridder_module, Only : charge_gridder
  Use charge_gridder_calculate_module, Only : charge_gridder_calculate
  Use charge_gridder_recurse_module, Only : charge_gridder_recurse

  Implicit None

  Integer, Dimension( 1:3 ), Parameter :: ng = [ 105, 107, 106 ]
!!$  Integer, Dimension( 1:3 ), Parameter :: ng = [ 205, 207, 206 ]
!!$  Integer, Dimension( 1:3 ), Parameter :: ng = [ 5, 7, 6 ]
  Integer, Dimension( 1:3 ), Parameter :: proc_coords = [ 1, 3, 2 ]
  Integer, Dimension( 1:3 ), Parameter :: gs = ng * proc_coords
  Integer, Dimension( 1:3 ), Parameter :: gf = gs + ng - 1
  Integer, Dimension( 1:3 ), Parameter :: origin = [ 0, 0, 0 ]

  Real( wp ), Parameter :: qtol = 1e-14_wp

  Class( charge_gridder ), Allocatable :: qgridder

  Real( wp ), Dimension( gs( 1 ):gf( 1 ), gs( 2 ):gf( 2 ), gs( 3 ):gf( 3 ) ) :: grid_ref  
  Real( wp ), Dimension( gs( 1 ):gf( 1 ), gs( 2 ):gf( 2 ), gs( 3 ):gf( 3 ) ) :: grid
  Real( wp ), Dimension( gs( 1 ):gf( 1 ), gs( 2 ):gf( 2 ), gs( 3 ):gf( 3 ) ) :: grid_calc
  Real( wp ), Dimension( gs( 1 ):gf( 1 ), gs( 2 ):gf( 2 ), gs( 3 ):gf( 3 ) ) :: grid_recurse

  Real( wp ), Dimension( 1:3, 1:3 ) :: gvecs_dir
  Real( wp ), Dimension( :, : ), Allocatable :: gvecs_inv

  Integer, Parameter :: nq = 1000
!!$  Integer, Parameter :: nq = 10
  Real( wp ), Dimension( 1:3, 1:nq ) :: rq
  Real( wp ), Dimension(      1:nq ) ::  q

  Real( wp ), Dimension( 1:3 ) :: r_frac
  Real( wp ), Dimension( 1:3 ) :: r

  Real( wp ) :: V
  Real( wp ) :: alpha = 0.6
  Real( wp ) :: norm
  Real( wp ) :: time

  Integer :: i1, i2, i3
  Integer :: i

  Integer :: start, finish, rate

  V = 0.0_wp
  Do While( Abs( V ) < 0.1_wp * 0.1_wp * 0.1_wp )
     Call Random_number( gvecs_dir )
     gvecs_dir = gvecs_dir - 0.5_wp
     gvecs_dir = gvecs_dir / 25.0_wp
     Do i = 1, 3
        gvecs_dir( i, i ) = gvecs_dir( i, i ) + 0.2_wp
     End Do
     Call givens_invert( gvecs_dir, gvecs_inv, V )
  End Do

  Write( *, '( "Volume of grid cell = ", f8.4 )' ) Abs( V )
  Write( *, '( "Direct  grid vectors", /, 3( 3( f8.4, 1x ) / ) )' ) Transpose( gvecs_dir )
  Write( *, '( "Inverse grid vectors", /, 3( 3( f8.4, 1x ) / ) )' ) Transpose( gvecs_inv )

  Do i = 1, nq
     Call random_number( r_frac )
     ! make possible for charge to be just outside the area of the grid so as to
     ! simulate a halo
     r_frac = 0.95_wp * ( gs + 1.1_wp * r_frac * ng )
     !!CHECK THIS!!
     r_frac = r_frac - origin
     rq( :, i ) = r_frac( 1 ) * gvecs_dir( :, 1 ) + &
          r_frac( 2 ) * gvecs_dir( :, 2 ) + &
          r_frac( 3 ) * gvecs_dir( :, 3 )
     q( i ) = Mod( i, 2 )
  End Do

  Write( *, * ) 'Done set pos'

  norm = ( alpha ** 3 ) / ( pi ** 1.5_wp )
  Call system_clock( start, rate )
  grid = 0.0_wp
  Do i = 1, nq
     Do i3 = gs( 3 ), gf( 3 )
        Do i2 = gs( 2 ), gf( 2 )
           Do i1 = gs( 1 ), gf( 1 )
              r = i1 * gvecs_dir( :, 1 ) + i2 * gvecs_dir( :, 2 ) + i3 * gvecs_dir( :, 3 )
              grid( i1, i2, i3 ) = grid( i1, i2, i3 ) + norm * q( i ) * &
                   Exp( - alpha * alpha * Dot_product( r - rq( :, i ), r - rq( :, i ) ) )
           End Do
        End Do
     End Do
  End Do
  Call system_clock( finish, rate )
  time = Real( finish - start ) / rate
  Write( *, '( a, t15, f8.4, 1x, f8.6, 1x, g13.4 )' ) 'naive    times ', time, time / nq, time / Size( grid ) / nq

  grid_ref = grid
  Write( 10 ) grid_ref

  Allocate( charge_gridder_calculate :: qgridder )

  Call system_clock( start, rate )
  !$omp parallel
  Call qgridder%grid_charge( ng, Lbound( grid ), origin, Lbound( grid ), Ubound( grid ), &
       gvecs_dir, gvecs_inv, qtol, alpha, rq, q, grid )
  Call system_clock( finish, rate )
  !$omp end parallel
  time = Real( finish - start ) / rate
  Write( *, '( a, t15, f8.4, 1x, f8.6, 1x, g13.4 )' ) 'calc     times ', time, time / nq, time / Size( grid ) / nq

  grid_calc = grid
  Write( 11 ) grid_calc

!!$  Do i3 = gs( 3 ), gf( 3 )
!!$     Do i2 = gs( 2 ), gf( 2 )
!!$        Do i1 = gs( 1 ), gf( 1 )
!!$           Write( *, * ) i1, i2, i3, grid_ref( i1, i2, i3 ), grid_calc( i1, i2, i3 ), &
!!$                grid_ref( i1, i2, i3 ) - grid_calc( i1, i2, i3 )
!!$        End Do
!!$     End Do
!!$  End Do

  Deallocate( qgridder )
  Allocate( charge_gridder_recurse :: qgridder )

  Call system_clock( start, rate )
  Call qgridder%reset_count
  !$omp parallel
  Call qgridder%grid_charge( ng, Lbound( grid ), origin, Lbound( grid ), Ubound( grid ), &
       gvecs_dir, gvecs_inv, qtol, alpha, rq, q, grid )
  Call system_clock( finish, rate )
  !$omp end parallel
  time = Real( finish - start ) / rate
  Write( *, '( a, t15, f8.4, 1x, f8.6, 1x, g13.4 )' ) 'recurse  times', time, time / nq, time / Size( grid ) / nq
  Write( *, '( "~Gflops for recursive method: ", 1x, f8.4 )' ) 4.0_wp * qgridder%get_count() / time / 1.0e9_wp

  grid_recurse = grid
  Write( 12 ) grid_recurse

  Write( *, * ) 'Max error for calculate ', Maxval( Abs( grid_calc    - grid_ref ) ) 
  Write( *, * ) 'Max error for recurse   ', Maxval( Abs( grid_recurse - grid_ref ) ) 

!!$  Do i3 = gs( 3 ), gf( 3 )
!!$     Do i2 = gs( 2 ), gf( 2 )
!!$        Do i1 = gs( 1 ), gf( 1 )
!!$           Write( *, * ) i1, i2, i3, grid_ref( i1, i2, i3 ), grid_recurse( i1, i2, i3 ), &
!!$                grid_ref( i1, i2, i3 ) - grid_recurse( i1, i2, i3 )
!!$        End Do
!!$     End Do
!!$  End Do

Contains

 Subroutine givens_invert( A, B, det )
    !!----------------------------------------------------
    !! Invert the matrix A, returning the inverse in B
    !! and the determinant of A in det
    !!
    !! NOT FOR USE ON LARGE MATRICES
    !!
    !! Really designed as a robust, not totally disastorous, way
    !! to find inverse lattice vectors and the volume of the unit cell
    !! ( i.e. the determinant )
    !!
    !! The method used is to QR factorise ( in fact LQ factorise ),
    !! invert the triangular matrix and the from that form the
    !! inverse of the original matrix. The determinant is simply the
    !! product of the diagonal elements of the triangular matrix.
    !! Method chosen as easy to implement and nicely numerically stable,
    !! especially for the small matrices need here.
    !!
    !! Written by I.J. Bush
    !!----------------------------------------------------

    Real( wp ), Dimension( 1:, 1: ),               Intent( In    )           :: A
    Real( wp ), Dimension(  :,  : ), Allocatable,  Intent(   Out )           :: B
    Real( wp )                     ,               Intent(   Out ), Optional :: det

    Real( wp ), Dimension( :, : ), Allocatable :: Q, L

    Real( wp ) :: theta, c, s
    Real( wp ) :: ann, anm
    Real( wp ) :: Lim, Lin
    Real( wp ) :: Qmi, Qni
    Real( wp ) :: mv

    Integer :: dim
    Integer :: m, n
    Integer :: i

    Allocate( B, Mold = A )

    ! Size of the problem
    dim = Size( A, dim = 1 )

    ! The factors of the matrix A
    Allocate( L( 1:dim, 1:dim ) )
    Allocate( Q( 1:dim, 1:dim ) )

    ! Factorise A by Givens rotations
    L = A
    Q = 0.0_wp
    Do i = 1, dim
      Q( i, i ) = 1.0_wp
    End Do
    LQ_factorise: Do m = 2, dim
      Do n = 1, m - 1
        ! For element mn generate the rotation which
        ! will zero that element
        ann = L( n, n )
        anm = L( n, m )
        theta = atan2( - anm, ann )
        c = Cos( theta )
        s = Sin( theta )
        ! And apply that rotation to the original matrix
        Do i = 1, Dim
          Lim = L( i, m )
          Lin = L( i, n )
          L( i, n ) = c * Lin - s * Lim
          L( i, m ) = c * Lim + s * Lin
        End Do
        ! And then apply the rotation to the orthognal matrix
        Do i = 1, Dim
          Qmi = Q( m, i )
          Qni = Q( n, i )
          Q( n, i ) = c * Qni - s * Qmi
          Q( m, i ) = c * Qmi + s * Qni
        End Do
      End Do
    End Do LQ_factorise

    ! Now invert the triangular matrix, returning the inverse in B
    B(  1, 1  ) = 1.0_wp / L( 1, 1 )
    B(  1, 2: ) = 0.0_wp
    Do m = 2, dim
      Do n = 1, m
        mv = 0.0_wp
        Do i = 1, m - 1
          mv = mv - L( m, i ) * B( i, n )
        End Do
        If( m == n ) Then
          mv = mv + 1.0_wp
        End If
        mv = mv / L( m, m )
        B( m, n ) = mv
      End Do
      B( m, n: ) = 0.0_wp
    End Do

    ! The inverse of A is now trivially formed
    B = Matmul( Transpose( Q ), B )

    ! And generate the determinant
    If( Present( det ) ) Then
      det = 1.0_wp
      Do i = 1, dim
        det = det * L( i, i )
      End Do
    End If

    Deallocate( Q )
    Deallocate( L )

  End Subroutine givens_invert

End Program testit
