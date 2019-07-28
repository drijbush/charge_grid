# Set up library
add_library(q_grid cache_module.f90 charge_gridder_calculate_module.f90 charge_gridder_module.f90 charge_gridder_recurse_module.f90 charge_ranger_cubic_module.f90 charge_ranger_module.f90 charge_ranger_radial_module.f90 numbers_module.f90 vec_operators_module.f90 )

add_executable(test_q_grid testit.f90)
target_link_libraries(test_q_grid q_grid m)
