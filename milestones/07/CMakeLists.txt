configure_file(cluster_147_0.xyz .)
configure_file(cluster_309_0.xyz .)
configure_file(cluster_561_0.xyz .)
configure_file(cluster_923_0.xyz .)
configure_file(cluster_2057_0.xyz .)
configure_file(cluster_3871_0.xyz .)
configure_file(cluster_5083_0.xyz .)
configure_file(cluster_10179_0.xyz .)

# Defining the executable target
add_executable(07 main.cpp)

# Linking against our MD code (propagates all dependencies)
target_link_libraries(07 PUBLIC my_md_lib cxxopts)
