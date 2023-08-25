configure_file(cluster_923.xyz .)
configure_file(cluster_3871.xyz .)
configure_file(equilibrate.sh .)
configure_file(heat_cap.sh .)

# Defining the executable target
add_executable(07 main.cpp)

# Linking against our MD code (propagates all dependencies)
target_link_libraries(07 PUBLIC my_md_lib cxxopts)
