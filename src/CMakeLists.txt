
# Find dependencies
find_package(Eigen3 REQUIRED)
find_package(Boost 1.65 )

# Source files
set(SRC
    mesh.hpp
    basis.hpp
    geometry.hpp
    functionspace.hpp
    error.hpp
    tests.cpp
    basis.cpp
)

# Add executable for tests
add_executable(fem_tests tests.cpp basis.cpp mesh.cpp functionspace.cpp geomap.cpp error.cpp)
target_include_directories(fem_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(fem_tests PRIVATE Eigen3::Eigen)

enable_testing()
add_test(NAME FEMTests COMMAND fem_tests)