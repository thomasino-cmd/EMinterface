# Install script for directory: C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/pkgs/eigen3_x64-windows")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "OFF")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Devel" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/AdolcForward"
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/AlignedVector3"
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/ArpackSupport"
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/AutoDiff"
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/BVH"
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/EulerAngles"
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/FFT"
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/IterativeSolvers"
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/KroneckerProduct"
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/LevenbergMarquardt"
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/MatrixFunctions"
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/MoreVectorization"
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/MPRealSupport"
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/NonLinearOptimization"
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/NumericalDiff"
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/OpenGLSupport"
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/Polynomials"
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/Skyline"
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/SparseExtra"
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/SpecialFunctions"
    "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/Splines"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Devel" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/x64-windows-rel/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()

