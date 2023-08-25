cmake_minimum_required(VERSION 3.12)  # Adjust the required version as needed
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
project(karst)  # Change to your project name

set(CMAKE_CXX_STANDARD 11)  # Set the C++ standard
set(CMAKE_VERBOSE_MAKEFILE ON)
add_compile_options(-Wall -pedantic -O2)

#if (UNIX AND NOT APPLE)
#    set(CC clang)
#    set(CXX clang++)
#    set(FL gfortran)
#    set(MDIR ~/.nix-profile/opt/mumps)
#    set(FLIBS -L~/.nix-profile/lib -lgfortran -lquadmath -lblas)
#else()
set(CC /opt/homebrew/bin/gcc-12)
set(CXX /opt/homebrew/bin/g++-12)
set(FL /opt/homebrew/bin/gfortran)
set(CMAKE_C_COMPILER /opt/homebrew/bin/gcc-12)
set(CMAKE_CXX_COMPILER /opt/homebrew/bin/g++-12)
set(MDIR /Users/agnieszka/Desktop/KARST/MUMPS_5.6.0)
set(FLIBS -L/opt/homebrew/lib -L/opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/12 -lblas -llapack -lesmumps -lscotch -lscotcherr -lmumps_common -lmetis -lgomp -lgfortran)
set(MLIBS  -L${MDIR} -L${MDIR}/lib -L${MDIR}/libseq)
#endif()


# Add the executable target
file(GLOB SOURCES ./src/*.cc)
add_executable(karst ${SOURCES})

# Set compiler flags
target_compile_options(karst PRIVATE
        -Wall -pedantic -O2
        -std=c++11 -Weffc++ -Wextra -Wfatal-errors -Wunused -pedantic -w -Winline -O2
        )

# Set include directories
target_include_directories(karst PRIVATE
        ${MDIR}/include
#        ${MDIR}/libseq
#        ${MDIR}/lib
#        /opt/homebrew/lib
#        /opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/12
#        /Users/agnieszka/Desktop/KARST/MUMPS_5.6.0/lib
#        /Users/agnieszka/Desktop/KARST/MUMPS_5.6.0/libseq
#        /Users/agnieszka/Desktop/KARST/MUMPS_5.6.0/include
        ./algorithms
        ./include
        )
#
## Add the directory containing the algorithms library to link directories
link_directories(./algorithms)


# Link libraries
target_link_libraries(karst PRIVATE
        ${FLIBS}
        ${MLIBS}
        -ldmumps
        -lmpiseq
        -lpthread -lpord -lstdc++
        -lm -L../algorithms -lalgorithms
        ${CMAKE_SOURCE_DIR}/algorithms/libalgorithms.a
#        /Users/agnieszka/Desktop/KARST/MUMPS_5.6.0/lib/libmumps_common.a
#        /Users/agnieszka/Desktop/KARST/MUMPS_5.6.0/lib/libdmumps.a
        )


