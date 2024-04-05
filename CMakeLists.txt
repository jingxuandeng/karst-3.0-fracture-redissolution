#cmake_minimum_required(VERSION 3.12)  # Adjust the required version as needed
cmake_minimum_required(VERSION 3.6)  # Adjust the required version as needed
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
project(karst)  # Change to your project name

set(CMAKE_CXX_STANDARD 11)  # Set the C++ standard
set(CMAKE_VERBOSE_MAKEFILE ON)
add_compile_options(-Wall -pedantic -O2)




if (UNIX AND NOT APPLE)
    message(STATUS "Running on a Linux.")
    set(CC clang)
    set(CXX clang++)
    set(FL gfortran)
    set(MDIR ~/.nix-profile/opt/mumps)
    set(FLIBS -L~/.nix-profile/lib -lgfortran -lquadmath -lblas)
else()
    # Execute 'sw_vers' command to get macOS version
    execute_process(
            COMMAND sw_vers -productVersion
            OUTPUT_VARIABLE MACOS_VERSION
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    # Check macOS version and set CXX compiler accordingly
    if (MACOS_VERSION VERSION_GREATER "14.1")
        message(STATUS "Running on a new MACOS")
        set(CC /opt/homebrew/bin/gcc-13)
        set(CXX /opt/homebrew/bin/g++-13)
        set(FL /opt/homebrew/bin/gfortran)
        set(CMAKE_C_COMPILER /opt/homebrew/bin/gcc-13)
        set(CMAKE_CXX_COMPILER /opt/homebrew/bin/g++-13)
        set(MDIR /Users/agnieszka/Desktop/KARST/MUMPS_5.6.0)
        set(FLIBS -L/opt/homebrew/lib -L/opt/homebrew/Cellar/gcc/13.2.0/lib/gcc/13 -lblas -llapack -lesmumps -lscotch -lscotcherr -lmumps_common -lmetis -lgomp -lgfortran)
        set(MLIBS  -L${MDIR} -L${MDIR}/lib -L${MDIR}/libseq)
    else (MACOS_VERSION GREATER_EQUAL "9.14")
        message(STATUS "Running on an old  MACOS")
        set(CMAKE_C_COMPILER "clang")
        set(CMAKE_CXX_COMPILER "clang++")
        set(CC clang)
        set(CXX clang++)
        set(FL gfortran)
        set(MDIR /Users/agnieszka/Desktop/KARST/MUMPS_4.7.3)
        set(FLIBS -lgfortran -lquadmath -lblas)
        set(MLIBS  -L${MDIR} -L${MDIR}/lib -L${MDIR}/libseq)
    endif()
endif()

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


