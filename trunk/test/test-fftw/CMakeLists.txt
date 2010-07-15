
#PROJECT(TEST)

#ADD_DEFINITIONS(-DD3)
#SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

# Make sure the compiler can find include files from our lib library.
INCLUDE_DIRECTORIES(/Users/mkv/local/include )

# Make sure the linker can find the Hello library once it is built.
LINK_DIRECTORIES(/Users/mkv/local/lib )

# Add executable called "helloDemo" that is built from the source files
# "demo.cxx" and "demo_b.cxx".  The extensions are automatically found.
#ADD_EXECUTABLE(test-fftw main_omp.cc)
ADD_EXECUTABLE(test-fftw test_mpi.cc)

# Link the executable to the Hello library.
#TARGET_LINK_LIBRARIES(test-fftw fourier  m  blitz fftw3  fftw3_threads pthread)

TARGET_LINK_LIBRARIES(test-fftw blitz fftw3_mpi fftw3  mpich m)