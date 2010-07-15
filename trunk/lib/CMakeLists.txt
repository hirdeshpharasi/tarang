
#SUBDIRS(basis_basicfn four scft fields IncFlow IncFluid)

SUBDIRS(IncFluid IncFlow fields plainfields four scft basis_basicfn)

SET(CMAKE_BUILD_TYPE distribution)
SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")
SET(CMAKE_C_FLAGS_DISTRIBUTION "-O3")
