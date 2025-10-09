set(CMAKE_CXX_STANDARD 17)

# User options
option(MAXIMUM "Enable flags for maximum clique search" OFF)
option(MAXIMAL "Enable flags for maximal clique enumeration" OFF)

# For maximum clique search
if(MAXIMUM)
  add_compile_options(-DCTD_MULTIROUND)
  add_compile_options(-DREDUCING)
endif()

# For maximal clique enumeration
if(MAXIMAL)
  add_compile_options(-DSKIP_LEAF)
  add_compile_options(-DINDUCING)
  add_compile_options(-DLISTING_MAXIMAL)
endif()

# For individual technique evaluation
add_compile_options(-DTWOHOP_REDUCTION)
add_compile_options(-DWORST_CASE_OPTIMAL)
add_compile_options(-DPIVOTING)
# add_compile_options(-DNO_PRACTICAL)
# add_compile_options(-DRANDOM_PIVOT)
# add_compile_options(-DONLY_PREPROCESSING)
# add_compile_options(-DNO_COLOR_REDUCTION)


# For debuging
# add_compile_options(-DNO_REDUCTION)
# add_compile_options(-DPRINT_STAT)
# add_compile_options(-DTIGHT_COLORING)
# add_compile_options(-DTIGHT_RECOLORING)
# add_compile_options(-fsanitize=address)
# add_link_options(-fsanitize=address)
# add_compile_options(-lprofiler)
# add_link_options(-lprofiler)

# preprocessing
# add_compile_options(-DCTD_MULTIROUND)
# add_compile_options(-DNO_COLOR_TRUSS)
# add_compile_options(-DPRINT_STAT)
# add_compile_options(-DREDUCING)

# options
add_compile_options(-DNO_COLOR_TRUSS)
add_compile_options(-DRECOLOR)
add_compile_options(-DNO_PIVOT_FROM_X)

add_executable(
  main.out
  main.cc
)
