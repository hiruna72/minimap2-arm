# Disable in-source builds to prevent source tree corruption.
if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "FATAL: In-source builds are not allowed.
       You should create a separate directory for build files.")
endif()

# choose the architecture
set(DEPLOY_PLATFORM "unknown" CACHE STRING "choose the instruction set architecture")
set_property(CACHE DEPLOY_PLATFORM PROPERTY STRINGS unknown x86 arm64-v8a armeabi-v7a)
IF (${DEPLOY_PLATFORM} STREQUAL "unknown")
    message(FATAL_ERROR "choose the DEPLOY_PLATFORM")
    return() # This is to stop proceeding further and to stop opencv getting set to the default ANDROID_ABI
ENDIF()


# Set minimum cmake version + cmake settings
cmake_minimum_required (VERSION 3.9)

# Specify the project, and set the use of c++ and c
project(MINIMAP2 LANGUAGES CXX C)

message(STATUS "Configure Minimap2")

# Set project name
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -fPIC -Wall -O2 -Wc++-compat -D_FILE_OFFSET_BITS=64 -fsigned-char" )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -O2 -std=c++11 -D_FILE_OFFSET_BITS=64 -fsigned-char -DHAVE_KALLOC")


include_directories(${PROJECT_SOURCE_DIR})

IF (${ANDROID_ARM_NEON})
	message(STATUS "NEON build...")
	include_directories(${PROJECT_SOURCE_DIR}/sse2neon)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DKSW_SSE2_ONLY -D__SSE2__" ) 
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DKSW_SSE2_ONLY -D__SSE2__")
ENDIF()


set (align align.c)
set (bseq bseq.c)
set (chain chain.c)
set (esterr esterr.c)
set (format format.c)
set (hit hit.c)
set (index index.c)
set (kalloc kalloc.c)
set (kthread kthread.c)
set (main main.c)
set (map map.c)
set (misc misc.c)
set (options options.c)
set (pe pe.c)
set (sdust sdust.c)
set (sketch sketch.c)
set (splitidx splitidx.c)

set (ksw2_extd2_sse ksw2_extd2_sse.c)
set (ksw2_exts2_sse ksw2_exts2_sse.c)
set (ksw2_extz2_sse ksw2_extz2_sse.c)
set (ksw2_ll_sse ksw2_ll_sse.c)


# Build a shared lib
add_library(minimapshared SHARED ${main} ${align} ${bseq} ${chain} ${esterr} ${format} ${hit} ${index} ${kalloc} ${ksw2_extd2_sse}
	${ksw2_exts2_sse} ${ksw2_extz2_sse} ${ksw2_ll_sse} ${kthread}  ${map} ${misc} ${options} ${pe} ${sdust} ${sketch} ${splitidx})

# link with dependencies
IF (${DEPLOY_PLATFORM} STREQUAL "x86")
    target_link_libraries(minimapshared -lz -ldl -lm -lpthread)
ELSE()
    target_link_libraries(minimapshared -lz -ldl -lm)
ENDIF()

#  compile init_minimap2_example and link to shared lib
# Create the executable
add_executable(minimap2_example examples/init_minimap2_example.c)

# Link the required libraries
target_link_libraries(minimap2_example minimapshared)
