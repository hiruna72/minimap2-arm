# Find all the files which need to be compiled
include_directories(${PROJECT_SOURCE_DIR}/3rdparty)

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
cmake_policy(SET CMP0048 OLD)

# Specify the project, and set the use of c
project(MINIMAP2 LANGUAGES C)

# Set project name
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wall -O2 -Wc++-compat -D_FILE_OFFSET_BITS=64 -fsigned-char")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -rdynamic -Wall -O2 -std=c++11 ") -static

include_directories(${CMAKE_CURRENT_BINARY_DIR})

set (ALIGN align.c minimap.h mmpriv.h bseq.h ksw2.h kalloc.h)
set (BSEQ bseq.c bseq.h kvec.h kalloc.h kseq.h)
set (CHAIN chain.c minimap.h mmpriv.h bseq.h kalloc.h)
set (ESTERR esterr.c mmpriv.h minimap.h bseq.h)
set (FORMAT format.c kalloc.h mmpriv.h minimap.h bseq.h)
set (HIT hit.c mmpriv.h minimap.h bseq.h kalloc.h khash.h)
set (INDEX index.c kthread.h bseq.h minimap.h mmpriv.h kvec.h kalloc.h khash.h)
set (KALLOC kalloc.c kalloc.h)
set (KTHREAD kthread.c kthread.h)
set (MAIN main.c bseq.h minimap.h mmpriv.h getopt.h)
set (MAP map.c kthread.h kvec.h kalloc.h sdust.h mmpriv.h minimap.h bseq.h khash.h ksort.h)
set (MISC misc.c mmpriv.h minimap.h bseq.h ksort.h)
set (OPTIONS options.c mmpriv.h minimap.h bseq.h)
set (PE pe.c mmpriv.h minimap.h bseq.h kvec.h kalloc.h ksort.h)
set (SDUST sdust.c kalloc.h kdq.h kvec.h sdust.h)
set (SKETCH sketch.c kvec.h kalloc.h mmpriv.h minimap.h bseq.h)
set (MERGE merge.c merge.h)

# Build a shared lib
add_library(minimapshared SHARED ${ALIGN} ${BSEQ} ${CHAIN} ${ESTERR} ${FORMAT} ${HIT} ${INDEX} ${KALLOC}
		${KTHREAD} ${MAIN} ${MAP} ${MISC} ${OPTIONS} ${PE} ${SDUST} ${SKETCH} ${MERGE})

# platform dependent libraries
set (KSW2_EXTD2_SSE "${PROJECT_SOURCE_DIR}/3rdparty/${DEPLOY_PLATFORM}/ksw2_extd2_neon.o")
set (KSW2_EXTS2_SSE "${PROJECT_SOURCE_DIR}/3rdparty/${DEPLOY_PLATFORM}/ksw2_exts2_neon.o")
set (KSW2_EXTZ2_SSE "${PROJECT_SOURCE_DIR}/3rdparty/${DEPLOY_PLATFORM}/ksw2_extz2_neon.o")
set (KSW2_LL_SSE "${PROJECT_SOURCE_DIR}/3rdparty/${DEPLOY_PLATFORM}/ksw2_ll_sse.o")

# link with dependencies
IF (${DEPLOY_PLATFORM} STREQUAL "x86")
	target_link_libraries(minimapshared -lz -lpthread -ldl)
	target_link_libraries(minimapshared ${KSW2_EXTD2_SSE} ${KSW2_EXTS2_SSE} ${KSW2_EXTZ2_SSE} ${KSW2_LL_SSE} -lz -ldl -lpthread)
ELSE()
	find_library(ANDROID_LOG_LIB log)
	target_link_libraries(minimapshared ${ANDROID_LOG_LIB} ${KSW2_EXTD2_SSE} ${KSW2_EXTS2_SSE} ${KSW2_EXTZ2_SSE} ${KSW2_LL_SSE} -lz -ldl)
ENDIF()

