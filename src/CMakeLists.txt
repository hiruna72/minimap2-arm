include_directories(${CMAKE_CURRENT_BINARY_DIR})

## until making f5c a C project use g++
#set_source_files_properties(align.c PROPERTIES LANGUAGE CXX)
#set_source_files_properties(bseq.c PROPERTIES LANGUAGE CXX)
#set_source_files_properties(chain.c PROPERTIES LANGUAGE CXX)
#set_source_files_properties(esterr.c PROPERTIES LANGUAGE CXX)
#set_source_files_properties(format.c PROPERTIES LANGUAGE CXX)
#set_source_files_properties(hit.c PROPERTIES LANGUAGE CXX)
#set_source_files_properties(index.c PROPERTIES LANGUAGE CXX)
#set_source_files_properties(kalloc.c PROPERTIES LANGUAGE CXX)
#set_source_files_properties(kthread.c PROPERTIES LANGUAGE CXX)
#set_source_files_properties(main.c PROPERTIES LANGUAGE CXX)
#set_source_files_properties(map.c PROPERTIES LANGUAGE CXX)
#set_source_files_properties(misc.c PROPERTIES LANGUAGE CXX)
#set_source_files_properties(options.c PROPERTIES LANGUAGE CXX)
#set_source_files_properties(pe.c PROPERTIES LANGUAGE CXX)
#set_source_files_properties(sdust.c PROPERTIES LANGUAGE CXX)
#set_source_files_properties(sketch.c PROPERTIES LANGUAGE CXX)
#set_source_files_properties(merge.c PROPERTIES LANGUAGE CXX)

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

# platform dependent libraries
IF (${DEPLOY_PLATFORM} STREQUAL "arm64-v8a")

	set (KSW2_EXTD2_SSE "${PROJECT_SOURCE_DIR}/3rdparty/${DEPLOY_PLATFORM}/ksw2_extd2_neon.o")
	set (KSW2_EXTS2_SSE "${PROJECT_SOURCE_DIR}/3rdparty/${DEPLOY_PLATFORM}/ksw2_exts2_neon.o")
	set (KSW2_EXTZ2_SSE "${PROJECT_SOURCE_DIR}/3rdparty/${DEPLOY_PLATFORM}/ksw2_extz2_neon.o")
	set (KSW2_LL_SSE "${PROJECT_SOURCE_DIR}/3rdparty/${DEPLOY_PLATFORM}/ksw2_ll_sse.o")

	# Build a shared lib
	add_library(minimapshared SHARED ${ALIGN} ${BSEQ} ${CHAIN} ${ESTERR} ${FORMAT} ${HIT} ${INDEX} ${KALLOC}
			${KTHREAD} ${MAIN} ${MAP} ${MISC} ${OPTIONS} ${PE} ${SDUST} ${SKETCH} ${MERGE})

ELSE()

	set (KSW2_EXTD2_SSE ksw2_extd2_sse.c ksw2.h kalloc.h)
	set (KSW2_EXTS2_SSE ksw2_exts2_sse.c ksw2.h kalloc.h)
	set (KSW2_EXTZ2_SSE ksw2_extz2_sse.c ksw2.h kalloc.h)
	set (KSW2_LL_SSE ksw2_ll_sse.c ksw2.h kalloc.h)

	# Build a shared lib
	add_library(minimapshared SHARED ${ALIGN} ${BSEQ} ${CHAIN} ${ESTERR} ${FORMAT} ${HIT} ${INDEX} ${KALLOC} ${KSW2_EXTD2_SSE} ${KSW2_EXTS2_SSE} ${KSW2_EXTZ2_SSE} ${KSW2_LL_SSE}
			${KTHREAD} ${MAIN} ${MAP} ${MISC} ${OPTIONS} ${PE} ${SDUST} ${SKETCH} ${MERGE})

	# to create an executatble
	#add_executable(minimapexe ${ALIGN} ${BSEQ} ${CHAIN} ${ESTERR} ${FORMAT} ${HIT} ${INDEX} ${KALLOC} ${KSW2_EXTD2_SSE} ${KSW2_EXTS2_SSE} ${KSW2_EXTZ2_SSE} ${KSW2_LL_SSE} ${KTHREAD} ${MAIN} ${MAP} ${MISC} ${OPTIONS} ${PE} ${SDUST} ${SKETCH} ${MERGE})
	#target_link_libraries(minimapexe -lz -lpthread -ldl)

ENDIF()




# link with dependencies
IF (${DEPLOY_PLATFORM} STREQUAL "x86")
	target_link_libraries(minimapshared -lz -lpthread -ldl)
ELSE()
	find_library(ANDROID_LOG_LIB log)
	target_link_libraries(minimapshared ${ANDROID_LOG_LIB} ${KSW2_EXTD2_SSE} ${KSW2_EXTS2_SSE} ${KSW2_EXTZ2_SSE} ${KSW2_LL_SSE} -lz -ldl)
ENDIF()