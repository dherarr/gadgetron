IF (WIN32)
	ADD_DEFINITIONS(-DTIXML_USE_STL)
	ADD_DEFINITIONS(-D__BUILD_GADGETRON_MRICORE__)
	ADD_DEFINITIONS(-D__BUILD_GADGETS__)
ENDIF (WIN32)

include_directories(   ${ACE_INCLUDE_DIR} 
				       ${Boost_INCLUDE_DIR}
				       ${FFTW3_INCLUDE_DIR}
					   ${CMAKE_SOURCE_DIR}/dependencies/tinyxml
				       ${CMAKE_SOURCE_DIR}/apps/gadgetron 
				       ${CMAKE_SOURCE_DIR}/toolboxes/ndarray
				       ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
				       ${CMAKE_SOURCE_DIR}/toolboxes/linalg
				       ${CMAKE_SOURCE_DIR}/gadgets/core
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpucore)

add_library(gadgetronmricore SHARED MRINoiseAdjustGadget.cpp)

target_link_libraries(gadgetronmricore tinyxml hondarray linalg ${FFTW3_LIBRARIES} optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

install (FILES MRINoiseAdjustGadget.h DESTINATION include)

install(TARGETS gadgetronmricore DESTINATION lib)
install(FILES default_optimized.xml DESTINATION config)