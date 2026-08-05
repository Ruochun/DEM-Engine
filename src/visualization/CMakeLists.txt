# The visualizer stays in a separate object target so the C++ library and Python extension can both include it without
# forcing raylib into DEME's solver internals.
add_library(visualization OBJECT DEMVisualizer.cpp)

target_include_directories(visualization
	PUBLIC ${ProjectIncludeSource}
	PUBLIC ${ProjectIncludeGenerated}
)

target_link_libraries(visualization PUBLIC raylib)

set_target_properties(visualization PROPERTIES
	POSITION_INDEPENDENT_CODE TRUE
	CXX_STANDARD ${CXXSTD_SUPPORTED}
)
