add_library(PolylineFactory
    src/polylinefactory.cpp
    src/casebuilder.cpp
    include/polylinefactory.hpp
    include/casebuilder.hpp
    include/testhelpers.hpp
    )

target_include_directories(PolylineFactory
    PUBLIC include/)

target_link_libraries(PolylineFactory
    PUBLIC ${CAVC_CPP_HEADER_ONLY_LIB})

target_link_libraries(PolylineFactory
    PUBLIC ${CAVC_C_API_LIB})
