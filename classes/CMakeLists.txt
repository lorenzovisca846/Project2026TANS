set(LIB_SOURCES
    Cylinder.cxx
    MyPoint.cxx
    Particle.cxx
    SimRandom.cxx
    Config.cxx
    Tracklet.cxx
    Hit.cxx
    VertexReconstructor.cxx
)

set(LIB_HEADERS
    Cylinder.h
    MyPoint.h
    Particle.h
    SimRandom.h
    Config.h
    Tracklet.h
    Hit.h
    VertexReconstructor.h
)

root_generate_dictionary(G__Classes ${LIB_HEADERS} LINKDEF LinkDef.h)

add_library(ProjectClasses SHARED ${LIB_SOURCES} G__Classes.cxx)

target_link_libraries(ProjectClasses PUBLIC ${ROOT_LIBRARIES})