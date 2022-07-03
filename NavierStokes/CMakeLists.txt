file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(Mymain main.cpp ${sources} ${headers})

add_executable(NavierStokes ${Mymain})

target_link_libraries(NavierStokes NeoPZ::pz)

