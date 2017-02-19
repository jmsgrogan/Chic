find_package(Chaste COMPONENTS ode)
chaste_do_project(Chic)

# Add extra dependencies
set(MUSCLE_DIR "/opt/muscle-svn" CACHE FILEPATH "Path to the muscle installation.")
list(APPEND Chaste_INCLUDES "${MUSCLE_DIR}/include")
message(STATUS "${MUSCLE_DIR}/include")
list(APPEND Chaste_LINK_LIBRARIES "${MUSCLE_DIR}/libmuscle2.so")

set( CMAKE_SKIP_BUILD_RPATH true PARENT_SCOPE)

target_include_directories(chaste_project_Chic PUBLIC ${MUSCLE_DIR}/include)
target_link_libraries(chaste_project_Chic PUBLIC "${MUSCLE_DIR}/lib/libmuscle2.so")
