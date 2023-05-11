cmake_minimum_required(VERSION 3.10)

set(VERSION "0.1.0")

project(AMPT-Genesis VERSION ${VERSION} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

#configure a header file to pass the version number to the source code:
#configure_file(TutorialConfig.h.in TutorialConfig.h)

# set default build type
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
    #set(CMAKE_BUILD_TYPE "Debug")
endif()

set( CMAKE_CXX_FLAGS "${SSE_FLAGS} -std=c++17")

find_package(Boost REQUIRED COMPONENTS program_options)

if(Boost_FOUND)
    message(STATUS "Boost_INCLUDE_DIR: ${Boost_INCLUDE_DIR}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

#Tells the source code that root is present
find_package(ROOT)
if(ROOT_FOUND)
    message(STATUS "ROOT_INCLUDE_DIRS: ${ROOT_INCLUDE_DIRS}")
    message(STATUS "ROOT_LIBRARIES: ${ROOT_LIBRARIES}")
    message(STATUS "ROOT_VERSION: ${ROOT_VERSION}")
endif()

include_directories(SYSTEM ${Boost_INCLUDE_DIR} ${ROOT_INCLUDE_DIRS})

#set names of the test executables
set(AMPT-Genesis.exe "AMPT-Genesis.exe")
add_subdirectory(src)


add_executable(AMPT-Genesis.exe src/main.cpp)
target_compile_options(AMPT-Genesis.exe PUBLIC -I${Boost_INCLUDE_DIRS})
target_link_libraries(AMPT-Genesis.exe PUBLIC ROOT::Eve)

add_dependencies(AMPT-Genesis.exe
                 AMPT_smearer
                 hydrodynamizer
                 parton_collision
                 AMPTGenesis
                 )

target_link_libraries("AMPT-Genesis.exe" PUBLIC Boost::program_options
                       ${ROOT_LIBRARIES}
                       )
