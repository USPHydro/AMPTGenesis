include_directories(SYSTEM ${Boost_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS})

add_library(AMPT_smearer AMPT_smearer.cpp)
add_library(hydrodynamizer hydrodynamizer.cpp)
add_library(parton_collision parton_collision.cpp)
add_library(AMPTGenesis AMPTGenesis.cpp)

target_compile_options(AMPT_smearer PUBLIC "-I${Boost_INCLUDE_DIRS}" )
target_compile_options(hydrodynamizer PUBLIC "-I${ROOT_INCLUDE_DIRS}" )
target_compile_options(AMPTGenesis PUBLIC "-I${Boost_INCLUDE_DIRS}"  )

#target_link_libraries(AMPTGenesis PUBLIC ROOT::Eve)
#target_link_libraries(AMPTGenesis PUBLIC  ${ROOT_LIBRARIES})

#set_target_properties (AMPT_smearer PROPERTIES COMPILE_FLAGS "${CompileFlags}")
#set_target_properties (hydrodynamizer PROPERTIES COMPILE_FLAGS "${CompileFlags}")
#set_target_properties (parton_collision PROPERTIES COMPILE_FLAGS "${CompileFlags}")image.png
#set_target_properties (progress_bar PROPERTIES COMPILE_FLAGS "${CompileFlags}")


