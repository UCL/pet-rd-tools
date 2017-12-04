add_executable(stub stub.cxx  )
target_link_libraries(stub 
      ${GDCM_LIBRARIES} 
      ${Boost_LIBRARIES} 
      glog 
    )

install(TARGETS stub DESTINATION bin)