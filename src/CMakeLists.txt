add_executable(stub stub.cpp  )
target_link_libraries(stub 
      ${GDCM_LIBRARIES} 
      ${Boost_LIBRARIES} 
      glog 
    )

install(TARGETS stub DESTINATION bin)

add_executable(nm_validate NMValidate.cpp  )
target_link_libraries(nm_validate 
      ${GDCM_LIBRARIES} 
      ${Boost_LIBRARIES} 
      glog 
    )

install(TARGETS nm_validate DESTINATION bin)