add_executable(nm_validate NMValidate.cpp  )
target_link_libraries(nm_validate 
      ${GDCM_LIBRARIES} 
      ${Boost_LIBRARIES} 
      glog 
    )

add_executable(nm_extract NMExtract.cpp  )
target_link_libraries(nm_extract 
      ${GDCM_LIBRARIES} 
      ${Boost_LIBRARIES} 
      glog 
    )

add_executable(nm_mrac2mu NMmrac2mu.cpp  )
target_link_libraries(nm_mrac2mu 
      ${GDCM_LIBRARIES} 
      ${Boost_LIBRARIES}
      ${ITK_LIBRARIES} 
      glog 
    )

install(TARGETS nm_validate DESTINATION bin)
install(TARGETS nm_extract DESTINATION bin)
install(TARGETS nm_mrac2mu DESTINATION bin)