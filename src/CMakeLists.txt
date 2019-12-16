add_executable(nm_validate NMValidate.cpp  )
target_link_libraries(nm_validate 
      ${ITK_LIBRARIES}  
      ${Boost_LIBRARIES}
      glog::glog 
    )

add_executable(nm_extract NMExtract.cpp  )
target_link_libraries(nm_extract 
      ${Boost_LIBRARIES} 
      ${ITK_LIBRARIES} 
      glog::glog 
    )

add_executable(nm_mrac2mu NMmrac2mu.cpp  )
target_link_libraries(nm_mrac2mu 
      ${Boost_LIBRARIES}
      ${ITK_LIBRARIES} 
      glog::glog 
    )

add_executable(nm_signa2mu signa.cpp  )
target_link_libraries(nm_signa2mu
        ${Boost_LIBRARIES}
        ${ITK_LIBRARIES}
        glog::glog
        )

install(TARGETS nm_validate DESTINATION bin)
install(TARGETS nm_extract DESTINATION bin)
install(TARGETS nm_mrac2mu DESTINATION bin)
install(TARGETS nm_signa2mu DESTINATION bin)