add_executable(nm_validate NMValidate.cpp  )
target_link_libraries(nm_validate PUBLIC
      ${ITK_LIBRARIES}  
      ${Boost_LIBRARIES}
      glog::glog
      nlohmann_json::nlohmann_json
    )
# possibly not necessary, but might prevent problems with linking boost static/shared libs
set_target_properties(nm_validate PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_executable(nm_extract NMExtract.cpp  )
target_link_libraries(nm_extract PUBLIC
      ${Boost_LIBRARIES} 
      ${ITK_LIBRARIES} 
      glog::glog
      nlohmann_json::nlohmann_json
    )
set_target_properties(nm_extract PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_executable(nm_mrac2mu NMmrac2mu.cpp  )
target_link_libraries(nm_mrac2mu PUBLIC
      ${Boost_LIBRARIES}
      ${ITK_LIBRARIES} 
      glog::glog 
      nlohmann_json::nlohmann_json
    )
set_target_properties(nm_mrac2mu PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_executable(nm_signa2mu signa.cpp  )
target_link_libraries(nm_signa2mu PUBLIC
        ${Boost_LIBRARIES}
        ${ITK_LIBRARIES}
        glog::glog
        nlohmann_json::nlohmann_json
        )
set_target_properties(nm_signa2mu PROPERTIES POSITION_INDEPENDENT_CODE ON)

install(TARGETS nm_validate DESTINATION bin)
install(TARGETS nm_extract DESTINATION bin)
install(TARGETS nm_mrac2mu DESTINATION bin)
install(TARGETS nm_signa2mu DESTINATION bin)