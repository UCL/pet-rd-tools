add_executable(nm_validate NMValidate.cpp)
target_link_libraries(nm_validate
      ${ITK_LIBRARIES}
      Boost::filesystem Boost::program_options Boost::regex
      glog::glog
      nlohmann_json::nlohmann_json
    )

add_executable(nm_extract NMExtract.cpp)
target_link_libraries(nm_extract
      ${ITK_LIBRARIES}
      Boost::filesystem Boost::program_options Boost::regex
      glog::glog
      nlohmann_json::nlohmann_json
    )

add_executable(nm_mrac2mu NMmrac2mu.cpp)
target_link_libraries(nm_mrac2mu
      ${ITK_LIBRARIES}
      Boost::filesystem Boost::program_options Boost::regex
      glog::glog
      nlohmann_json::nlohmann_json
    )

add_executable(nm_signa2mu signa.cpp)
target_link_libraries(nm_signa2mu
      ${ITK_LIBRARIES}
      Boost::filesystem Boost::program_options Boost::regex
      glog::glog
      nlohmann_json::nlohmann_json
    )

install(TARGETS nm_validate nm_extract nm_mrac2mu nm_signa2mu DESTINATION bin)
