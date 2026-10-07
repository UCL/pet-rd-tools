# cmake>=4.1: no FindBoost (CMP0167)
# boost>=1.70: own config (preferred)
find_package(Boost 1.54.0 CONFIG QUIET
  COMPONENTS
    chrono
    date_time
    filesystem
    program_options
    regex
    system
    thread
)
if (NOT Boost_FOUND)
  find_package(Boost 1.54.0 MODULE
    COMPONENTS
      chrono
      date_time
      filesystem
      program_options
      regex
      system
      thread
    REQUIRED
  )
endif()

include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

if (WIN32)
   add_definitions(-DBOOST_ALL_NO_LIB)
   add_definitions(-DBOOST_ALL_DYN_LINK)
endif()

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
#if (NOT ITKReview_LOADED)
#	message(FATAL_ERROR "ITK should be built with the Module_ITKReview enabled.")
#endif()

find_package(glog REQUIRED)

find_package(nlohmann_json 3.2.0 CONFIG)

