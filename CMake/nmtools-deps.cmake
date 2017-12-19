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

include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(glog CONFIG REQUIRED)
include_directories(${GLOG_INCLUDE_DIRS}
