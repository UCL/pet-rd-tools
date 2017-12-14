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

find_package(GDCM REQUIRED)
include(${GDCM_USE_FILE})

set(GDCM_LIBRARIES
  gdcmCommon
  gdcmDICT
  gdcmDSED
  gdcmexpat
  gdcmIOD
  gdcmjpeg12
  gdcmjpeg16
  gdcmjpeg8
  gdcmMEXD
  gdcmMSFF
  gdcmopenjpeg
  gdcmzlib
  socketxx
)

find_package(glog REQUIRED)
