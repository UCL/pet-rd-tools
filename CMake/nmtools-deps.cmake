# CMake >= 4 removed FindBoost, so rely on Boost's own per-library CMake configs
# (Boost >= 1.82, or distro packages that install them). Linking the imported
# Boost::<lib> targets also picks up the per-library *_DYN_LINK definitions
# Boost needs on Windows.
find_package(Boost 1.82 REQUIRED COMPONENTS filesystem program_options regex)

# Requires the itk::SpatialOrientationEnums API (ITK >= 5.3).
find_package(ITK 5.3 REQUIRED)
include(${ITK_USE_FILE})
#if (NOT ITKReview_LOADED)
#	message(FATAL_ERROR "ITK should be built with the Module_ITKReview enabled.")
#endif()

find_package(glog REQUIRED)

find_package(nlohmann_json 3.2.0 CONFIG REQUIRED)

