/*
   NMmrac2mu.cpp

   Author:      Benjamin A. Thomas

   Copyright 2017 Institute of Nuclear Medicine, University College London.
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   This program generates a mu-map from Siemens mMR MRAC for PET reconstruction.
 */

#include <gdcmReader.h>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <glog/logging.h>
#include <memory>

#include "nmtools/MRAC.hpp"
#include "EnvironmentInfo.h"

int main(int argc, char **argv)
{

  const char* APP_NAME = "nm_mrac2mu";

  std::string inputDirPath;
  std::string outputFilePath = "";
  std::string coordOrientation = "RAI";

  //Set-up command line options
  namespace po = boost::program_options;
  namespace fs = boost::filesystem;
  namespace nm = nmtools;

  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "Print help information")
    ("version","Print version number")
    //("verbose,v", "Be verbose")
    ("input,i", po::value<std::string>(&inputDirPath)->required(), "Input directory")
    ("output,o", po::value<std::string>(&outputFilePath)->required(), "Output file")
    ("orient", po::value<std::string>(&coordOrientation), "Output orientation: RAI, RAS or LPS")
    ("log,l", "Write log file");

  //Evaluate command line options
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc),
      vm); // can throw

    /** --help option
    */
    if (vm.count("help")) {
      std::cout << APP_NAME << std::endl
        << desc << std::endl;
      return EXIT_SUCCESS;
    }

    if (vm.count("version") ) {
      std::cout << APP_NAME << " : v" << VERSION_NO << std::endl;
      return EXIT_SUCCESS;
    }

    po::notify(vm); // throws on error

  } catch (po::error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
    std::cerr << desc << std::endl;
    return EXIT_FAILURE;
  }

  //Configure logging
  fs::path log_path = fs::complete(fs::current_path());
  log_path /= APP_NAME;
  log_path += "-";

  //Pretty coloured logging (if supported)
  FLAGS_colorlogtostderr = 1;

  if (vm.count("log")){
    FLAGS_alsologtostderr = 1;
  }
  else {
    FLAGS_logtostderr = 1;
  }

  google::InitGoogleLogging(argv[0]);
  google::SetLogDestination(google::INFO, log_path.string().c_str());

  std::time_t startTime = std::time( 0 ) ;

  //Application starts here
  LOG(INFO) << "Started: " << std::asctime(std::localtime(&startTime));
  LOG(INFO) << "Running '" << APP_NAME << "' version: " << VERSION_NO;

  fs::path srcPath = inputDirPath;

  //Check if input file even exists!
  if (!fs::exists(srcPath)) {
    LOG(ERROR) << "Input path " << srcPath << " does not exist!";
    return EXIT_FAILURE;
  }

  //Check if the input is a file.
  if (!fs::is_directory(srcPath)) {
    LOG(ERROR) << srcPath.native() << " does not appear to be a  directory!";
    return EXIT_FAILURE;
  }

  std::unique_ptr<nm::MRAC2MU> mrac;

  try {
    mrac.reset(new nm::MRAC2MU(srcPath, coordOrientation));
  } catch (bool){
    LOG(ERROR) << "Failed to create MRAC converter!";
    return EXIT_FAILURE;
  }

  //mrac->SetDesiredCoordinateOrientation(coordOrientation);
 
  if (mrac->Update()){
    LOG(INFO) << "Scaling and reslicing complete";
  } else {
    LOG(ERROR) << "Failed to scale and reslice";
    return EXIT_FAILURE;
  }

  if (mrac->Write(outputFilePath)){
    LOG(INFO) << "Writing complete";
  } else {
    LOG(ERROR) << "Failed to write output file!";
    return EXIT_FAILURE;
  }

  //Print total execution time
  std::time_t stopTime = std::time( 0 ) ;
  unsigned int totalTime = stopTime - startTime;
  LOG(INFO) << "Time taken: " << totalTime << " seconds";
  LOG(INFO) << "Ended: " << std::asctime(std::localtime(&stopTime));

  return EXIT_SUCCESS;
}