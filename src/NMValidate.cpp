/*
   NMValidate.cpp

   Author:      Benjamin A. Thomas
   Author:      Kris Thielemans

   Copyright 2017, 2020 Institute of Nuclear Medicine, University College London.
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   This program confirms the validity of mMR raw data files.
 */

#include <gdcmReader.h>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <glog/logging.h>

#include "nmtools/MMR.hpp"
#include "nmtools/GEPET.hpp"
#include "EnvironmentInfo.h"

int main(int argc, char **argv)
{

  const char* APP_NAME = "nm_validate";

  std::string inputFilePath;
  std::string outputFileName;

  //Set-up command line options
  namespace po = boost::program_options;
  namespace fs = boost::filesystem;
  namespace nm = nmtools;

  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "Print help information")
    ("version","Print version number")
    //("verbose,v", "Be verbose")
    ("input,i", po::value<std::string>(&inputFilePath)->required(), "Input file")
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
  fs::path log_path = fs::absolute(fs::current_path());
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

  fs::path srcPath = inputFilePath;
  
  //Check if input file path exists
  if (! fs::exists( srcPath ) )
  {
    LOG(ERROR) << "Input path: " << srcPath << " does not exist!";
    return 1;
  }

  //Check if it is a file.
  if (! fs::is_regular_file( srcPath ) )
  {
    LOG(ERROR) << srcPath.native() << " does not appear to be a file!";
    return EXIT_FAILURE;
  }

  //Try and read input file
  std::unique_ptr<nm::IDicomExtractor> reader;
  try {
    std::unique_ptr<nm::SiemensPETFactory> Siemens_factory(new nm::SiemensPETFactory);
    reader = Siemens_factory->Create(srcPath);
  }
  catch (...)
  {
    LOG(INFO) << "Not a Siemens file. Trying GE (only minimal checks will be performed though).";
  }
  
  if (reader == nullptr){
    std::unique_ptr<nm::GEPETFactory> GE_factory(new nm::GEPETFactory);
    reader = GE_factory->Create(srcPath);
  }

  //If no valid reader found, exit.
  if (reader == nullptr){
    LOG(ERROR) << "File appears to be INVALID (not recognised)";
    return EXIT_FAILURE;
  }

  //Check if the file is correct for the identified type.
  if (!reader->IsValid()) {
    LOG(ERROR) << "File appears to be INVALID";
    return EXIT_FAILURE;    
  }

  LOG(INFO) << "File appears to be VALID";

  //Print total execution time
  std::time_t stopTime = std::time( 0 ) ;
  unsigned int totalTime = stopTime - startTime;
  LOG(INFO) << "Time taken: " << totalTime << " seconds";
  LOG(INFO) << "Ended: " << std::asctime(std::localtime(&stopTime));

  return EXIT_SUCCESS;
}
