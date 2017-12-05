/*
   NMExtract.cpp

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

   This program extracts raw mMR data.
 */

#include <gdcmReader.h>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <glog/logging.h>

#include "MMR.hpp"
#include "EnvironmentInfo.h"

int main(int argc, char **argv)
{

  const char* APP_NAME = "nm_extract";

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


  return 0;
}