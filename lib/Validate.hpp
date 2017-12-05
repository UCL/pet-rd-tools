/*
   Validate.hpp

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

 */

#ifndef VALIDATE_HPP
#define VALIDATE_HPP

#include <memory>

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

#include <glog/logging.h>

#include "Common.hpp"

namespace nmtools {

FileStatusCode CheckForSiemensBFFile( boost::filesystem::path src, uint64_t numOfWords  ) {
    //Takes a filepath src, switches the extension to .bf, then checks if:
    //    - The file exists
    //    - The total length = numOfWords

  using namespace nmtools;

  FILE *inFile;

  boost::filesystem::path bfPath = src;
  bfPath = boost::filesystem::change_extension(bfPath, ".bf").string();

  if (! (inFile = fopen( bfPath.string().c_str(), "rb" ))) {
      LOG(INFO) << "Cannot open " << bfPath.native();
      return FileStatusCode::EIOERROR;
  }

  fseeko64( inFile, 0L, SEEK_END );
  int64_t endOfFile = ftello64(inFile);

  fclose( inFile );

  LOG(INFO) << "File size in bytes: " << endOfFile;
  uint64_t actualNoWords = endOfFile / 4;

  LOG(INFO) << endOfFile << " / 4 = " << actualNoWords << " words";

  if ( endOfFile != numOfWords*4 )
  {
    LOG(INFO) << "Expected no. of LM words does not equal no. read!";
    return FileStatusCode::EBAD;
  }

  return FileStatusCode::EGOOD;

}

FileStatusCode ReadAsSiemensDICOM( boost::filesystem::path src ){
  //Try and open as a DICOM file and assume this is from a Siemens machine.

  using namespace nmtools;

  std::unique_ptr<gdcm::Reader> dicomReader(new gdcm::Reader);
  dicomReader->SetFileName( src.string().c_str() );

  if ( ! dicomReader->Read() ) {
    LOG(INFO) << "Unable to read as DICOM file";
    return FileStatusCode::EIOERROR;
  }

  //Get DICOM data
  const gdcm::DataSet& ds = dicomReader->GetFile().GetDataSet();

  //Get manufacturer
  const gdcm::Tag manufacturer(0x008,0x1090);
  const gdcm::DataElement& pdde= ds.GetDataElement(manufacturer);
  std::stringstream ss;
  ss << (char*) pdde.GetByteValue()->GetPointer();
  LOG(INFO) << "Manufacturer: " << ss.str() ;

  //Get Image Type
  const gdcm::Tag imageType(0x0008, 0x0008);
  const gdcm::DataElement& pdde2 = ds.GetDataElement(imageType);
  std::stringstream ss2;
  ss2 << (char*)pdde2.GetByteValue()->GetPointer();
  LOG(INFO) << "Image type: " << ss2.str();

  //Get Interfile header
  const gdcm::Tag headerTag(0x029, 0x1010);
  const gdcm::DataElement &headerData = ds.GetDataElement(headerTag);
  std::stringstream headerString;
  std::stringstream headerStringTmp;
  headerStringTmp << (char *) headerData.GetByteValue()->GetPointer();

  //If this is actually a SMS-MI v 3.2 file then get header from 0029,1110.
  if (headerStringTmp.str().find("SV10") != std::string::npos) {
    gdcm::Tag altHeaderTag(0x029, 0x1110);
    const gdcm::DataElement &altHeaderData = ds.GetDataElement(altHeaderTag);
    std::stringstream tmpstring;
    tmpstring << (char *) altHeaderData.GetByteValue()->GetPointer();
    headerString << tmpstring.str();
  }
  else {
    headerString << headerStringTmp.str();
  }

  //Look for expected no. of list mode words
  std::string target = "%total listmode word counts";
  if ( headerString.str().find( target ) == std::string::npos ) {
    LOG(INFO) << "No word count found in Interfile header";
    return FileStatusCode::EBAD;
  }

  std::string lastline = headerString.str().substr(headerString.str().find(target),headerString.str().length());
  lastline = lastline.substr(0,lastline.find("\n"));

  boost::regex pattern("[0-9]+");
  boost::smatch result;

  uint64_t expectedNoWords = 0;
  if (boost::regex_search(lastline, result, pattern)) {
    expectedNoWords = boost::lexical_cast<unsigned long>(result);
  }
  else {
    LOG(INFO) << "No word count found in Interfile header";
    return FileStatusCode::EBAD;
  }

  LOG(INFO) << "Expected number of LM words: " << expectedNoWords;

    //Check length of listmode data stored in tag (7fe1,1010).
  const gdcm::Tag lmDataTag(0x7fe1, 0x1010);
  const gdcm::DataElement& lmData= ds.GetDataElement(lmDataTag);
  const gdcm::ByteValue *bv = lmData.GetByteValue();

  uint64_t lmLength = bv->GetLength();
  LOG(INFO) << lmLength << " bytes in LM field";

  uint64_t actualNoWords = lmLength / 4;

  LOG(INFO) << lmLength << " / 4 = " << actualNoWords << " words";

    //If the length of data in field (7fe1,1010) != number of words expected in
    //the interfile header, check for the existence of an associated .bf file
    //and see if the length of the .bf is correct.
  if (lmLength != expectedNoWords * 4) {
    LOG(INFO) << "Expected no. of LM words does not equal no. read!";
    LOG(INFO) << "Looking for BF file...";
    return CheckForSiemensBFFile( src, expectedNoWords );
  }

  return FileStatusCode::EGOOD;

}

FileStatusCode ReadAsSiemensPTD( boost::filesystem::path src ){
    //Checks if a given raw data file is in .ptd format.
    //A .ptd file starts with the raw listmode data and then has a DICOM header
    //placed at the end.
    //
    //This function moves to the end of the file and scans backwards until it
    //matches 'DICM' for the start of the DICOM header. The header is then
    //scanned to find the INTERFILE header. The expected number of list mode
    //words is extracted from the header and compared to the position of the
    //start of the DICOM header (less 128 bytes).

    //Open file
  FILE *inFile;

  if (! (inFile = fopen( src.string().c_str(), "rb" ))) {
    LOG(INFO) << "Cannot open " << src.native();
    return FileStatusCode::EIOERROR;
  }

    //Skip to the end
  fseeko64( inFile, 0L, SEEK_END );
  int64_t endOfFile = ftello64(inFile);

  LOG(INFO) << "File size in bytes: " << endOfFile;

  std::string target = "DICM";
  std::vector<char> headerData;
  std::string tmpStr;

  int64_t dicomHeaderPos = -1;

  //Start at the end of the file and search backwards for DICOM header start
  for (int64_t c = endOfFile; c > endOfFile - 50000; c--) {
        //TODO: probably makes more sense to step back in 132 byte chuncks.
    fseeko64( inFile, c, SEEK_SET );
    headerData.insert( headerData.begin(), fgetc(inFile) );

    tmpStr.clear();
    tmpStr= std::string( headerData.begin(), headerData.begin()+target.length() );

    if (tmpStr == target) {
      LOG(INFO) << "Found DICOM header at: " << c << " bytes";
      dicomHeaderPos = c;
      break;
    }
  }

    //Don't need input file any more.
  fclose(inFile);

  if ( dicomHeaderPos == -1 ) {
    LOG(INFO) << "No DICOM header found";
    return FileStatusCode::EBAD;
  }

    //Try and find Interfile header
  std::string interfileHeader = std::string( headerData.begin(), headerData.end() );
  target.clear();
  target = "!INTERFILE";

  if ( interfileHeader.find( target ) == std::string::npos ) {
    LOG(INFO) << "No Interfile header found";
    return FileStatusCode::EBAD;
  }

  std::size_t inPoint = interfileHeader.find(target);

  target.clear();
  target = "%comment";

  if ( interfileHeader.find( target ) == std::string::npos ) {
    LOG(INFO) << "No end of Interfile header found";
    return FileStatusCode::EBAD;
  }

    //Find the last line of the Interfile header, then dump everything else.
  std::string lastline = interfileHeader.substr( interfileHeader.find(target), interfileHeader.length());
  lastline = lastline.substr(0,lastline.find("\n"));
  interfileHeader = interfileHeader.substr(0,interfileHeader.find(target))+lastline;
  interfileHeader = interfileHeader.substr( inPoint, interfileHeader.length() );

    //Find total no. of words expected in listmode.
  target = "%total listmode word counts";
  if ( interfileHeader.find( target ) == std::string::npos ) {
    LOG(INFO) << "No word count found in Interfile header";
    return FileStatusCode::EBAD;
  }

  lastline = interfileHeader.substr(interfileHeader.find(target),interfileHeader.length());
  lastline = lastline.substr(0,lastline.find("\n"));

  boost::regex pattern("[0-9]+");
  boost::smatch result;

  uint64_t expectedNoWords = 0;
  if (boost::regex_search(lastline, result, pattern)) {
    expectedNoWords = boost::lexical_cast<unsigned long>(result);
  }
  else {
    LOG(INFO) << "No word count found in Interfile header";
    return FileStatusCode::EBAD;
  }

  LOG(INFO) << "Expected number of LM words: " << expectedNoWords;

    //DICM is prefixed with 128 bytes of 0x00
    //Length of the (32-bit) listmode data must always be divisible by 4!
  if ( (dicomHeaderPos - 128) % 4 != 0 ) {
    LOG(INFO) << (dicomHeaderPos -128) / 4 << " words found";
    LOG(INFO) << "Incorrect number of bytes";
    return FileStatusCode::EBAD;
  }

  uint64_t actualNoWords = (dicomHeaderPos - 128) / 4;
  LOG(INFO) << actualNoWords << " LM words found";

  if ( actualNoWords != expectedNoWords ) {
    LOG(INFO) << "Expected no. of LM words does not equal no. read!";
    return FileStatusCode::EBAD;
  }

  return FileStatusCode::EGOOD;
}

}// namespace nmtools
#endif
