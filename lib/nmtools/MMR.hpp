/*
   MMR.hpp

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

   Classes for reading and modifying Siemens mMR raw data.

 */

#ifndef MMR_HPP
#define MMR_HPP

#include <memory>

#include <gdcmReader.h>
#include <boost/filesystem.hpp>
#include <glog/logging.h>
#include <boost/regex.hpp>

#include "Common.hpp"

namespace nmtools {

class IMMR: public IDicomExtractor {

//Base class for mMR. Splits Interfile header and raw data
public:

  explicit IMMR(boost::filesystem::path src);

    //FileType GetFileType( boost::filesystem::path src );

  virtual bool ExtractHeader( const boost::filesystem::path dst );
  bool ModifyHeader( const boost::filesystem::path src, const boost::filesystem::path dataFile);
  //Deal with EOF in norm header.
  std::string CleanUpLineEncoding( std::string );

  virtual ~IMMR(){};

protected:

  bool ReadHeader();

  FileStatusCode CheckForSiemensBFFile(boost::filesystem::path src, uint64_t numOfWords);
  std::string _headerString;

};

class MMR32BitList : public IMMR {
//Derived class for handling 32-bit list mode data.

  using IMMR::IMMR;
public:
  bool IsValid();
  bool ExtractData( const boost::filesystem::path dst );
  boost::filesystem::path GetStdFileName( boost::filesystem::path srcFile, ContentType ctype);

};

class MMRSino : public IMMR {
//Derived class for handling sinogram data.

  using IMMR::IMMR;
public:
  bool IsValid();
  bool ExtractData( const boost::filesystem::path dst );
  boost::filesystem::path GetStdFileName( boost::filesystem::path srcFile, ContentType ctype);
};

class MMRNorm : public IMMR {
//Derived class for handling normalisation files.

  using IMMR::IMMR;
public:
  //MMR Norm file byte length
  //({344,127}+{9,344}+{504,64}+{837}+{64}+{64}+{9}+{837}) * 4
  const uint32_t MMRNORMBYTELENGTH = 323404;

  bool IsValid();
  bool ExtractData( const boost::filesystem::path dst );
  bool ModifyHeader( const boost::filesystem::path src, const boost::filesystem::path dataFile);
  boost::filesystem::path GetStdFileName( boost::filesystem::path srcFile, ContentType ctype);
protected:
};

class SiemensPETFactory : public IRawDataFactory{
public:
  enum class FileType { EMMRSINO, EMMRLIST, EMMRNORM,
                        EUNKNOWN, EERROR };

  FileType GetFileType( boost::filesystem::path src){

    //Extracts information via DICOM to determine what kind of mMR
    //raw data type we're dealing with.
    //
    //Will check for list mode, sinograms and norms.
    //
    //TODO: Support physio files? (11/12/2017)

    FileType foundFileType = FileType::EUNKNOWN;

    if (!Open(src)) {
      return FileType::EERROR;
    }

    const gdcm::File &file = dicomReader->GetFile();

    if (manufacturerName.find("SIEMENS") != std::string::npos) {
      DLOG(INFO) << "Manufacturer = SIEMENS";

      //Get image tpye description
      const gdcm::Tag imageType(0x0008, 0x0008);
      std::string imageTypeValue;
      if (!GetTagInfo(file,imageType,imageTypeValue)){
        LOG(ERROR) << "Unable to image type!";
        return FileType::EERROR;  
      }
      LOG(INFO) << "Image type: " << imageTypeValue;

      if (modelName.find("Biograph_mMR") != std::string::npos) {
        DLOG(INFO) << "Scanner = MMR";

        if (imageTypeValue.find("ORIGINAL\\PRIMARY\\PET_LISTMODE") != std::string::npos)
          foundFileType = FileType::EMMRLIST;
        if (imageTypeValue.find("ORIGINAL\\PRIMARY\\PET_EM_SINO") != std::string::npos)
          foundFileType = FileType::EMMRSINO;
        if (imageTypeValue.find("ORIGINAL\\PRIMARY\\PET_NORM") != std::string::npos)
          foundFileType = FileType::EMMRNORM;
      }
    }
    return foundFileType;
  }  
private:
  IMMR* Create_ptr( boost::filesystem::path inFile ) {

    FileType fType = GetFileType( inFile );

    if (fType == FileType::EMMRLIST){
      IMMR* instance(new MMR32BitList(inFile));
      return instance;
    }

    if (fType == FileType::EMMRSINO){
      IMMR* instance(new MMRSino(inFile)); 
      return instance;
    }

    if (fType == FileType::EMMRNORM){
      IMMR* instance(new MMRNorm(inFile));  
      return instance;    
    }

    if (fType == FileType::EUNKNOWN){
      LOG(ERROR) << "Unsupported file type (only handling list/sino/norm)";
    }
    return nullptr;
  }
};

//Constructor with path.
IMMR::IMMR(boost::filesystem::path src)
  : IDicomExtractor(src)
{}

//Header extraction
bool IMMR::ReadHeader() {

  if (!_dicomReader) {
      LOG(ERROR) << "DICOM reader not initialised. Internal error.";
      return false;
  }

  const gdcm::File &file = _dicomReader->GetFile();

  const gdcm::Tag headerTag(0x029, 0x1010);

  std::string headerString;
  std::string headerStringTmp;

  if (!GetTagInfo(file,headerTag,headerStringTmp)){
      LOG(WARNING) << "Unable to read header from " << headerTag;
    //return false;  
  }
  
  //If this is actually a SMS-MI v 3.2 file then get header from 0029,1110.
  if ( (headerStringTmp.find("SV10") != std::string::npos) || (headerStringTmp.size() == 0) ) {
    gdcm::Tag altHeaderTag(0x029, 0x1110);
    if (!GetTagInfo(file,altHeaderTag,headerString)){
      LOG(ERROR) << "Unable to read header (SV10)";
      return false;  
    }
  }
  else {
    headerString = headerStringTmp;
  }

  _headerString = headerString;

  if (_headerString.size() > 0)
    return true;

  return false;
}

//Header extraction and writing to file dst.
bool IMMR::ExtractHeader( const boost::filesystem::path dst ){

  bool bStatus = false;

  const gdcm::File &file = _dicomReader->GetFile();
  //const gdcm::DataSet &ds = _dicomReader->GetFile().GetDataSet();

  const gdcm::Tag headerTag(0x029, 0x1010);

  std::string headerString;
  std::string headerStringTmp;

  if (!GetTagInfo(file,headerTag,headerStringTmp)){
    LOG(ERROR) << "Unable to header from " << headerTag;
    return false;  
  }
  //If this is actually a SMS-MI v 3.2 file then get header from 0029,1110.
  if (headerStringTmp.find("SV10") != std::string::npos) {
    gdcm::Tag altHeaderTag(0x029, 0x1110);
    if (!GetTagInfo(file,altHeaderTag,headerString)){
      LOG(ERROR) << "Unable to header (SV10) from " << altHeaderTag;
      return false;  
    }
  }
  else {
    headerString = headerStringTmp;
  }

  if (boost::filesystem::exists(dst)) {
      LOG(ERROR) << "Header already exists at destination!";
      LOG(ERROR) << "Refusing to over-write!";
      return false;
  }

  std::ofstream outfile(dst.string().c_str(), std::ios::out | std::ios::binary);
  if (!outfile.is_open()) {
    LOG(ERROR) << "Unable to write header to " << dst;
    return false;
  }
  else {
    outfile << headerString;
    outfile.close();
    bStatus = true;
  }

  if (bStatus == true)
    LOG(INFO) << "Successfully extracted raw header.";
  else
    LOG(ERROR) << "Failed to extract raw header!";

  return bStatus;
}

FileStatusCode IMMR::CheckForSiemensBFFile(boost::filesystem::path src, uint64_t numOfBytes) {
  //Test for existence of associated bf file and if the length is correct
  //according to the Interfile header.

  FILE *inFile;

  boost::filesystem::path bfPath = src;
  bfPath.replace_extension(".bf");

  if (!(inFile = fopen(bfPath.string().c_str(), "rb"))) {
      LOG(INFO) << "Cannot open " << bfPath.native();
      return FileStatusCode::EIOERROR;
  }

  fseeko64(inFile, 0L, SEEK_END);
  int64_t endOfFile = ftello64(inFile);

  fclose(inFile);

  LOG(INFO) << ".bf file size in bytes: " << endOfFile;
  //uint64_t actualNoWords = endOfFile / 4;
  //BOOST_LOG_TRIVIAL(info) << endOfFile << " / 4 = " << actualNoWords << " words";

  if (endOfFile != numOfBytes) {
    LOG(INFO) << "Expected no. of bytes does not equal no. read!";
    return FileStatusCode::EBAD;
  }

  LOG(INFO) << bfPath << " is valid raw data file for this header.";
  return FileStatusCode::EGOOD;

}

//Extract raw data and write to dst.
bool MMR32BitList::ExtractData( const boost::filesystem::path dst ){

  namespace fs = boost::filesystem;

  bool bStatus = false;

  if (!this->ReadHeader()){
    LOG(ERROR) << "Unable to read header!";
    return false;
  }

  if (boost::filesystem::exists(dst)) {
    LOG(ERROR) << "The data file already exists!";
    LOG(ERROR) << "Refusing to over-write!";
    return false;
  }

  const gdcm::DataSet &ds = _dicomReader->GetFile().GetDataSet();

  std::string target = "%total listmode word counts";
  if (_headerString.find(target) == std::string::npos) {
    LOG(INFO) << "No word count tag found in Interfile header";
    return false;
  }

  std::string lastline = _headerString.substr(_headerString.find(target), _headerString.length());
  lastline = lastline.substr(0, lastline.find("\n"));

  boost::regex pattern("[0-9]+");
  boost::smatch result;

  uint64_t expectedNoWords = 0;
  if (boost::regex_search(lastline, result, pattern)) {
    expectedNoWords = boost::lexical_cast<unsigned long>(result);
  } else {
    LOG(INFO) << "No word count number found in Interfile header";
    return false;
  }

  LOG(INFO) << "Expected number of LM words: " << expectedNoWords;

  const gdcm::Tag lmDataTag(0x7fe1, 0x1010);
  const gdcm::DataElement &lmData = ds.GetDataElement(lmDataTag);
  const gdcm::ByteValue *bv = lmData.GetByteValue();

  uint64_t lmLength = bv->GetLength();
  LOG(INFO) << lmLength << " bytes in LM field";

  uint64_t actualNoWords = lmLength / 4;

  LOG(INFO) << lmLength << " / 4 = " << actualNoWords << " words";

  if (lmLength != expectedNoWords * 4) {
    LOG(INFO) << "Expected no. of LM words does not equal no. read!";
    LOG(INFO) << "Looking for BF file...";

    DLOG(INFO) << "SRC: " << this->_srcPath;
    FileStatusCode bfStatus = CheckForSiemensBFFile(this->_srcPath, expectedNoWords*4);

    if ( bfStatus == FileStatusCode::EGOOD ) {
      
      boost::filesystem::path bfPath = _srcPath;
      bfPath.replace_extension(".bf").string();
      try {
        if (boost::filesystem::exists(dst)) {
          LOG(ERROR) << "The .bf file already exists!";
          return false;
        }
        
        boost::filesystem::copy(bfPath, dst);
        bStatus = true;
      }
      catch(boost::filesystem::filesystem_error const &e){
        LOG(ERROR) << "Unable to copy listmode from .bf file!";
        LOG(ERROR) << e.what();
        return false;
      }

    }
    else {
      LOG(ERROR) << "No listmode data found in either header or .bf file!";
      return false;
    }
  } else {
    std::ofstream outfile(dst.string().c_str(), std::ios::out | std::ios::binary);
    if (!outfile.is_open()) {
      LOG(ERROR) << "Unable to write listmode to " << dst;
      return false;
    }
    bv->WriteBuffer(outfile);
    outfile.close();
    bStatus = true;
  }

  return bStatus;
}

//Check if mMR list mode file is valid.
bool MMR32BitList::IsValid(){

  namespace fs = boost::filesystem;

  bool bStatus = false;

  if (!this->ReadHeader()){
    LOG(ERROR) << "Unable to read header!";
    return false;
  }

  const gdcm::DataSet &ds = _dicomReader->GetFile().GetDataSet();

  std::string target = "%total listmode word counts";
  if (_headerString.find(target) == std::string::npos) {
    LOG(INFO) << "No word count tag found in Interfile header";
    return false;
  }

  std::string lastline = _headerString.substr(_headerString.find(target), _headerString.length());
  lastline = lastline.substr(0, lastline.find("\n"));

  boost::regex pattern("[0-9]+");
  boost::smatch result;

  uint64_t expectedNoWords = 0;
  if (boost::regex_search(lastline, result, pattern)) {
    expectedNoWords = boost::lexical_cast<unsigned long>(result);
  } else {
    LOG(INFO) << "No word count number found in Interfile header";
    return false;
  }

  LOG(INFO) << "Expected number of LM words: " << expectedNoWords;

  const gdcm::Tag lmDataTag(0x7fe1, 0x1010);
  const gdcm::DataElement &lmData = ds.GetDataElement(lmDataTag);
  const gdcm::ByteValue *bv = lmData.GetByteValue();

  uint64_t lmLength = bv->GetLength();
  LOG(INFO) << lmLength << " bytes in LM field";

  uint64_t actualNoWords = lmLength / 4;

  LOG(INFO) << lmLength << " / 4 = " << actualNoWords << " words";

  if (lmLength != expectedNoWords * 4) {
    LOG(INFO) << "Expected no. of LM words does not equal no. read!";
    LOG(INFO) << "Looking for BF file...";

    DLOG(INFO) << "SRC: " << this->_srcPath;
    FileStatusCode bfStatus = CheckForSiemensBFFile(this->_srcPath, expectedNoWords*4);

    if ( bfStatus == FileStatusCode::EGOOD ) {
      return true;
    }
    else {
      LOG(ERROR) << "No listmode data found in either header or .bf file!";
      return false;
    }
  } 
  else {
    bStatus = true;
  }

  return bStatus;
}

//Extract sinogram data.
bool MMRSino::ExtractData( const boost::filesystem::path dst ){

  namespace fs = boost::filesystem;

  bool bStatus = false;

  if (!this->ReadHeader()){
    LOG(ERROR) << "Unable to read header!";
    return false;
  }

  const gdcm::DataSet &ds = _dicomReader->GetFile().GetDataSet();

  const gdcm::Tag lmDataTag(0x7fe1, 0x1010);
  const gdcm::DataElement &lmData = ds.GetDataElement(lmDataTag);
  const gdcm::ByteValue *bv = lmData.GetByteValue();

  uint64_t lmLength = bv->GetLength();
  LOG(INFO) << lmLength << " bytes in data field (0x7fe1, 0x1010)";

  DLOG(INFO) << "SRC: " << this->_srcPath;

  if (boost::filesystem::exists(dst)) {
    LOG(ERROR) << "The data file already exists!";
    LOG(ERROR) << "Refusing to over-write!";
    return false;
  }

  boost::filesystem::path bfPath = _srcPath;
  bfPath.replace_extension(".bf");

  if ( boost::filesystem::exists(bfPath) ){

    boost::filesystem::path bfPath = _srcPath;
    bfPath.replace_extension(".bf").string();
    try {
      boost::filesystem::copy(bfPath, dst);
      bStatus = true;
    }
    catch(boost::filesystem::filesystem_error const &e){
      LOG(ERROR) << "Unable to copy sinogram from .bf file!";
      LOG(ERROR) << e.what();
      return false;
    }
  }
  else {
    std::ofstream outfile(dst.string().c_str(), std::ios::out | std::ios::binary);
    if (!outfile.is_open()) {
      LOG(ERROR) << "Unable to write sinogram to " << dst;
      return false;
    }
    bv->WriteBuffer(outfile);
    outfile.close();
    bStatus = true;
  }

  return bStatus;
}

//Check if sinogram is valid.
bool MMRSino::IsValid(){

  namespace fs = boost::filesystem;

  bool bStatus = false;

  if (!this->ReadHeader()){
    LOG(ERROR) << "Unable to read header!";
    return false;
  }

  //Due to compression, this isn't a real check.
  //TODO: check for span-1 at least (11/12/2017).
  LOG(WARNING) << "Cannot check sinogram length due to compression.";

  const gdcm::DataSet &ds = _dicomReader->GetFile().GetDataSet();

  const gdcm::Tag lmDataTag(0x7fe1, 0x1010);
  const gdcm::DataElement &lmData = ds.GetDataElement(lmDataTag);
  const gdcm::ByteValue *bv = lmData.GetByteValue();

  uint64_t lmLength = bv->GetLength();
  LOG(INFO) << lmLength << " bytes in data field (0x7fe1, 0x1010)";
  DLOG(INFO) << "SRC: " << this->_srcPath;

  boost::filesystem::path bfPath = _srcPath;
  bfPath.replace_extension(".bf").string();

  if ( boost::filesystem::exists(bfPath) ){
    LOG(INFO) << ".bf file exists.";
    return true;
  }
  else {
    if (lmLength == 0)
      bStatus = false;
    else
      bStatus = true;
  }

  return bStatus;
}

//Extract norm raw data.
bool MMRNorm::ExtractData( const boost::filesystem::path dst ){

  namespace fs = boost::filesystem;

  bool bStatus = false;

  if (!this->ReadHeader()){
    LOG(ERROR) << "Unable to read header!";
    return false;
  }

  const gdcm::DataSet &ds = _dicomReader->GetFile().GetDataSet();

  LOG(INFO) << "Expected number of bytes: " << MMRNORMBYTELENGTH;

  const gdcm::Tag lmDataTag(0x7fe1, 0x1010);
  const gdcm::DataElement &lmData = ds.GetDataElement(lmDataTag);
  const gdcm::ByteValue *bv = lmData.GetByteValue();

  uint64_t lmLength = bv->GetLength();
  LOG(INFO) << lmLength << " bytes in data field (0x7fe1, 0x1010)";

  if (boost::filesystem::exists(dst)) {
    LOG(ERROR) << "The data file already exists!";
    LOG(ERROR) << "Refusing to over-write!";
    return false;
  }

  if (lmLength != MMRNORMBYTELENGTH) {
    LOG(INFO) << "Expected no. of bytes does not equal no. read!";
    LOG(INFO) << "Looking for BF file...";

    DLOG(INFO) << "SRC: " << this->_srcPath;
    FileStatusCode bfStatus = CheckForSiemensBFFile(this->_srcPath, MMRNORMBYTELENGTH);

    if ( bfStatus == FileStatusCode::EGOOD ) {

      boost::filesystem::path bfPath = _srcPath;
      bfPath.replace_extension(".bf").string();
      try {
        boost::filesystem::copy(bfPath, dst);
        bStatus = true;
      }
      catch(boost::filesystem::filesystem_error const &e){
        LOG(ERROR) << "Unable to copy norm from .bf file!";
        LOG(ERROR) << e.what();
        return false;
      }
    }
    else {
      LOG(ERROR) << "No norm data found in either header or .bf file!";
      return false;
    }
  } else {
    std::ofstream outfile(dst.string().c_str(), std::ios::out | std::ios::binary);

    if (!outfile.is_open()) {
      LOG(ERROR) << "Unable to write norm to " << dst;
      return false;
    }
    bv->WriteBuffer(outfile);
    outfile.close();
    bStatus = true;
  }

  return bStatus;
}

//Check if norm is valid.
bool MMRNorm::IsValid(){

  namespace fs = boost::filesystem;

  bool bStatus = false;

  if (!this->ReadHeader()){
    LOG(ERROR) << "Unable to read header!";
    return false;
  }

  const gdcm::DataSet &ds = _dicomReader->GetFile().GetDataSet();

  LOG(INFO) << "Expected number of bytes: " << MMRNORMBYTELENGTH;

  const gdcm::Tag lmDataTag(0x7fe1, 0x1010);
  const gdcm::DataElement &lmData = ds.GetDataElement(lmDataTag);
  const gdcm::ByteValue *bv = lmData.GetByteValue();

  uint64_t lmLength = bv->GetLength();
  LOG(INFO) << lmLength << " bytes in data field (0x7fe1, 0x1010)";

  if (lmLength != MMRNORMBYTELENGTH) {
    LOG(INFO) << "Expected no. of bytes does not equal no. read!";
    LOG(INFO) << "Looking for BF file...";

    DLOG(INFO) << "SRC: " << this->_srcPath;
    FileStatusCode bfStatus = CheckForSiemensBFFile(this->_srcPath, MMRNORMBYTELENGTH);

    if ( bfStatus == FileStatusCode::EGOOD ) {
      return true;
    }
    else {
      LOG(ERROR) << "No norm data found in either header or .bf file!";
      return false;
    }
  } 
  else {
    bStatus = true;
  }

  return bStatus;
}

//Re-write list mode data file location in header.
bool IMMR::ModifyHeader(const boost::filesystem::path src, const boost::filesystem::path dataFile){

  std::ifstream headerFile( boost::filesystem::canonical( src ).string().c_str() );
  std::stringstream buffer;
  buffer << headerFile.rdbuf();

  headerFile.close();

  DLOG(INFO) << "Read " << src;

  std::string headerInfo = buffer.str();

  std::string::size_type pos = 0;
  std::string target = "name of data file";

  pos = headerInfo.find(target);
  std::string line = headerInfo.substr(pos,headerInfo.length());
  std::string::size_type end = std::min(line.find("\n"), line.find("\r"));
  line = line.substr(0, end);

  std::string newLine = "name of data file:=" + dataFile.filename().string();
  headerInfo.replace(pos,line.length(), newLine);

  std::ofstream outfile( src.string().c_str(), std::ios::out | std::ios::binary);

  if ( ! outfile.is_open() ) {
      LOG(ERROR) << "Unable to update header in " << src;
      return false;
  }

  //Find /r/r/n and replace with /r/n.
  headerInfo = CleanUpLineEncoding( headerInfo );

  outfile << headerInfo;
  outfile.close();

  return true;

}

//Re-write norm data file location in header.
bool MMRNorm::ModifyHeader(const boost::filesystem::path src, const boost::filesystem::path dataFile){

  std::ifstream headerFile( boost::filesystem::canonical( src ).string().c_str() );
  std::stringstream buffer;
  buffer << headerFile.rdbuf();
  headerFile.close();

  std::string headerInfo = buffer.str();

  std::string::size_type pos = 0;
  std::string target = "%data set [1]:={0,,";

  pos = headerInfo.find(target);
  std::string line = headerInfo.substr(pos,headerInfo.length());
  std::string::size_type end = std::min(line.find("\n"), line.find("\r"));
  line = line.substr(0, end);

  std::string newLine = "%data set [1]:={0,," +  dataFile.filename().string() + "}";
  headerInfo.replace(pos,line.length(), newLine);

  std::ofstream outfile( src.string().c_str(), std::ios::out | std::ios::binary);

  if ( ! outfile.is_open() ) {
      LOG(ERROR) << "Unable to update norm header in " << src;
      return false;
  }

  //Find /r/r/n and replace with /r/n.
  headerInfo = CleanUpLineEncoding( headerInfo );

  outfile << headerInfo;
  outfile.close();

  // Also need to run the superclass to update name of data file and run cleanup
  IMMR::ModifyHeader(src, dataFile);

  return true;

}

//Removes \r\r\n line endings in mMR norm header and replaces it with \r\n
std::string IMMR::CleanUpLineEncoding( const std::string origStr ){

  std::stringstream ss;
  ss << origStr;
  std::string outstr;

  std::string line;
  std::string::size_type pos;

  while (std::getline(ss, line)) {
      // \r\r -> \r
      pos = line.find("\r\r");
      if (pos != std::string::npos) {
            DLOG(INFO) << "fixed: " << line;
            line.replace(pos,1,"");
      } else {
            DLOG(INFO) << "left:  " << line;
      }

      // no \r -> \r
      pos = line.find("\r");
      if (pos == std::string::npos) {
            outstr += line + "\r\n";
      } else {
            outstr += line + "\n";
      }
  }

  //Add carriage return at EOF
  outstr += "\r\n";

  return outstr;
}

//Create destination filename for list mode.
boost::filesystem::path MMR32BitList::GetStdFileName( boost::filesystem::path srcFile, ContentType ctype){

  boost::filesystem::path outputPath = srcFile.filename().stem();
  outputPath += ".l";

  if (ctype == ContentType::EHEADER)
    outputPath += ".hdr";

  DLOG(INFO) << "Created filename: " << outputPath;
  return outputPath;
}

//Create destination filename for sinogram.
boost::filesystem::path MMRSino::GetStdFileName( boost::filesystem::path srcFile, ContentType ctype){

  boost::filesystem::path outputPath = srcFile.filename().stem();
  outputPath += ".s";

  if (ctype == ContentType::EHEADER)
    outputPath += ".hdr";

  DLOG(INFO) << "Created filename: " << outputPath;
  return outputPath;
}

//Create destination filename for norm.
boost::filesystem::path MMRNorm::GetStdFileName( boost::filesystem::path srcFile, ContentType ctype){

  boost::filesystem::path outputPath = srcFile.filename().stem();
  outputPath += ".n";

  if (ctype == ContentType::EHEADER)
    outputPath += ".hdr";

  DLOG(INFO) << "Created filename: " << outputPath;
  return outputPath;
}

} // namespace nmtools

#endif 
