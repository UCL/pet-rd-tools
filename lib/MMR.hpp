/*
   MMR.hpp

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

   Classes for reading and modifying Siemens mMR raw data.

 */

// http://ikeptwalking.com/factory-method-pattern/

#ifndef MMR_HPP
#define MMR_HPP

#include <memory>

#include <gdcmReader.h>
#include <boost/filesystem.hpp>
#include <glog/logging.h>
#include <boost/regex.hpp>

#include "Common.hpp"

namespace nmtools {

class IMMR {

public:
  bool SetInputFile ( boost::filesystem::path src );
  FileType GetFileType( boost::filesystem::path src );

  virtual bool ExtractHeader( const boost::filesystem::path dst );
  virtual bool ExtractData( const boost::filesystem::path dst )=0;
  virtual bool ModifyHeader( const boost::filesystem::path src, const boost::filesystem::path dataFile)=0;

protected:

  bool ReadHeader();

  virtual boost::filesystem::path GetStdFileName( boost::filesystem::path srcFile, ContentType ctype) = 0;
  FileStatusCode CheckForSiemensBFFile(boost::filesystem::path src, uint64_t numOfWords);
  gdcm::Reader * _dicomReader;
  std::string _headerString;

  boost::filesystem::path _srcPath;

};

class MMR32BitList : public IMMR {
public:
  bool ExtractData( const boost::filesystem::path dst );
protected:
  boost::filesystem::path GetStdFileName( boost::filesystem::path srcFile, ContentType ctype);
};

class MMRSino : public IMMR {

protected:
  boost::filesystem::path GetStdFileName( boost::filesystem::path srcFile, ContentType ctype);
};

class MMRNorm : public IMMR {


protected:
  boost::filesystem::path GetStdFileName( boost::filesystem::path srcFile, ContentType ctype);
  std::string cleanUpLineEncoding( std::string );
};

class IRawDataFactory {
public:
  virtual IMMR* Create( boost::filesystem::path inFile )=0;
};

class MMRFactory : public IRawDataFactory{
public:
  IMMR* Create( boost::filesystem::path inFile ){

    std::unique_ptr<IMMR> rawData(new IMMR);

    FileType fType = rawData->GetFileType( inFile );

    if (fType == FileType::EMMRLIST)
      return new MMR32BitList();

    if (fType == FileType::EMMRSINO)
      return new MMRSino(); 

    if (fType == FileType::EMMRNORM)
      return new MMRNorm();      

    return nullptr;
  }
};


FileType IMMR::GetFileType( boost::filesystem::path src){

  FileType foundFileType = FileType::EUNKNOWN;

  std::unique_ptr<gdcm::Reader> dicomReader(new gdcm::Reader);
  dicomReader->SetFileName(src.string().c_str());

  if (!dicomReader->Read()) {
    LOG(INFO) << "Unable to read as DICOM file";
    return FileType::EERROR;
  }

  const gdcm::DataSet &ds = dicomReader->GetFile().GetDataSet();

  const gdcm::Tag manufacturer(0x008, 0x0070);
  const gdcm::DataElement &pdde = ds.GetDataElement(manufacturer);

  std::stringstream manufacturerName;
  manufacturerName << (char *) pdde.GetByteValue()->GetPointer();

  LOG(INFO) << "Manufacturer: " << manufacturerName.str();

  const gdcm::Tag model(0x008, 0x1090);
  const gdcm::DataElement &pdde2 = ds.GetDataElement(model);

  std::stringstream modelName;
  modelName << (char *) pdde2.GetByteValue()->GetPointer();

  LOG(INFO) << "Model name: " << modelName.str();

  const gdcm::Tag imageType(0x0008, 0x0008);
  const gdcm::DataElement &pdde3 = ds.GetDataElement(imageType);

  std::stringstream imageTypeValue;
  imageTypeValue << (char *) pdde3.GetByteValue()->GetPointer();

  LOG(INFO) << "Image type: " << imageTypeValue.str();

  if (manufacturerName.str().find("SIEMENS") != std::string::npos) {
    DLOG(INFO) << "Manufacturer = SIEMENS";

    if (modelName.str().find("Biograph_mMR") != std::string::npos) {
      DLOG(INFO) << "Scanner = MMR";

      if (imageTypeValue.str().find("ORIGINAL\\PRIMARY\\PET_LISTMODE") != std::string::npos)
        foundFileType = FileType::EMMRLIST;
      if (imageTypeValue.str().find("ORIGINAL\\PRIMARY\\PET_EM_SINO") != std::string::npos)
        foundFileType = FileType::EMMRSINO;
      if (imageTypeValue.str().find("ORIGINAL\\PRIMARY\\PET_NORM") != std::string::npos)
        foundFileType = FileType::EMMRNORM;
    }
  }

  return foundFileType;
}

bool IMMR::ReadHeader() {

  const gdcm::DataSet &ds = _dicomReader->GetFile().GetDataSet();

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

  _headerString = headerString.str();

  if (_headerString.size() > 0)
    return true;

  return false;
}

bool IMMR::ExtractHeader( const boost::filesystem::path dst ){

  bool bStatus = false;

  const gdcm::DataSet &ds = _dicomReader->GetFile().GetDataSet();

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

  std::ofstream outfile(dst.string().c_str(), std::ios::out | std::ios::binary);
  if (!outfile.is_open()) {
    LOG(ERROR) << "Unable to write header to " << dst;
    return false;
  }
  else {
    outfile << headerString.str();
    outfile.close();
    bStatus = true;
  }

  if (bStatus == true)
    LOG(INFO) << "Successfully extracted raw header.";
  else
    LOG(ERROR) << "Failed to extract raw header!";

  return bStatus;
}

bool MMR32BitList::ExtractData( const boost::filesystem::path dst ){

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
      
      fs::path bfPath = _srcPath;
      bfPath = boost::filesystem::change_extension(bfPath, ".bf").string();
      try {
        if (fs::exists(dst)) {
          DLOG(ERROR) << "The .bf file already exists!";
          return false;
        }
        
        fs::copy(bfPath, dst);
        bStatus = true;
      }
      catch(fs::filesystem_error const &e){
        DLOG(ERROR) << "Unable to copy listmode from .bf file!";
        DLOG(ERROR) << e.what();
        return false;
      }

    }
    else {
      DLOG(ERROR) << "No listmode data found in either header or .bf file!";
      return false;
    }
  } else {
    std::ofstream outfile(dst.string().c_str(), std::ios::out | std::ios::binary);
    if (!outfile.is_open()) {
      DLOG(ERROR) << "Unable to write listmode to " << dst;
      return false;
    }
    bv->WriteBuffer(outfile);
    outfile.close();
    bStatus = true;
  }

  return bStatus;
}

bool MMRNorm::ExtractData( const boost::filesystem::path dst ){

  namespace fs = boost::filesystem;

  bool bStatus = false;

  if (!this->ReadHeader()){
    LOG(ERROR) << "Unable to read header!";
    return false;
  }

  const gdcm::DataSet &ds = _dicomReader->GetFile().GetDataSet();


  return bStatus;
}


boost::filesystem::path MMR32BitList::GetStdFileName( boost::filesystem::path srcFile, ContentType ctype){

  boost::filesystem::path outputPath = srcFile.filename().stem();
  outputPath += ".l";

  if (ctype == ContentType::EHEADER)
    outputPath += ".hdr";

  DLOG(INFO) << "Created filename: " << outputPath;
  return outputPath;
}

boost::filesystem::path MMRSino::GetStdFileName( boost::filesystem::path srcFile, ContentType ctype){

  boost::filesystem::path outputPath = srcFile.filename().stem();
  outputPath += ".s";

  if (ctype == ContentType::EHEADER)
    outputPath += ".hdr";

  DLOG(INFO) << "Created filename: " << outputPath;
  return outputPath;
}

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
