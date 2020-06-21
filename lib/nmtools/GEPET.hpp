/*
   GEPET.hpp

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

#ifndef GEPET_HPP
#define GEPET_HPP

#include <memory>

#include <gdcmReader.h>
#include <boost/filesystem.hpp>
#include <glog/logging.h>
#include <boost/regex.hpp>

#include "Common.hpp"

namespace nmtools {

class IGEPET : public IDicomExtractor {


public:

  explicit IGEPET(boost::filesystem::path src);

  //FileType GetFileType( boost::filesystem::path src );

  virtual bool IsValid() { return true; }

  virtual bool ExtractHeader( const boost::filesystem::path dst )
  { return ExtractRDF( dst ); }
  virtual bool ExtractData( const boost::filesystem::path dst )
  { return true; }
  virtual bool ModifyHeader( const boost::filesystem::path src, const boost::filesystem::path dataFile)
  { return true; }
  virtual boost::filesystem::path GetStdFileName( boost::filesystem::path srcFile, ContentType ctype)
  {
    // GE doesn't have header/raw data, it's all just one RDF file
    if (ctype == ContentType::EHEADER)
      return GetStdFileName(srcFile);
    else
      return boost::filesystem::path();
  }
  // GE specific
  virtual bool ExtractRDF( const boost::filesystem::path dst );
  virtual boost::filesystem::path GetStdFileName( boost::filesystem::path srcFile) = 0;

protected:

  //Extract raw data.
  bool ExtractBlob( const boost::filesystem::path dst, const gdcm::Tag DataTag );
};

class GEPETList : public IGEPET {
//Derived class for handling list mode data.

  using IGEPET::IGEPET;
public:
  //bool IsValid();
  boost::filesystem::path GetStdFileName( boost::filesystem::path srcFile);

};

class GEPETSino : public IGEPET {
//Derived class for handling sinogram data.

  using IGEPET::IGEPET;
public:
  //bool IsValid();
  boost::filesystem::path GetStdFileName( boost::filesystem::path srcFile);
};

class GEPETNorm : public IGEPET {
//Derived class for handling normalisation files.

  using IGEPET::IGEPET;
public:
  //bool IsValid();
  boost::filesystem::path GetStdFileName( boost::filesystem::path srcFile);
};

class GEPETGeo : public IGEPET {
//Derived class for handling normalisation files.

  using IGEPET::IGEPET;
public:
  //bool IsValid();
  boost::filesystem::path GetStdFileName( boost::filesystem::path srcFile);
};

//Constructor with path.
IGEPET::IGEPET(boost::filesystem::path src)
  : IDicomExtractor(src)
{}

//Extract raw data.
bool IGEPET::ExtractBlob( const boost::filesystem::path dst, const gdcm::Tag DataTag ){

  namespace fs = boost::filesystem;

  bool bStatus = false;

  const gdcm::DataSet &ds = _dicomReader->GetFile().GetDataSet();

  const gdcm::DataElement &Data = ds.GetDataElement(DataTag);
  const gdcm::ByteValue *bv = Data.GetByteValue();

  uint64_t lmLength = bv->GetLength();
  LOG(INFO) << lmLength << " bytes in data field";

  if (boost::filesystem::exists(dst)) {
    LOG(ERROR) << "The data file already exists!";
    LOG(ERROR) << "Refusing to over-write!";
    return false;
  }

 {
    std::ofstream outfile(dst.string().c_str(), std::ios::out | std::ios::binary);

    if (!outfile.is_open()) {
      LOG(ERROR) << "Unable to write to " << dst;
      return false;
    }
    bv->WriteBuffer(outfile);
    outfile.close();
    bStatus = true;
  }

  return bStatus;
}

//Extract raw data and write to dst.
bool IGEPET::ExtractRDF( const boost::filesystem::path dst ){

  const gdcm::Tag DataTag(0x0023, 0x1002);
  return ExtractBlob(dst, DataTag);
}


//Create destination filename for list mode.
boost::filesystem::path GEPETList::GetStdFileName( boost::filesystem::path srcFile){

  boost::filesystem::path outputPath = srcFile.filename().stem();
  outputPath += ".BLF";

  DLOG(INFO) << "Created filename: " << outputPath;
  return outputPath;
}

//Create destination filename for sinogram.
boost::filesystem::path GEPETSino::GetStdFileName( boost::filesystem::path srcFile){

  boost::filesystem::path outputPath = srcFile.filename().stem();
  outputPath += ".sino.rdf";

  DLOG(INFO) << "Created filename: " << outputPath;
  return outputPath;
}

//Create destination filename for norm.
boost::filesystem::path GEPETGeo::GetStdFileName( boost::filesystem::path srcFile){

  boost::filesystem::path outputPath = srcFile.filename().stem();
  outputPath += "geo.rdf";

  DLOG(INFO) << "Created filename: " << outputPath;
  return outputPath;
}

//Create destination filename for norm.
boost::filesystem::path GEPETNorm::GetStdFileName( boost::filesystem::path srcFile){

  boost::filesystem::path outputPath = srcFile.filename().stem();
  outputPath += ".norm.rdf";

  DLOG(INFO) << "Created filename: " << outputPath;
  return outputPath;
}

class GEPETFactory : public IRawDataFactory {
public:
  enum class FileType { EGEPETCTAC, EGEPETSINO, EGEPETLIST, EGEPETNORM2D, EGEPETNORM3D, EGEPETWCC, EGEPETGEO,
                        EUNKNOWN, EERROR };

  FileType GetFileType( boost::filesystem::path src) {

    //Extracts information via DICOM to determine what kind of raw data type we're dealing with.
    FileType foundFileType = FileType::EUNKNOWN;

    if (!Open(src)) {
        return FileType::EERROR;
    }

    const gdcm::File &file = dicomReader->GetFile();
    
    if (manufacturerName.find("GE MEDICAL SYSTEMS") != std::string::npos) {
        DLOG(INFO) << "Manufacturer = GE";

        const gdcm::Tag rawDataType(0x0021, 0x1001);
        std::string rawDataTypeValue;
        if (!GetTagInfo(file,rawDataType,rawDataTypeValue)){
          LOG(ERROR) << "Unable to type of raw data!";
          return FileType::EERROR;  
        }
        LOG(INFO) << "type of raw data: " << rawDataTypeValue;

        if (rawDataTypeValue.find("3") != std::string::npos) {
          // sino or CTAC in sino format
          const gdcm::Tag sinoType(0x0009, 0x1019);
          std::string sinoTypeValue;
          if (!GetTagInfo(file,sinoType,sinoTypeValue)){
            LOG(ERROR) << "Unable to type of sino data!";
            return FileType::EERROR;  
          }
          LOG(INFO) << "type of sino data: " << sinoTypeValue;
          if (sinoTypeValue.find("0") != std::string::npos) {
            foundFileType = FileType::EGEPETSINO;
          }
          else if (sinoTypeValue.find("5") != std::string::npos) {
            foundFileType = FileType::EGEPETCTAC;
          }
        }
        else if (rawDataTypeValue.find("4") != std::string::npos) {
          // norm file
          const gdcm::Tag calType(0x0017, 0x1006);
          std::string calTypeValue;
          if (!GetTagInfo(file,calType,calTypeValue)){
            LOG(ERROR) << "Unable to find type of normalisation data";
            return FileType::EERROR;  
          }
          LOG(INFO) << "type of normalisation data: " << calTypeValue;
          if (calTypeValue.find("0") != std::string::npos) {
            foundFileType = FileType::EGEPETNORM2D;
          }
          else if (calTypeValue.find("2") != std::string::npos) {
            foundFileType = FileType::EGEPETNORM3D;
          }
        }
        else if (rawDataTypeValue.find("5") != std::string::npos) {
          // geometric calibration file
          const gdcm::Tag calType(0x0017, 0x1006);
          std::string calTypeValue;
          if (!GetTagInfo(file,calType,calTypeValue)){
            LOG(ERROR) << "Unable to find type of calibration data";
            return FileType::EERROR;  
          }
          LOG(INFO) << "type of geo data: " << calTypeValue;
          if (calTypeValue.find("3") != std::string::npos) {
            foundFileType = FileType::EGEPETGEO;
          }
        }
        else if (rawDataTypeValue.find("7") != std::string::npos) {
          // WCC file
          LOG(ERROR) << "pet-rd-tools does not support GE Well-counter-calibration (WCC) files yet";
          // foundFileType = FileType::EGEPETGEO;
        }

    }
    return foundFileType;
  }

  private:
      IGEPET* Create_ptr( boost::filesystem::path inFile ) {

        FileType fType = GetFileType( inFile );
    
        if (fType == FileType::EGEPETLIST){
          IGEPET* instance(new GEPETList(inFile));
          return instance;
        }

        if (fType == FileType::EGEPETSINO){
          IGEPET* instance(new GEPETSino(inFile)); 
          return instance;
        }

        if (fType == FileType::EGEPETNORM3D || fType == FileType::EGEPETNORM2D){
          IGEPET* instance(new GEPETNorm(inFile));  
          return instance;    
        }

        if (fType == FileType::EGEPETGEO){
          IGEPET* instance(new GEPETGeo(inFile));
          return instance;    
        }

        return nullptr;
      }
  };

} // namespace nmtools

#endif 
