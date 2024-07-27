/*
   Common.hpp

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

  Common utils. DICOM reading.
  
 */

#ifndef COMMON_HPP
#define COMMON_HPP

#include <itkImage.h>
#include <gdcmStringFilter.h>
#include <exception>
#include <sstream>

namespace nmtools {

#ifdef __APPLE__
    #define fseeko64 fseeko
    #define ftello64 ftello
#endif

#ifdef WIN32
    #define fseeko64 _fseeki64
    #define ftello64 _ftelli64
#endif

enum class ContentType { EHEADER, ERAWDATA };
enum class FileStatusCode { EGOOD, EBAD, EIOERROR };

bool GetTagInfo(const gdcm::File &file, const gdcm::Tag tag, std::string &dst){

  //Extracts information for a given DICOM tag from a gdcm file.
  //Tag contents are returned as a string in dst variable.

  //TODO: Do actual check for valid content.

  //Tries to read the element associated with the tag. If the read fails, the
  //DataElement should have a ByteValue of NULL.

  std::stringstream inStream;
  
  const gdcm::DataSet &ds = file.GetDataSet();
  gdcm::DataElement element = ds.GetDataElement(tag);

  dst = "";

  try {

    gdcm::StringFilter sf;
    sf.SetFile(file);

    if (element.GetByteValue() != NULL) {
      dst = sf.ToString(tag);
    }
    //else return false;

  } catch (std::bad_alloc){
    LOG(ERROR) << "GetTagInfo : Cannot read!";
    return false;
  }

  if (dst.size() == 0){
    LOG(WARNING) << "GetTagInfo : Empty field - " << tag;
    LOG(WARNING) << "Reverting to using element.GetValue()";
    inStream << std::fixed << element.GetValue();
    std::string tmpDst = inStream.str();

    std::string smatch = "Loaded:";

    if (tmpDst.compare(0, smatch.length(), smatch ) != 0){
      LOG(INFO) << "Value found via element.GetValue()";
      dst = tmpDst;
    }
  }

  return true;
}

class IDicomExtractor {

//Base class for extracting headers etc from a (probably DICOM) file
public:

  IDicomExtractor();
  explicit IDicomExtractor(boost::filesystem::path src);

  virtual bool SetInputFile ( boost::filesystem::path src );
  //FileType GetFileType( boost::filesystem::path src );

  virtual bool IsValid()=0;

  virtual bool ExtractHeader( const boost::filesystem::path dst ) = 0;
  virtual bool ExtractData( const boost::filesystem::path dst ) = 0;
  virtual boost::filesystem::path GetStdFileName( boost::filesystem::path srcFile, ContentType ctype) = 0;
  virtual bool ModifyHeader( const boost::filesystem::path src, const boost::filesystem::path dataFile) = 0;

  virtual ~IDicomExtractor(){};

protected:
  std::unique_ptr<gdcm::Reader> _dicomReader = nullptr;

  boost::filesystem::path _srcPath;

};

IDicomExtractor::IDicomExtractor() {
    std::unique_ptr<gdcm::Reader> dcm(new gdcm::Reader);
  _dicomReader = std::move(dcm);
}

//Constructor with path.
IDicomExtractor::IDicomExtractor(boost::filesystem::path src){

  std::unique_ptr<gdcm::Reader> dcm(new gdcm::Reader);
  _dicomReader = std::move(dcm);

  if (! this->SetInputFile(src) ) {
    LOG(ERROR) << "Unable to read data in: " << src;
    throw std::invalid_argument("Unable to read \"" + src.string() + "\" as DICOM");
  }

}

//On set input, try and read.
bool IDicomExtractor::SetInputFile(boost::filesystem::path src){

    _dicomReader->SetFileName(src.string().c_str());

    if (!_dicomReader->Read()) {
        LOG(ERROR) << "Unable to read as DICOM file";
        return false;
    }

    _srcPath = src;

    return true;
}

class IRawDataFactory {
//Factory that returns suitable child for given data.

public:
  virtual std::unique_ptr<IDicomExtractor> Create( boost::filesystem::path inFile )
  {
    return std::unique_ptr<IDicomExtractor>(Create_ptr( inFile ));
  }
protected:
  std::unique_ptr<gdcm::Reader> dicomReader;
  std::string manufacturerName;
  std::string modelName;

  bool Open(boost::filesystem::path inFile);
  virtual IDicomExtractor* Create_ptr( boost::filesystem::path inFile ) = 0;
};

bool IRawDataFactory::Open(boost::filesystem::path inFile) {
  dicomReader = std::unique_ptr<gdcm::Reader>(new gdcm::Reader);
  dicomReader->SetFileName(inFile.string().c_str());

  if (!dicomReader->Read()) {
    LOG(ERROR) << "Unable to read '" << inFile.string() << "' as DICOM file";
    return false;
  }

  //Get dataset via GDCM
  const gdcm::File &file = dicomReader->GetFile();
  //const gdcm::DataSet &ds = dicomReader->GetFile().GetDataSet();

  //Read manufacturer name.
  const gdcm::Tag manufacturer(0x008, 0x0070);
  if (!GetTagInfo(file,manufacturer,manufacturerName)){
    LOG(ERROR) << "Unable to manufacturer name";
    return false;
  }

  LOG(INFO) << "Manufacturer: " << manufacturerName;

  //Read model of scanner
  const gdcm::Tag model(0x008, 0x1090);
  if (!GetTagInfo(file,model,modelName)){
    LOG(ERROR) << "Unable to scanner model name";
    return false;
  }
  LOG(INFO) << "Model name: " << modelName;

  return true;
}

itk::SpatialOrientation::CoordinateTerms GetOrientationCode(char &c){

  c = toupper(c);

  const std::string validVals = "RLPAIS";

  if (validVals.find(c) == std::string::npos){
    LOG(ERROR) << c << " is not a valid orientation code value!";
    return itk::SpatialOrientation::ITK_COORDINATE_UNKNOWN;
  }

  if (c == 'R')
    return itk::SpatialOrientation::ITK_COORDINATE_Right;

  if (c == 'L')
    return itk::SpatialOrientation::ITK_COORDINATE_Left;

  if (c == 'P')
    return itk::SpatialOrientation::ITK_COORDINATE_Posterior;

  if (c == 'A')
    return itk::SpatialOrientation::ITK_COORDINATE_Anterior;

  if (c == 'I')
    return itk::SpatialOrientation::ITK_COORDINATE_Inferior;

  if (c == 'S')
    return itk::SpatialOrientation::ITK_COORDINATE_Superior;

  return itk::SpatialOrientation::ITK_COORDINATE_UNKNOWN;

}

bool SetDesiredCoordinateOrientation(const std::string &target,
                                     itk::SpatialOrientation::ValidCoordinateOrientationFlags &finalOrientation){

  std::vector<int> coordVals(3);

  std::string orient = target;
  //Check we have three letter code.
  if (orient.size() != 3){
    LOG(ERROR) << "Expected three letter orientation code. Read: " << orient;
    return false;
  }

  //Check they are all valid identifiers
  for (int i = 0; i < 3; i++) {
    coordVals[i] = static_cast<int>(GetOrientationCode(orient[i]));
    if (coordVals[i] == 0){
      LOG(ERROR) << "Unknown coordinate: " << orient[i];
      return false;
    }
  }

  //See itkSpatialOrientation.h
  itk::SpatialOrientation::ValidCoordinateOrientationFlags o =
      (itk::SpatialOrientation::ValidCoordinateOrientationFlags)(
          ( coordVals[0] << static_cast<int>(itk::SpatialOrientation::ITK_COORDINATE_PrimaryMinor )) +
          ( coordVals[1] << static_cast<int>(itk::SpatialOrientation::ITK_COORDINATE_SecondaryMinor )) +
          ( coordVals[2] << static_cast<int>(itk::SpatialOrientation::ITK_COORDINATE_TertiaryMinor )));

  //Check we don't have an duplicates.
  std::sort(coordVals.begin(), coordVals.end());
  auto last = std::unique(coordVals.begin(), coordVals.end());
  coordVals.erase(last, coordVals.end());

  if (coordVals.size() != 3){
    LOG(ERROR) << "Duplicate coordinate codes found: " << orient;
    return false;
  }

  LOG(INFO) << "Using orientation code: " << orient;
  finalOrientation = o;

  return true;

}

}

#endif
