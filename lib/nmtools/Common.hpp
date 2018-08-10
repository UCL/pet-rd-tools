/*
   Common.hpp

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

  Common utils. DICOM reading.
  
 */

#ifndef COMMON_HPP
#define COMMON_HPP

#include <itkImage.h>

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
enum class FileType { EMMRSINO, EMMRLIST, EMMRNORM, EUNKNOWN, EERROR };
enum class FileStatusCode { EGOOD, EBAD, EIOERROR };

bool GetTagInfo(const gdcm::DataSet &ds, const gdcm::Tag tag, std::string &dst){

  //Extracts information for a given DICOM tag from a gdcm dataset.
  //Tag contents are returned as a string in dst variable.

  //TODO: Do actual check for valid content.

  //Tries to read the element associated with the tag. If the read fails, the
  //DataElement should have a ByteValue of NULL.

  try {
    std::stringstream inStream;
    inStream.exceptions(std::ios::badbit);
    const gdcm::DataElement& element = ds.GetDataElement(tag);
    if (element.GetByteValue() != NULL) {
      inStream << std::fixed << element.GetValue();
      dst = inStream.str();
    }
    else return false;

  } catch (std::bad_alloc){
    LOG(ERROR) << "GetTagInfo : Cannot read!";
    return false;
  }

  return true;
}

FileType GetFileType( boost::filesystem::path src){

  //Extracts information via DICOM to determine what kind of mMR
  //raw data type we're dealing with.
  //
  //Will check for list mode, sinograms and norms.
  //
  //TODO: Support physio files? (11/12/2017)

  FileType foundFileType = FileType::EUNKNOWN;

  std::unique_ptr<gdcm::Reader> dicomReader(new gdcm::Reader);
  dicomReader->SetFileName(src.string().c_str());

  if (!dicomReader->Read()) {
    LOG(ERROR) << "Unable to read as DICOM file";
    return FileType::EERROR;
  }

  //Get dataset via GDCM
  const gdcm::DataSet &ds = dicomReader->GetFile().GetDataSet();

  //Read manufacturer name.
  const gdcm::Tag manufacturer(0x008, 0x0070);
  std::string manufacturerName;
  if (!GetTagInfo(ds,manufacturer,manufacturerName)){
    LOG(ERROR) << "Unable to manufacturer name";
    return FileType::EERROR;  
  }

  LOG(INFO) << "Manufacturer: " << manufacturerName;

  //Read model of scanner
  const gdcm::Tag model(0x008, 0x1090);
  std::string modelName;
  if (!GetTagInfo(ds,model,modelName)){
    LOG(ERROR) << "Unable to scanner model name";
    return FileType::EERROR;  
  }
  LOG(INFO) << "Model name: " << modelName;

  //Get image tpye description
  const gdcm::Tag imageType(0x0008, 0x0008);
  std::string imageTypeValue;
  if (!GetTagInfo(ds,imageType,imageTypeValue)){
    LOG(ERROR) << "Unable to image type!";
    return FileType::EERROR;  
  }
  LOG(INFO) << "Image type: " << imageTypeValue;

  if (manufacturerName.find("SIEMENS") != std::string::npos) {
    DLOG(INFO) << "Manufacturer = SIEMENS";

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
    coordVals[i] = GetOrientationCode(orient[i]);
    if (coordVals[i] == 0){
      LOG(ERROR) << "Unknown coordinate: " << orient[i];
      return false;
    }
  }

  //See itkSpatialOrientation.h
  itk::SpatialOrientation::ValidCoordinateOrientationFlags o =
      (itk::SpatialOrientation::ValidCoordinateOrientationFlags)(
          ( coordVals[0] << itk::SpatialOrientation::ITK_COORDINATE_PrimaryMinor ) +
          ( coordVals[1] << itk::SpatialOrientation::ITK_COORDINATE_SecondaryMinor ) +
          ( coordVals[2] << itk::SpatialOrientation::ITK_COORDINATE_TertiaryMinor ));

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