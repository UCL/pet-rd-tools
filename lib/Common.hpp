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

 */

#ifndef COMMON_HPP
#define COMMON_HPP

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

FileType GetFileType( boost::filesystem::path src){

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

}

#endif