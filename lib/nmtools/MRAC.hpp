/*
   MRAC.hpp

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

   Generates a mu-map from Siemens mMR MRAC for PET reconstruction.
 */

#ifndef MRAC_HPP
#define MRAC_HPP

#include <iostream>
#include <string>

#include <itkImage.h>
#include <itkImageSeriesReader.h>
#include <itkImageFileWriter.h>
#include <itkImageDuplicator.h>
#include <itkGDCMImageIO.h>
#include <itkGDCMSeriesFileNames.h>
#include <itkResampleImageFilter.h>
#include <itkIdentityTransform.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkConstantPadImageFilter.h>
#include <itkCropImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkOrientImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>

#include <glog/logging.h>

#include <boost/filesystem.hpp>
#include <boost/any.hpp>
#include <boost/lexical_cast.hpp>

#include "nmtools/Common.hpp"
#include "json/json.hpp"

namespace nmtools {

//All images are 3D 32-bit float ITK images.
typedef typename itk::Image<float, 3 >  MuMapImageType;

//Default parameters for reslicing
//FOV = 700mm; voxel size = [2.09,2.09,2.03];
//Matrix size: [344,344,127].
const nlohmann::json resliceDefaultParams = R"(
{
  "FOV": 700.0,
  "px": 2.08626,
  "py": 2.08626,
  "pz": 2.03125,
  "sx": 344,
  "sy": 344,
  "sz": 127
}
)"_json;

class MRAC2MU {
  //Class for converting from mMR MRAC to mu values.

  typedef itk::GDCMImageIO ImageIOType;

public:

  //Either just empty constructor, with input directory or
  //with input directory and user-specified json params.
  MRAC2MU(){};
  explicit MRAC2MU(boost::filesystem::path src, std::string orientationCode);
  MRAC2MU(boost::filesystem::path src, nlohmann::json params, std::string orientationCode);

  //Set input file and attempt to read.
  bool SetInput(boost::filesystem::path src);

  //Accept alternative reslicing parameters.
  void SetParams(nlohmann::json params);

  //Set output orientation
  bool SetDesiredCoordinateOrientation(const std::string &target);

  void SetIsHead(bool bStatus){ _isHead = bStatus; };

  //Trigger exexcution
  bool Update();

  //Return final image.
  const typename MuMapImageType::Pointer GetOutput();

  //Get manufactured Interfile header. 
  std::string GetInterfileHdr() const;

  //Write file(s) to dst.
  bool Write(boost::filesystem::path dst);

protected:

  //File reading
  bool Read();

  //Do reslicing etc.
  bool Scale();
  bool ScaleAndResliceHead();

  //Write interfile case.
  bool WriteToInterFile(boost::filesystem::path dst);

  //Modify header 
  bool UpdateInterfile(const std::string &key, const boost::any info);

  //Grab info from DICOM data.
  bool GetStudyDate(std::string &studyDate);
  bool GetStudyTime(std::string &studyTime);

  //Original image
  typename MuMapImageType::Pointer _inputImage;

  //Output image
  typename MuMapImageType::Pointer _muImage;

  //Interfile header
  std::string _header;

  //DICOM data
  typename ImageIOType::Pointer _pDicomInfo;

  //Source path
  boost::filesystem::path _srcPath;

  //JSON params for reslicing.
  nlohmann::json _params = resliceDefaultParams;

  //Default image orientation is RAI
  itk::SpatialOrientation::ValidCoordinateOrientationFlags _outputOrientation 
    = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI;

  //Reslice and crop into 344x344 matrix for brain. Off by default.
  bool _isHead = false;

};

//Construct object from source directory.
MRAC2MU::MRAC2MU(boost::filesystem::path src, std::string orientationCode = "RAI"){

  if (!SetDesiredCoordinateOrientation(orientationCode)){
    throw false;
  }

  if (!SetInput(src)) {
    LOG(ERROR) << "Input path not set!";
    throw false;
  }

  DLOG(INFO) << "SRC = " << _srcPath;

}

//Construct object from source directory and use user-specified
//params.
MRAC2MU::MRAC2MU(boost::filesystem::path src, nlohmann::json params, std::string orientationCode = "RAI"){

  if (!SetDesiredCoordinateOrientation(orientationCode)){
    throw false;
  }

  if (!SetInput(src)) {
    LOG(ERROR) << "Input path not set!";
    throw false;
  }

  _params = params;
  DLOG(INFO) << "JSON = " << std::setw(4) << _params;
}

//Check if the source dir is what is expected.
bool MRAC2MU::SetInput(boost::filesystem::path src){
  //TODO: Should probably check whether src contains
  //DICOM data. (11/12/2017)

  //Check if input file even exists!
  if (!boost::filesystem::exists(src)) {
    LOG(ERROR) << "Input path " << src << " does not exist!";
    return false;
  }

  //Check if the input is a file.
  if (!boost::filesystem::is_directory(src)) {
    LOG(ERROR) << src.native() << " does not appear to be a  directory!";
    return false;
  }

  _srcPath = src;

  DLOG(INFO) << "SRC = " << _srcPath;

  return Read();
}

//Use user-specified reslicing parameters.
void MRAC2MU::SetParams(nlohmann::json params){
  //TODO: validate JSON inputs (11/12/2017)

  _params = params;
}

//Run pipeline
bool MRAC2MU::Update(){
  if (_isHead)
    return ScaleAndResliceHead();

  return Scale();
}

//Get final image.
const typename MuMapImageType::Pointer MRAC2MU::GetOutput(){
  return _muImage;
}

//Return header.
std::string MRAC2MU::GetInterfileHdr() const {
  return _header;
}

//- Reads input directory with GDCM (via ITK).
//- Creates image from slices and orients to LPS.
//- Creates Interfile skeleton.
bool MRAC2MU::Read(){

  //Generic DICOM to ITK image read.
  //Will only convert first series in a folder.

  DLOG(INFO) << "Reading DICOMDIR";

  if (!boost::filesystem::exists(_srcPath))
  {
    LOG(ERROR) << "Input path " << _srcPath << " does not exist!";
    return false;
  }

  _pDicomInfo = ImageIOType::New();

  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

  nameGenerator->SetUseSeriesDetails(true);
  nameGenerator->AddSeriesRestriction("0008|0021"); //Series date
  nameGenerator->AddSeriesRestriction("0020|0037"); //Patient orientation
  nameGenerator->SetLoadSequences(true);
  nameGenerator->SetLoadPrivateTags(true);

  _pDicomInfo->SetMaxSizeLoadEntry(0xffffffffffffffff);
  _pDicomInfo->SetLoadPrivateTags(true);
  _pDicomInfo->SetLoadSequences(true);
  //_pDicomInfo->SetLoadPrivateTagsDefault(true);

  //Create ITK ImageSeriesReader
  typedef typename itk::ImageSeriesReader<MuMapImageType> ReaderType;
  typename ReaderType::Pointer dicomReader = ReaderType::New();

  try
  {
    DLOG(INFO) << "DICOMDIR directory: " << _srcPath;

    std::string pathToSeries = _srcPath.string();
    //Load directory
    nameGenerator->SetDirectory(pathToSeries);

    //TODO: Should probably check Series UIDs here 14/12/2017
    //typedef std::vector<std::string> SeriesIdContainer;
    //const SeriesIdContainer &seriesUID = nameGenerator->GetSeriesUIDs();

    const typename ReaderType::FileNamesContainer &fileNames = nameGenerator->GetInputFileNames();
    nameGenerator->Update();

    //Set reader to update _pDicomInfo with header information.
    dicomReader->SetImageIO(_pDicomInfo);
    dicomReader->SetFileNames(fileNames);

    if (fileNames.size() == 0)
    {
      LOG(ERROR) << "No valid DICOM series found";
      return false;
    }
  }
  catch (itk::ExceptionObject &ex)
  {
    //std::cout << ex << std::endl;
    LOG(ERROR) << "Cannot read DICOM directory!";
    return false;
  }

  try
  {
    //Execute pipeline
    dicomReader->Update();

    //Re-orient
    typedef typename itk::OrientImageFilter<MuMapImageType,MuMapImageType> OrienterType;
    typename OrienterType::Pointer orienter = OrienterType::New();
  
    orienter->UseImageDirectionOn();
    orienter->SetDesiredCoordinateOrientation(_outputOrientation);
    orienter->SetInput(dicomReader->GetOutput());
    orienter->Update();

    //Duplicate contents of reader into _inputImage.
    typedef typename itk::ImageDuplicator<MuMapImageType> DuplicatorType;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(orienter->GetOutput());
    duplicator->Update();
    _inputImage = duplicator->GetOutput();

    DLOG(INFO) << "DICOM Origin: " << _inputImage->GetOrigin();
  }
  catch (itk::ExceptionObject &ex)
  {
    //std::cout << ex << std::endl;
    LOG(ERROR) << "Unable to get image from DICOM series";
    return false;
  }

  DLOG(INFO) << "Reading complete";

  //TODO: Finish filling Interfile header
  _header.clear();

  std::stringstream ss;

  ss << "!INTERFILE:=" << std::endl;
  ss << "%comment:=created with nm_mrac2mu for mMR data" << std::endl;
  ss << "!originating system:=2008" << std::endl;

  ss << std::endl << "!GENERAL DATA:=" << std::endl;
  ss << "!name of data file:=<%%DATAFILE%%>";
  
  ss << std::endl << "!GENERAL IMAGE DATA:=" << std::endl;
  ss << "!type of data := PET" << std::endl;

  ss << std::endl << "%study date (yyyy:mm:dd):=<%%STUDYDATE%%>"  << std::endl;
  ss << "%study time (hh:mm:ss GMT+00:00):=<%%STUDYTIME%%>" << std::endl;
  ss << "imagedata byte order:=LITTLEENDIAN" << std::endl;
  ss << "%patient orientation:=HFS" << std::endl;
  ss << "!PET data type:=image" << std::endl;
  ss << "number format:=float" << std::endl;
  ss << "!number of bytes per pixel:=4" << std::endl;
  ss << "number of dimensions:=3" << std::endl;
  ss << "matrix axis label[1]:=x" << std::endl;
  ss << "matrix axis label[2]:=y" << std::endl;
  ss << "matrix axis label[3]:=z" << std::endl;
  ss << "matrix size[1]:=<%%NX%%>" << std::endl;
  ss << "matrix size[2]:=<%%NY%%>" << std::endl;
  ss << "matrix size[3]:=<%%NZ%%>" << std::endl;
  ss << "scaling factor (mm/pixel) [1]:=<%%SX%%>" << std::endl;
  ss << "scaling factor (mm/pixel) [2]:=<%%SY%%>" << std::endl;
  ss << "scaling factor (mm/pixel) [3]:=<%%SZ%%>" << std::endl;
  ss << "start horizontal bed position (mm):=-10" << std::endl;
  ss << "end horizontal bed position (mm):=-10" << std::endl;
  ss << "start vertical bed position (mm):=0.0" << std::endl;

  ss << std::endl << "!IMAGE DATA DESCRIPTION:=" << std::endl;
  ss << "!total number of data sets:=1" << std::endl;
  ss << "number of time frames:=1" << std::endl;
  ss << "!image duration (sec)[1]:=0" << std::endl;
  ss << "!image relative start time (sec)[1]:=0" << std::endl;

  ss << std::endl << "%SUPPLEMENTARY ATTRIBUTES:=" << std::endl;
  ss << "quantification units:=1/cm" << std::endl;
  ss << "slice orientation:=Transverse" << std::endl;
  ss << "%image zoom:=1" << std::endl;
  ss << "%x-offset (mm):=0.0" << std::endl;
  ss << "%y-offset (mm):=0.0" << std::endl;
  ss << "%image slope:=1" << std::endl;
  ss << "%image intercept:=0.0" << std::endl;
  ss << "maximum pixel count:=<%%MAXVAL%%>" << std::endl;
  ss << "minimum pixel count:=<%%MINVAL%%>" << std::endl;

  ss << "!END OF INTERFILE :=" << std::endl;

  _header = ss.str();

  return true;

}

//Insert info. into Interfile header.
//Casting might not be totally safe.
bool MRAC2MU::UpdateInterfile(const std::string &key, const boost::any info){

  std::string updateStr;

  bool bStatus = false;

  try {
    updateStr = boost::any_cast<std::string>(info);
    bStatus=true;
  }
  catch (boost::bad_any_cast &e) {

  }

  if (!bStatus){
    try {
      updateStr = boost::lexical_cast<std::string>(boost::any_cast<const char *>(info));
      bStatus=true;
    }
    catch (boost::bad_any_cast &e) {

    }
  }

  if (!bStatus){
    try {
      updateStr = boost::lexical_cast<std::string>(boost::any_cast<float>(info));
      bStatus=true;
    }
    catch (boost::bad_any_cast &e) {
      
    }
  }

  if (!bStatus){
    try {
      updateStr = boost::lexical_cast<std::string>(boost::any_cast<int>(info));
      bStatus=true;
    }
    catch (boost::bad_any_cast &e) {
      
    }
  }

  if (!bStatus){
    LOG(WARNING) << "Unable to find conversion for Interfile header update";
    return false;
  }

  std::string target = "<%%" + key + "%%>";
  std::string::size_type n;

  n = _header.find(target);

  if (n == std::string::npos){
    LOG(WARNING) << "Interfile replacement key: " << target << " not found!";
    return false; 
  }

  _header.replace(n,target.length(),updateStr);

  return bStatus;
}

//Get study date from DICOM and convert from 'YYYYMMDD' to 'YYYY:MM:DD'.
bool MRAC2MU::GetStudyDate(std::string &studyDate){

  std::string infoStr;

  bool bStatus = _pDicomInfo->GetValueFromTag("0008|0020",infoStr);

  if (!bStatus){
    studyDate = "";
    return false;
  }

  LOG(INFO) << "Study date: " << infoStr;

  studyDate = infoStr.substr(0,4) + ":" + infoStr.substr(4,2) + ":" + infoStr.substr(6,2);

  return true;
}

//Get study time from DICOM and convert from 'HHMMSS.mmmmmm' to 'HH:MM:SS'.
bool MRAC2MU::GetStudyTime(std::string &studyTime){

  std::string infoStr;

  bool bStatus = _pDicomInfo->GetValueFromTag("0008|0030",infoStr);

  if (!bStatus){
    studyTime = "";
    return false;
  }

  LOG(INFO) << "Study time: " << infoStr;
  studyTime = infoStr.substr(0,2) + ":" + infoStr.substr(2,2) + ":" + infoStr.substr(4,2);

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

};

bool MRAC2MU::SetDesiredCoordinateOrientation(const std::string &target){ 

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
  itk::SpatialOrientation::ValidCoordinateOrientationFlags o = (itk::SpatialOrientation::ValidCoordinateOrientationFlags)(
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
  _outputOrientation = o; 

  return true;

};

//Divide by 10000 to get mu-values (cm-1).
//Interpolate and reslice according to JSON params.
bool MRAC2MU::Scale(){

  typedef typename itk::DivideImageFilter<MuMapImageType, MuMapImageType, MuMapImageType> DivideFilterType;

  //Scale to mu-values
  DivideFilterType::Pointer divide = DivideFilterType::New();
  divide->SetInput( _inputImage );
  divide->SetConstant( 10000.0 );

  try {
    divide->Update();
  } catch (itk::ExceptionObject &ex){
    //std::cout << ex << std::endl;
    LOG(ERROR) << "Unable to scale!";
    return false;
  }

  try {
    //Duplicate contents or reader into _muImage.
    typedef typename itk::ImageDuplicator<MuMapImageType> DuplicatorType;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(divide->GetOutput());
    duplicator->Update();
    _muImage = duplicator->GetOutput();
  } catch (itk::ExceptionObject &ex){
    //std::cout << ex << std::endl;
    LOG(ERROR) << "Unable to scale to mu!";
    return false;
  }

  //Get max and min voxel values to insert into Interfile header.
  typedef typename itk::MinimumMaximumImageCalculator <MuMapImageType> ImageCalculatorFilterType;
  typename ImageCalculatorFilterType::Pointer minmax = ImageCalculatorFilterType::New();
  minmax->SetImage(_muImage);

  try {
    minmax->Compute();
  } catch (itk::ExceptionObject &ex){
    //std::cout << ex << std::endl;
    LOG(ERROR) << "Unable to calculate min/max!";
    return false;
  }

  //Update the Interfile header with new sizes etc.
  const MuMapImageType::SizeType& size = _muImage->GetLargestPossibleRegion().GetSize();
  this->UpdateInterfile("NX", int(size[0]));
  this->UpdateInterfile("NY", int(size[1]));
  this->UpdateInterfile("NZ", int(size[2]));

  const MuMapImageType::SpacingType& voxSize = _muImage->GetSpacing();
  this->UpdateInterfile("SX", float(voxSize[0]));
  this->UpdateInterfile("SY", float(voxSize[1]));
  this->UpdateInterfile("SZ", float(voxSize[2]));  

  this->UpdateInterfile("MAXVAL", float(minmax->GetMaximum()));
  this->UpdateInterfile("MINVAL", float(minmax->GetMinimum()));

  std::string studyDate;
  if (GetStudyDate(studyDate))
    this->UpdateInterfile("STUDYDATE", studyDate);

  std::string studyTime;
  if (GetStudyTime(studyTime))
    this->UpdateInterfile("STUDYTIME", studyTime);

  return true;
}

//Divide by 10000 to get mu-values (cm-1).
//Interpolate and reslice according to JSON params.
bool MRAC2MU::ScaleAndResliceHead(){

  typedef typename itk::DivideImageFilter<MuMapImageType, MuMapImageType, MuMapImageType> DivideFilterType;

  typedef typename itk::IdentityTransform< double, 3 > TransformType;
  typedef typename itk::LinearInterpolateImageFunction< MuMapImageType, double > InterpolatorType;

  //Linear interpolation and identity transform.
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();

  //Grab original voxel and matrix size.
  const MuMapImageType::SpacingType& inputSpacing =
    _inputImage->GetSpacing();
  const MuMapImageType::RegionType& inputRegion =
    _inputImage->GetLargestPossibleRegion();
  const MuMapImageType::SizeType& inputSize =
    inputRegion.GetSize();

  //Get new voxel size from JSON params.
  //Unsafe.
  MuMapImageType::SpacingType outputSpacing;
  outputSpacing[0] = _params["px"];
  outputSpacing[1] = _params["py"];
  outputSpacing[2] = _params["pz"];

  //Reslice to new matrix size.
  MuMapImageType::SizeType   outputSize;
  typedef MuMapImageType::SizeType::SizeValueType SizeValueType;
  outputSize[0] = static_cast<SizeValueType>(inputSize[0] * inputSpacing[0] / outputSpacing[0] + .5);
  outputSize[1] = static_cast<SizeValueType>(inputSize[1] * inputSpacing[1] / outputSpacing[1] + .5);
  outputSize[2] = static_cast<SizeValueType>(inputSize[2] * inputSpacing[2] / outputSpacing[2] + .5);

  if ((outputSize[0] % 2 == 1) or (outputSize[1] % 2 == 1)){
    LOG(ERROR) << "Input x or y size is odd. Unsure how to resample!";
    return false;
  }

  typedef typename itk::ResampleImageFilter< MuMapImageType, MuMapImageType > ResampleFilterType;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetInput( _inputImage );
  resampler->SetTransform( transform );
  resampler->SetInterpolator( interpolator );
  resampler->SetOutputOrigin ( _inputImage->GetOrigin());
  resampler->SetOutputSpacing ( outputSpacing );
  resampler->SetOutputDirection ( _inputImage->GetDirection());
  resampler->SetSize ( outputSize );

  try { 
  resampler->Update ();
  }
  catch (itk::ExceptionObject &ex){
    //std::cout << ex << std::endl;
    LOG(ERROR) << "Unable to resample!";
    return false;
  }

  //Scale to mu-values
  DivideFilterType::Pointer divide = DivideFilterType::New();
  divide->SetInput( resampler->GetOutput() );
  divide->SetConstant( 10000.0 );

  //Pad x-y
  MuMapImageType::SizeType lowerExtendRegion;
  int pad_x = (_params["sx"].get<int>() - outputSize[0]) / 2;
  if (pad_x < 0) {
    pad_x = 0;
  }
  int pad_y = (_params["sy"].get<int>() - outputSize[1]) / 2;
  if (pad_y < 0) {
    pad_y = 0;
  }
  lowerExtendRegion[0] = pad_x;
  lowerExtendRegion[1] = pad_y;
  lowerExtendRegion[2] = 0;

  MuMapImageType::SizeType upperExtendRegion;
  upperExtendRegion[0] = pad_x;
  upperExtendRegion[1] = pad_y;
  upperExtendRegion[2] = 0;

  MuMapImageType::PixelType constantPixel = 0;

  typedef typename itk::ConstantPadImageFilter <MuMapImageType, MuMapImageType> ConstantPadImageFilterType;
  ConstantPadImageFilterType::Pointer padFilter
    = ConstantPadImageFilterType::New();
  padFilter->SetInput( divide->GetOutput() );
  padFilter->SetPadLowerBound(lowerExtendRegion);
  padFilter->SetPadUpperBound(upperExtendRegion);
  padFilter->SetConstant(constantPixel);

  try {
    padFilter->Update();
  } catch (itk::ExceptionObject &ex){
    //std::cout << ex << std::endl;
    LOG(ERROR) << "Unable to pad!";
    return false;
  }

  // magic numbers
  int z_lcrop = 11;
  int z_ucrop = 10;
  //crop in z direction.
  MuMapImageType::SizeType lcropSize;
  lcropSize[0] = 0;
  lcropSize[1] = 0;
  lcropSize[2] = z_lcrop;

  MuMapImageType::SizeType ucropSize;
  ucropSize[0] = 0;
  ucropSize[1] = 0;
  ucropSize[2] = z_ucrop;

  typedef typename itk::CropImageFilter <MuMapImageType, MuMapImageType> CropImageFilterType;
  CropImageFilterType::Pointer cropFilter = CropImageFilterType::New();
  cropFilter->SetInput( padFilter->GetOutput() );
  cropFilter->SetLowerBoundaryCropSize(lcropSize);
  cropFilter->SetUpperBoundaryCropSize(ucropSize);

  try {
    cropFilter->Update();

    //Duplicate contents or reader into _muImage.
    typedef typename itk::ImageDuplicator<MuMapImageType> DuplicatorType;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(cropFilter->GetOutput());
    duplicator->Update();
    _muImage = duplicator->GetOutput();
  } catch (itk::ExceptionObject &ex){
    //std::cout << ex << std::endl;
    LOG(ERROR) << "Unable to scale to mu!";
    return false;
  }

  //Get max and min voxel values to insert into Interfile header.
  typedef typename itk::MinimumMaximumImageCalculator <MuMapImageType> ImageCalculatorFilterType;
  typename ImageCalculatorFilterType::Pointer minmax = ImageCalculatorFilterType::New();
  minmax->SetImage(_muImage);

  try {
    minmax->Compute();
  } catch (itk::ExceptionObject &ex){
    //std::cout << ex << std::endl;
    LOG(ERROR) << "Unable to calculate min/max!";
    return false;
  }

  //Update the Interfile header with new sizes etc.
  const MuMapImageType::SizeType& size = _muImage->GetLargestPossibleRegion().GetSize();
  this->UpdateInterfile("NX", int(size[0]));
  this->UpdateInterfile("NY", int(size[1]));
  this->UpdateInterfile("NZ", int(size[2]));

  const MuMapImageType::SpacingType& voxSize = _muImage->GetSpacing();
  this->UpdateInterfile("SX", float(voxSize[0]));
  this->UpdateInterfile("SY", float(voxSize[1]));
  this->UpdateInterfile("SZ", float(voxSize[2]));  

  this->UpdateInterfile("MAXVAL", float(minmax->GetMaximum()));
  this->UpdateInterfile("MINVAL", float(minmax->GetMinimum()));

  std::string studyDate;
  if (GetStudyDate(studyDate))
    this->UpdateInterfile("STUDYDATE", studyDate);

  std::string studyTime;
  if (GetStudyTime(studyTime))
    this->UpdateInterfile("STUDYTIME", studyTime);

  return true;
}

//Dump image (and header if applicable) to disk.
bool MRAC2MU::Write(boost::filesystem::path dst) {

  if (dst.extension() == ".hv"){
   return WriteToInterFile(dst);
  }

  //Write output file
  typedef typename itk::ImageFileWriter<MuMapImageType> WriterType; 

  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( dst.string().c_str() );
  writer->SetInput( _muImage );

  try {
    writer->Update();
  }
  catch (itk::ExceptionObject &ex){
    LOG(ERROR) << " Could not write output file!";
    return false;    
  }

  return true;

}

//Write Interfile header and image pair to disk.
bool MRAC2MU::WriteToInterFile(boost::filesystem::path dst){

  boost::filesystem::path altPath = boost::filesystem::change_extension(dst, ".mhd");
  //Write output file
  typedef typename itk::ImageFileWriter<MuMapImageType> WriterType; 

  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( altPath.string().c_str() );
  writer->SetInput( _muImage );

  try {
    writer->Update();
  }
  catch (itk::ExceptionObject &ex){
    LOG(ERROR) << " Could not write output data!";
    return false;    
  }

  //Put new data file in header.
  boost::filesystem::path dataFile =  dst.replace_extension(".raw");
  this->UpdateInterfile("DATAFILE", dataFile.filename().string());

  std::string outputHeader = GetInterfileHdr();

  altPath = boost::filesystem::change_extension(dst, ".hv");

  std::ofstream infoStream;

  if ( infoStream ) {
    infoStream.open(altPath.string().c_str());
    infoStream << outputHeader;
    infoStream.close();
    LOG(INFO) << "Wrote Interfile header to " << altPath;
  } else {
    LOG(INFO) << "Could not write Interfile header to " << altPath;
    return false;
  }

  return true;
}

} //namespace nmtools

#endif