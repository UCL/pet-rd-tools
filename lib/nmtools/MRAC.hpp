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

#include <glog/logging.h>

#include <boost/filesystem.hpp>

#include "json/json.hpp"

namespace nmtools {

typedef typename itk::Image<float, 3 >  MuMapImageType;

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

  typedef itk::GDCMImageIO ImageIOType;

public:

  MRAC2MU(){};
  MRAC2MU(boost::filesystem::path src);
  MRAC2MU(boost::filesystem::path src, nlohmann::json params);

  bool SetInput(boost::filesystem::path src);
  void SetParams(nlohmann::json params);

  bool Update();

  const typename MuMapImageType::Pointer GetOutput();
  std::string GetInterfileHdr() const;

  bool Write(boost::filesystem::path dst) const;

protected:

  bool Read();
  bool ScaleAndReslice();

  bool WriteToInterFile(boost::filesystem::path dst) const;

  typename MuMapImageType::Pointer _inputImage;
  typename MuMapImageType::Pointer _muImage;
  std::string _header;

  typename ImageIOType::Pointer _pDicomInfo;

  boost::filesystem::path _srcPath;
  nlohmann::json _params = resliceDefaultParams;

};

MRAC2MU::MRAC2MU(boost::filesystem::path src){

  if (!SetInput(src))
    LOG(ERROR) << "Input path not set!";

  DLOG(INFO) << "SRC = " << _srcPath;

}

MRAC2MU::MRAC2MU(boost::filesystem::path src, nlohmann::json params){

  if (!SetInput(src))
    LOG(ERROR) << "Input path not set!";

  _params = params;
  DLOG(INFO) << "JSON = " << std::setw(4) << _params;
}

bool MRAC2MU::SetInput(boost::filesystem::path src){

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

void MRAC2MU::SetParams(nlohmann::json params){
  //TODO: validate JSON inputs

  _params = params;
}

bool MRAC2MU::Update(){
  return ScaleAndReslice();
}

const typename MuMapImageType::Pointer MRAC2MU::GetOutput(){
  return _muImage;
}

std::string MRAC2MU::GetInterfileHdr() const {
  return _header;
}

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

    typedef std::vector<std::string> SeriesIdContainer;
    const SeriesIdContainer &seriesUID = nameGenerator->GetSeriesUIDs();
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

    //Orient to LPS
    typedef typename itk::OrientImageFilter<MuMapImageType,MuMapImageType> OrienterType;
    typename OrienterType::Pointer orienter = OrienterType::New();
  
    orienter->UseImageDirectionOn();
    orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPS);
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
  _header += "!INTERFILE";

  return true;

}

bool MRAC2MU::ScaleAndReslice(){

  typedef typename itk::DivideImageFilter<MuMapImageType, MuMapImageType, MuMapImageType> DivideFilterType;

  typedef typename itk::IdentityTransform< double, 3 > TransformType;
  typedef typename itk::LinearInterpolateImageFunction< MuMapImageType, double > InterpolatorType;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();

  const MuMapImageType::SpacingType& inputSpacing =
    _inputImage->GetSpacing();
  const MuMapImageType::RegionType& inputRegion =
    _inputImage->GetLargestPossibleRegion();
  const MuMapImageType::SizeType& inputSize =
    inputRegion.GetSize();

  MuMapImageType::SpacingType outputSpacing;
  outputSpacing[0] = _params["px"];
  outputSpacing[1] = _params["py"];
  outputSpacing[2] = _params["pz"];

  MuMapImageType::SizeType   outputSize;
  typedef MuMapImageType::SizeType::SizeValueType SizeValueType;
  outputSize[0] = static_cast<SizeValueType>(inputSize[0] * inputSpacing[0] / outputSpacing[0] + .5);
  outputSize[1] = static_cast<SizeValueType>(inputSize[1] * inputSpacing[1] / outputSpacing[1] + .5);
  outputSize[2] = static_cast<SizeValueType>(inputSize[2] * inputSpacing[2] / outputSpacing[2] + .5);

  typedef typename itk::ResampleImageFilter< MuMapImageType, MuMapImageType > ResampleFilterType;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetInput( _inputImage );
  resampler->SetTransform( transform );
  resampler->SetInterpolator( interpolator );
  resampler->SetOutputOrigin ( _inputImage->GetOrigin());
  resampler->SetOutputSpacing ( outputSpacing );
  resampler->SetOutputDirection ( _inputImage->GetDirection());
  resampler->SetSize ( outputSize );

  //Z-flip ?

  try { 
  resampler->Update ();
  }
  catch (itk::ExceptionObject &ex){
    //std::cout << ex << std::endl;
    LOG(ERROR) << "Unable to resample!";
    return false;
  }

  DivideFilterType::Pointer divide = DivideFilterType::New();
  divide->SetInput( resampler->GetOutput() );
  divide->SetConstant( 10000.0 );

  //pad xy
  MuMapImageType::SizeType lowerExtendRegion;
  lowerExtendRegion[0] = 100;
  lowerExtendRegion[1] = 100;
  lowerExtendRegion[2] = 0;

  MuMapImageType::SizeType upperExtendRegion;
  upperExtendRegion[0] = 100;
  upperExtendRegion[1] = 100;
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

  //crop in z
  MuMapImageType::SizeType lcropSize;
  lcropSize[0] = 0;
  lcropSize[1] = 0;
  lcropSize[2] = 11;

  MuMapImageType::SizeType ucropSize;
  ucropSize[0] = 0;
  ucropSize[1] = 0;
  ucropSize[2] = 10;

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

  return true;
}

bool MRAC2MU::Write(boost::filesystem::path dst) const {

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

bool MRAC2MU::WriteToInterFile(boost::filesystem::path dst) const {

  boost::filesystem:: path altPath = boost::filesystem::change_extension(dst, ".mhd");
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

  std::string outputHeader = GetInterfileHdr();

  altPath = boost::filesystem::change_extension(dst, ".hv");
  
  std::ofstream infoStream;

  if ( infoStream ) {
    infoStream.open(altPath.string().c_str());
    infoStream << outputHeader << std::endl;
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