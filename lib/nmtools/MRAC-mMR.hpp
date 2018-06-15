/*
   MRAC-mMR.hpp

   Author:      Benjamin A. Thomas

   Copyright 2017-2018 Institute of Nuclear Medicine, University College London.
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Generating mu-maps from mMR MRAC for PET reconstruction.
 */

#ifndef MRAC_MMR_HPP
#define MRAC_MMR_HPP

#include <iostream>
#include <string>

#include "MRAC.hpp"

namespace nmtools {

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

class MMRMRAC : public MRAC2MU {
  using MRAC2MU::MRAC2MU;

public:
  bool Update();

protected:
  
  bool ScaleAndResliceHead();

};

//Run pipeline
bool MMRMRAC::Update(){
  if (_isHead) {
    LOG(INFO) << "Performing requested mMR head reslicing.";
    return ScaleAndResliceHead();
  }

  return Scale();
}

//Divide by 10000 to get mu-values (cm-1).
//Interpolate and reslice according to JSON params.
bool MMRMRAC::ScaleAndResliceHead(){

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

  //JSON params for reslicing.
  if (_params.empty())
    _params = resliceDefaultParams;

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

} //end namespace nmtools


#endif
