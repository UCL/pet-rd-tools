# petmr-rd-tools

[![Build Status](https://travis-ci.org/UCL/petmr-rd-tools.svg?branch=master)](https://travis-ci.org/UCL/petmr-rd-tools) [![Codacy Badge](https://api.codacy.com/project/badge/Grade/d71cdf9cba3d4f9f9f973f371624bfe7)](https://www.codacy.com/app/bathomas/petmr-rd-tools?utm_source=github.com&utm_medium=referral&utm_content=UCL/petmr-rd-tools&utm_campaign=badger)

Command line tools for PET-MR (pre)-processing.

Tools for validating and extracting raw PET data, and associated files, for the purposes of image reconstruction via [STIR](https://github.com/UCL/STIR) and [SIRF](https://github.com/CCPPETMR/SIRF). Currently these tools are specifically for the [Siemens mMR system](https://www.healthcare.siemens.com/magnetic-resonance-imaging/mr-pet-scanner/biograph-mmr) although will likely expand in the future.

## Requirements for building

- CMake (>= 3.7)
- ITK (>= 4.7)
- Boost (>= 1.55)
- GLOG ([https://github.com/google/glog](https://github.com/google/glog))

---
## Running the applications
### `nm_validate`

The purpose of `nm_validate` is to confirm if a raw data file (or file pair) contains all the expected data.

#### Usage:

```bash
nm_validate -i <DICOM file>
```

Example output for a valid list mode file:
```
I1205 17:15:37.418709 3333379008 NMValidate.cpp:99] Started: Tue Dec  5 17:15:37 2017
I1205 17:15:37.419257 3333379008 NMValidate.cpp:100] Running 'nm_validate' version: 0.1.0
I1205 17:15:39.136943 3333379008 Common.hpp:89] Manufacturer: SIEMENS 
I1205 17:15:39.136970 3333379008 Common.hpp:97] Model name: Biograph_mMR
I1205 17:15:39.136979 3333379008 Common.hpp:105] Image type: ORIGINAL\PRIMARY\PET_LISTMODE 
I1205 17:15:41.942068 3333379008 MMR.hpp:399] Expected number of LM words: 74308087
I1205 17:15:41.942098 3333379008 MMR.hpp:406] 297232348 bytes in LM field
I1205 17:15:41.942103 3333379008 MMR.hpp:410] 297232348 / 4 = 74308087 words
I1205 17:15:41.942114 3333379008 NMValidate.cpp:134] File appears to be VALID
I1205 17:15:41.942119 3333379008 NMValidate.cpp:139] Time taken: 4 seconds
I1205 17:15:41.942123 3333379008 NMValidate.cpp:140] Ended: Tue Dec  5 17:15:41 2017
```

`nm_validate` will check list mode, sinogram and normalisation (norm) files. The size of the anticipated raw data is checked, but not the actual contents. Due to compression of the sinogram data, only the existence of files are tested.

### `nm_extract`

Raw PET data from the mMR scanner can be in one of two forms: a single DICOM file or a pair of files (one DICOM header and a raw binary file). `nm_extract` reads the DICOM data and extracts the Interfile header and the raw data for either of the two forms. Once extracted, the Interfile header can be used for image reconstruction with STIR. 

#### Usage:

```bash
nm_extract -i <DICOM file> [-o <OUTPUTDIR> -p <PREFIX> --noupdate ]
```
where `<DICOM file>` is the input file for extraction, `<OUTPUTDIR>` is the target output directory and `<PREFIX>` is the desired filename prefix for the output files. `--noupdate` will extract the raw Interfile without modification (mainly for debugging). If the `<OUTPUTDIR>` does not exist, `nm_validate` will attempt to create it. If `<OUTPUTDIR>` is not specified, the output will be written to the same directory as the input. 

#### Output extensions

- List mode files will be extracted with `.l` extensions for the list mode data and `.l.hdr` for the associated Interfile header.
- Sinograms files will have the extension `.s` for the sinogram data and `.s.hdr` for the associated Interfile header.
- Normalisation files will be extracted with `.n` and `.n.hdr` extensions.

### `nm_mrac2mu`

`nm_mrac2mu` extracts the patient mu-map from the DICOM data, scales to linear attenuation coefficients (LACs) and reslices into a full-size matrix (344 x 344 x 127) for PET reconstruction. The mu-map is oriented in LPS. 

#### Usage: 

```bash
nm_mrac2mu -i <DICOMDIR> -o <OUTPUT file>
```

where `<DICOMDIR>` is the path to the MRAC DICOM folder and `<OUTPUT file>` is the destination file.

#### Output extensions

- The output file type is determined by the extension of the destination file. To produce a compressed NiFTi file, specify the extension `.nii.gz` e.g. `myoutput.nii.gz`. 
- If the extension `.hv` is given, an Interfile header is generated in addition to the volume.