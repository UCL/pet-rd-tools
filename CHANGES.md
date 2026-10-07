# ChangeLog

## v2.1.1

* add pre-commit
* add GitHub Actions checks

## v2.1.0

* fix win build: `path::string()`, `google::GLOG_INFO`

## v2.0.2

* more Siemens debug info
* fix ITK enum change
* remove deprecated Boost functions
* CMake `glog` fixes

## v2.0.1

* fix reading of Siemens data

## v2.0.0

* add capability to extract GE PET raw data. Just use `nm_extract`
* `nm_extract` write to input path by default

## v1.1.0

* Add default.nix for Nix package manager. (#30)
* Add `nm_signa2mu`
* `nm_mrac2mu`: adjust x-y padding (with --head) to meet default size.

## v1.0.0

* Support Siemens mMR
