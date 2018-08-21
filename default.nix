/*
  petmr-rd-tools package for the Nix package manager.
  From the root directory of the project, run:

    nix-build

  This will build to ./result.

  This file is distributed as part of the petmr-rd-tools.

  Author: Ashley Gillman

  Copyright 2018 Commonwealth Scientific and Industrial Research Organisation
                 Australian eHealth Research Centre

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

let pkgs = import <nixpkgs> {}; in

({ stdenv
, fetchFromGitHub
, cmake
, boost
, itk
, glog
, ...
}:

stdenv.mkDerivation rec {
  name = "petmr-rd-tools";

  # src = fetchFromGitHub {
  #   owner = "UCL";
  #   repo = "petmr-rd-tools";
  #   rev = "b88281f";  # 20180307
  #   sha256 = "1pxa9j4hy5vwdarcilfxfn7kd51gcfp2jgh7m4h5bnvxdfb8i1lq";
  # };
  src = builtins.fetchGit { url = ./.; };  # respect .gitignore

  enableParallelBuilding = true;

  buildInputs = [ cmake boost itk glog ];
  cmakeBuildType = "Debug";

  meta = with stdenv.lib; {
    description = "Command line tools for PET-MR (pre)-processing";
    homepage = https://github.com/UCL/petmr-rd-tools;
    license = licenses.asl20;  # Apache 2.0
  };
}) pkgs
