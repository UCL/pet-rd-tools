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
  name = "petmr-rd-tools-master";

  # src = fetchFromGitHub {
  #   owner = "UCL";
  #   repo = "petmr-rd-tools";
  #   rev = "b88281f";  # 20180307
  #   sha256 = "1pxa9j4hy5vwdarcilfxfn7kd51gcfp2jgh7m4h5bnvxdfb8i1lq";
  # };
  src = ./.;

  enableParallelBuilding = true;

  buildInputs = [ cmake boost itk glog ];
  cmakeBuildType = "Debug";

  meta = {
    description = "Command line tools for PET-MR (pre)-processing";
    homepage = https://github.com/UCL/petmr-rd-tools;
  };
}) pkgs
