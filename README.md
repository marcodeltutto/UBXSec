# UBXSec

[![Build Status](https://travis-ci.org/marcodeltutto/UBXSec.svg?branch=master)](https://travis-ci.org/marcodeltutto/UBXSec)

Wiki Page: https://github.com/marcodeltutto/UBXSec/wiki

Doxygen Documentation: https://marcodeltutto.github.io/UBXSec/documentation/html/annotated.html 

Clone with `--recursive` option.

Works with MCC8.6.

## How to install

On uboonegpvm, from a fresh terminal session:

```
export UBCODE_RELEASE=v06_26_01_10
export QUAL1=e10
export QUAL2=prof

source /grid/fermiapp/products/uboone/setup_uboone.sh

setup mrb
export MRB_PROJECT=larsoft

mkdir mylarsoft
setup uboonecode $UBCODE_RELEASE -q ${QUAL1}:${QUAL2}
cd mylarsoft
mrb newDev
source localProducts_larsoft_${UBCODE_RELEASE}_${QUAL1}_${QUAL2}/setup
cd srcs
mrb g -t $UBCODE_RELEASE uboonecode
cd uboonecode/uboone/
git clone --recursive git@github.com:marcodeltutto/UBXSec.git
[add UBXSec to the CMakeLists.txt]
cd ${MRB_TOP}/build_slf6.x86_64
mrbsetenv
mrb i
cd ..
mrbslp
```

## How to run

To run over a BNB+Cosmic file:
`lar -c run_ubxsec_mc_bnbcosmic.fcl -s art_root_file.root`

To run over a BNBON data file:
`lar -c run_ubxsec_data_bnbon.fcl -s art_root_file.root`

To run over a EXTBNB data file:
`lar -c run_ubxsec_data_extbnb.fcl -s art_root_file.root`


