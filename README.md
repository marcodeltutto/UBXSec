# UBXSec

[![Build Status](https://travis-ci.org/marcodeltutto/UBXSec.svg?branch=develop)](https://travis-ci.org/marcodeltutto/UBXSec)

Doxygen Documentation: https://marcodeltutto.github.io/UBXSec/documentation/html/index.html

## How it works

On uboonegpvm, from a fresh terminal session:

```
export UBCODE_RELEASE=v06_26_01_01
export QUAL1=e10
export QUAL2=prof

source /grid/fermiapp/products/uboone/setup_uboone.sh

setup git
setup gitflow
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
git clone git@github.com:marcodeltutto/UBXSec.git
[add UBXSec to the CMakeLists.txt]
cd ${MRB_TOP}/build_slf6.x86_64
mrbsetenv
mrb i
cd ..
mrbslp
```


