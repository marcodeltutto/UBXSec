#ifndef COSMICTAGBYHITINTEGRAL_CXX
#define COSMICTAGBYHITINTEGRAL_CXX

#include "CosmicTagByHitIntegral.h"
#include <iostream>

namespace ubana {

  bool CosmicTagByHitIntegral::IsCosmic(SimpleHitVector simple_hit_v) {

    std::cout << ">>>>>>>>>>>>>>>>> Hi! Test is " << _test << std::endl;

    for (auto s_hit : simple_hit_v){
      //std::cout << "Hit Wire " << s_hit.wire << ", Time " << s_hit.time << std::endl;
      _csvfile << _counter << "," << s_hit.wire << "," << s_hit.integral << std::endl;
    }

    _counter++;
    return true;
  }

}


#endif
