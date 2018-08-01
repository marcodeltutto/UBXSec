/**
 * \class CosmicTagToolInterface
 *
 * \ingroup UBXSec
 *
 * \brief Interface for art module
 * 
 *
 * \author $Author: Marco Del Tutto<marco.deltutto@physics.ox.ac.uk> $
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2017/08/18 $
 *
 * Contact: marco.deltutto@physics.ox.ac.uk
 *
 * Created on: Friday, August 18, 2017 at 08:53:50
 *
 */


#ifndef COSMICTAGTOOLINTERFACE_H
#define COSMICTAGTOOLINTERFACE_H

#include <iostream>
#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
//#include "uboone/UBXSec/DataTypes/UBXSecFMWKInterface.h"

namespace ubana {
  using SimpleHit_t = struct SimpleHit {
    double time;
    unsigned int wire;
    double integral;
  };

  using SimpleHitVector = std::vector<ubana::SimpleHit_t>;
}


namespace ubana {

  class CosmicTagToolInterface {
  
  public:

    /// Default destructor
    virtual ~CosmicTagToolInterface() noexcept = default;

    /// Performs the matching
    virtual bool IsCosmic(SimpleHitVector) = 0;

  };
}


#endif //  COSMICTAGTOOLINTERFACE_H
