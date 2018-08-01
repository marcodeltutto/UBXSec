/**
 * \class RecoTrueMatchInterface
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


#ifndef RECOTRUEMATCHINTERFACE_H
#define RECOTRUEMATCHINTERFACE_H

#include <iostream>
#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
//#include "uboone/UBXSec/DataTypes/UBXSecFMWKInterface.h"


namespace ubxsec {

  class RecoTrueMatchInterface {
  
  public:

    /// Default destructor
    virtual ~RecoTrueMatchInterface() noexcept = default;

    /// Performs the matching
    virtual void GetRecoToTrueMatches(lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits) = 0;
}

#endif //  RECOTRUEMATCHINTERFACE_H
