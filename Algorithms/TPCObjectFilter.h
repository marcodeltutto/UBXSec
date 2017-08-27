/**
 * \file TPCObjectFilter.h
 *
 * \ingroup UBXSec
 * 
 * \brief Class def header for a class TPCObjectFilter
 *
 * @author Marco Del Tutto
 */

/** \addtogroup UBXSec

    @{*/
#ifndef TPCOBJECTFILTER_H
#define TPCOBJECTFILTER_H

#include <iostream>
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

int factorial(int n) {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

namespace ubana{
  
  /**
   \class TPCObjectFilter
   User defined class TPCObjectFilter ... these comments are used to generate
   doxygen documentation!
 */

  class TPCObjectFilter {
    
  public:
    
    /// Default constructor
    TPCObjectFilter();

    /// Default destructor
    ~TPCObjectFilter(){}

    /// Configure function parameters
    void Configure(fhicl::ParameterSet const& p);

    /// Printd the current configuration
    void PrintConfig();

    /// Returns a filtered vector of PFParticles: the new TPCObject 
    lar_pandora::PFParticleVector Filter(lar_pandora::PFParticleVector pfp_v, lar_pandora::PFParticlesToTracks pfp_to_tracks, lar_pandora::PFParticlesToShowers pfp_to_showers, lar_pandora::PFParticlesToVertices pfp_to_vertices);

  protected:

    double _tolerance; ///< Tolerance for reclustering
    bool _debug;       ///< Enabels debug printouts
  };
}

#endif
/** @} */ // end of doxygen group 

