/**
 * \file VertexCheck.h
 *
 * \ingroup UBXSec
 * 
 * \brief Class def header for a class VertexCheck
 *
 * @author Marco Del Tutto
 */

/** \addtogroup UBXSec

    @{*/
#ifndef VERTEXCHECK_H
#define VERTEXCHECK_H

#include <iostream>
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

namespace ubxsec{
  
  /**
   \class VertexCheck
   User defined class VertexCheck ... these comments are used to generate
   doxygen documentation!
 */

  class VertexCheck {
    
  public:
    
    /// Default constructor
    VertexCheck();
    
    /// Default destructor
    ~VertexCheck(){}

    /// Configure function parameters
    void Configure(fhicl::ParameterSet const& p);

    /// Printd the current configuration
    void PrintConfig();

    /// Sets the TPC object
    void SetTPCObj(lar_pandora::TrackVector);
    
    /// Sets the TPC object
    void SetVtx(recob::Vertex);

    /// Returns the distance between two points
    double GetDistance(double *, TVector3);

    /// Returns the angle between the two longest tracks in the TPC obj
    double AngleBetweenLongestTracks();

  protected:

    lar_pandora::TrackVector _track_v; ///< the TPC object
    recob::Vertex            _vtx;     ///< the vertex
    double _max_distance;
    bool _tpcObjIsSet, _vtxIsSet;
  };
}

#endif
/** @} */ // end of doxygen group 

