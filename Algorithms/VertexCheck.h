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

    /// Construct setting TPC object and Vertex
    VertexCheck(lar_pandora::TrackVector, recob::Vertex);

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

    /// Restores flags
    void Clear();

    /// Returns the distance between two points (the first a 3d array, the second a TVector3...I know...)
    double GetDistance(double *, TVector3);

    /// Returns the angle between the two longest tracks in the TPC obj
    double AngleBetweenLongestTracks();

  protected:

    lar_pandora::TrackVector _track_v; ///< the TPC object
    recob::Vertex            _vtx;     ///< the vertex
    double _max_distance;
    bool _tpcObjIsSet, _vtxIsSet;


    double xyz[3]; ///< will store vertex location
    TVector3 start, end; ///< will store track start and end
    double dist; ///< will store track distance
    std::vector<recob::Track> keep_track; ///< will store tracks to be kept for calculation

  };
}

#endif
/** @} */ // end of doxygen group 

