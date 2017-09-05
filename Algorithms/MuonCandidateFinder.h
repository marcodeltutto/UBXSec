/**
 * \file MuonCandidateFinder.h
 *
 * \ingroup ubana
 * 
 * \brief Class def header for a class MuonCandidateFinder
 *
 * @author Marco Del Tutto
 */

/** \addtogroup ubana

    @{*/
#ifndef MUONCANDIDATEFINDER_H
#define MUONCANDIDATEFINDER_H

#include <iostream>
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"

namespace ubana{
  
  /**
   \class MuonCandidateFinder
   User defined class MuonCandidateFinder ... these comments are used to generate
   doxygen documentation!
 */

  class MuonCandidateFinder {
    
  public:
    
    /// Default constructor
    MuonCandidateFinder();

    /// Default destructor
    ~MuonCandidateFinder(){}

    /// Configure function parameters
    void Configure(fhicl::ParameterSet const& p);

    /// Print the current configuration
    void PrintConfig();

    /// Set the TPC object
    void SetTPCObject(ubana::TPCObject theobject) {_tpc_object = theobject; _tpcobject_is_set = true;}

    /// Set map from tracks to PID for the TPC object
    void SetTrackToPIDMap(std::map<recob::Track,anab::ParticleID> themap) {_track_to_pid = themap; _tracktopidmap_is_set = true;} 

    /// Sets the TPC object
    bool GetCandidateTrack(recob::Track &);
    
  protected:

    ubana::TPCObject _tpc_object;
    std::map<recob::Track,anab::ParticleID> _track_to_pid;

    double _max_distance;
    bool _tpcobject_is_set;
    bool _tracktopidmap_is_set;
  };
}

#endif
/** @} */ // end of doxygen group 

