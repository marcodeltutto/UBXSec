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
#include <math.h>
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

    /// Reset flags
    void Reset();

    /// Set the tracks 
    void SetTracks(std::vector<art::Ptr<recob::Track>> thetracks) {_tracks = thetracks; _tracks_are_set = true;}

    /// Set map from tracks to PID for the TPC object
    void SetTrackToPIDMap(std::map<art::Ptr<recob::Track>,art::Ptr<anab::ParticleID>> themap) {_track_to_pid = themap; _tracktopidmap_is_set = true;} 

    /// Sets the TPC object
    bool GetCandidateTrack(art::Ptr<recob::Track> &);

    /// Check track is consistent with muon hypothesis
    bool MIPConsistency(double dqds, double length);

  protected:

    std::vector<art::Ptr<recob::Track>> _tracks;
    std::map<art::Ptr<recob::Track>,art::Ptr<anab::ParticleID>> _track_to_pid;

    bool _use_pida_cut;

    bool _tracks_are_set;
    bool _tracktopidmap_is_set;

    std::vector<double> _svm_x;

  };
}

#endif
/** @} */ // end of doxygen group 

