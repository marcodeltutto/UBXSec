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

    /// Check track is consistent with muon hypothesis
    bool SVMPredict(double dqds, double length);

  protected:

    /// SVM Polynomial Kernel
    double _SVM_kernel(std::pair<double, double> support_vector, std::pair<double, double> user_vector);

    /// Feature Scaling
    std::pair<double, double>  _SVM_feature_scaling(double dqds, double length);

    std::vector<art::Ptr<recob::Track>> _tracks;
    std::map<art::Ptr<recob::Track>,art::Ptr<anab::ParticleID>> _track_to_pid;

    bool _use_pida_cut;

    double _gain_calib;

    bool _tracks_are_set;
    bool _tracktopidmap_is_set;

    std::vector<double> _svm_x;
    std::vector<double> _alpha;     
    std::vector<double> _support_vectors_x;
    std::vector<double> _support_vectors_y;
    double _rho;
    double _loop;
    double _gamma;
    double _degree;
    double _r;
    double _x1_mean;
    double _x2_mean;
    double _x1_range;
    double _x2_range;

  };
}

#endif
/** @} */ // end of doxygen group 

