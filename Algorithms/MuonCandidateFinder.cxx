#ifndef MUONCANDIDATEFINDER_CXX
#define MUONCANDIDATEFINDER_CXX

#include "MuonCandidateFinder.h"
#include <iostream>

namespace ubana {

  MuonCandidateFinder::MuonCandidateFinder()
  {
    _use_pida_cut   = false;
    _tracks_are_set = false;
    _tracktopidmap_is_set = false;
  }

  void MuonCandidateFinder::Configure(fhicl::ParameterSet const& pset)
  {
    _use_pida_cut   = pset.get<bool> ( "UsePIDACut", false );
    _svm_x          = pset.get<std::vector<double>> ("SVM_X");

    if (_svm_x.size() != 1000) {
      std::cout << "[MuonCandidateFinder] _svm_x size error, is " << _svm_x.size() << ", should be 1000" << std::endl;
      throw std::exception();
    }
  }

  void MuonCandidateFinder::PrintConfig() {

    std::cout << "--- MuonCandidateFinder configuration:" << std::endl;
    std::cout << "---   _use_pida_cut  = " << _use_pida_cut << std::endl;

  }

  void MuonCandidateFinder::Reset() 
  {
    _tracks_are_set = false;
    _tracktopidmap_is_set = false;    
  }

  bool MuonCandidateFinder::GetCandidateTrack(art::Ptr<recob::Track> & out_track) {
   
    if (!_tracktopidmap_is_set || !_tracks_are_set) {
      std::cerr << "[MuonCandidateFinder] Call to GetCandidateTrack without having set input Tracks and map to PIDs. Exiting now." << std::endl;
      throw std::exception();
    }

    if (_tracks.size() == 0)
      return false;

    std::vector<art::Ptr<recob::Track>> possible_tracks;

    if (_use_pida_cut) {
      for (auto track : _tracks) {

        double length, pida;

        length = track->Length();

        auto iter = _track_to_pid.find(track);
        if(iter != _track_to_pid.end()) {
          pida = iter->second->PIDA();
        } else {
          possible_tracks.emplace_back(track);
          continue;
        }

        if (length < 90. && pida > 13) {
          continue;
        } else {
          possible_tracks.emplace_back(track);
        }
      }

      if (possible_tracks.size() == 0) {
        return false;
      }else if (possible_tracks.size() == 1) {
        out_track = possible_tracks.at(0);
        return true;
      }
    }
    else {
      possible_tracks = _tracks;
    }

    // Sort remaining tracks by length
    std::sort(possible_tracks.begin(), possible_tracks.end(),
              [](art::Ptr<recob::Track> a, art::Ptr<recob::Track> b) -> bool
              {
                return a->Length() > b->Length();
              });

    out_track = possible_tracks.at(0);

    return true;
  }



  bool MuonCandidateFinder::MIPConsistency(double dqds, double length) {

    if (length > 1000)
      return true;

    if (length < 0) {
      std::cout << "[MuonCandidateFinder] Track length < 0?!" << std::endl;
      return false;
    }

    int l = std::round(length);
    double dqds_cut = _svm_x.at(l);

    std::cout << "[MuonCandidateFinder] Track length is " << length << ", dqds_cut is " << dqds_cut << std::endl;

    if (dqds <= dqds_cut)
      return true;
  

    return false;

  }

  /*
  bool MuonCandidateFinder::SVMPredict(double dqds, double length) {

    double gamma = 0.1;
    double degree = 2.;
    double r = 0.;

    double rho = 14.96135081;

    int loop = 374;

  }
  */
}


#endif
