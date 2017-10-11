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
    _use_pida_cut   = pset.get< bool > ( "UsePIDACut", false );
  }

  void MuonCandidateFinder::PrintConfig() {

    std::cout << "--- MuonCandidateFinder configuration:" << std::endl;
    std::cout << "---   _use_pida_cut  = " << _use_pida_cut << std::endl;

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

}


#endif
