#ifndef MUONCANDIDATEFINDER_CXX
#define MUONCANDIDATEFINDER_CXX

#include "MuonCandidateFinder.h"
#include <iostream>

namespace ubana {

  MuonCandidateFinder::MuonCandidateFinder()
  {
    _max_distance   = 10.;
    _tpcobject_is_set = false;
    _tracktopidmap_is_set = false;
  }

  void MuonCandidateFinder::Configure(fhicl::ParameterSet const& pset)
  {
    _max_distance   = pset.get< double > ( "MaxDistance" );
  }

  void MuonCandidateFinder::PrintConfig() {

    std::cout << "--- MuonCandidateFinder configuration:" << std::endl;
    std::cout << "---   _max_distance  = " << _max_distance << std::endl;

  }

  bool MuonCandidateFinder::GetCandidateTrack(recob::Track & out_track) {
   
    if (!_tracktopidmap_is_set || !_tpcobject_is_set) {
      std::cerr << "[MuonCandidateFinder] Call to GetCandidateTrack without having set TPCObject and map to PIDs. Exiting now." << std::endl;
      throw std::exception();
    }

    std::vector<recob::Track> tracks = _tpc_object.GetTracks();

    if (tracks.size() == 0)
      return false;

    std::vector<recob::Track> possible_tracks;

    for (auto track : tracks) {

      double length, pida;

      length = track.Length();

      auto iter = _track_to_pid.find(track);
      if(iter != _track_to_pid.end()) {
        pida = iter->second.PIDA();
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
    
    // Otherwise sort remaining tracks by length
    std::sort(possible_tracks.begin(), possible_tracks.end(),
              [](recob::Track a, recob::Track b) -> bool
              {
                return a.Length() > b.Length();
              });

    out_track = possible_tracks.at(0);

    return true;
  }

}


#endif
