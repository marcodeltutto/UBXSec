#ifndef MUONCANDIDATEFINDER_CXX
#define MUONCANDIDATEFINDER_CXX

#include "MuonCandidateFinder.h"
#include <iostream>

namespace ubxsec {

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

    return true;
  }

}


#endif
