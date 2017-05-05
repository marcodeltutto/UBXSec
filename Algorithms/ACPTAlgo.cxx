#ifndef ACPTALGO_CXX
#define ACPTALGO_CXX

#include "ACPTAlgo.h"
#include <iostream>

namespace ubxsec {

  ACPTAlgo::ACPTAlgo()
  {
    _tpcObjIsSet    = false;
    _vtxIsSet       = false;
    _max_distance   = 10.;
  }

  ACPTAlgo::ACPTAlgo(lar_pandora::TrackVector tpcObj, recob::Vertex vtx) {
    _tpcObjIsSet = _vtxIsSet = true;
    _track_v = tpcObj;
    _vtx = vtx;
    _max_distance   = 10.;
  }

  void ACPTAlgo::Configure(fhicl::ParameterSet const& pset)
  {
    _max_distance   = pset.get< double > ( "MaxDistance" );
  }

  void ACPTAlgo::PrintConfig() {

    std::cout << "--- ACPTAlgo configuration:" << std::endl;
    std::cout << "---   _max_distance  = " << _max_distance << std::endl;

  }

  void ACPTAlgo::SetTPCObj(lar_pandora::TrackVector tpcObj) {
    _tpcObjIsSet = true;
    _track_v = tpcObj;
  }

  void ACPTAlgo::SetVtx(recob::Vertex vtx){
    _vtxIsSet = true;
    _vtx = vtx;
  }

  void ACPTAlgo::Clear(){
    _vtxIsSet = _tpcObjIsSet = false;
  }

  double ACPTAlgo::AngleBetweenLongestTracks() {

    if (!_tpcObjIsSet || !_vtxIsSet){
      std::cerr << "The TPC object or the vertex was not set. Exiting now." << std::endl;
      exit(0);
    }

    // Get vertex location
    _vtx.XYZ(xyz);

    keep_track.clear();
    track_dir.clear();

    // Loop over tracks 
    for (size_t t = 0; t < _track_v.size(); t++) {
      start = _track_v[t]->Vertex();
      end   = _track_v[t]->End();
      
      dist = GetDistance(xyz, start);
      if (dist < _max_distance){
        keep_track.emplace_back(*(_track_v[t]));
        track_dir.emplace_back(_track_v[t]->VertexDirection());
        continue;
      }

      dist = GetDistance(xyz, end);
      if (dist < _max_distance){
        keep_track.emplace_back(*(_track_v[t]));
        track_dir.emplace_back(-(_track_v[t]->EndDirection()));
      }  
    }

    // Sort tracks by length
    std::sort(keep_track.begin(), keep_track.end(),
              [](recob::Track a, recob::Track b) -> bool
              {
                return a.Length() > b.Length();
              });

    if (keep_track.size() < 2) return -9999;

    // Calculate angle between two longest tracks
    return track_dir[0].Angle(track_dir[1]);  
  }




  double ACPTAlgo::GetDistance(double *x1, TVector3 x2) {

    TVector3 x1v(x1[0], x1[1], x1[2]);   
    return (x1v - x2).Mag();

  }

}


#endif
