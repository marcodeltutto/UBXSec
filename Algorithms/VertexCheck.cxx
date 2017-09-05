#ifndef VERTEXCHECK_CXX
#define VERTEXCHECK_CXX

#include "VertexCheck.h"
#include <iostream>

namespace ubxsec {

  VertexCheck::VertexCheck()
  {
    _tpcObjIsSet    = false;
    _vtxIsSet       = false;
    _max_distance   = 10.;
  }

  VertexCheck::VertexCheck(lar_pandora::TrackVector tpcObj, recob::Vertex vtx) {
    _tpcObjIsSet = _vtxIsSet = true;
    _track_v = tpcObj;
    _vtx = vtx;
    _max_distance   = 10.;
  }

  void VertexCheck::Configure(fhicl::ParameterSet const& pset)
  {
    _max_distance   = pset.get< double > ( "MaxDistance" );
  }

  void VertexCheck::PrintConfig() {

    std::cout << "--- VertexCheck configuration:" << std::endl;
    std::cout << "---   _max_distance  = " << _max_distance << std::endl;

  }

  void VertexCheck::SetTPCObj(lar_pandora::TrackVector tpcObj) {
    _tpcObjIsSet = true;
    _track_v = tpcObj;
  }

  void VertexCheck::SetVtx(recob::Vertex vtx){
    _vtxIsSet = true;
    _vtx = vtx;
  }

  void VertexCheck::Clear(){
    _vtxIsSet = _tpcObjIsSet = false;
  }

  double VertexCheck::AngleBetweenLongestTracks() {

    if (!_tpcObjIsSet || !_vtxIsSet){
      std::cerr << "[VertexCheck] The TPC object or the vertex was not set. Exiting now." << std::endl;
      throw std::exception();
    }

    // Get vertex location
    _vtx.XYZ(xyz);

    keep_track.clear();

    // Loop over tracks 
    for (size_t t = 0; t < _track_v.size(); t++) {
      start = _track_v[t]->Vertex();
      end   = _track_v[t]->End();
      
      dist = GetDistance(xyz, start);
      if (dist < _max_distance){
        keep_track.emplace_back(*(_track_v[t]));
        continue;
      }

      dist = GetDistance(xyz, end);
      if (dist < _max_distance){
        keep_track.emplace_back(*(_track_v[t]));
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
    TVector3 dir1 = keep_track.at(0).VertexDirection();
    TVector3 dir2 = keep_track.at(1).VertexDirection();

    return dir1.Angle(dir2);
  }




  double VertexCheck::GetDistance(double *x1, TVector3 x2) {

    TVector3 x1v(x1[0], x1[1], x1[2]);   
    return (x1v - x2).Mag();

  }

}


#endif
