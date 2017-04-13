#ifndef VERTEXCHECK_CXX
#define VERTEXCHECK_CXX

#include "VertexCheck.h"
#include <iostream>

namespace ubxsec {

  VertexCheck::VertexCheck()
  {
    _tpcObjIsSet    = false;
    _vtxIsSet       = false;
    _max_distance    = 10.;
  }

  void VertexCheck::Configure(fhicl::ParameterSet const& pset)
  {
    _max_distance   = pset.get< double > ( "MaxDistance"   );
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

  double VertexCheck::AngleBetweenLongestTracks() {

    if (!_tpcObjIsSet || !_vtxIsSet){
      std::cerr << "The TPC object or the vertex was not set. Exiting now." << std::endl;
      exit(0);
    }

    // Get vertex location
    double xyz[3];
    _vtx.XYZ(xyz);

    std::cout << "VERTEX LOCATION x=" << xyz[0] << " y=" << xyz[1] << " z=" << xyz[2] << std::endl;
    TVector3 start, end;
    double dist;
    std::vector<recob::Track> keep_track;
    keep_track.clear();
    std::vector<TVector3> track_dir;
    track_dir.clear();

    std::cout << "Number of tracks in this TPC object: " << _track_v.size() << std::endl;

    // Loop over tracks 
    for (size_t t = 0; t < _track_v.size(); t++) {
      std::cout << "  Looping over track " << t << std::endl;
      start = _track_v[t]->Vertex();
      end   = _track_v[t]->End();
      std::cout << "  TRACK START x=" << start.X() << " y=" << start.Y() << " z=" << start.Z() << std::endl;
      std::cout << "  TRACK END x=" << end.X() << " y=" << end.Y() << " z=" << end.Z() << std::endl;
      
      dist = GetDistance(xyz, start);
      std::cout << "    Distance from track start is " << dist << std::endl;
      if (dist < _max_distance){
        keep_track.emplace_back(*(_track_v[t]));
        track_dir.emplace_back(_track_v[t]->VertexDirection());
        continue;
      }

      dist = GetDistance(xyz, end);
      std::cout << "    Distance from track end is " << dist << std::endl;    
      if (dist < _max_distance){
        keep_track.emplace_back(*(_track_v[t]));
        track_dir.emplace_back(-(_track_v[t]->EndDirection()));
      }  
    }

    for (size_t i = 0; i < keep_track.size(); i++) {
      std::cout << std::endl << "  Track length id " << keep_track[i].Length() << std::endl; 
    }

    // Sort tracks by length
    std::sort(keep_track.begin(), keep_track.end(),
              [](recob::Track a, recob::Track b) -> bool
              {
                return a.Length() > b.Length();
              });

    for (size_t i = 0; i < keep_track.size(); i++) {
      std::cout << std::endl << "  Track length id " << keep_track[i].Length() << std::endl;
    }

    if (keep_track.size() < 2) return -9999;

    // Calculate angle between two longest tracks
    return track_dir[0].Angle(track_dir[1]);  
  }

  double VertexCheck::GetDistance(double *x1, TVector3 x2) {

    TVector3 x1v(x1[0], x1[1], x1[2]);
    
    return (x1v - x2).Mag();
  }
}


#endif
