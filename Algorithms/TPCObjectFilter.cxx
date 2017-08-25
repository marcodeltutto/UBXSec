#ifndef TPCOBJECTFILTER_CXX
#define TPCOBJECTFILTER_CXX

#include "TPCObjectFilter.h"
#include <iostream>

namespace ubana {

  TPCObjectFilter::TPCObjectFilter()
  {
    _tolerance = 50.;
    _debug = false;
  }

  void TPCObjectFilter::Configure(fhicl::ParameterSet const& pset)
  {
    _tolerance   = pset.get< double > ( "Tolerance", 50. );
    _debug       = pset.get< bool   > ( "DebugMode", false );
  }

  void TPCObjectFilter::PrintConfig() {

    std::cout << "--- TPCObjectFilter configuration:" << std::endl;
    std::cout << "---   _tolerance  = " << _tolerance << std::endl;

  }

  struct two_points {
    TVector3 pt1;
    TVector3 pt2;
    bool good;
    bool is_nu;
  };

  lar_pandora::PFParticleVector TPCObjectFilter::Filter(lar_pandora::PFParticleVector pfp_v, 
                                                        lar_pandora::PFParticlesToTracks pfp_to_tracks, 
                                                        lar_pandora::PFParticlesToShowers pfp_to_showers,
                                                        lar_pandora::PFParticlesToVertices  pfp_to_vertices) {
    

    // Convert PFP to a vector of two points, and find the nu vertex
    std::vector<two_points> tp_v;
    TVector3 nu_vtx;
    for (auto pfp : pfp_v) {

      if (_debug) std::cout << "[TPCObjectFilter] Preparing, pfp " << pfp->Self() << std::endl;

      two_points tp;
      tp.good = false;
      tp.is_nu = false;

      // Special case for neutrino pfp
      if (lar_pandora::LArPandoraHelper::IsNeutrino(pfp)) {
        if (_debug) std::cout << "[TPCObjectFilter] \t Is neutrino" << std::endl;
        tp.is_nu = true;
        tp_v.emplace_back(tp);
        auto iter = pfp_to_vertices.find(pfp);
        if (iter != pfp_to_vertices.end()) {
          double xyz[3];
          iter->second[0]->XYZ(xyz);
          nu_vtx.SetX(xyz[0]); nu_vtx.SetY(xyz[1]); nu_vtx.SetZ(xyz[2]);
        } else {
          std::cerr << "[TPCObjectFilter] Warning, can't find neutrino vertex." << std::endl;
          return pfp_v;
        }
        continue;
      }

      // Tracks
      auto iter1 = pfp_to_tracks.find(pfp);
      if (iter1 != pfp_to_tracks.end()) {
        if (_debug) std::cout << "[TPCObjectFilter] \t Found track for this pfp" << std::endl;
        tp.pt1 = (*iter1).second[0]->Vertex();
        tp.pt2 = (*iter1).second[0]->End(); 
        if (_debug) std::cout << "[TPCObjectFilter] \t tp.pt1 = " << tp.pt1.X() << " " << tp.pt1.Y() << " " << tp.pt1.Z() << ", tp.pt2 = " << tp.pt2.X() << " " << tp.pt2.Y() << " " << tp.pt2.Z() << std::endl;
        tp.good = true;
      }

      // Showers
      auto iter2 = pfp_to_showers.find(pfp);
      if (iter2 != pfp_to_showers.end()) {

        if (_debug) std::cout << "[TPCObjectFilter] \t Found shower for this pfp" << std::endl;
        TVector3 s_start = (*iter2).second[0]->ShowerStart();
        double s_len = (*iter2).second[0]->Length();
        TVector3 s_dir = (*iter2).second[0]->Direction();
        s_dir = s_dir.Unit();
        TVector3 s_end = s_start + s_dir * s_len;

        tp.pt1 = s_start;
        tp.pt2 = s_end;
        if (_debug) std::cout << "[TPCObjectFilter] \t tp.pt1 = " << tp.pt1.X() << " " << tp.pt1.Y() << " " << tp.pt1.Z() << ", tp.pt2 = " << tp.pt2.X() << " " << tp.pt2.Y() << " "
           << tp.pt2.Z() << std::endl;
        tp.good = true;
      }

      tp_v.emplace_back(tp);
    }


    std::vector<int> remove_v; // 0: remove, 1: don't remove, 2: todo
    remove_v.resize(tp_v.size());
    for (auto & rm : remove_v) rm = 2;

    std::vector<two_points> tp_close_to_nuvtx;


    // Now try to understand if there are separate PFPs
    bool do_loop = true;
    int cycle_counter = 0;

    while (do_loop) {

      for (size_t i = 0; i < tp_v.size(); i++) {

        // Skip the pfp if already categorized
        if (remove_v.at(i) != 2) continue;

        if (_debug) std::cout << "[TPCObjectFilter] Looking at tp " << i << std::endl;  
          
        if(tp_v.at(i).is_nu) {
          remove_v.at(i) = 1;
          if (_debug) std::cout << "[TPCObjectFilter] \t Is neutrino, will keep " << std::endl;
          continue;
        }

        if(!tp_v.at(i).good) {
          remove_v.at(i) = 0;
          if (_debug) std::cout << "[TPCObjectFilter] \t Is not good, will be removed " << std::endl;
          continue;
        }

        // Check if this pfp is close to the neutrino vertex
        if ( (tp_v.at(i).pt1 - nu_vtx).Mag() < _tolerance
          || (tp_v.at(i).pt2 - nu_vtx).Mag() < _tolerance ) {
          if (_debug) std::cout << "[TPCObjectFilter] \t Is close to nu vtx" << std::endl;
          tp_close_to_nuvtx.emplace_back(tp_v.at(i));
          remove_v.at(i) = 1;
          continue;
        }

        // Otherwise, check if it's close to one of the pfp close to the nu vtx
        for (auto tp : tp_close_to_nuvtx) {

          if ( (tp_v.at(i).pt1 - tp.pt1).Mag() < _tolerance
            || (tp_v.at(i).pt1 - tp.pt2).Mag() < _tolerance
            || (tp_v.at(i).pt2 - tp.pt1).Mag() < _tolerance
            || (tp_v.at(i).pt2 - tp.pt2).Mag() < _tolerance ) {

            if (_debug) std::cout << "[TPCObjectFilter] \t Is close to a pfp close to nu vtx " << std::endl;
            tp_close_to_nuvtx.emplace_back(tp_v.at(i));
            remove_v.at(i) = 1;
          }
        }
      }

      // Check if there are pfp not tagged
      do_loop = false;
      for (auto rm : remove_v) {
        if (rm == 2) do_loop = true;
      }
      cycle_counter++;

      if (cycle_counter > factorial(tp_v.size()) ) do_loop = false;
    } // end while loop

/*
      // Otherwise check is it's closer to anothe pfp
      if (!is_close_to_another) {
        for (size_t j = 0; j < tp_v.size() && i != j; j++) {
          std::cout << "here 2" << std::endl;
          std::cout << "tolerance = " << _tolerance << std::endl;
          std::cout << "1) " << (tp_v.at(i).pt1 - tp_v.at(j).pt1).Mag() << std::endl;
          std::cout << "2) " << (tp_v.at(i).pt1 - tp_v.at(j).pt2).Mag() << std::endl;
          std::cout << "3) " << (tp_v.at(i).pt2 - tp_v.at(j).pt1).Mag() << std::endl;
          std::cout << "4)" << (tp_v.at(i).pt2 - tp_v.at(j).pt2).Mag() << std::endl;
          if ( (tp_v.at(i).pt1 - tp_v.at(j).pt1).Mag() < _tolerance  
            || (tp_v.at(i).pt1 - tp_v.at(j).pt2).Mag() < _tolerance 
            || (tp_v.at(i).pt2 - tp_v.at(j).pt1).Mag() < _tolerance 
            || (tp_v.at(i).pt2 - tp_v.at(j).pt2).Mag() < _tolerance ) {

            is_close_to_another = true;
            std::cout << "[TPCObjectFilter] \t Is close to another, tolerance = " << _tolerance << std::endl;
            break; 
          }
        }
      }

      if (!is_close_to_another)
        remove_v.at(i) = true;
    }
*/

    // Instantiate the output
    lar_pandora::PFParticleVector out_pfp_v;

    for (size_t i = 0; i < pfp_v.size(); i++) {

      if (_debug) std::cout << "[TPCObjectFilter] Looking at pfp " << pfp_v.at(i)->Self() << std::endl;
      if (remove_v.at(i) == 1) {
        if (_debug) std::cout << "[TPCObjectFilter] \t It made it" << std::endl;
        out_pfp_v.emplace_back(pfp_v.at(i));
      }
    }

    return out_pfp_v;

  }

}


#endif
