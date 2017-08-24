#ifndef TPCOBJECTFILTER_CXX
#define TPCOBJECTFILTER_CXX

#include "TPCObjectFilter.h"
#include <iostream>

namespace ubana {

  TPCObjectFilter::TPCObjectFilter()
  {
    _tolerance = 10.;
  }

  void TPCObjectFilter::Configure(fhicl::ParameterSet const& pset)
  {
    _tolerance   = pset.get< double > ( "Tolerance", 10. );
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
                                                        lar_pandora::PFParticlesToShowers pfp_to_showers) {


    // Convert PFP to a vector of two points
    std::vector<two_points> tp_v;
    for (auto pfp : pfp_v) {

      two_points tp;
      tp.good = false;
      
      if (lar_pandora::LArPandoraHelper::IsNeutrino(pfp)) {
        tp.is_nu = true;
        tp_v.emplace_back(tp);
        continue;
      }

      // Tracks
      auto iter1 = pfp_to_tracks.find(pfp);
      if (iter1 != pfp_to_tracks.end()) {
        tp.pt1 = (*iter1).second[0]->Vertex();
        tp.pt2 = (*iter1).second[0]->End(); 
        tp.good = true;
      }

      // Showers
      auto iter2 = pfp_to_showers.find(pfp);
      if (iter2 != pfp_to_showers.end()) {

        TVector3 s_start = (*iter2).second[0]->ShowerStart();
        double s_len = (*iter2).second[0]->Length();
        TVector3 s_dir = (*iter2).second[0]->Direction();
        s_dir = s_dir.Unit();
        TVector3 s_end = s_start + s_dir * s_len;

        tp.pt1 = s_start;
        tp.pt2 = s_end;
        tp.good = true;
      }

      tp_v.emplace_back(tp);
    }


    std::vector<bool> remove_v;
    remove_v.resize(tp_v.size());

    for (size_t i = 0; i < tp_v.size(); i++) {

      if(tp_v.at(i).is_nu) {
        remove_v.at(i) = false;
        continue;
      }

      if(!tp_v.at(i).good) {
        remove_v.at(i) = true;
        continue;
      }

      bool is_close_to_another = false;
      for (size_t j = 0; j < tp_v.size() && i != j; j++) {

        if ( (tp_v.at(i).pt1 - tp_v.at(j).pt1).Mag() < _tolerance  
          || (tp_v.at(i).pt1 - tp_v.at(j).pt2).Mag() < _tolerance 
          || (tp_v.at(i).pt2 - tp_v.at(j).pt1).Mag() < _tolerance 
          || (tp_v.at(i).pt2 - tp_v.at(j).pt2).Mag() < _tolerance ) {

          is_close_to_another = true;
          break; 
        }
      }

      if (!is_close_to_another)
        remove_v.at(i) = true;
    }

    lar_pandora::PFParticleVector out_pfp_v;
    for (size_t i = 0; i < pfp_v.size(); i++) {

      if (!remove_v.at(i)) {
        out_pfp_v.emplace_back(pfp_v.at(i));
      }
    }

    return out_pfp_v;

  }

}


#endif
