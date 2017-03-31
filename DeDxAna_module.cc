////////////////////////////////////////////////////////////////////////
// Class:       DeDxAna
// Plugin Type: analyzer (art v2_05_00)
// File:        DeDxAna_module.cc
//
// Generated at Tue Mar 14 16:04:35 2017 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" 
#include "lardata/Utilities/AssociationUtil.h"
#include "larcore/Geometry/Geometry.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "uboone/UBXSec/UBXSecHelper.h"

#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "TString.h"
#include "TTree.h"

class DeDxAna;


class DeDxAna : public art::EDAnalyzer {
public:
  explicit DeDxAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DeDxAna(DeDxAna const &) = delete;
  DeDxAna(DeDxAna &&) = delete;
  DeDxAna & operator = (DeDxAna const &) = delete;
  DeDxAna & operator = (DeDxAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  TTree* _tree1;
  int _run, _subrun, _event;
  std::vector<std::vector<double>> _res_range, _dedx, _dqdx;
  std::vector<int> _trk_number;
  std::vector<std::vector<int>> _plane;

  std::vector<int> _hit_ch, _hit_localind, _hit_mult;
  std::vector<double> _hit_int, _hit_goodoffit, _hit_sigmapeaktime, _hit_rms;

  // Declare member data here.

};


DeDxAna::DeDxAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{

    art::ServiceHandle<art::TFileService> fs;
    _tree1 = fs->make<TTree>("tree","");
    _tree1->Branch("run",                &_run,                "run/I");
    _tree1->Branch("subrun",             &_subrun,             "subrun/I");
    _tree1->Branch("event",              &_event,              "event/I");
    _tree1->Branch("res_range",          "std::vector<std::vector<double>>",          &_res_range);
    _tree1->Branch("dedx",               "std::vector<std::vector<double>>",               &_dedx);
    _tree1->Branch("dqdx",               "std::vector<std::vector<double>>",               &_dqdx);
    _tree1->Branch("plane",              "std::vector<std::vector<int>>",                 &_plane);
    _tree1->Branch("trk_number",         "std::vector<int>",                         &_trk_number);

    _tree1->Branch("hit_ch", "std::vector<int>", & _hit_ch);
    _tree1->Branch("hit_int", "std::vector<double>", & _hit_int);
    _tree1->Branch("hit_mult", "std::vector<int>", & _hit_mult);
    _tree1->Branch("hit_localind", "std::vector<int>", & _hit_localind);
    _tree1->Branch("hit_goodoffit", "std::vector<double>", & _hit_goodoffit);
    _tree1->Branch("hit_sigmapeaktime", "std::vector<double>", & _hit_sigmapeaktime);
    _tree1->Branch("hit_rms", "std::vector<double>", & _hit_rms);



}

void DeDxAna::analyze(art::Event const & e) {

  art::Handle<std::vector<recob::Track>> track_h;
  e.getByLabel("pandoraNu",track_h);

  if(!track_h.isValid() || track_h->empty()) {
    std::cout << "Track handle is not valid or empty" << std::endl;
    return;
  }

  art::FindMany<anab::Calorimetry> calo_track_ass(track_h, e, "pandoraNucalo");

  if (!calo_track_ass.isValid() ){
    std::cout << "Calo is not valid." << std::endl;
    return;
  }
  std::cout << "calo has size " << calo_track_ass.size() << std::endl;

  art::FindMany<recob::Hit> track_hit_ass (track_h, e, "pandoraNu");
  if (!track_hit_ass.isValid() ){
    std::cout << "Hit track ass is not valid." << std::endl;
    return;
  } 
  std::cout << "track hit ass has size " << track_hit_ass.size() << std::endl;

  _res_range.resize(track_h->size()); 
  _dedx.resize(track_h->size());  
  _dqdx.resize(track_h->size()) ;
  _trk_number.resize(track_h->size()) ;
  _plane.resize(track_h->size()) ;

  ::art::ServiceHandle<geo::Geometry> geo;

  // Wire check
  double point_in_tpc[3];
  for (double zpos = 0; zpos < 1000; zpos += 0.3){

    point_in_tpc[0] = 100.;
    point_in_tpc[1] = 0.;
    point_in_tpc[2] = zpos;
    raw::ChannelID_t ch = geo->NearestChannel(point_in_tpc, 2);
    const lariov::ChannelStatusProvider& chanFilt = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    std::cout << ">>> this is ch " << ch << "  with status status " << chanFilt.Status(ch) << std::endl;
  }



  // Track loop
  for (unsigned int trk = 0; trk < track_h->size(); trk++) {

    std::cout << "***** Track " << trk << std::endl;
    //_trk_number(trk) = (int) trk;

    auto const& track = (*track_h)[trk];

    if (track.Length() < 20) continue;

    int vtx_ok;
    if(!UBXSecHelper::IsCrossingBoundary(track, vtx_ok)) continue;
    std::cout << "vtx_ok " << vtx_ok << std::endl;
    // Understand dead region
    double point_in_tpc[3];
    if (vtx_ok == 0) {
      point_in_tpc[0] = track.Vertex().X();
      point_in_tpc[1] = track.Vertex().Y();
      point_in_tpc[2] = track.Vertex().Z();
    } else if (vtx_ok == 1) {
      point_in_tpc[0] = track.End().X();
      point_in_tpc[1] = track.End().Y();
      point_in_tpc[2] = track.End().Z();
    }
    // Get nearest channel
    raw::ChannelID_t ch = geo->NearestChannel(point_in_tpc, 2);
    std::cout << "nearest channel is " << ch << std::endl;
    
    // Get channel status
    const lariov::ChannelStatusProvider& chanFilt = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    std::cout << "ch status " << chanFilt.Status(ch) << std::endl;
    if( chanFilt.Status(ch) < 3) continue;

    // Now check close wires
    double new_point_in_tpc[3];
    bool bad = false;
    for (int z_off = -5; z_off < 5; z_off += 2.5) {
      new_point_in_tpc[0] = point_in_tpc[0];
      new_point_in_tpc[1] = point_in_tpc[1];
      new_point_in_tpc[2] = point_in_tpc[2] + z_off;
      std::cout << "trying with point " << new_point_in_tpc[0] << " " << new_point_in_tpc[1] << " " << new_point_in_tpc[2] << std::endl;
      
      raw::ChannelID_t ch = geo->NearestChannel(new_point_in_tpc, 2);
      if( chanFilt.Status(ch) < 3) bad = true;
    }

    if (bad) continue;

    std::vector<const anab::Calorimetry*> calos = calo_track_ass.at(trk);

    for (size_t ical = 0; ical<calos.size(); ++ical){
      if (!calos[ical]) continue;
      if (!calos[ical]->PlaneID().isValid) continue;
      int planenum = calos[ical]->PlaneID().Plane;
      if (planenum<0||planenum>2) continue;
      std::cout << "plane " << planenum << std::endl;
      std::cout << " ke " << calos[ical]->KineticEnergy() << std::endl;
      std::cout << " range " << calos[ical]->Range() << std::endl;

      const size_t NHits = calos[ical] -> dEdx().size();

      for(size_t iTrkHit = 0; iTrkHit < NHits; ++iTrkHit) {

        _dedx.at(trk).emplace_back((calos[ical]->dEdx())[iTrkHit]);
        _dqdx.at(trk).emplace_back((calos[ical]->dQdx())[iTrkHit]);
        _res_range.at(trk).emplace_back((calos[ical]->ResidualRange())[iTrkHit]);
        _plane.at(trk).emplace_back(planenum);
      }
    }


    // Other method

    std::vector<const recob::Hit*> hit_v = track_hit_ass.at(trk);
    for (unsigned int hit = 0; hit < hit_v.size(); hit++) {

      if (hit_v[hit]->View() == 2) {
        std::cout << "this is hit " << hit << " with channel " << hit_v[hit]->Channel() << " adc is " << hit_v[hit]->Integral() << " Multiplicity " << hit_v[hit]->Multiplicity() << " LocalIndex " << hit_v[hit]->LocalIndex() << " GoodnessOfFit " << hit_v[hit]->GoodnessOfFit() << " SigmaPeakTime " << hit_v[hit]->SigmaPeakTime() << std::endl;

        if (trk == 0){
        _hit_ch.emplace_back(hit_v[hit]->Channel());
        _hit_int.emplace_back(hit_v[hit]->Integral());
        _hit_mult.emplace_back(hit_v[hit]->Multiplicity());
        _hit_localind.emplace_back(hit_v[hit]->LocalIndex());
        _hit_goodoffit.emplace_back(hit_v[hit]->GoodnessOfFit());
        _hit_sigmapeaktime.emplace_back(hit_v[hit]->SigmaPeakTime());
        _hit_rms.emplace_back(hit_v[hit]->RMS());
        }
        }
    }
  }

  _tree1->Fill();


}

DEFINE_ART_MODULE(DeDxAna)
