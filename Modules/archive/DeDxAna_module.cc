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

#include "uboone/UBXSec/Algorithms/UBXSecHelper.h"
#include "uboone/UBXSec/Algorithms/FindDeadRegions.h"

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

  FindDeadRegions deadRegionsFinder;

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

  /* Wire check
  double point_in_tpc[3];
  for (double zpos = 0; zpos < 1000; zpos += 0.3){

    point_in_tpc[0] = 100.;
    point_in_tpc[1] = 0.;
    point_in_tpc[2] = zpos;
    raw::ChannelID_t ch = geo->NearestChannel(point_in_tpc, 2);
    const lariov::ChannelStatusProvider& chanFilt = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    std::cout << ">>> this is ch " << ch << "  with status status " << chanFilt.Status(ch) << std::endl;
  }*/



  // Track loop
  for (unsigned int trk = 0; trk < track_h->size(); trk++) {

    std::cout << "***** Track " << trk << std::endl;
    //_trk_number(trk) = (int) trk;

    auto const& track = (*track_h)[trk];

    double track_length = track.Length();
    if (track_length < 50) {
      std::cout << "Track length is less than a 50cm. Continue." << std::endl;       
      continue;
    }

    // Get only tracks that cross the boundary
    int vtx_ok;
    if(!UBXSecHelper::IsCrossingBoundary(track, vtx_ok)) {
      std::cout << "Track is not crossing boundaries. Continue." << std::endl;
      continue;
    }
    std::cout << "vtx_ok " << vtx_ok << std::endl;

    // Understand dead region
    bool isCloseToDeadRegion = false;
    double point_in_tpc[3];
    if (vtx_ok == 0) {
      point_in_tpc[0] = track.Vertex().X();
      point_in_tpc[1] = track.Vertex().Y();
      point_in_tpc[2] = track.Vertex().Z();
    } else {
      point_in_tpc[0] = track.End().X();
      point_in_tpc[1] = track.End().Y();
      point_in_tpc[2] = track.End().Z();
    }
    isCloseToDeadRegion = deadRegionsFinder.NearDeadReg2P(point_in_tpc[1], point_in_tpc[2], 0.6);

    if (isCloseToDeadRegion) {
      std::cout << "Point is close to a dead region. Continue." << std::endl;
      continue;
    }


    // Look at dedx along the track

    std::vector<const anab::Calorimetry*> calos = calo_track_ass.at(trk);

    for (size_t ical = 0; ical<calos.size(); ++ical){
      if (!calos[ical]) continue;
      if (!calos[ical]->PlaneID().isValid) continue;
      int planenum = calos[ical]->PlaneID().Plane;
      if (planenum<0||planenum>2) continue;
      if (planenum != 2) continue;

      // Understand if the calo module flipped the track
      double dqdx_start = (calos[ical]->dQdx())[0] + (calos[ical]->dQdx())[1] + (calos[ical]->dQdx())[2];
      double dqdx_end   = (calos[ical]->dQdx())[calos[ical]->dQdx().size()-1] + (calos[ical]->dQdx())[calos[ical]->dQdx().size()-2] + (calos[ical]->dQdx())[calos[ical]->dQdx().size()-3];
      bool caloFlippedTrack = dqdx_start < dqdx_end;

      // We want two vectors: one with the residual range and the other with the dedx
      std::vector<double> res_range, dedx;
      res_range.resize(calos[ical]->dEdx().size(), -9999);
      dedx.resize(calos[ical]->dEdx().size(), -9999);

      for(size_t iTrkHit = 0; iTrkHit < calos[ical]->dEdx().size(); ++iTrkHit) {
        res_range[iTrkHit] = (calos[ical]->ResidualRange())[iTrkHit];
        dedx[iTrkHit]      = (calos[ical]->dEdx())[iTrkHit];
      }

      // Now apply the smoothing
      std::cout << "Applying smoothing now." << std::endl;
      int smoothing_steps = std::round(0.0733*track_length); //40;
      std::cout << "Number of smoothing_steps: " << smoothing_steps << std::endl;
      std::vector<double> res_range_smooth, dedx_smooth;
      std::vector<double> dedx_temp;
      std::vector<double> res_range_temp;
      for (unsigned int irange = 0; irange < res_range.size(); irange++){
        dedx_temp.emplace_back(dedx[irange]);
        res_range_temp.emplace_back(res_range[irange]);

        if (irange % smoothing_steps == 0) {
         // Median
         double median;
         size_t size = dedx_temp.size();
         std::sort(dedx_temp.begin(), dedx_temp.end());
         if (size % 2 == 0){
          median = (dedx_temp[size/2 - 1] + dedx_temp[size/2]) / 2;
         }
         else{
           median = dedx_temp[size/2];
         }
         dedx_smooth.emplace_back(median);
      
         double sum = std::accumulate(res_range_temp.begin(), res_range_temp.end(), 0.0);
         double mean = sum / res_range_temp.size();
         res_range_smooth.emplace_back(mean);
      
         dedx_temp.clear();
         res_range_temp.clear(); 

        }
      }
      // smoothing ended

      int nBins = std::round(0.0133 * track_length);
      double _dedx_threshold_bragg = 1.9;
      std::cout << "nBins used is: " << nBins << std::endl;
      int nhigher = 0;
      if ( (vtx_ok == 0 && !caloFlippedTrack) || (vtx_ok == 1 && caloFlippedTrack) ) {
        // The track start is in the FV, look around this point to see if there is a Bragg peak
        // This means looking at the region where the residual range is bigger
        for (size_t range_bin = 0; range_bin < 4; range_bin++) {
          if (dedx_smooth[range_bin] > _dedx_threshold_bragg) nhigher ++;
        }
      } else {
        // The end start is in the FV, look around this point to see if there is a Bragg peak
        // This means looking at the region where the residual range is smaller
        for (size_t range_bin = res_range_smooth.size()-1; range_bin > res_range_smooth.size()-1-4; range_bin--) {
          if (dedx_smooth[range_bin] > _dedx_threshold_bragg) nhigher ++;
        }
      }

      bool isStoppingMuon = false;
      if (nhigher >= 3) isStoppingMuon = true; 
      
      if (isStoppingMuon) std::cout << "This is a cosmic stopping muon." << std::endl;
      else std::cout << "This is NOT a cosmic stopping muon." << std::endl;

      //std::cout << "plane " << planenum << std::endl;
      //std::cout << " ke " << calos[ical]->KineticEnergy() << std::endl;
      //std::cout << " range " << calos[ical]->Range() << std::endl;

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
        //std::cout << "this is hit " << hit << " with channel " << hit_v[hit]->Channel() << " adc is " << hit_v[hit]->Integral() << " Multiplicity " << hit_v[hit]->Multiplicity() << " LocalIndex " << hit_v[hit]->LocalIndex() << " GoodnessOfFit " << hit_v[hit]->GoodnessOfFit() << " SigmaPeakTime " << hit_v[hit]->SigmaPeakTime() << std::endl;

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
