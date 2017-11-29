////////////////////////////////////////////////////////////////////////
// Class:       StoppingMuonTagger
// Plugin Type: producer (art v2_05_00)
// File:        StoppingMuonTagger_module.cc
//
// Generated at Fri Oct  6 15:55:20 2017 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
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

#include <memory>

// Data products include
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "uboone/UBXSec/DataTypes/FlashMatch.h"
#include "uboone/UBXSec/DataTypes/MCGhost.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"

// LArSoft include
#include "uboone/UBFlashFinder/PECalib.h"
#include "larsim/MCCheater/BackTracker.h"

//Algorithms include
#include "uboone/UBXSec/Algorithms/FiducialVolume.h"

// Root include
#include "TString.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"


class StoppingMuonTagger;


class StoppingMuonTagger : public art::EDProducer {
public:
  explicit StoppingMuonTagger(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  StoppingMuonTagger(StoppingMuonTagger const &) = delete;
  StoppingMuonTagger(StoppingMuonTagger &&) = delete;
  StoppingMuonTagger & operator = (StoppingMuonTagger const &) = delete;
  StoppingMuonTagger & operator = (StoppingMuonTagger &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  //::art::ServiceHandle<cheat::BackTracker> bt;
  ::art::ServiceHandle<geo::Geometry> geo;
  ::ubana::FiducialVolume _fiducial_volume;

  std::string _tpcobject_producer;


  TH1D * _h_nstopmu;
  TH1D * _h_stopmu_type;
};


StoppingMuonTagger::StoppingMuonTagger(fhicl::ParameterSet const & p) {


  _fiducial_volume.Configure(p.get<fhicl::ParameterSet>("FiducialVolumeSettings"),
                             geo->DetHalfHeight(),
                             2.*geo->DetHalfWidth(),
                             geo->DetLength());

  _fiducial_volume.PrintConfig();

  _tpcobject_producer             = p.get<std::string>("TPCObjectProducer", "TPCObjectMaker");


  art::ServiceHandle<art::TFileService> fs;
  _h_nstopmu = fs->make<TH1D>("h_nstopmu", ";Stopping Muons Per Event;", 10, 0, 10);
  _h_stopmu_type = fs->make<TH1D>("h_stopmu_type", "Stopping Muon Chategories;;", 10, 0, 10);
}

void StoppingMuonTagger::produce(art::Event & e) {


  if (e.isRealData())
    return;

  std::cout <<"[StoppingMuonTagger] Starts." << std::endl;

  // Get TPCObjects from the Event
  art::Handle<std::vector<ubana::TPCObject>> tpcobj_h;
  e.getByLabel(_tpcobject_producer, tpcobj_h);
  if (!tpcobj_h.isValid()) {
    std::cout << "[StoppingMuonTagger] Cannote locate ubana::TPCObject." << std::endl;
  }
  //art::FindManyP<ubana::FlashMatch> tpcobjToFlashMatchAssns(tpcobj_h, e, _neutrino_flash_match_producer);
  art::FindManyP<recob::Track>      tpcobjToTrackAssns(tpcobj_h, e, _tpcobject_producer);
  art::FindManyP<recob::PFParticle> tpcobjToPFPAssns(tpcobj_h, e, _tpcobject_producer);
  art::FindManyP<recob::Shower>     tpcobjToShowerAssns(tpcobj_h, e, _tpcobject_producer);
  //art::FindManyP<anab::CosmicTag>   tpcobjToCosmicTagAssns(tpcobj_h, e, _geocosmictag_producer);

  std::vector<art::Ptr<ubana::TPCObject>> tpcobj_v;
  art::fill_ptr_vector(tpcobj_v, tpcobj_h);

  int n_stopmu = 0;

  for (size_t i = 0; i < tpcobj_v.size(); i++) {

    art::Ptr<ubana::TPCObject> tpcobj = tpcobj_v.at(i);

    if (tpcobj->GetOrigin() != ubana::TPCObjectOrigin::kCosmicRay) continue;

    if (tpcobj->GetOriginExtra() != ubana::TPCObjectOriginExtra::kStoppingMuon) continue;

    std::cout << "[StoppingMuonTagger] Found TPCObject representing a stopping cosmic muon." << std::endl;
    n_stopmu++;

    std::vector<art::Ptr<recob::Track>>      tracks  = tpcobjToTrackAssns.at(tpcobj.key());
    std::vector<art::Ptr<recob::PFParticle>> pfps    = tpcobjToPFPAssns.at(tpcobj.key());
    std::vector<art::Ptr<recob::Shower>>     showers = tpcobjToShowerAssns.at(tpcobj.key());

    if (tracks.size() == 1 && showers.size() == 0)          // One track - bin 0
      _h_stopmu_type->Fill(0);
    else if (tracks.size() == 0 && showers.size() == 1)     // One shower - bin 1
      _h_stopmu_type->Fill(1);
    else if (tracks.size() == 1 && showers.size() == 1)     // One track One shower - bin 2
      _h_stopmu_type->Fill(2);
    else if (tracks.size() == 2 && showers.size() == 0)     // Two track Zero shower - bin 3
      _h_stopmu_type->Fill(3);
    else if (tracks.size() == 0 && showers.size() == 2)     // Zero track Two shower - bin 4
      _h_stopmu_type->Fill(4);
    else                                                    // Others
      _h_stopmu_type->Fill(5);


    std::cout <<"[StoppingMuonTagger]\t Number of tracks for this TPCObject:  " << tracks.size()  << std::endl;
    std::cout <<"[StoppingMuonTagger]\t Number of showers for this TPCObject: " << showers.size() << std::endl;
    std::cout <<"[StoppingMuonTagger]\t Number of PFPs for this TPCObject:    " << pfps.size()    << std::endl;

  }

  _h_nstopmu->Fill(n_stopmu);

  std::cout <<"[StoppingMuonTagger] Ends." << std::endl;

}

DEFINE_ART_MODULE(StoppingMuonTagger)
