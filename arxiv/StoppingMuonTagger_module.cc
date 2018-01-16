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

#include <memory>

// Data products include
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Track.h"
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

// LArSoft include
#include "uboone/UBFlashFinder/PECalib.h"
#include "larsim/MCCheater/BackTracker.h"



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

};


StoppingMuonTagger::StoppingMuonTagger(fhicl::ParameterSet const & p) {


  _fiducial_volume.Configure(p.get<fhicl::ParameterSet>("FiducialVolumeSettings"),
                             geo->DetHalfHeight(),
                             2.*geo->DetHalfWidth(),
                             geo->DetLength());

  _fiducial_volume.PrintConfig();
}

void StoppingMuonTagger::produce(art::Event & e) {


  // Getting PFP from the event
  art::Handle<std::vector<recob::PFParticle> > pfp_h;
  e.getByLabel(_pfp_producer,pfp_h);
  if(!pfp_h.isValid()){
    std::cout << "[StoppingMuonTagger] PFP product " << _pfp_producer << " not found..." << std::endl;
    //throw std::exception();
  }
  if(pfp_h->empty()) {
    std::cout << "[StoppingMuonTagger] PFP " << _pfp_producer << " is empty." << std::endl;
  }
  art::FindManyP<recob::Track> tracks_from_pfp(pfp_h, e, _pfp_producer);

  // Getting MCGhost from the event and their association
  // to PFP and MCP
  art::Handle<std::vector<ubana::MCGhost> > ghost_h;
  e.getByLabel(_mc_ghost_producer,ghost_h);
  if(!ghost_h.isValid()){
    std::cout << "[StoppingMuonTagger] MCGhost product " << _mc_ghost_producer << " not found..." << std::endl;
    //throw std::exception();
  }
  art::FindManyP<ubana::MCGhost>   mcghost_from_pfp   (pfp_h,   e, _mc_ghost_producer);
  art::FindManyP<simb::MCParticle> mcpar_from_mcghost (ghost_h, e, _mc_ghost_producer);


  // PFP loop
  for (unsigned int p = 0; p < pfp_h->size(); p++) {

    art::Ptr<PFParticle> pfp = (*pfp_h)[p];

    // Getting the MCGhost related to this PFP
    std::vector<art::Ptr<ubana::MCGhost>> mcghosts = mcghost_from_pfp.at(pfp.key());
    if (mcghosts.size() == 0 || mcghosts.size() > 1 ) {
      std::cout << "[StoppingMuonTagger] \t\t mcghosts is ether 0 or >1" << std::endl;
      continue;
    }

    // Getting the MCParticle related to this MCGhost
    art::Ptr<simb::MCParticle>mcpar = mcpar_from_mcghost.at(mcghosts[0].key())[0];
    pdg = mcpar->PdgCode();
    std::cout << "[StoppingMuonTagger] \t\t MCPar has pdg " << pdg << std::endl;
    //const auto mc_truth = bt->TrackIDToMCTruth(mcpar->TrackId());
    if (!mc_truth) {
      std::cerr << "[StoppingMuonTagger] Problem with MCTruth pointer." << std::endl;
      continue;
    }

    //if (!mc_truth->Origin() == simb::kCosmicRay) continue;                             // A cosmic ray
    //if (!(mcpar->PdgCode() == 13 || mcpar->PdgCode() == -13)) continue;                // A muon
    //if (!_fiducial_volume.InFV(mcpar->EndX(), mcpar->EndY(), mcpar->EndZ())) continue; // That truly stops in the TPC


  }
}

DEFINE_ART_MODULE(StoppingMuonTagger)
