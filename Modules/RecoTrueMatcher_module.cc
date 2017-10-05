////////////////////////////////////////////////////////////////////////
// Class:       RecoTrueMatcher
// Plugin Type: producer (art v2_05_00)
// File:        RecoTrueMatcher_module.cc
//
// Generated at Fri Aug 18 08:43:34 2017 by Marco Del Tutto using cetskelgen
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
#include "lardata/Utilities/AssociationUtil.h"

// Data product include
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Wire.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/BackTracker.h"
#include "uboone/UBXSec/DataTypes/MCGhost.h"

#include <memory>

// Algorithms include
#include "uboone/UBXSec/Algorithms/McPfpMatch.h"


class RecoTrueMatcher;


class RecoTrueMatcher : public art::EDProducer {
public:
  explicit RecoTrueMatcher(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  RecoTrueMatcher(RecoTrueMatcher const &) = delete;
  RecoTrueMatcher(RecoTrueMatcher &&) = delete;
  RecoTrueMatcher & operator = (RecoTrueMatcher const &) = delete;
  RecoTrueMatcher & operator = (RecoTrueMatcher &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  ubxsec::McPfpMatch mcpfpMatcher;

  std::string _pfp_producer;
  std::string _spacepointLabel;
  std::string _hitfinderLabel;
  std::string _geantModuleLabel;

  bool _is_data;
  bool _debug;
};


RecoTrueMatcher::RecoTrueMatcher(fhicl::ParameterSet const & p) {

  _pfp_producer                   = p.get<std::string>("PFParticleProducer");
  _hitfinderLabel                 = p.get<std::string>("HitProducer");
  _geantModuleLabel               = p.get<std::string>("GeantModule");
  _spacepointLabel                = p.get<std::string>("SpacePointProducer");

  _debug                          = p.get<bool>("DebugMode");

  produces< std::vector<ubana::MCGhost>>();
  produces< art::Assns<simb::MCParticle, ubana::MCGhost>>();
  produces< art::Assns<recob::PFParticle, ubana::MCGhost>>();
}

void RecoTrueMatcher::produce(art::Event & e)
{

  if(_debug) std::cout << "[RecoTrueMatcher] Starts" << std::endl;
  if(_debug) std::cout << "[RecoTrueMatcher] event: " << e.id().event() << std::endl;

  // Instantiate the output
  std::unique_ptr< std::vector< ubana::MCGhost > >                mcGhostVector   (new std::vector<ubana::MCGhost>);
  std::unique_ptr< art::Assns<simb::MCParticle, ubana::MCGhost>>  assnOutGhostMCP (new art::Assns<simb::MCParticle, ubana::MCGhost>);
  std::unique_ptr< art::Assns<recob::PFParticle, ubana::MCGhost>> assnOutGhostPFP (new art::Assns<recob::PFParticle, ubana::MCGhost>);


  _is_data = e.isRealData();

  if (_is_data) {
    std::cout << "[RecoTrueMatcher] Running on a real data file. No MC-PFP matching will be attempted." << std::endl;
    e.put(std::move(mcGhostVector));
    e.put(std::move(assnOutGhostMCP));
    e.put(std::move(assnOutGhostPFP));
    return;
  } 
    
  mcpfpMatcher.Configure(e, _pfp_producer, _spacepointLabel, _hitfinderLabel, _geantModuleLabel);

  lar_pandora::MCParticlesToPFParticles matchedMCToPFParticles;    // This is a map: MCParticle to matched PFParticle
  lar_pandora::MCParticlesToHits        matchedParticleHits;

  mcpfpMatcher.GetRecoToTrueMatches(matchedMCToPFParticles, matchedParticleHits);

  std::cout << "[RecoTrueMatcher] Generating " << matchedMCToPFParticles.size() << " MCGhosts." << std::endl;

  for (auto const& iter : matchedMCToPFParticles) {

    art::Ptr<simb::MCParticle>  mc_par = iter.first;   // The MCParticle 
    art::Ptr<recob::PFParticle> pf_par = iter.second;  // The matched PFParticle 

     std::cout << "[RecoTrueMatching]\t PFP with ID " << pf_par->Self() << ", and PDG " << pf_par->PdgCode() << std::endl;
     std::cout << "[RecoTrueMatching]\t\t ...matched to MCPAR with PDG " << mc_par->PdgCode() << std::endl;

    ubana::MCGhost mcGhost;
    mcGhost.SetMode("depEnergy");

    mcGhostVector->emplace_back(mcGhost);
    util::CreateAssn(*this, e, *mcGhostVector, pf_par, *assnOutGhostPFP);
    util::CreateAssn(*this, e, *mcGhostVector, mc_par, *assnOutGhostMCP);
  }

  e.put(std::move(mcGhostVector));
  e.put(std::move(assnOutGhostMCP));
  e.put(std::move(assnOutGhostPFP));
}

DEFINE_ART_MODULE(RecoTrueMatcher)
