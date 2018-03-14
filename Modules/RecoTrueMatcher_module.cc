////////////////////////////////////////////////////////////////////////
// Class:       RecoTrueMatcher
// Plugin Type: producer (art v2_05_00)
// File:        RecoTrueMatcher_module.cc
//
// Generated at Fri Aug 18 08:43:34 2017 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

/**
 * \class RecoTrueMatcher
 *
 * \ingroup UBXSec
 *
 * \brief Art producer module that performs the PFParticle-MCParticle matching
 * 
 *
 * \author Marco Del Tutto <marco.deltutto@physics.ox.ac.uk>
 *
 * \version producer (art v2_05_00)
 *
 * \date 2017/03/10
 *
 * Contact: marco.deltutto@physics.ox.ac.uk
 *
 * Created on: Fri Aug 18 08:43:34 2017
 *
 */

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

#include "larcore/Geometry/Geometry.h"

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
#include "uboone/UBXSec/Algorithms/FiducialVolume.h"

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

  ubana::McPfpMatch _mcpfpMatcher;
  ::ubana::FiducialVolume _fiducial_volume;

  std::string _pfp_producer;
  std::string _spacepointLabel;
  std::string _hitfinderLabel;
  std::string _geantModuleLabel;
  std::string _mcpHitAssLabel;

  bool _is_data;
  bool _debug;
  bool _verbose;
  bool _use_premade_ass;

  void PrintInfo(lar_pandora::PFParticlesToMCParticles matched_pfp_to_mcp_map);
};


RecoTrueMatcher::RecoTrueMatcher(fhicl::ParameterSet const & p) {

  std::cout << "[RecoTrueMatcher] Initialize" << std::endl;


  ::art::ServiceHandle<geo::Geometry> geo;

  _pfp_producer                   = p.get<std::string>("PFParticleProducer");
  _hitfinderLabel                 = p.get<std::string>("HitProducer");
  _geantModuleLabel               = p.get<std::string>("GeantModule");
  _spacepointLabel                = p.get<std::string>("SpacePointProducer");
  _mcpHitAssLabel                 = p.get<std::string>("MCPHitAssProducer", "pandoraCosmicHitRemoval");

  _use_premade_ass                = p.get<bool>("UsePremadeMCPHitAss");

  _debug                          = p.get<bool>("DebugMode");
  _verbose                        = p.get<bool>("Verbose");

  _fiducial_volume.Configure(p.get<fhicl::ParameterSet>("FiducialVolumeSettings"),
                             geo->DetHalfHeight(),
                             2.*geo->DetHalfWidth(),
                             geo->DetLength());

  produces< std::vector<ubana::MCGhost>>();
  produces< art::Assns<simb::MCParticle, ubana::MCGhost>>();
  produces< art::Assns<recob::PFParticle, ubana::MCGhost>>();

  std::cout << "[RecoTrueMatcher] End Initialize" << std::endl;
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
 
  if (_use_premade_ass)   
    _mcpfpMatcher.Configure(e, _pfp_producer, _spacepointLabel, _hitfinderLabel, _geantModuleLabel, _mcpHitAssLabel, lar_pandora::LArPandoraHelper::kAddDaughters);
  else 
    _mcpfpMatcher.Configure(e, _pfp_producer, _spacepointLabel, _hitfinderLabel, _geantModuleLabel);

  lar_pandora::PFParticlesToMCParticles matched_pfp_to_mcp_map;    // This is a map: PFParticle to matched MCParticle

  _mcpfpMatcher.GetRecoToTrueMatches(matched_pfp_to_mcp_map);

  //if (_verbose)
    //this->PrintInfo(matched_pfp_to_mcp_map);

  if(_debug) std::cout << "[RecoTrueMatcher] Generating " << matched_pfp_to_mcp_map.size() << " MCGhosts." << std::endl;

  for (auto const& iter : matched_pfp_to_mcp_map) {

    art::Ptr<recob::PFParticle> pf_par = iter.first;    // The PFParticle 
    art::Ptr<simb::MCParticle>  mc_par = iter.second;   // The matched MCParticle 

    if(_debug) {
      std::cout << "[RecoTrueMatcher]\t PFP with ID " << pf_par->Self() << ", and PDG " << pf_par->PdgCode() << std::endl;
      std::cout << "[RecoTrueMatcher]\t\t ...matched to MCPAR with PDG " << mc_par->PdgCode()
                << " and Vx, Vy, Vz = " << mc_par->Vx() << ", " << mc_par->Vy() << ", " << mc_par->Vz() << std::endl;
    }

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






void RecoTrueMatcher::PrintInfo(lar_pandora::PFParticlesToMCParticles matched_pfp_to_mcp_map) {

  std::cout << "[RecoTrueMatcher] ~~~~~~~~~~ Printing Additional Info" << std::endl;

  ::art::ServiceHandle<cheat::BackTracker> bt;

  for (auto iter : matched_pfp_to_mcp_map) {

    art::Ptr<recob::PFParticle> pf_par = iter.first;    // The PFParticle 
    art::Ptr<simb::MCParticle>  mc_par = iter.second;   // The matched MCParticle

    const art::Ptr<simb::MCTruth> mc_truth = bt->TrackIDToMCTruth(mc_par->TrackId());

    if (!mc_truth) {
      std::cerr << "[RecoTrueMatcher] Problem with MCTruth pointer." << std::endl;
      continue;
    }

    if (mc_truth->Origin() == simb::kCosmicRay) {

      //cosmicOriginPFP.emplace_back(pf_par);

      double end[3];
      end[0] = mc_par->EndX();
      end[1] = mc_par->EndY();
      end[2] = mc_par->EndZ();
      if ( (mc_par->PdgCode() == 13 || mc_par->PdgCode() == -13) && _fiducial_volume.InFV(end) ){
        std::cout << "[RecoTrueMatcher] Is cosmic stopping muon" << std::endl;

        /*
        lar_pandora::VertexVector          vertexVector;
        lar_pandora::PFParticlesToVertices particlesToVertices;
        lar_pandora::LArPandoraHelper::CollectVertices(e, _pfp_producer, vertexVector, particlesToVertices);

        auto iter = particlesToVertices.find(pf_par);
        if (iter != particlesToVertices.end()) {
          lar_pandora::VertexVector vertex_v = particlesToVertices.find(pf_par)->second;
          double xyz[3];
          vertex_v[0]->XYZ(xyz);
          std::cout << "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM The PFP has vtx x="<<xyz[0]<<" y="<<xyz[1]<<" z="<<xyz[2] << std::endl;
        } else {
          std::cout << "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Can't find ass vertex" << std::endl;
        }*/
      }
    } // cosmic origin



    if (mc_truth->Origin() == simb::kBeamNeutrino) {
       
      std::cout << "[RecoTrueMatcher] Neutrino related track found." << std::endl;
      std::cout << "[RecoTrueMatcher] Process (0==CC, 1==NC) " << mc_truth->GetNeutrino().CCNC()         << std::endl;
      std::cout << "[RecoTrueMatcher] Neutrino PDG           " << mc_truth->GetNeutrino().Nu().PdgCode() << std::endl;
      std::cout << "[RecoTrueMatcher] Particle PDG  " << mc_par->PdgCode() << std::endl;
      std::cout << "[RecoTrueMatcher] Particle Mass " << mc_par->Mass()    << std::endl;
      std::cout << "[RecoTrueMatcher] Particle Proc " << mc_par->Process() << std::endl;
      std::cout << "[RecoTrueMatcher] Particle Vx   " << mc_par->Vx()      << std::endl;
      std::cout << "[RecoTrueMatcher] Particle Vy   " << mc_par->Vy()      << std::endl;
      std::cout << "[RecoTrueMatcher] Particle Vz   " << mc_par->Vz()      << std::endl;
      std::cout << "[RecoTrueMatcher] Particle T    " << mc_par->T()       << std::endl;
      double timeCorrection = 343.75;
      std::cout << "[RecoTrueMatcher] Remeber a time correction of " << timeCorrection << std::endl;
       //std::cout << "Related hits: " << (iter->second).size() << std::endl;
       
       
      std::cout << "[RecoTrueMatcher] The related PFP: " << std::endl;
      std::cout << "[RecoTrueMatcher] has ID: " << pf_par->Self() << std::endl;
       

       //neutrinoOriginPFP.emplace_back(pf_par);
/*
       // If we matched a muon
       if (mc_par->PdgCode() == 13 && mc_par->Mother() == 0) {
         muonMCParticle = mc_par;
         muonPFP = pf_par;
         _muon_is_reco = 1;

         // Muon track puritity and efficiency
         _muon_reco_pur = _muon_reco_eff = -9999;
         auto iter = recoParticlesToHits.find(pf_par);
         if (iter != recoParticlesToHits.end()) {
           UBXSecHelper::GetTrackPurityAndEfficiency((*iter).second, _muon_reco_pur, _muon_reco_eff);
         }
         _true_muon_mom_matched = mc_par->P();
         //std::cout << "YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY efficiency: " << eff << "  purity "  << pur << std::endl;

        
         lar_pandora::PFParticlesToTracks::const_iterator it =  pfParticleToTrackMap.find(pf_par);
         if (it != pfParticleToTrackMap.end()) {

           lar_pandora::TrackVector trk_v = it->second;
           std::cout << "Track vector length is: " << trk_v.size() << std::endl;
           art::Ptr<recob::Track>   trk   = trk_v[0];

           std::cout << "Recon muon start x " << trk->Vertex().X() << std::endl;
           std::cout << "Recon muon start y " << trk->Vertex().Y() << std::endl;
           std::cout << "Recon muon start z " << trk->Vertex().Z() << std::endl;

           _recon_muon_start_x = trk->Vertex().X();
           _recon_muon_start_y = trk->Vertex().Y();
           _recon_muon_start_z = trk->Vertex().Z();

           std::cout << "Recon muon end x " << trk->End().X() << std::endl;
           std::cout << "Recon muon end y " << trk->End().Y() << std::endl;
           std::cout << "Recon muon end z " << trk->End().Z() << std::endl;

           _recon_muon_end_x = trk->End().X();
           _recon_muon_end_y = trk->End().Y();
           _recon_muon_end_z = trk->End().Z();

           std::cout << "MC muon start x " << mc_par->Vx() << std::endl;
           std::cout << "MC muon start y " << mc_par->Vy() << std::endl;
           std::cout << "MC muon start z " << mc_par->Vz() << std::endl;

           _mc_muon_start_x = mc_par->Vx();
           _mc_muon_start_y = mc_par->Vy();
           _mc_muon_start_z = mc_par->Vz();

           std::cout << "MC muon end x " << mc_par->EndX() << std::endl;
           std::cout << "MC muon end y " << mc_par->EndY() << std::endl;
           std::cout << "MC muon end z " << mc_par->EndZ() << std::endl;

           _mc_muon_end_x = mc_par->EndX();
           _mc_muon_end_y = mc_par->EndY();
           _mc_muon_end_z = mc_par->EndZ();
         }

         double start[3] = {_mc_muon_start_x, _mc_muon_start_y, _mc_muon_start_z};
         double stop[3]  = {_mc_muon_end_x,   _mc_muon_end_y,   _mc_muon_end_z};

        if (_fiducial_volume.InFV(start) && _fiducial_volume.InFV(stop))
          _mc_muon_contained = 1;
        
      }*/
    } // neutrino origin
  } // loop over matched particles


  std::cout << "[RecoTrueMatcher] ~~~~~~~~~~ Printing Additional Info Ends" << std::endl;

  return;

}

DEFINE_ART_MODULE(RecoTrueMatcher)
