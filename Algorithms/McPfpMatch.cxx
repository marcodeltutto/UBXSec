#ifndef MCPFPMATCH_CXX
#define MCPFPMATCH_CXX

#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

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
//#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "McPfpMatch.h"
#include "UBXSecHelper.h"


namespace ubana {

  McPfpMatch::McPfpMatch(){

    _configured = false;

  }

  void McPfpMatch::_Configure_(const Config_t &pset){

    _debug             = pset.get<bool>("Debug", false);
    _verbose           = pset.get<bool>("Verbose", false);

  }

   
  void McPfpMatch::Configure(art::Event const & e, 
                             std::string _pfp_producer, 
                             std::string _spacepoint_producer, 
                             std::string _hitfinder_producer, 
                             std::string _geant_producer) {

    // Collect hits
    lar_pandora::HitVector hitVector;
    lar_pandora::LArPandoraHelper::CollectHits(e, _hitfinder_producer, hitVector);

    // Collect PFParticles and match Reco Particles to Hits
    lar_pandora::PFParticleVector  recoParticleVector;
    lar_pandora::PFParticleVector  recoNeutrinoVector;
    lar_pandora::PFParticlesToHits pfp_to_hits_map;
    lar_pandora::HitsToPFParticles recoHitsToParticles;

    lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, recoParticleVector);
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(recoParticleVector, recoNeutrinoVector);
    lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e, 
                                                          _pfp_producer, 
                                                          _spacepoint_producer, 
                                                          pfp_to_hits_map, 
                                                          recoHitsToParticles, 
                                                          lar_pandora::LArPandoraHelper::kAddDaughters,
                                                          true); // Use clusters to go from pfp to hits

    if (_verbose) {
      std::cout << "[McPfpMatch] RecoNeutrinos: " << recoNeutrinoVector.size() << std::endl;
      std::cout << "[McPfpMatch] RecoParticles: " << recoParticleVector.size() << std::endl;
    }

    // Collect MCParticles and match True Particles to Hits
    lar_pandora::MCParticleVector     trueParticleVector;
    lar_pandora::MCTruthToMCParticles truthToParticles;
    lar_pandora::MCParticlesToMCTruth particlesToTruth;
    lar_pandora::MCParticlesToHits    trueParticlesToHits;
    lar_pandora::HitsToMCParticles    hit_to_mcps_map;

    if (!e.isRealData()) {
      lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geant_producer, trueParticleVector);
      lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geant_producer, truthToParticles, particlesToTruth);
      lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(e, 
                                                            _geant_producer, 
                                                            hitVector, 
                                                            trueParticlesToHits, 
                                                            hit_to_mcps_map, 
                                                            lar_pandora::LArPandoraHelper::kAddDaughters); // Consider daughters as independent mcps
    }

    if (_verbose) {
      std::cout << "[McPfpMatch] TrueParticles: " << particlesToTruth.size() << std::endl;
      std::cout << "[McPfpMatch] TrueEvents: " << truthToParticles.size() << std::endl;
    }  

    // Now set the things we need for the future
    _hit_to_mcps_map = hit_to_mcps_map;
    _pfp_to_hits_map = pfp_to_hits_map;

    if (_debug) { 
      std::cout << "[McPfpMatch] This is event " << e.id().event() << std::endl;
      art::ServiceHandle<cheat::BackTracker> bt;
      std::cout << "[McPfpMatch] Number of MCParticles matched to hits: " << trueParticlesToHits.size() << std::endl;
      for (const auto & iter : trueParticlesToHits) {
        const art::Ptr<simb::MCTruth> mc_truth = bt->TrackIDToMCTruth((iter.first)->TrackId());
        std::cout << "[McPfpMatch] MCParticle with pdg " << (iter.first)->PdgCode()
                  << " and origin " << (mc_truth->Origin() == 1 ? "neutrino" : "cosmic")  
                  << " has " << (iter.second).size() << " hits ass." << std::endl;
        if (mc_truth->Origin() == 1) {
          lar_pandora::HitVector hits = (iter.second);
        }        
      }
    }

    _configured = true;
  }





  void McPfpMatch::Configure(art::Event const & e,
                             std::string _pfp_producer,
                             std::string _spacepoint_producer,
                             std::string _hitfinder_producer,
                             std::string _geant_producer,
                             std::string _hit_mcp_producer,
                             lar_pandora::LArPandoraHelper::DaughterMode daughterMode) {

    // Collect hits
    lar_pandora::HitVector hitVector;
    //lar_pandora::LArPandoraHelper::CollectHits(e, _hitfinder_producer, hitVector);
    art::Handle<std::vector<recob::Hit>> hit_h;
    e.getByLabel(_hitfinder_producer, hit_h);
    if (!hit_h.isValid()) {
      std::cout << "[McPfpMatch] Hit Handle is not valid." << std::endl;
      throw std::exception();
    }
    art::fill_ptr_vector(hitVector, hit_h);
    art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> mcps_from_hit(hit_h, e, _hit_mcp_producer);

    // Collect PFParticles and match Reco Particles to Hits
    lar_pandora::PFParticleVector  recoParticleVector;
    lar_pandora::PFParticleVector  recoNeutrinoVector;
    lar_pandora::PFParticlesToHits pfp_to_hits_map;
    lar_pandora::HitsToPFParticles recoHitsToParticles;

    lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, recoParticleVector);
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(recoParticleVector, recoNeutrinoVector);
    lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e,
                                                          _pfp_producer,
                                                          _spacepoint_producer,
                                                          pfp_to_hits_map,
                                                          recoHitsToParticles,
                                                          daughterMode,
                                                          true); // Use clusters to go from pfp to hits

    if (_verbose) {
      std::cout << "[McPfpMatch] RecoNeutrinos: " << recoNeutrinoVector.size() << std::endl;
      std::cout << "[McPfpMatch] RecoParticles: " << recoParticleVector.size() << std::endl;
    }

    // Collect MCParticles and match True Particles to Hits
    lar_pandora::MCParticleVector     trueParticleVector;
    lar_pandora::MCTruthToMCParticles truthToParticles;
    lar_pandora::MCParticlesToMCTruth particlesToTruth;
    lar_pandora::MCParticlesToHits    trueParticlesToHits;
    lar_pandora::HitsToMCParticles    hit_to_mcps_map;

    if (!e.isRealData()) {
      lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geant_producer, trueParticleVector);
      lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geant_producer, truthToParticles, particlesToTruth);

      // Construct a Particle Map (trackID to MCParticle)
      lar_pandora::MCParticleMap particleMap;

      for (lar_pandora::MCTruthToMCParticles::const_iterator iter1 = truthToParticles.begin(), iterEnd1 = truthToParticles.end(); iter1 != iterEnd1; ++iter1)
      {
          const lar_pandora::MCParticleVector &particleVector = iter1->second;
          for (lar_pandora::MCParticleVector::const_iterator iter2 = particleVector.begin(), iterEnd2 = particleVector.end(); iter2 != iterEnd2; ++iter2)
          {
              const art::Ptr<simb::MCParticle> particle = *iter2;
              particleMap[particle->TrackId()] = particle;
          }
      }

      // Loop over the hits, get the ass MCP, and then tru to link 

      std::vector<art::Ptr<simb::MCParticle>> mcp_v;
      std::vector<anab::BackTrackerHitMatchingData const*> match_v;

      for (auto hit : hitVector) {

        mcp_v.clear(); match_v.clear();
        mcps_from_hit.get(hit.key(), mcp_v, match_v);

        double max_energy = -1;
        int best_match_id = -1;

        for (size_t m = 0; m < match_v.size(); m++) {
          double this_energy = match_v[m]->energy;
          if (this_energy > max_energy) {
            best_match_id = m;
            max_energy = this_energy;
          }
        }
        
        if (best_match_id > -1) {

          try {

            const art::Ptr<simb::MCParticle> thisParticle = mcp_v.at(best_match_id);
            const art::Ptr<simb::MCParticle> primaryParticle(lar_pandora::LArPandoraHelper::GetFinalStateMCParticle(particleMap, thisParticle));
            const art::Ptr<simb::MCParticle> selectedParticle((lar_pandora::LArPandoraHelper::kAddDaughters == daughterMode) ? primaryParticle : thisParticle);

            if ((lar_pandora::LArPandoraHelper::kIgnoreDaughters == daughterMode) && (selectedParticle != primaryParticle))
              continue;

            if (!(lar_pandora::LArPandoraHelper::IsVisible(selectedParticle)))
              continue;

            hit_to_mcps_map[hit] = selectedParticle;

          } catch (cet::exception &e) {
          }
          //const auto mct = UBXSecHelper::TrackIDToMCTruth(e, "largeant", selectedParticle->TrackId());
          //if (mct->Origin() == simb::kBeamNeutrino && selectedParticle->PdgCode() == 13 && selectedParticle->Mother() == 0) {
            //std::cout << "Muon from neutrino ass to hit " << hit->PeakTime() << ", "<< hit->WireID () << std::endl;
          //}
        }

      }
    }

    // Now set the things we need for the future
    _hit_to_mcps_map = hit_to_mcps_map;
    _pfp_to_hits_map = pfp_to_hits_map;

    if (_debug) std::cout << "hit_to_mcps_map size " << hit_to_mcps_map.size() << std::endl;

    _configured = true;
  }

  //___________________________________________________________________________________________________
  void McPfpMatch::GetRecoToTrueMatches(lar_pandora::PFParticlesToMCParticles &matchedParticles) 
  {

    if (!_configured) {
      std::cout << "Call to " << __PRETTY_FUNCTION__ << " whitout having done configuration. Abort." << std::endl;
      throw std::exception();
    }
      
    // Loop over the reco particles
    for (auto iter1 : _pfp_to_hits_map) {

      // The PFParticle
      const art::Ptr<recob::PFParticle> recoParticle = iter1.first;
 
      if (_debug) std::cout << "[McPfpMatch::GetRecoToTrueMatches] Looking at PFP with ID " << recoParticle->Self() << std::endl;

      // The PFParticle's hits
      const lar_pandora::HitVector &hitVector = iter1.second;

      if (_debug) std::cout << "[McPfpMatch::GetRecoToTrueMatches] \t This PFP has " << hitVector.size() << " hits." << std::endl;
  
      lar_pandora::MCParticlesToHits truthContributionMap;
  
      // Loop over all the hits associated to this reco particle
      for (auto hit : hitVector) {
  
        // Find the MCParticle that share this same hit (if any)
        auto iter3 = _hit_to_mcps_map.find(hit);
        if (_hit_to_mcps_map.end() == iter3)
          continue;
  
        // If exists, get the MCParticle
        const art::Ptr<simb::MCParticle> trueParticle = iter3->second;

        if (_debug) std::cout << "[McPfpMatch::GetRecoToTrueMatches] \t Found a hit shared with MCParticle with PDG " << trueParticle->PdgCode() << std::endl;

        // This map will contain all the true particles that match some or all of the hits of the reco particle
        truthContributionMap[trueParticle].push_back(hit);
      }
  
      // Now we want to find the true particle that has more hits in common with this reco particle than the others
      lar_pandora::MCParticlesToHits::const_iterator mIter = truthContributionMap.end();
  
      for (lar_pandora::MCParticlesToHits::const_iterator iter4 = truthContributionMap.begin(), iterEnd4 = truthContributionMap.end();
           iter4 != iterEnd4; ++iter4) 
      {
        if ((truthContributionMap.end() == mIter) || (iter4->second.size() > mIter->second.size())) 
        {
            mIter = iter4;
        }
      }
 
      if (truthContributionMap.end() != mIter)
      {
        const art::Ptr<simb::MCParticle> trueParticle = mIter->first;
  
        if (_debug) std::cout << "[McPfpMatch::GetRecoToTrueMatches] \t >>> Match found with MCParticle with PDG " << trueParticle->PdgCode() << std::endl;

        // Emplace into the output map
        matchedParticles[recoParticle] = trueParticle;
      }
      
    } // _pfp_to_hits_map loop ends
  
  }
  
} // namespace


#endif
