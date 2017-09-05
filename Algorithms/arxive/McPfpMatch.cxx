#ifndef MCPFPMATCH_CXX
#define MCPFPMATCH_CXX

#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Handle.h"
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
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/BackTracker.h"

#include "McPfpMatch.h"



namespace ubxsec {

  McPfpMatch::McPfpMatch(){

  }

  void McPfpMatch::_Configure_(const Config_t &pset){

    _debug             = pset.get<bool>("Debug", false);
    _verbose           = pset.get<bool>("Verbose", false);
    _recursiveMatching = pset.get<bool>("RecursiveMatching", false);

  }

  void McPfpMatch::Configure(art::Event const & e, 
                             std::string _pfp_producer, 
                             std::string _spacepointLabel, 
                             std::string _hitfinderLabel, 
                             std::string _geantModuleLabel) {

    // Collect hits
    lar_pandora::HitVector hitVector;
    lar_pandora::LArPandoraHelper::CollectHits(e, _hitfinderLabel, hitVector);

    // Collect PFParticles and match Reco Particles to Hits
    lar_pandora::PFParticleVector  recoParticleVector;
    lar_pandora::PFParticleVector  recoNeutrinoVector;
    lar_pandora::PFParticlesToHits recoParticlesToHits;
    lar_pandora::HitsToPFParticles recoHitsToParticles;

    lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, recoParticleVector);
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(recoParticleVector, recoNeutrinoVector);
    lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e, _pfp_producer, _spacepointLabel, recoParticlesToHits, recoHitsToParticles, lar_pandora::LArPandoraHelper::kAddDaughters);

    if (_verbose) {
      std::cout << "[McPfpMatch] RecoNeutrinos: " << recoNeutrinoVector.size() << std::endl;
      std::cout << "[McPfpMatch] RecoParticles: " << recoParticleVector.size() << std::endl;
    }

    // Collect MCParticles and match True Particles to Hits
    lar_pandora::MCParticleVector     trueParticleVector;
    lar_pandora::MCTruthToMCParticles truthToParticles;
    lar_pandora::MCParticlesToMCTruth particlesToTruth;
    lar_pandora::MCParticlesToHits    trueParticlesToHits;
    lar_pandora::HitsToMCParticles    trueHitsToParticles;

    if (!e.isRealData()) {
      lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geantModuleLabel, trueParticleVector);
      lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geantModuleLabel, truthToParticles, particlesToTruth);
      lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(e, _geantModuleLabel, hitVector, trueParticlesToHits, trueHitsToParticles, lar_pandora::LArPandoraHelper::kAddDaughters);
    }

    if (_verbose) {
      std::cout << "[McPfpMatch] TrueParticles: " << particlesToTruth.size() << std::endl;
      std::cout << "[McPfpMatch] TrueEvents: " << truthToParticles.size() << std::endl;
    }  

    // Now set the things we need for the future
    _trueHitsToParticles = trueHitsToParticles;
    _recoParticlesToHits = recoParticlesToHits;

    if (_debug) { // yes, don't do it
      std::cout << "[McPfpMatch] This is event " << e.id().run() << std::endl;
      art::ServiceHandle<cheat::BackTracker> bt;
      std::cout << "[McPfpMatch] Number of MCParticles matched to hits: " << trueParticlesToHits.size() << std::endl;
      for (const auto & iter : trueParticlesToHits) {
        const art::Ptr<simb::MCTruth> mc_truth = bt->TrackIDToMCTruth((iter.first)->TrackId());
        std::cout << "[McPfpMatch] MCParticle with pdg " << (iter.first)->PdgCode()
                  << " and origin " << (mc_truth->Origin() == 1 ? "neutrino" : "cosmic")  
                  << " has " << (iter.second).size() << " hits ass." << std::endl;
        if (mc_truth->Origin() == 1) {
          lar_pandora::HitVector hits = (iter.second);
          for (const auto & hit : hits){
            std::cout << "[McPfpMatch]   > Hit on plane " << hit->View() 
                      << " on wire " << hit->WireID() 
                      << " with time " << hit->PeakTime() << std::endl;
          }     
        }        
      }
    }
  }

  //___________________________________________________________________________________________________
  void McPfpMatch::GetRecoToTrueMatches(lar_pandora::MCParticlesToPFParticles &matchedParticles,
                                        lar_pandora::MCParticlesToHits &matchedParticleHits) {
  
 
    GetRecoToTrueMatches(_recoParticlesToHits,
                         _trueHitsToParticles,
                         matchedParticles,
                         matchedParticleHits);
  }
  
  
  //___________________________________________________________________________________________________
  void McPfpMatch::GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits,
                                        const lar_pandora::HitsToMCParticles &trueHitsToParticles,
                                        lar_pandora::MCParticlesToPFParticles &matchedParticles,
                                        lar_pandora::MCParticlesToHits &matchedHits) {

      PFParticleSet recoVeto; MCParticleSet trueVeto;
  
      GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles, matchedParticles, matchedHits, recoVeto, trueVeto, _recursiveMatching);
  }
  
  //___________________________________________________________________________________________________
  void McPfpMatch::GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits,
                                             const lar_pandora::HitsToMCParticles &trueHitsToParticles,
                                             lar_pandora::MCParticlesToPFParticles &matchedParticles,
                                             lar_pandora::MCParticlesToHits &matchedHits,
                                             PFParticleSet &vetoReco,
                                             MCParticleSet &vetoTrue,
                                             bool _recursiveMatching) 
  {
      bool foundMatches(false);
  
      // Loop over the reco particles
      for (lar_pandora::PFParticlesToHits::const_iterator iter1 = recoParticlesToHits.begin(), iterEnd1 = recoParticlesToHits.end();
          iter1 != iterEnd1; ++iter1)
      {
          const art::Ptr<recob::PFParticle> recoParticle = iter1->first;
          if (vetoReco.count(recoParticle) > 0)
              continue;
  
          const lar_pandora::HitVector &hitVector = iter1->second;
  
          lar_pandora::MCParticlesToHits truthContributionMap;
  
          // Loop over all the hits associated to this reco particle
          for (lar_pandora::HitVector::const_iterator iter2 = hitVector.begin(), iterEnd2 = hitVector.end(); iter2 != iterEnd2; ++iter2)
          {
              const art::Ptr<recob::Hit> hit = *iter2;
  
              lar_pandora::HitsToMCParticles::const_iterator iter3 = trueHitsToParticles.find(hit);
              if (trueHitsToParticles.end() == iter3)
                  continue;
  
              const art::Ptr<simb::MCParticle> trueParticle = iter3->second;
              if (vetoTrue.count(trueParticle) > 0)
                  continue;
  
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
  
              lar_pandora::MCParticlesToHits::const_iterator iter5 = matchedHits.find(trueParticle);
  
              if ((matchedHits.end() == iter5) || (mIter->second.size() > iter5->second.size()))
              {
                  matchedParticles[trueParticle] = recoParticle;
                  matchedHits[trueParticle] = mIter->second;
                  foundMatches = true;
              }
          }
      } // recoParticlesToHits loop ends
  
      if (!foundMatches)
          return;
  
      for (lar_pandora::MCParticlesToPFParticles::const_iterator pIter = matchedParticles.begin(), pIterEnd = matchedParticles.end();
          pIter != pIterEnd; ++pIter)
      {
          vetoTrue.insert(pIter->first);
          vetoReco.insert(pIter->second);
      }
  
      if (_recursiveMatching)
          GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles, matchedParticles, matchedHits, vetoReco, vetoTrue, _recursiveMatching);
  }
  
} // namespace


#endif
