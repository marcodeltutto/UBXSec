////////////////////////////////////////////////////////////////////////
// Class:       CosmicTaggerAna
// Plugin Type: analyzer (art v2_05_00)
// File:        CosmicTaggerAna_module.cc
//
// Generated at Fri Dec  9 09:44:39 2016 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

/**
 * \class CosmicTaggerAna
 *
 * \ingroup UBXSec
 *
 * \brief Art analyzer module that analyzes cosmic tags created for cosmic removal
 * 
 *
 * \author Marco Del Tutto <marco.deltutto@physics.ox.ac.uk>
 *
 * \version analyzer (art v2_05_00)
 *
 * \date 2017/03/10
 *
 * Contact: marco.deltutto@physics.ox.ac.uk
 *
 * Created on: Fri Dec  9 09:44:39 2016
 *
 */

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

#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "TString.h"
#include "TTree.h"

class CosmicTaggerAna;


class CosmicTaggerAna : public art::EDAnalyzer {
public:
  explicit CosmicTaggerAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CosmicTaggerAna(CosmicTaggerAna const &) = delete;
  CosmicTaggerAna(CosmicTaggerAna &&) = delete;
  CosmicTaggerAna & operator = (CosmicTaggerAna const &) = delete;
  CosmicTaggerAna & operator = (CosmicTaggerAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  std::string _hitfinderLabel;
  std::string _pfp_producer;
  std::string _geantModuleLabel;
  std::string _spacepointLabel;
  std::string _cosmic_flash_tag_producer;
  std::string _cosmic_geo_tag_producer;
  std::string _cosmic_acpt_tag_producer;
  std::string _cosmic_stopmu_tag_producer;
  double _cosmic_flash_tag_score_cut; ///< Score cut used in the analysis to consider the PFP as cosmic (applied to flash tagger)
  double _cosmic_geo_tag_score_cut;   ///< Score cut used in the analysis to consider the PFP as cosmic (applied to geo tagger)
  double _cosmic_acpt_tag_score_cut;  ///< Score cut used in the analysis to consider the PFP as cosmic (applied to acpt tagger)
  double _cosmic_stopmu_tag_score_cut;  ///< Score cut used in the analysis to consider the PFP as cosmic (applied to stopmu tagger)

  bool _recursiveMatching = true;
  bool _debug = false;

  const simb::Origin_t NEUTRINO_ORIGIN = simb::kBeamNeutrino;

  /// Maps used for PFParticle truth matching
  typedef std::map< art::Ptr<recob::PFParticle>, unsigned int > RecoParticleToNMatchedHits;
  typedef std::map< art::Ptr<simb::MCParticle>,  RecoParticleToNMatchedHits > ParticleMatchingMap;
  typedef std::set< art::Ptr<recob::PFParticle> > PFParticleSet;
  typedef std::set< art::Ptr<simb::MCParticle> >  MCParticleSet;

  /**
   *  @brief Perform matching between true and reconstructed particles
   *
   *  @param recoParticlesToHits the mapping from reconstructed particles to hits
   *  @param trueHitsToParticles the mapping from hits to true particles
   *  @param matchedParticles the output matches between reconstructed and true particles
   *  @param matchedHits the output matches between reconstructed particles and hits
   */
  void GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, const lar_pandora::HitsToMCParticles &trueHitsToParticles,
       lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits) const;
   /**
   *  @brief Perform matching between true and reconstructed particles
   *
   *  @param recoParticlesToHits the mapping from reconstructed particles to hits
   *  @param trueHitsToParticles the mapping from hits to true particles
   *  @param matchedParticles the output matches between reconstructed and true particles
   *  @param matchedHits the output matches between reconstructed particles and hits
   *  @param recoVeto the veto list for reconstructed particles
   *  @param trueVeto the veto list for true particles
   */
  void GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, const lar_pandora::HitsToMCParticles &trueHitsToParticles,
               lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits, PFParticleSet &recoVeto, MCParticleSet &trueVeto) const;


  void GetTaggedPFP(art::Event const & e, std::string cosmictag_producer, double score_cut, lar_pandora::PFParticleVector & pfpTaggedOut, std::vector<int> & tagid_v);
  int PFPInCommon(lar_pandora::PFParticleVector taggedPFP, lar_pandora::PFParticleVector pfpGeoTagged);
  bool InFV(double * nu_vertex_xyz);
 
  TTree* _tree1;
  int _run, _subrun, _event;
  int _ccnc, _pdg, _fv;
  double _nu_e;
  int _n_pfp;                ///< Number of PFP
  int _n_pfp_primary;        ///< Number of primary PFP
  int _n_pfp_flash_tagged;   ///< Number of PFP tagged by the Flash Tagger algo
  int _nu_pfp_flash_tagged;  ///< Number of neutrino origin PFP tagged by the Flash Tagger algo 
  int _n_pfp_geo_tagged;     ///< Number of PFP tagged by the Geo Tagger algo
  int _nu_pfp_geo_tagged;    ///< Number of neutrino origin PFP tagged by the Geo Tagger algo
  int _n_pfp_acpt_tagged;    ///< Number of PFP tagged by the ACPT Tagger algo
  int _nu_pfp_acpt_tagged;   ///< Number of neutrino origin PFP tagged by the ACPT Tagger algo
  int _n_pfp_stopmu_tagged;  ///< Number of PFP tagged by the ACPT Tagger algo
  int _nu_pfp_stopmu_tagged; ///< Number of neutrino origin PFP tagged by the ACPT Tagger algo
  int _nu_pfp_tagged_total;  ///< Number of neutrino origin PFP tagged in total
  int _geo_flash_incommon;   ///< Number of tagged PFP in common between Geo and Flash Tagger algo
  int _acpt_flash_incommon;  ///< Number of tagged PFP in common between ACPT and Flash Tagger algo
  int _acpt_geo_incommon;    ///< Number of tagged PFP in common between ACPT and Geo Tagger algo
  bool _nu_pfp_is_reco;      ///< Is true if the (a)muon or (a)nue from the nu interaction is recon.
};


CosmicTaggerAna::CosmicTaggerAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p) 
{
  _pfp_producer               = p.get<std::string>("PFParticleProducer");
  _hitfinderLabel             = p.get<std::string>("HitProducer");
  _geantModuleLabel           = p.get<std::string>("GeantModule");
  _spacepointLabel            = p.get<std::string>("SpacePointProducer");
  _cosmic_flash_tag_producer  = p.get<std::string>("CosmicFlashTagProducer");
  _cosmic_geo_tag_producer    = p.get<std::string>("CosmicGeoTagProducer");
  _cosmic_acpt_tag_producer   = p.get<std::string>("CosmicACPTTagProducer");
  _cosmic_stopmu_tag_producer = p.get<std::string>("CosmicStopMuTagProducer");
  _cosmic_flash_tag_score_cut = p.get<double>("CosmicFlashTagScoreCut");
  _cosmic_geo_tag_score_cut   = p.get<double>("CosmicGeoTagScoreCut");
  _cosmic_acpt_tag_score_cut  = p.get<double>("CosmicACPTTagScoreCut");
  _cosmic_stopmu_tag_score_cut= p.get<double>("CosmicStopMuTagScoreCut");
  _debug                      = p.get<bool>("DebugMode"); 
  
  art::ServiceHandle<art::TFileService> fs;
  _tree1 = fs->make<TTree>("tree","");
  _tree1->Branch("run",                  &_run,                  "run/I");
  _tree1->Branch("subrun",               &_subrun,               "subrun/I");
  _tree1->Branch("event",                &_event,                "event/I");
  _tree1->Branch("ccnc",                 &_ccnc,                 "ccnc/I");
  _tree1->Branch("pdg",                  &_pdg,                  "pdg/I");
  _tree1->Branch("nu_e",                 &_nu_e,                 "nu_e/D");
  _tree1->Branch("fv",                   &_fv,                   "fv/I");
  _tree1->Branch("n_pfp",                &_n_pfp,                "n_pfp/I");
  _tree1->Branch("n_pfp_primary",        &_n_pfp_primary,        "n_pfp_primary/I");
  _tree1->Branch("n_pfp_flash_tagged",   &_n_pfp_flash_tagged,   "n_pfp_flash_tagged/I");
  _tree1->Branch("n_pfp_geo_tagged",     &_n_pfp_geo_tagged,     "n_pfp_geo_tagged/I");
  _tree1->Branch("nu_pfp_flash_tagged",  &_nu_pfp_flash_tagged,  "nu_pfp_flash_tagged/I");
  _tree1->Branch("nu_pfp_geo_tagged",    &_nu_pfp_geo_tagged,    "nu_pfp_geo_tagged/I");
  _tree1->Branch("nu_pfp_tagged_total",  &_nu_pfp_tagged_total,  "nu_pfp_tagged_total/I");
  _tree1->Branch("geo_flash_incommon",   &_geo_flash_incommon,   "geo_flash_incommon/I");
  _tree1->Branch("acpt_flash_incommon",  &_acpt_flash_incommon,  "acpt_flash_incommon/I");
  _tree1->Branch("acpt_geo_incommon",    &_acpt_geo_incommon,    "acpt_geo_incommon/I");
  _tree1->Branch("n_pfp_acpt_tagged",    &_n_pfp_acpt_tagged,    "n_pfp_acpt_tagged/I");
  _tree1->Branch("nu_pfp_acpt_tagged",   &_nu_pfp_acpt_tagged,   "nu_pfp_acpt_tagged/I");
  _tree1->Branch("n_pfp_stopmu_tagged",  &_n_pfp_stopmu_tagged,  "n_pfp_stopmu_tagged/I");
  _tree1->Branch("nu_pfp_stopmu_tagged", &_nu_pfp_stopmu_tagged, "nu_pfp_stopmu_tagged/I");
  _tree1->Branch("nu_pfp_is_reco",       &_nu_pfp_is_reco,       "nu_pfp_is_reco/O");
}

void CosmicTaggerAna::analyze(art::Event const & e)
{

  if (_debug) std::cout << "CosmicTaggerAna starts" << std::endl;
  if (_debug) std::cout << "Event " << e.id().event() << std::endl;

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();

  _ccnc   = -1;
  _pdg    = -1;

  art::ServiceHandle<cheat::BackTracker> bt;

  // *******************
  // Pandora MCParticle to PFParticle matching
  // *******************

  // --- Collect hits
  lar_pandora::HitVector hitVector;
  lar_pandora::LArPandoraHelper::CollectHits(e, _hitfinderLabel, hitVector);

  // --- Collect PFParticles and match Reco Particles to Hits
  lar_pandora::PFParticleVector  recoParticleVector;
  lar_pandora::PFParticleVector  recoNeutrinoVector;
  lar_pandora::PFParticlesToHits recoParticlesToHits;
  lar_pandora::HitsToPFParticles recoHitsToParticles;

  lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, recoParticleVector);
  lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(recoParticleVector, recoNeutrinoVector);
  lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e, _pfp_producer, _spacepointLabel, recoParticlesToHits, recoHitsToParticles, lar_pandora::LArPandoraHelper::kAddDaughters);

  if (_debug)
    std::cout << "  RecoNeutrinos: " << recoNeutrinoVector.size() << std::endl;

  if (_debug)
    std::cout << "  RecoParticles: " << recoParticleVector.size() << std::endl;

  _n_pfp = recoParticleVector.size();

  _n_pfp_primary = 0;
  for (int i = 0; i < _n_pfp; i++) {
    if (!recoParticleVector.at(i)->IsPrimary()) continue;
    _n_pfp_primary++;
  }

  // --- Collect MCParticles and match True Particles to Hits
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

  if (_debug)
    std::cout << "  TrueParticles: " << particlesToTruth.size() << std::endl;

  if (_debug)
    std::cout << "  TrueEvents: " << truthToParticles.size() << std::endl;

  lar_pandora::MCParticlesToPFParticles matchedParticles;    // This is a map: MCParticle to matched PFParticle
  lar_pandora::MCParticlesToHits        matchedParticleHits;

  // --- Do the matching
  this->GetRecoToTrueMatches(recoParticlesToHits, 
                             trueHitsToParticles, 
                             matchedParticles, 
                             matchedParticleHits);


  // *******************
  // Save PFP with neutrino origin
  // *******************

  std::vector<art::Ptr<recob::PFParticle>> taggedPFP;
  std::vector<art::Ptr<recob::PFParticle>> neutrinoOriginPFP;

  _nu_pfp_is_reco = false;

  // Loop over true particle and find the neutrino related ones
  for (lar_pandora::MCParticlesToPFParticles::const_iterator iter1 = matchedParticles.begin(), iterEnd1 = matchedParticles.end();
             iter1 != iterEnd1; ++iter1) {

     art::Ptr<simb::MCParticle>  mc_par = iter1->first;   // The MCParticle 
     art::Ptr<recob::PFParticle> pf_par = iter1->second;  // The matched PFParticle

     const art::Ptr<simb::MCTruth> mc_truth = bt->TrackIDToMCTruth(mc_par->TrackId());
     if (!mc_truth) continue;  
     if (mc_truth->Origin() == NEUTRINO_ORIGIN) {
       if (_debug) {
         std::cout << "Neutrino related track found." << std::endl;
         std::cout << "Process (0==CC, 1==NC) " << mc_truth->GetNeutrino().CCNC()         << std::endl;
         std::cout << "Neutrino PDG           " << mc_truth->GetNeutrino().Nu().PdgCode() << std::endl;
         std::cout << "PDG  " << mc_par->PdgCode() << std::endl;
         std::cout << "Mass " << mc_par->Mass()    << std::endl;
         std::cout << "Proc " << mc_par->Process() << std::endl;
         std::cout << "Vx   " << mc_par->Vx()      << std::endl;
         std::cout << "Vy   " << mc_par->Vy()      << std::endl;
         std::cout << "Vz   " << mc_par->Vz()      << std::endl;
         std::cout << "T    " << mc_par->T()       << std::endl;
         double timeCorrection = 343.75;
         std::cout << "Remeber a time correction of " << timeCorrection << std::endl;
       }    
       if (_debug) {
         std::cout << "  The related PFP: " << std::endl;
         std::cout << "  has ID: " << pf_par->Self() << std::endl;
         //std::cout << "  Vx " << mc_par->Vx() << std::endl;
         //std::cout << "  Vy " << mc_par->Vy() << std::endl;
         //std::cout << "  Vz " << mc_par->Vz() << std::endl;
       }

       // esclude neutrons, as they may travel away and be tagged as out-of-time, 
       // but we don't really care about those
       //if (mc_par->PdgCode() != 2112) neutrinoOriginPFP.emplace_back(pf_par);
       // actually only consider muons and electrons
       if (mc_par->PdgCode() == 13 || mc_par->PdgCode() == -13 || 
           mc_par->PdgCode() == 12 || mc_par->PdgCode() == -12) {
         neutrinoOriginPFP.emplace_back(pf_par);
         _nu_pfp_is_reco = true;
       }

       _ccnc = mc_truth->GetNeutrino().CCNC();
       _pdg  = mc_truth->GetNeutrino().Nu().PdgCode();
       _nu_e = mc_truth->GetNeutrino().Nu().E();
       double pos[3] = {mc_truth->GetNeutrino().Nu().Vx(), mc_truth->GetNeutrino().Nu().Vy(), mc_truth->GetNeutrino().Nu().Vz()};
       _fv   = (this->InFV(pos) ? 1 : 0);
     }
  }
  if (_debug) std::cout << "Neutrino related PFPs in this event: " << neutrinoOriginPFP.size() << std::endl;


  _nu_pfp_tagged_total = 0;


  // *******************
  // GEO
  // *******************

  lar_pandora::PFParticleVector pfpFlashTagged;
  std::vector<int> tagid_v;
  this->GetTaggedPFP(e, _cosmic_flash_tag_producer, _cosmic_flash_tag_score_cut, pfpFlashTagged,  tagid_v);
  if (_debug) std::cout << pfpFlashTagged.size() << " PFP have been flash tagged." << std::endl;
  _n_pfp_flash_tagged = pfpFlashTagged.size();

   // geo Loop through the taggedPFP and see if there is a neutrino related one
  _nu_pfp_flash_tagged = 0;
  for (unsigned int i = 0; i < pfpFlashTagged.size(); i++) {
    for (unsigned int j = 0; j < neutrinoOriginPFP.size(); j++) {
      if(pfpFlashTagged[i] == neutrinoOriginPFP[j]) {
        if (_debug) std::cout << ">>>>>>>>>>>>>>>>> A neutrino related PFP (with ID " << neutrinoOriginPFP[j]->Self() << ") was tagged by the flash tagger with tag " << tagid_v[i] << std::endl;
        _nu_pfp_flash_tagged = 1;
        _nu_pfp_tagged_total = 1;
      }
    }
  }



  // ****
  // Geo/Time
  // ****

  lar_pandora::PFParticleVector pfpGeoTagged;
  this->GetTaggedPFP(e, _cosmic_geo_tag_producer, 0.6, pfpGeoTagged,  tagid_v);
  if (_debug) std::cout << pfpGeoTagged.size() << " PFP have been geo tagged." << std::endl;
  _n_pfp_geo_tagged = pfpGeoTagged.size();

   // Loop through the taggedPFP and see if there is a neutrino related one
  _nu_pfp_geo_tagged = 0;
  for (unsigned int i = 0; i < pfpGeoTagged.size(); i++) {
    for (unsigned int j = 0; j < neutrinoOriginPFP.size(); j++) {
      if(pfpGeoTagged[i] == neutrinoOriginPFP[j]) {
        if (_debug) std::cout << ">>>>>>>>>>>>>>>>> A neutrino related PFP (with ID " << neutrinoOriginPFP[j]->Self() << ") was tagged by the geo tagger with tag " << tagid_v[i] << std::endl;
        _nu_pfp_geo_tagged = 1;
        _nu_pfp_tagged_total = 1;
      }
    }
  }

  _geo_flash_incommon = this->PFPInCommon(pfpFlashTagged, pfpGeoTagged);
  if (_debug) std::cout << _geo_flash_incommon << " PFP are in common between geo and flash tagger" << std::endl;

  // ****
  // ACPT
  // ****

  lar_pandora::PFParticleVector pfpACPTTagged;
  this->GetTaggedPFP(e, _cosmic_acpt_tag_producer, 0.05, pfpACPTTagged, tagid_v);
  if (_debug) std::cout << pfpACPTTagged.size() << " PFP have been ACPT tagged." << std::endl;
  _n_pfp_acpt_tagged = pfpACPTTagged.size();

   // Loop through the taggedPFP and see if there is a neutrino related one
  _nu_pfp_acpt_tagged = 0;
  for (unsigned int i = 0; i < pfpACPTTagged.size(); i++) {
    for (unsigned int j = 0; j < neutrinoOriginPFP.size(); j++) {
      if(pfpACPTTagged[i] == neutrinoOriginPFP[j]) {
        if (_debug) std::cout << ">>>>>>>>>>>>>>>>> ACPT A neutrino related PFP (with ID " << neutrinoOriginPFP[j]->Self() << ") was tagged by the acpt tagger with tag " << tagid_v[i] << std::endl;
        _nu_pfp_acpt_tagged = 1;
        _nu_pfp_tagged_total = 1;
      }
    }
  }

  _acpt_flash_incommon = this->PFPInCommon(pfpFlashTagged, pfpACPTTagged);
  if (_debug) std::cout << _acpt_flash_incommon << " PFP are in common between acpt and flash tagger" << std::endl;

  _acpt_geo_incommon = this->PFPInCommon(pfpGeoTagged, pfpACPTTagged);
  if (_debug) std::cout << _acpt_geo_incommon << " PFP are in common between acpt and geo tagger" << std::endl;



  // ****
  // StopMu
  // ****

  lar_pandora::PFParticleVector pfpStopMuTagged;
  this->GetTaggedPFP(e, _cosmic_stopmu_tag_producer, 0.05, pfpStopMuTagged, tagid_v);
  if (_debug) std::cout << pfpStopMuTagged.size() << " PFP have been StopMu tagged." << std::endl;
  _n_pfp_stopmu_tagged = pfpStopMuTagged.size();

   // Loop through the taggedPFP and see if there is a neutrino related one
  _nu_pfp_stopmu_tagged = 0;
  for (unsigned int i = 0; i < pfpStopMuTagged.size(); i++) {
    for (unsigned int j = 0; j < neutrinoOriginPFP.size(); j++) {
      if(pfpStopMuTagged[i] == neutrinoOriginPFP[j]) {
        if (_debug) std::cout << ">>>>>>>>>>>>>>>>> STOPMU A neutrino related PFP (with ID " << neutrinoOriginPFP[j]->Self() << ") was tagged by the stopmu tagger with tag " << tagid_v[i] << std::endl;
        _nu_pfp_stopmu_tagged = 1;
        _nu_pfp_tagged_total = 1;
      }
    }
  }

  _tree1->Fill();
}



//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicTaggerAna::GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, 
                                                const lar_pandora::HitsToMCParticles &trueHitsToParticles,
                                                lar_pandora::MCParticlesToPFParticles &matchedParticles, 
                                                lar_pandora::MCParticlesToHits &matchedHits) const
{   
  PFParticleSet recoVeto; MCParticleSet trueVeto;
    
  this->GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles, matchedParticles, matchedHits, recoVeto, trueVeto);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicTaggerAna::GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, 
                                                const lar_pandora::HitsToMCParticles &trueHitsToParticles,
                                                lar_pandora::MCParticlesToPFParticles &matchedParticles, 
                                                lar_pandora::MCParticlesToHits &matchedHits, 
                                                PFParticleSet &vetoReco, 
                                                MCParticleSet &vetoTrue) const
{
    bool foundMatches(false);

    for (lar_pandora::PFParticlesToHits::const_iterator iter1 = recoParticlesToHits.begin(), iterEnd1 = recoParticlesToHits.end();
        iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> recoParticle = iter1->first;
        if (vetoReco.count(recoParticle) > 0)
            continue;

        const lar_pandora::HitVector &hitVector = iter1->second;

        lar_pandora::MCParticlesToHits truthContributionMap;

        for (lar_pandora::HitVector::const_iterator iter2 = hitVector.begin(), iterEnd2 = hitVector.end(); iter2 != iterEnd2; ++iter2)
        {
            const art::Ptr<recob::Hit> hit = *iter2;

            lar_pandora::HitsToMCParticles::const_iterator iter3 = trueHitsToParticles.find(hit);
            if (trueHitsToParticles.end() == iter3)
                continue;

            const art::Ptr<simb::MCParticle> trueParticle = iter3->second;
            if (vetoTrue.count(trueParticle) > 0)
                continue;

            truthContributionMap[trueParticle].push_back(hit);
        }

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
    }

    if (!foundMatches)
        return;

    for (lar_pandora::MCParticlesToPFParticles::const_iterator pIter = matchedParticles.begin(), pIterEnd = matchedParticles.end();
        pIter != pIterEnd; ++pIter)
    {
        vetoTrue.insert(pIter->first);
        vetoReco.insert(pIter->second);
    }

    if (_recursiveMatching)
        this->GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles, matchedParticles, matchedHits, vetoReco, vetoTrue);
}






//____________________________________________________________________________________________
void CosmicTaggerAna::GetTaggedPFP(art::Event const & e, std::string cosmictag_producer, double score_cut, lar_pandora::PFParticleVector & pfpTaggedOut,std::vector<int> & tagid_v){

  pfpTaggedOut.clear();
  tagid_v.clear();

  if (_debug) std::cout << "Getting cosmic tags from " << cosmictag_producer << std::endl;

  // Get the CosmicTag from the ART event
  art::Handle<std::vector<anab::CosmicTag>> cosmicTagHandle;
  e.getByLabel(cosmictag_producer, cosmicTagHandle);

  if (!cosmicTagHandle.isValid() || cosmicTagHandle->empty()){
    std::cerr << "Cosmic tag " << cosmictag_producer << " is not valid or empty." << std::endl;
    return;
  }

  // Look up the associations to PFPs
  art::FindManyP<recob::PFParticle> cosmicPFPAssns(cosmicTagHandle, e, cosmictag_producer);

  //if (_debug) std::cout << " cosmicPFPAssns.size(): " << cosmicPFPAssns.size() << std::endl;

  // Loop over the cosmic tags
  for (unsigned int ct = 0; ct < cosmicPFPAssns.size(); ct++) {

    // Get the cosmic tag
    art::Ptr<anab::CosmicTag> cosmicTag(cosmicTagHandle, ct);
    //if(_debug) std::cout << "This cosmic tag (" << ct << ") has type: " << cosmicTag->CosmicType() << " and score: " << cosmicTag->CosmicScore() << std::endl;

    // Get the PFP associated with this CT
    std::vector<art::Ptr<recob::PFParticle>> cosmicTagToPFP_v = cosmicPFPAssns.at(cosmicTag.key());
    //if(_debug) std::cout << "Number of PFP associated with this Cosmic Tag: " << cosmicTagToPFP_v.size() << std::endl;

    if (score_cut < 0) {
      pfpTaggedOut.emplace_back(cosmicTagToPFP_v.at(0));
      tagid_v.emplace_back(cosmicTag->CosmicType());
    } else {
      if (cosmicTag->CosmicScore() > score_cut) {
        pfpTaggedOut.emplace_back(cosmicTagToPFP_v.at(0));
        tagid_v.emplace_back(cosmicTag->CosmicType());
      }
    }
  }

}
//____________________________________________________________________________________________
int CosmicTaggerAna::PFPInCommon(lar_pandora::PFParticleVector taggedPFP, lar_pandora::PFParticleVector pfpGeoTagged){

  int nInCommon = 0;

  for (unsigned int flash = 0; flash < taggedPFP.size(); flash++){

    for (unsigned int geo = 0; geo < pfpGeoTagged.size(); geo++){
 
      if(taggedPFP[flash] == pfpGeoTagged[geo]) {

        nInCommon++;

        //std::cout << "In common found, flash is  " << flash << " and geo is " << geo << std::endl; 

      }
    }

  }
  return nInCommon;

}


//______________________________________________________________________________________________________________________________________
bool CosmicTaggerAna::InFV(double * nu_vertex_xyz){

  double x = nu_vertex_xyz[0];
  double y = nu_vertex_xyz[1];
  double z = nu_vertex_xyz[2];

  //This defines our current settings for the fiducial volume
  double FVx = 256.35;
  double FVy = 233;
  double FVz = 1036.8;
  double borderx = 10.;
  double bordery = 20.;
  double borderz = 10.;
  //double cryoradius = 191.61;
  //double cryoz = 1086.49 + 2*67.63;

  if(x < (FVx - borderx) && (x > borderx) && (y < (FVy/2. - bordery)) && (y > (-FVy/2. + bordery)) && (z < (FVz - borderz)) && (z > borderz)) return true;
  return false;

}
DEFINE_ART_MODULE(CosmicTaggerAna)
