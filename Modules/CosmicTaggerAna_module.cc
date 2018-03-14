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
//#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "uboone/UBXSec/DataTypes/MCGhost.h"
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "uboone/UBXSec/Algorithms/UBXSecHelper.h"
#include "uboone/UBXSec/Algorithms/FiducialVolume.h"

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

  ::ubana::FiducialVolume _fiducial_volume;

  std::string _hitfinderLabel;
  std::string _pfp_producer;
  std::string _geantModuleLabel;
  std::string _spacepointLabel;
  std::string _cosmic_flash_tag_producer;
  std::string _cosmic_geo_tag_producer;
  std::string _cosmic_acpt_tag_producer;
  std::string _cosmic_stopmu_tag_producer;
  std::string _mc_ghost_producer;
  double _cosmic_flash_tag_score_cut; ///< Score cut used in the analysis to consider the PFP as cosmic (applied to flash tagger)
  double _cosmic_geo_tag_score_cut;   ///< Score cut used in the analysis to consider the PFP as cosmic (applied to geo tagger)
  double _cosmic_acpt_tag_score_cut;  ///< Score cut used in the analysis to consider the PFP as cosmic (applied to acpt tagger)
  double _cosmic_stopmu_tag_score_cut;  ///< Score cut used in the analysis to consider the PFP as cosmic (applied to stopmu tagger)

  bool _recursiveMatching = true;
  bool _debug = false;

  const simb::Origin_t NEUTRINO_ORIGIN = simb::kBeamNeutrino;

  void GetTaggedPFP(art::Event const & e, std::string cosmictag_producer, double score_cut, lar_pandora::PFParticleVector & pfpTaggedOut, std::vector<int> & tagid_v);
  int PFPInCommon(lar_pandora::PFParticleVector taggedPFP, lar_pandora::PFParticleVector pfpGeoTagged);
 
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

  ::art::ServiceHandle<geo::Geometry> geo;

  _pfp_producer               = p.get<std::string>("PFParticleProducer");
  _hitfinderLabel             = p.get<std::string>("HitProducer");
  _geantModuleLabel           = p.get<std::string>("GeantModule");
  _spacepointLabel            = p.get<std::string>("SpacePointProducer");
  _cosmic_flash_tag_producer  = p.get<std::string>("CosmicFlashTagProducer");
  _cosmic_geo_tag_producer    = p.get<std::string>("CosmicGeoTagProducer");
  _cosmic_acpt_tag_producer   = p.get<std::string>("CosmicACPTTagProducer");
  _cosmic_stopmu_tag_producer = p.get<std::string>("CosmicStopMuTagProducer");
  _mc_ghost_producer           = p.get<std::string>("MCGhostProducer");
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

  _fiducial_volume.Configure(p.get<fhicl::ParameterSet>("FiducialVolumeSettings"), 
                             geo->DetHalfHeight(),
                             2.*geo->DetHalfWidth(),
                             geo->DetLength());

  if (_debug) _fiducial_volume.PrintConfig();
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

  //art::ServiceHandle<cheat::BackTracker> bt;


/*
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
*/

  // Get PFP
  art::Handle<std::vector<recob::PFParticle> > pfp_h;
  e.getByLabel(_pfp_producer,pfp_h);
  if(!pfp_h.isValid()){
    std::cout << "[UBXSec] PFP product " << _pfp_producer << " not found..." << std::endl;
    //throw std::exception();
  }
  if(pfp_h->empty()) {
    std::cout << "[UBXSec] PFP " << _pfp_producer << " is empty." << std::endl;
  }
  std::vector<art::Ptr<recob::PFParticle>> pfp_v;
  art::fill_ptr_vector(pfp_v, pfp_h);


  // Get Ghosts
  art::Handle<std::vector<ubana::MCGhost> > ghost_h;
  e.getByLabel(_mc_ghost_producer,ghost_h);
  if(!ghost_h.isValid()){
    std::cout << "[UBXSec] MCGhost product " << _mc_ghost_producer << " not found..." << std::endl;
    //throw std::exception();
  }
  art::FindManyP<ubana::MCGhost>   mcghost_from_pfp   (pfp_h,   e, _mc_ghost_producer);
  art::FindManyP<simb::MCParticle> mcpar_from_mcghost (ghost_h, e, _mc_ghost_producer); 


  // *******************
  // Save PFP with neutrino origin
  // *******************

  //std::vector<art::Ptr<recob::PFParticle>> taggedPFP;
  //std::vector<art::Ptr<recob::PFParticle>> neutrinoOriginPFP;

  int n_neutrinoOriginPFP = 0;

  _nu_pfp_is_reco = false;

  // Loop over true particle and find the neutrino related ones
  for (auto p : pfp_v) {

    art::Ptr<simb::MCParticle>  mc_par;      // The MCParticle 
    art::Ptr<recob::PFParticle> pf_par = p;  // The matched PFParticle

    auto mcghosts = mcghost_from_pfp.at(p.key());
    if (mcghosts.size() == 0) 
      continue;
    mc_par = mcpar_from_mcghost.at(mcghosts.at(0).key()).at(0);

    const art::Ptr<simb::MCTruth> mc_truth = UBXSecHelper::TrackIDToMCTruth(e, "largeant", mc_par->TrackId());
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
        n_neutrinoOriginPFP++;
        _nu_pfp_is_reco = true;
      }

      _ccnc = mc_truth->GetNeutrino().CCNC();
      _pdg  = mc_truth->GetNeutrino().Nu().PdgCode();
      _nu_e = mc_truth->GetNeutrino().Nu().E();
      double pos[3] = {mc_truth->GetNeutrino().Nu().Vx(), mc_truth->GetNeutrino().Nu().Vy(), mc_truth->GetNeutrino().Nu().Vz()};
      _fv   = _fiducial_volume.InFV(pos);
    }
  }
  if (_debug) std::cout << "Neutrino related PFPs in this event: " << n_neutrinoOriginPFP << std::endl;


  _nu_pfp_tagged_total = 0;


  // *******************
  // Flash
  // *******************

  lar_pandora::PFParticleVector pfpFlashTagged;
  std::vector<int> tagid_v;
  this->GetTaggedPFP(e, _cosmic_flash_tag_producer, _cosmic_flash_tag_score_cut, pfpFlashTagged,  tagid_v);
  if (_debug) std::cout << pfpFlashTagged.size() << " PFP have been flash tagged." << std::endl;
  _n_pfp_flash_tagged = pfpFlashTagged.size();

  // geo Loop through the taggedPFP and see if there is a neutrino related one
  _nu_pfp_flash_tagged = 0;
  for (unsigned int i = 0; i < pfpFlashTagged.size(); i++) {
    auto mcghosts = mcghost_from_pfp.at(pfpFlashTagged.at(i).key());
    if (mcghosts.size() > 0) {
      art::Ptr<simb::MCParticle> mcpar = mcpar_from_mcghost.at(mcghosts.at(0).key()).at(0);
      const auto mc_truth = UBXSecHelper::TrackIDToMCTruth(e, "largeant", mcpar->TrackId());
      if (mc_truth) {
        if (mc_truth->Origin() == simb::kBeamNeutrino &&
           (mcpar->PdgCode() == 13 || mcpar->PdgCode() == -13 || 
            mcpar->PdgCode() == 12 || mcpar->PdgCode() == -12)) {
          if (_debug) std::cout << ">>>>>>>>>>>>>>>>> A neutrino related PFP (with ID " << pfpFlashTagged.at(i)->Self() << ") was tagged by the flash tagger with tag " << tagid_v[i] << std::endl;
          _nu_pfp_flash_tagged = 1;
          _nu_pfp_tagged_total = 1;
        }
      }
    }
  }




  // ****
  // Geometry
  // ****

  lar_pandora::PFParticleVector pfpGeoTagged;
  this->GetTaggedPFP(e, _cosmic_geo_tag_producer, 0.6, pfpGeoTagged,  tagid_v);
  if (_debug) std::cout << pfpGeoTagged.size() << " PFP have been geo tagged." << std::endl;
  _n_pfp_geo_tagged = pfpGeoTagged.size();

   // Loop through the taggedPFP and see if there is a neutrino related one
  _nu_pfp_geo_tagged = 0;
  for (unsigned int i = 0; i < pfpGeoTagged.size(); i++) {
    auto mcghosts = mcghost_from_pfp.at(pfpGeoTagged.at(i).key());
    if (mcghosts.size() > 0) {
      art::Ptr<simb::MCParticle> mcpar = mcpar_from_mcghost.at(mcghosts.at(0).key()).at(0);
      const auto mc_truth = UBXSecHelper::TrackIDToMCTruth(e, "largeant", mcpar->TrackId());
      if (mc_truth) {
        if (mc_truth->Origin() == simb::kBeamNeutrino &&
           (mcpar->PdgCode() == 13 || mcpar->PdgCode() == -13 || 
            mcpar->PdgCode() == 12 || mcpar->PdgCode() == -12)) { 
          if (_debug) std::cout << ">>>>>>>>>>>>>>>>> A neutrino related PFP (with ID " << pfpGeoTagged.at(i)->Self() << ") was tagged by the geo tagger with tag " << tagid_v[i] << std::endl;
          _nu_pfp_geo_tagged = 1;
          _nu_pfp_tagged_total = 1;
        }
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
    auto mcghosts = mcghost_from_pfp.at(pfpACPTTagged.at(i).key());
    if (mcghosts.size() > 0) {
      art::Ptr<simb::MCParticle> mcpar = mcpar_from_mcghost.at(mcghosts.at(0).key()).at(0);
      const auto mc_truth = UBXSecHelper::TrackIDToMCTruth(e, "largeant", mcpar->TrackId());
      if (mc_truth) {
        if (mc_truth->Origin() == simb::kBeamNeutrino &&
           (mcpar->PdgCode() == 13 || mcpar->PdgCode() == -13 || 
            mcpar->PdgCode() == 12 || mcpar->PdgCode() == -12)) { 
          if (_debug) std::cout << ">>>>>>>>>>>>>>>>> ACPT A neutrino related PFP (with ID " << pfpACPTTagged.at(i)->Self() << ") was tagged by the acpt tagger with tag " << tagid_v[i] << std::endl;
          _nu_pfp_acpt_tagged = 1;
          _nu_pfp_tagged_total = 1;
        }
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
    auto mcghosts = mcghost_from_pfp.at(pfpStopMuTagged.at(i).key());
    if (mcghosts.size() > 0) {
      art::Ptr<simb::MCParticle> mcpar = mcpar_from_mcghost.at(mcghosts.at(0).key()).at(0);
      const auto mc_truth = UBXSecHelper::TrackIDToMCTruth(e, "largeant", mcpar->TrackId());
      if (mc_truth) {
        if (mc_truth->Origin() == simb::kBeamNeutrino &&
           (mcpar->PdgCode() == 13 || mcpar->PdgCode() == -13 || 
            mcpar->PdgCode() == 12 || mcpar->PdgCode() == -12)) { 
          if (_debug) std::cout << ">>>>>>>>>>>>>>>>> STOPMU A neutrino related PFP (with ID " << pfpStopMuTagged.at(i)->Self() << ") was tagged by the stopmu tagger with tag " << tagid_v[i] << std::endl;
          _nu_pfp_stopmu_tagged = 1;
          _nu_pfp_tagged_total = 1;
        }
      }
    }
  }

  _tree1->Fill();
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
int CosmicTaggerAna::PFPInCommon(lar_pandora::PFParticleVector first, lar_pandora::PFParticleVector second){

  int nInCommon = 0;


  for (unsigned int f = 0; f < first.size(); f++){

    for (unsigned int s = 0; s < second.size(); s++){
 
      if(first.at(f) == second.at(s)) {

        nInCommon++;

        //std::cout << "In common found, flash is  " << flash << " and geo is " << geo << std::endl; 

      }
    }

  }
  return nInCommon;

}



DEFINE_ART_MODULE(CosmicTaggerAna)
