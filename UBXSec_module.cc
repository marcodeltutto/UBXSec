////////////////////////////////////////////////////////////////////////
// Class:       UBXSec
// Plugin Type: analyzer (art v2_05_00)
// File:        UBXSec_module.cc
//
// Generated at Fri Jan  27 09:44:39 2017 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

/**
 * \class UBXSec
 *
 * \ingroup UBXSec
 *
 * \brief Art analyzer module
 * 
 *
 * \author $Author: Marco Del Tutto<marco.deltutto@physics.ox.ac.uk> $
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2017/03/10 $
 *
 * Contact: marco.deltutto@physics.ox.ac.uk
 *
 * Created on: Friday, March 10, 2017 at 12:32:31
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
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "uboone/UBXSec/FlashMatch.h" // new!
#include "lardata/Utilities/AssociationUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "uboone/UBXSec/UBXSecHelper.h"

#include "TString.h"
#include "TTree.h"

class UBXSec;


class UBXSec : public art::EDAnalyzer {
public:
  explicit UBXSec(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  UBXSec(UBXSec const &) = delete;
  UBXSec(UBXSec &&) = delete;
  UBXSec & operator = (UBXSec const &) = delete;
  UBXSec & operator = (UBXSec &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  std::string _hitfinderLabel;
  std::string _pfp_producer;
  std::string _geantModuleLabel;
  std::string _spacepointLabel;
  std::string _cosmic_tag_producer;
  std::string _neutrino_flash_match_producer;
  std::string _cosmic_flash_match_producer;
  std::string _opflash_producer_beam;
  bool _recursiveMatching = false;
  bool _debug = true;

  const simb::Origin_t NEUTRINO_ORIGIN = simb::kBeamNeutrino;
  const simb::Origin_t COSMIC_ORIGIN   = simb::kCosmicRay;

  /// Maps used for PFParticle truth matching
  typedef std::map< art::Ptr<recob::PFParticle>, unsigned int > RecoParticleToNMatchedHits;
  typedef std::map< art::Ptr<simb::MCParticle>,  RecoParticleToNMatchedHits > ParticleMatchingMap;
  typedef std::set< art::Ptr<recob::PFParticle> > PFParticleSet;
  typedef std::set< art::Ptr<simb::MCParticle> >  MCParticleSet;


  TTree* _tree1;
  int _run, _subrun, _event;
  int _muon_is_reco;
  int _nPFPtagged, _muon_is_flash_tagged;
  double _muon_tag_score;
  double _fm_score;
  int _fv, _ccnc, _nupdg;
  double _nu_e;
  double _recon_muon_start_x, _recon_muon_start_y, _recon_muon_start_z;
  double _recon_muon_end_x, _recon_muon_end_y, _recon_muon_end_z;
  double _mc_muon_start_x, _mc_muon_start_y, _mc_muon_start_z;
  double _mc_muon_end_x, _mc_muon_end_y, _mc_muon_end_z;
  int _mc_muon_contained;

  int _nslices;
  std::vector<double> _slc_flsmatch_score, _slc_flsmatch_xfixed_chi2, _slc_flsmatch_xfixed_ll;
  std::vector<double> _slc_flsmatch_cosmic_score, _slc_flsmatch_cosmic_t0;
  std::vector<double> _slc_nuvtx_x, _slc_nuvtx_y, _slc_nuvtx_z;
  std::vector<int> _slc_nuvtx_fv;
  std::vector<int> _slc_origin;
  std::vector<int> _slc_nhits_u, _slc_nhits_v, _slc_nhits_w;
  std::vector<double> _slc_longesttrack_length;
  std::vector<int> _slc_acpt_outoftime;
  std::vector<int> _slc_crosses_top_boundary;
  std::vector<int> _slc_nuvtx_closetodeadregion;

  int _nbeamfls;
  std::vector<double> _beamfls_time, _beamfls_pe;
  std::vector<std::vector<double>> _beamfls_spec, _slc_flshypo_xfixed_spec;
  std::vector<double> _numc_flash_spec;
  int _nsignal;
  int _is_swtriggered;

  TTree* _tree2;
  int _total_matches, _nmatch;
  std::vector<double> _hypo_spec, _beam_spec, _fixx_spec;
  double _score;
  int _is_muon;
};


UBXSec::UBXSec(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p) 
{
  _pfp_producer                   = p.get<std::string>("PFParticleProducer");
  _hitfinderLabel                 = p.get<std::string>("HitProducer");
  _geantModuleLabel               = p.get<std::string>("GeantModule");
  _spacepointLabel                = p.get<std::string>("SpacePointProducer");
  _neutrino_flash_match_producer  = p.get<std::string>("NeutrinoFlashMatchProducer");
  _cosmic_flash_match_producer    = p.get<std::string>("CosmicFlashMatchProducer");
  _opflash_producer_beam          = p.get<std::string>("OpFlashBeamProducer");

  art::ServiceHandle<art::TFileService> fs;
  _tree1 = fs->make<TTree>("tree","");
  _tree1->Branch("run",                &_run,                "run/I");
  _tree1->Branch("subrun",             &_subrun,             "subrun/I");
  _tree1->Branch("event",              &_event,              "event/I");
  _tree1->Branch("muon_is_reco",       &_muon_is_reco,       "muon_is_reco/I");
  _tree1->Branch("nPFPtagged",         &_nPFPtagged,         "nPFPtagged/I");
  _tree1->Branch("muon_is_flash_tagged",      &_muon_is_flash_tagged,      "muon_is_flash_tagged/I");
  _tree1->Branch("muon_tag_score",     &_muon_tag_score,     "muon_tag_score/D");
  _tree1->Branch("fm_score",           &_fm_score,           "fm_score/D");
  _tree1->Branch("fv",                 &_fv,                 "fv/I");
  _tree1->Branch("ccnc",               &_ccnc,               "ccnc/I");
  _tree1->Branch("nupdg",              &_nupdg,              "nupdg/I");
  _tree1->Branch("nu_e",               &_nu_e,               "nu_e/D");
  _tree1->Branch("recon_muon_start_x", &_recon_muon_start_x, "recon_muon_start_x/D");
  _tree1->Branch("recon_muon_start_y", &_recon_muon_start_y, "recon_muon_start_y/D");
  _tree1->Branch("recon_muon_start_z", &_recon_muon_start_z, "recon_muon_start_z/D");
  _tree1->Branch("recon_muon_end_x",   &_recon_muon_end_x,   "recon_muon_end_x/D");
  _tree1->Branch("recon_muon_end_y",   &_recon_muon_end_y,   "recon_muon_end_y/D");
  _tree1->Branch("recon_muon_end_z",   &_recon_muon_end_z,   "recon_muon_end_z/D");
  _tree1->Branch("mc_muon_start_x",    &_mc_muon_start_x,    "mc_muon_start_x/D");
  _tree1->Branch("mc_muon_start_y",    &_mc_muon_start_y,    "mc_muon_start_y/D");
  _tree1->Branch("mc_muon_start_z",    &_mc_muon_start_z,    "mc_muon_start_z/D");
  _tree1->Branch("mc_muon_end_x",      &_mc_muon_end_x,      "mc_muon_end_x/D");
  _tree1->Branch("mc_muon_end_y",      &_mc_muon_end_y,      "mc_muon_end_y/D");
  _tree1->Branch("mc_muon_end_z",      &_mc_muon_end_z,      "mc_muon_end_z/D");
  _tree1->Branch("mc_muon_contained",  &_mc_muon_contained,  "mc_muon_contained/I");
  _tree1->Branch("is_swtriggered",     &_is_swtriggered,     "is_swtriggered/I");

  _tree1->Branch("nslices",            &_nslices,            "nslices/I");
  _tree1->Branch("slc_flsmatch_score", "std::vector<double>", &_slc_flsmatch_score);
  _tree1->Branch("slc_flsmatch_xfixed_chi2", "std::vector<double>", &_slc_flsmatch_xfixed_chi2);
  _tree1->Branch("slc_flsmatch_xfixed_ll", "std::vector<double>", &_slc_flsmatch_xfixed_ll);
  _tree1->Branch("slc_flsmatch_cosmic_score", "std::vector<double>", &_slc_flsmatch_cosmic_score);
  _tree1->Branch("slc_flsmatch_cosmic_t0", "std::vector<double>", &_slc_flsmatch_cosmic_t0);
  _tree1->Branch("slc_nuvtx_x",        "std::vector<double>", &_slc_nuvtx_x);
  _tree1->Branch("slc_nuvtx_y",        "std::vector<double>", &_slc_nuvtx_y);
  _tree1->Branch("slc_nuvtx_z",        "std::vector<double>", &_slc_nuvtx_z);
  _tree1->Branch("slc_nuvtx_fv",       "std::vector<int>",    &_slc_nuvtx_fv);
  _tree1->Branch("slc_origin",         "std::vector<int>",    &_slc_origin);
  _tree1->Branch("slc_nhits_u",        "std::vector<int>",    &_slc_nhits_u);
  _tree1->Branch("slc_nhits_v",        "std::vector<int>",    &_slc_nhits_v);
  _tree1->Branch("slc_nhits_w",        "std::vector<int>",    &_slc_nhits_w);
  _tree1->Branch("slc_longesttrack_length",        "std::vector<double>",    &_slc_longesttrack_length);
  _tree1->Branch("slc_acpt_outoftime", "std::vector<int>",    &_slc_acpt_outoftime);
  _tree1->Branch("slc_crosses_top_boundary", "std::vector<int>",    &_slc_crosses_top_boundary);
  _tree1->Branch("slc_nuvtx_closetodeadregion",  "std::vector<int>",    &_slc_nuvtx_closetodeadregion);
  _tree1->Branch("nbeamfls",           &_nbeamfls,            "nbeamfls/I");
  _tree1->Branch("beamfls_time",       "std::vector<double>", &_beamfls_time);
  _tree1->Branch("beamfls_pe",         "std::vector<double>", &_beamfls_pe);
  _tree1->Branch("beamfls_spec",       "std::vector<std::vector<double>>", &_beamfls_spec);
  _tree1->Branch("numc_flash_spec",    "std::vector<double>", &_numc_flash_spec);
  _tree1->Branch("slc_flshypo_xfixed_spec", "std::vector<std::vector<double>>", &_slc_flshypo_xfixed_spec);
  _tree1->Branch("nsignal",            &_nsignal,             "nsignal/I");


  _tree2 = fs->make<TTree>("matchtree","");
  _tree2->Branch("run",                &_run,                "run/I");
  _tree2->Branch("subrun",             &_subrun,             "subrun/I");
  _tree2->Branch("event",              &_event,              "event/I");
  _tree2->Branch("total_matches",      &_total_matches,      "total_matches/I");
  _tree2->Branch("nmatch",             &_nmatch,             "nmatch/I");
  _tree2->Branch("score",              &_score,              "score/D");
  _tree2->Branch("hypo_spec",          "std::vector<double>", &_hypo_spec);
  _tree2->Branch("beam_spec",          "std::vector<double>", &_beam_spec);
  _tree2->Branch("fixx_spec",          "std::vector<double>", &_fixx_spec);
  _tree2->Branch("is_muon",            &_is_muon,            "is_muon/I");
  _tree2->Branch("muon_is_reco",       &_muon_is_reco,       "muon_is_reco/I");

}

void UBXSec::analyze(art::Event const & e)
{

  if(_debug) std::cout << "********** UBXSec starts" << std::endl;
  if(_debug) std::cout << "event: " << e.id().event() << std::endl;

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();


  art::ServiceHandle<cheat::BackTracker> bt;

  // *******************
  // Pandora MCParticle to PFParticle matching
  // *******************

  // --- Collect hits
  lar_pandora::HitVector hitVector;
  lar_pandora::LArPandoraHelper::CollectHits(e, _hitfinderLabel, hitVector);

  // --- Collect tracks
  lar_pandora::TrackVector            allPfParticleTracks;
  lar_pandora::PFParticlesToTracks    pfParticleToTrackMap;
  lar_pandora::LArPandoraHelper::CollectTracks(e, _pfp_producer, allPfParticleTracks, pfParticleToTrackMap);

  // --- Collect PFParticles and match Reco Particles to Hits
  lar_pandora::PFParticleVector  recoParticleVector;
  lar_pandora::PFParticleVector  recoNeutrinoVector;
  lar_pandora::PFParticlesToHits recoParticlesToHits;
  lar_pandora::HitsToPFParticles recoHitsToParticles;

  lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, recoParticleVector);
  lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(recoParticleVector, recoNeutrinoVector);
  lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e, _pfp_producer, _spacepointLabel, recoParticlesToHits, recoHitsToParticles, lar_pandora::LArPandoraHelper::kAddDaughters);

  if (_debug) {
    std::cout << "  RecoNeutrinos: " << recoNeutrinoVector.size() << std::endl;
    std::cout << "  RecoParticles: " << recoParticleVector.size() << std::endl;
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

  if (_debug) {
    std::cout << "  TrueParticles: " << particlesToTruth.size() << std::endl;
    std::cout << "  TrueEvents: " << truthToParticles.size() << std::endl;
  }

  lar_pandora::MCParticlesToPFParticles matchedParticles;    // This is a map: MCParticle to matched PFParticle
  lar_pandora::MCParticlesToHits        matchedParticleHits;

  // --- Do the matching
  UBXSecHelper::GetRecoToTrueMatches(recoParticlesToHits, 
                                        trueHitsToParticles, 
                                        matchedParticles, 
                                        matchedParticleHits);


  // *******************
  // Analysis
  // *******************

  std::vector<art::Ptr<recob::PFParticle>> taggedPFP;
  std::vector<double>                      taggedPFPscore;
  std::vector<art::Ptr<recob::PFParticle>> neutrinoOriginPFP;
  art::Ptr<recob::PFParticle>              muonPFP;

  _muon_is_reco = 0;
  _mc_muon_contained = 0;


  // Loop over true particle and find the cosmic related ones
  for (lar_pandora::MCParticlesToPFParticles::const_iterator iter1 = matchedParticles.begin(), iterEnd1 = matchedParticles.end();
      iter1 != iterEnd1; ++iter1) {

    art::Ptr<simb::MCParticle>  mc_par = iter1->first;   // The MCParticle 
    //art::Ptr<recob::PFParticle> pf_par = iter1->second;  // The matched PFParticle

    const art::Ptr<simb::MCTruth> mc_truth = bt->TrackIDToMCTruth(mc_par->TrackId());

    if (mc_truth->Origin() == COSMIC_ORIGIN) {

      double end[3];
      end[0] = mc_par->EndX();
      end[1] = mc_par->EndY();
      end[2] = mc_par->EndZ();
      if ( (mc_par->PdgCode() == 13 || mc_par->PdgCode() == -13) && UBXSecHelper::InFV(end) ){
        std::cout << "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Is stopping muon" << std::endl;
      }
         
    }
  }





  // Loop over true particle and find the neutrino related ones
  for (lar_pandora::MCParticlesToPFParticles::const_iterator iter1 = matchedParticles.begin(), iterEnd1 = matchedParticles.end();
             iter1 != iterEnd1; ++iter1) {

     art::Ptr<simb::MCParticle>  mc_par = iter1->first;   // The MCParticle 
     art::Ptr<recob::PFParticle> pf_par = iter1->second;  // The matched PFParticle

     const art::Ptr<simb::MCTruth> mc_truth = bt->TrackIDToMCTruth(mc_par->TrackId());
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
       }

       neutrinoOriginPFP.emplace_back(pf_par);

       // If we matched a muon
       if (mc_par->PdgCode() == 13) {
         muonPFP = pf_par;
         _muon_is_reco = 1;

         std::cout << "Here we are" << std::endl;

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

         if (UBXSecHelper::InFV(start) && UBXSecHelper::InFV(stop))
           _mc_muon_contained = 1;
       }  
     }
  }
  if (_debug) std::cout << "Neutrino related PFPs in this event: " << neutrinoOriginPFP.size() << std::endl;


 
  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  e.getByLabel(_pfp_producer, pfpHandle);
  art::FindManyP<ubana::FlashMatch> pfpToNeutrinoFlashMatchAssns(pfpHandle, e, _neutrino_flash_match_producer);
  art::FindManyP<ubana::FlashMatch> pfpToCosmicFlashMatchAssns(pfpHandle, e, _cosmic_flash_match_producer);

  // Get the FlashMatch tag from the ART event
  art::Handle<std::vector<ubana::FlashMatch>> flashMatchHandle;
  e.getByLabel(_cosmic_flash_match_producer, flashMatchHandle);

  if (!flashMatchHandle.isValid() || flashMatchHandle->empty()){
    std::cerr << "Cosmic tag is not valid or empty." << std::endl;
    //return;
  }

  // Look up the associations to pfparticle
  art::FindManyP<recob::PFParticle> flashMatchToPFPAssns(flashMatchHandle, e, _cosmic_flash_match_producer);

  if (_debug) std::cout << "Number of FlashMatch products in this event: " << flashMatchToPFPAssns.size() << std::endl;

  if (flashMatchToPFPAssns.size() == 0){
    std::cerr << "No flash match tags in this event." << std::endl;
    //return;
  }

  _total_matches = flashMatchToPFPAssns.size();

  // Loop over the flash match tags
  _nmatch = 0;
  for (unsigned int ct = 0; ct < flashMatchToPFPAssns.size(); ct++) {

    _nmatch++;

    // Get the flash match tag 
    art::Ptr<ubana::FlashMatch> flashMatch(flashMatchHandle, ct);
    if(_debug) std::cout << "This flash match product (" << ct << ") has score: " << flashMatch->GetScore() << std::endl;

    _fm_score = flashMatch->GetScore();
    _score = flashMatch->GetScore();
    _hypo_spec = flashMatch->GetHypoFlashSpec();
    _beam_spec = flashMatch->GetRecoFlashSpec();
    _fixx_spec = flashMatch->GetXFixedHypoFlashSpec();

    // Get the PFPs associated with this FM
    std::vector<art::Ptr<recob::PFParticle>> flashMatchToPFP_v = flashMatchToPFPAssns.at(flashMatch.key());
    if(_debug) std::cout << "Number of PFP associated with this Flash Match: " << flashMatchToPFP_v.size() << std::endl;

    _is_muon = 0;
    for (unsigned int pfp = 0; pfp < flashMatchToPFP_v.size(); pfp++){
      taggedPFP.emplace_back(flashMatchToPFP_v.at(pfp));
      taggedPFPscore.emplace_back(_fm_score);

      if(flashMatchToPFP_v.at(pfp) == muonPFP) {
        _is_muon = 1;
      }
    }

    _tree2->Fill();

  }

  if (_debug) std::cout << "Flash matched PFPs in this event: " << taggedPFP.size() << std::endl;
  _nPFPtagged = -1;
  _nPFPtagged = taggedPFP.size();

  
  // Loop through the taggedPFP and see if there is we flash matched the muon
  _muon_is_flash_tagged = 0;
  _muon_tag_score = -999;
  for (unsigned int i = 0; i < taggedPFP.size(); i++) {
    std::cout << "taggedPFP[" << i << "] has ID " << taggedPFP[i]->Self() << std::endl;
    for (unsigned int j = 0; j < neutrinoOriginPFP.size(); j++) {
      if(taggedPFP[i] == neutrinoOriginPFP[j]) {
        std::cout << ">>>>>>>>>>>>>>>>> A neutrino related PFP (with ID " << neutrinoOriginPFP[j]->Self() << ") was flash matched with score " << taggedPFPscore[i] << std::endl;
        if(taggedPFP[i] == muonPFP){
          std::cout << ">>>>>>>>>>>>>>>>> The muon recon PFP (with ID " << neutrinoOriginPFP[j]->Self() << ") was flash matched with score " << taggedPFPscore[i] << std::endl;
          _muon_is_flash_tagged = 1;
          _muon_tag_score = taggedPFPscore[i];
        }
      }
    }
  }
  


  // ACPT
  art::Handle<std::vector<anab::T0> > t0_h;
  e.getByLabel(Form("T0TrackTaggerCosmic%s",_pfp_producer.c_str()),t0_h);
  if(!t0_h.isValid()) {
    std::cout << "T0 product not found..." << std::endl;
    throw std::exception();
  }
  if(t0_h->empty()) {
    std::cout << "t0 is empty." << std::endl;
  }
    
  art::FindManyP<recob::Track> track_ptr_coll_v(t0_h, e, Form("T0TrackTaggerCosmic%s",_pfp_producer.c_str())); 

  art::Handle<std::vector<recob::Track>> track_h;
  e.getByLabel(_pfp_producer,track_h);
  if (!track_h.isValid() || track_h->empty()) {
    std::cout << "Track handle is not valid or empty." << std::endl;
    //throw std::exception();
  }

  art::FindManyP<recob::OpFlash> opfls_ptr_coll_v(track_h, e, Form("T0TrackTaggerCosmic%s",_pfp_producer.c_str()));





  // Check if truth nu in is FV
  // Collecting GENIE particles
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if (e.getByLabel("generator",mctruthListHandle))
    art::fill_ptr_vector(mclist, mctruthListHandle);

  int iList = 0; // 1 nu int per spill
  double truth_nu_vtx[3] = {mclist[iList]->GetNeutrino().Nu().Vx(),mclist[iList]->GetNeutrino().Nu().Vy(),mclist[iList]->GetNeutrino().Nu().Vz()};
  if (UBXSecHelper::InFV(truth_nu_vtx)) _fv = 1;
  else _fv = 0;

  _ccnc    = mclist[iList]->GetNeutrino().CCNC();
  _nupdg   = mclist[iList]->GetNeutrino().Nu().PdgCode();
  _nu_e    = mclist[iList]->GetNeutrino().Nu().E();

  _nsignal = 0;
  if(_nupdg==14 && _ccnc==0 && _fv==1) _nsignal=1; 

  // Save the number of slices in this event
  std::vector<lar_pandora::TrackVector     > track_v_v;
  std::vector<lar_pandora::PFParticleVector> pfp_v_v;

  UBXSecHelper::GetTPCObjects(e, _pfp_producer, pfp_v_v, track_v_v);

  _nslices = pfp_v_v.size();
  _slc_flsmatch_score.resize(_nslices, -9999);
  _slc_flsmatch_xfixed_chi2.resize(_nslices, -9999);
  _slc_flsmatch_xfixed_ll.resize(_nslices, -9999);
  _slc_nuvtx_x.resize(_nslices);
  _slc_nuvtx_y.resize(_nslices);
  _slc_nuvtx_z.resize(_nslices);
  _slc_nuvtx_fv.resize(_nslices);
  _slc_origin.resize(_nslices);
  _slc_flshypo_xfixed_spec.resize(_nslices);
  _slc_nhits_u.resize(_nslices, -9999);
  _slc_nhits_v.resize(_nslices, -9999);
  _slc_nhits_w.resize(_nslices, -9999);
  _slc_flsmatch_cosmic_score.resize(_nslices, -9999);
  _slc_flsmatch_cosmic_t0.resize(_nslices, -9999);
  _slc_longesttrack_length.resize(_nslices, -9999);
  _slc_acpt_outoftime.resize(_nslices, -9999);
  _slc_crosses_top_boundary.resize(_nslices, -9999);
  _slc_nuvtx_closetodeadregion.resize(_nslices, -9999);
    
  std::cout << "UBXSec - SAVING INFORMATION" << std::endl;
  for (unsigned int slice = 0; slice < pfp_v_v.size(); slice++){
    std::cout << ">>> SLICE" << slice << std::endl;

    // Slice origin (0 is neutrino, 1 is cosmic)
    _slc_origin[slice] = UBXSecHelper::GetSliceOrigin(neutrinoOriginPFP, pfp_v_v[slice]);

    // Reco vertex
    double reco_nu_vtx[3];
    UBXSecHelper::GetNuVertexFromTPCObject(e, _pfp_producer, pfp_v_v[slice], reco_nu_vtx);
    _slc_nuvtx_x[slice] = reco_nu_vtx[0];
    _slc_nuvtx_y[slice] = reco_nu_vtx[1];
    _slc_nuvtx_z[slice] = reco_nu_vtx[2];
    _slc_nuvtx_fv[slice] = (UBXSecHelper::InFV(reco_nu_vtx) ? 1 : 0);
    std::cout << "    Reco vertex saved" << std::endl;

    // Neutrino Flash match
    _slc_flsmatch_score[slice] = -9999;
    art::Ptr<recob::PFParticle> NuPFP = UBXSecHelper::GetNuPFP(pfp_v_v[slice]);
    //std::cout << "NuPFP has id " << NuPFP->Self() << std::endl;
    std::vector<art::Ptr<ubana::FlashMatch>> pfpToFlashMatch_v = pfpToNeutrinoFlashMatchAssns.at(NuPFP.key());
    if (pfpToFlashMatch_v.size() > 1) {
      std::cout << "    More than one flash match per nu pfp!" << std::endl;
      continue;
    } else if (pfpToFlashMatch_v.size() == 0){
    //  continue;
    } else {
      _slc_flsmatch_score[slice]       = pfpToFlashMatch_v[0]->GetScore(); 
      _slc_flsmatch_xfixed_chi2[slice] = pfpToFlashMatch_v[0]->GetXFixedChi2();
      _slc_flsmatch_xfixed_ll[slice]   = pfpToFlashMatch_v[0]->GetXFixedLl();
      _slc_flshypo_xfixed_spec[slice]  = pfpToFlashMatch_v[0]->GetXFixedHypoFlashSpec();
    }

    // Cosmic Flash Match
    _slc_flsmatch_cosmic_score[slice] = -9999;
    std::vector<art::Ptr<ubana::FlashMatch>> pfpToCosmicFlashMatch_v = pfpToCosmicFlashMatchAssns.at(NuPFP.key());
    if (pfpToCosmicFlashMatch_v.size() > 1) {
      std::cout << "    More than one flash match per nu pfp!" << std::endl;
      continue;
    } else if (pfpToCosmicFlashMatch_v.size() == 0){
      std::cout << "    PFP to flash match ass for cosmic is zero." << std::endl;
      //continue;
    } else if (pfpToCosmicFlashMatch_v.size() == 1){
      //std::cout << "pfpToCosmicFlashMatch_v[0]->GetScore() is " << pfpToCosmicFlashMatch_v[0]->GetScore() << std::endl;
      //std::cout << "pfpToCosmicFlashMatch_v[0]->GetT0() is " << pfpToCosmicFlashMatch_v[0]->GetT0() << std::endl;
      _slc_flsmatch_cosmic_score[slice] = pfpToCosmicFlashMatch_v[0]->GetScore();
      _slc_flsmatch_cosmic_t0[slice]    = pfpToCosmicFlashMatch_v[0]->GetT0();
    } else {
      std::cout << "    I don't know what fucking case this is." << std::endl;
    }

    // Hits
    int nhits_u, nhits_v, nhits_w;
    UBXSecHelper::GetNumberOfHitsPerPlane(e, _pfp_producer, track_v_v[slice], nhits_u, nhits_v, nhits_w);
    _slc_nhits_u[slice] = nhits_u;
    _slc_nhits_v[slice] = nhits_v;
    _slc_nhits_w[slice] = nhits_w;

    // Longest track and check boundary
    recob::Track lt;
    if (UBXSecHelper::GetLongestTrackFromTPCObj(track_v_v[slice], lt)){
      _slc_longesttrack_length[slice] = lt.Length();
      int vtx_ok;
      _slc_crosses_top_boundary[slice] = (UBXSecHelper::IsCrossingTopBoundary(lt, vtx_ok) ? 1 : 0);
    } else {
      _slc_longesttrack_length[slice] = -9999;
    }

    // ACPT
    _slc_acpt_outoftime[slice] = 0;
    for (unsigned int t = 0; t < track_v_v[slice].size(); t++) {
      if(opfls_ptr_coll_v.at(track_v_v[slice][t].key()).size()>1) {
        std::cout << "    More than 1 association found (ACPT)!" << std::endl;
        throw std::exception();
      } else if (opfls_ptr_coll_v.at(track_v_v[slice][t].key()).size()==0){
        continue;
      } else {
        art::Ptr<recob::OpFlash> flash_ptr = opfls_ptr_coll_v.at(track_v_v[slice][t].key()).at(0);
        if (flash_ptr->Time() < 3 || flash_ptr->Time() > 5)
        _slc_acpt_outoftime[slice] = 1;
      }
      
    }

    // Channel status
    _slc_nuvtx_closetodeadregion[slice] = (UBXSecHelper::PointIsCloseToDeadRegion(reco_nu_vtx, 2) ? 1 : 0);

    /*
    std::cout << "NEW--------------------- pfpToFlashMatch_v.size() " << pfpToFlashMatch_v.size() << std::endl;
    if ( std::find(taggedPFP.begin(), taggedPFP.end(), NuPFP) != taggedPFP.end() ){
      std::cout << "Slice " << slice << " was flash tagged" << std::endl;
      _slc_flsmatch_score[slice] = _fm_score;
    }
    */
    std::cout << "UBXSec - INFORMATION SAVED" << std::endl;
  }



  // Flashes
  ::art::Handle<std::vector<recob::OpFlash>> beamflash_h;
  e.getByLabel(_opflash_producer_beam,beamflash_h);
  if( !beamflash_h.isValid() || beamflash_h->empty() ) {
    std::cerr << "Don't have good flashes." << std::endl;
  }

  _nbeamfls = beamflash_h->size();
  _beamfls_pe.resize(_nbeamfls);
  _beamfls_time.resize(_nbeamfls);
  _beamfls_spec.resize(_nbeamfls);

  ::art::ServiceHandle<geo::Geometry> geo;

  for (size_t n = 0; n < beamflash_h->size(); n++) {
    auto const& flash = (*beamflash_h)[n];
    _beamfls_pe[n]   = flash.TotalPE();
    _beamfls_time[n] = flash.Time();

    _beamfls_spec[n].resize(32);
    for (unsigned int i = 0; i < 32; i++) {
      unsigned int opdet = geo->OpDetFromOpChannel(i);
      _beamfls_spec[n][opdet] = flash.PE(i);
    }
  }
  

  // SW Trigger
  art::Handle<raw::ubdaqSoftwareTriggerData> softwareTriggerHandle;
  e.getByLabel("swtrigger", softwareTriggerHandle);

  if (softwareTriggerHandle->getNumberOfAlgorithms() == 1) {
    std::vector<std::string> algoNames = softwareTriggerHandle->getListOfAlgorithms();
    std::cout << "SW trigger name: " << algoNames[0] << std::endl;
    _is_swtriggered = (softwareTriggerHandle->passedAlgo(algoNames[0]) ? 1 : 0);
  }


  // MC Flash
  ::art::Handle<std::vector<recob::OpFlash> > nuMcflash_h;
  e.getByLabel("NeutrinoMCFlash",nuMcflash_h);
  if( !nuMcflash_h.isValid() || nuMcflash_h->empty() ) {
    std::cerr << "Don't have neutrino MC flashes." << std::endl;
    return;
  }

  auto const& flash = (*nuMcflash_h)[0];
  _numc_flash_spec.resize(geo->NOpDets());
  for (unsigned int i = 0; i < geo->NOpDets(); i++) {
    unsigned int opdet = geo->OpDetFromOpChannel(i);
    _numc_flash_spec[opdet] = flash.PE(i);
  }

  _tree1->Fill();

  if(_debug) std::cout << "********** UBXSec ends" << std::endl;
}



DEFINE_ART_MODULE(UBXSec)
