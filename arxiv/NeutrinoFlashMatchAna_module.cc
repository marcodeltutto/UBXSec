////////////////////////////////////////////////////////////////////////
// Class:       NeutrinoFlashMatchAna
// Plugin Type: analyzer (art v2_05_00)
// File:        NeutrinoFlashMatchAna_module.cc
//
// Generated at Fri Jan  27 09:44:39 2017 by Marco Del Tutto using cetskelgen
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

#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "uboone/UBXSec/FlashMatch.h" // new!
#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "TString.h"
#include "TTree.h"

class NeutrinoFlashMatchAna;


class NeutrinoFlashMatchAna : public art::EDAnalyzer {
public:
  explicit NeutrinoFlashMatchAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NeutrinoFlashMatchAna(NeutrinoFlashMatchAna const &) = delete;
  NeutrinoFlashMatchAna(NeutrinoFlashMatchAna &&) = delete;
  NeutrinoFlashMatchAna & operator = (NeutrinoFlashMatchAna const &) = delete;
  NeutrinoFlashMatchAna & operator = (NeutrinoFlashMatchAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  std::string _hitfinderLabel;
  std::string _pfp_producer;
  std::string _geantModuleLabel;
  std::string _spacepointLabel;
  std::string _cosmic_tag_producer;
  std::string _flash_match_producer;
  bool _recursiveMatching = false;
  bool _debug = true;

  const simb::Origin_t NEUTRINO_ORIGIN = simb::kBeamNeutrino;

  /// Maps used for PFParticle truth matching
  typedef std::map< art::Ptr<recob::PFParticle>, unsigned int > RecoParticleToNMatchedHits;
  typedef std::map< art::Ptr<simb::MCParticle>,  RecoParticleToNMatchedHits > ParticleMatchingMap;
  typedef std::set< art::Ptr<recob::PFParticle> > PFParticleSet;
  typedef std::set< art::Ptr<simb::MCParticle> >  MCParticleSet;

  bool InFV(double * nu_vertex_xyz);

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

  TTree* _tree1;
  int _run, _subrun, _event;
  int _muon_is_reco;
  int _nPFPtagged, _muonWasTagged;
  double _muon_tag_score;
  double _fm_score;
  int _fv, _ccnc, _nupdg;
  double _recon_muon_start_x, _recon_muon_start_y, _recon_muon_start_z;
  double _recon_muon_end_x, _recon_muon_end_y, _recon_muon_end_z;
  double _mc_muon_start_x, _mc_muon_start_y, _mc_muon_start_z;
  double _mc_muon_end_x, _mc_muon_end_y, _mc_muon_end_z;
  int _mc_muon_contained;

  TTree* _tree2;
  int _total_matches, _nmatch;
  std::vector<double> _hypo_spec, _beam_spec, _fixx_spec;
  double _score;
  int _is_muon;
};


NeutrinoFlashMatchAna::NeutrinoFlashMatchAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p) 
{
  _pfp_producer            = p.get<std::string>("PFParticleProducer");
  _hitfinderLabel          = p.get<std::string>("HitProducer");
  _geantModuleLabel        = p.get<std::string>("GeantModule");
  _spacepointLabel         = p.get<std::string>("SpacePointProducer");
  _flash_match_producer    = p.get<std::string>("FlashMatchProducer");

  art::ServiceHandle<art::TFileService> fs;
  _tree1 = fs->make<TTree>("tree","");
  _tree1->Branch("run",                &_run,                "run/I");
  _tree1->Branch("subrun",             &_subrun,             "subrun/I");
  _tree1->Branch("event",              &_event,              "event/I");
  _tree1->Branch("muon_is_reco",       &_muon_is_reco,       "muon_is_reco/I");
  _tree1->Branch("nPFPtagged",         &_nPFPtagged,         "nPFPtagged/I");
  _tree1->Branch("muonWasTagged",      &_muonWasTagged,      "muonWasTagged/I");
  _tree1->Branch("muon_tag_score",     &_muon_tag_score,     "muon_tag_score/D");
  _tree1->Branch("fm_score",           &_fm_score,           "fm_score/D");
  _tree1->Branch("fv",                 &_fv,                 "fv/I");
  _tree1->Branch("ccnc",               &_ccnc,               "ccnc/I");
  _tree1->Branch("nupdg",              &_nupdg,              "nupdg/I");
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

void NeutrinoFlashMatchAna::analyze(art::Event const & e)
{

  if(_debug) std::cout << "********** NeutrinoFlashMatchAna starts" << std::endl;
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
  this->GetRecoToTrueMatches(recoParticlesToHits, 
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

         if (this->InFV(start) && this->InFV(stop))
           _mc_muon_contained = 1;
       }  
     }
  }
  if (_debug) std::cout << "Neutrino related PFPs in this event: " << neutrinoOriginPFP.size() << std::endl;


  
  // Get the FlashMatch tag from the ART event
  art::Handle<std::vector<ubana::FlashMatch>> flashMatchHandle;
  e.getByLabel(_flash_match_producer, flashMatchHandle);

  if (!flashMatchHandle.isValid() || flashMatchHandle->empty()){
    std::cerr << "Cosmic tag is not valid or empty." << std::endl;
    return;
  }

  // Look up the associations to pfparticle
  art::FindManyP<recob::PFParticle> flashMatchToPFPAssns(flashMatchHandle, e, _flash_match_producer);

  if (_debug) std::cout << "Number of FlashMatch products in this event: " << flashMatchToPFPAssns.size() << std::endl;

  if (flashMatchToPFPAssns.size() == 0){
    std::cerr << "No flash match tags in this event." << std::endl;
    return;
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
  _muonWasTagged = 0;
  _muon_tag_score = -999;
  for (unsigned int i = 0; i < taggedPFP.size(); i++) {
    std::cout << "taggedPFP[" << i << "] has ID " << taggedPFP[i]->Self() << std::endl;
    for (unsigned int j = 0; j < neutrinoOriginPFP.size(); j++) {
      if(taggedPFP[i] == neutrinoOriginPFP[j]) {
        std::cout << ">>>>>>>>>>>>>>>>> A neutrino related PFP (with ID " << neutrinoOriginPFP[j]->Self() << ") was flash matched with score " << taggedPFPscore[i] << std::endl;
        if(taggedPFP[i] == muonPFP){
          std::cout << ">>>>>>>>>>>>>>>>> The muon recon PFP (with ID " << neutrinoOriginPFP[j]->Self() << ") was flash matched with score " << taggedPFPscore[i] << std::endl;
          _muonWasTagged = 1;
          _muon_tag_score = taggedPFPscore[i];
        }
      }
    }
  }
  


  



  // Check if truth nu in is FV
  // Collecting GENIE particles
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if (e.getByLabel("generator",mctruthListHandle))
    art::fill_ptr_vector(mclist, mctruthListHandle);

  int iList = 0; // 1 nu int per spill
  double truth_nu_vtx[3] = {mclist[iList]->GetNeutrino().Nu().Vx(),mclist[iList]->GetNeutrino().Nu().Vy(),mclist[iList]->GetNeutrino().Nu().Vz()};
  if (this->InFV(truth_nu_vtx)) _fv = 1;
  else _fv = 0;

  _ccnc    = mclist[iList]->GetNeutrino().CCNC();
  _nupdg   = mclist[iList]->GetNeutrino().Nu().PdgCode();
  



  _tree1->Fill();

  if(_debug) std::cout << "********** NeutrinoFlashMatchAna ends" << std::endl;
}


//______________________________________________________________________________________________________________________________________
bool NeutrinoFlashMatchAna::InFV(double * nu_vertex_xyz){

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

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoFlashMatchAna::GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, 
                                                const lar_pandora::HitsToMCParticles &trueHitsToParticles,
                                                lar_pandora::MCParticlesToPFParticles &matchedParticles, 
                                                lar_pandora::MCParticlesToHits &matchedHits) const
{   
  PFParticleSet recoVeto; MCParticleSet trueVeto;
    
  this->GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles, matchedParticles, matchedHits, recoVeto, trueVeto);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoFlashMatchAna::GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, 
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
        this->GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles, matchedParticles, matchedHits, vetoReco, vetoTrue);
}


DEFINE_ART_MODULE(NeutrinoFlashMatchAna)
