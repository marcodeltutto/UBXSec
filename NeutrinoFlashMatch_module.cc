////////////////////////////////////////////////////////////////////////
// Class:       NeutrinoFlashMatch
// Plugin Type: producer (art v2_05_00)
// File:        NeutrinoFlashMatch_module.cc
//
// Generated at Wed Jan 25 10:00:40 2017 by Marco Del Tutto using cetskelgen
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
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/OpDetGeo.h"
#include "uboone/Geometry/UBOpReadoutMap.h"

#include "uboone/LLSelectionTool/OpT0Finder/Base/OpT0FinderTypes.h"
#include "uboone/LLBasicTool/GeoAlgo/GeoTrajectory.h"
#include "uboone/LLSelectionTool/OpT0Finder/Base/FlashMatchManager.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/LightPath.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/PhotonLibHypothesis.h"

#include "uboone/UBXSec/FlashMatch.h"
#include "Algorithms/UBXSecHelper.h"

#include "TTree.h"

#include <memory>

class NeutrinoFlashMatch;


class NeutrinoFlashMatch : public art::EDProducer {
public:
  explicit NeutrinoFlashMatch(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NeutrinoFlashMatch(NeutrinoFlashMatch const &) = delete;
  NeutrinoFlashMatch(NeutrinoFlashMatch &&) = delete;
  NeutrinoFlashMatch & operator = (NeutrinoFlashMatch const &) = delete;
  NeutrinoFlashMatch & operator = (NeutrinoFlashMatch &&) = delete;


  /**
   *  @brief Takes a vector of recob::Tracks and returns a QCluster made out from all the tracks in the vector
   *
   */
  flashana::QCluster_t GetQCluster(std::vector<art::Ptr<recob::Track>>);

  /**
   *  @brief Takes a vector of recob::PFParticle and two maps and returns a QCluster made out from all the tracks in the vector (test)
   *
   */
  flashana::QCluster_t GetQCluster(std::vector<art::Ptr<recob::PFParticle>>, lar_pandora::PFParticlesToSpacePoints pfp_to_spacept, lar_pandora::SpacePointsToHits spacept_to_hits);

  /**
   *  @brief Test method to calculate the flash-TPCobject compatibilty at a fixed x position
   *
   */
  flashana::Flash_t Trial(std::vector<art::Ptr<recob::Track>> track_v, flashana::Flash_t flashBeam, double & _chi2, double & _ll); 

  // Required functions.
  void produce(art::Event & e) override;

private:
  std::string _pfp_producer;           ///<
  std::string _track_producer;         ///<
  std::string _opflash_producer_beam;  ///<
  std::string _nuMcFlash_producer;     ///<
  double _flash_trange_start;          ///<
  double _flash_trange_end;            ///<
  bool _debug;                         ///<
  bool _use_genie_info;                ///<

  std::vector<::flashana::Flash_t>    beam_flashes;

  ::flashana::FlashMatchManager       _mgr;
  std::vector<flashana::FlashMatch_t> _result;

  std::vector<double>    _xfixed_hypo_spec;
  double _xfixed_chi2, _xfixed_ll;

  TTree* _tree1;
  int _run, _subrun, _event, _matchid, _flashid;
  std::vector<double>               _score, _t0;
  std::vector<double>               _qll_xmin, _tpc_xmin;
  std::vector<double>              _beam_flash_spec;
  std::vector<std::vector<double>> _hypo_flash_spec;
  std::vector<double>              _numc_flash_spec;
  int _fv, _ccnc, _nupdg;
};


NeutrinoFlashMatch::NeutrinoFlashMatch(fhicl::ParameterSet const & p)
{
  _pfp_producer            = p.get<std::string>("PFParticleModule",      "pandoraNu");
  _track_producer          = p.get<std::string>("TrackModule",           "pandoraNu");
  _nuMcFlash_producer      = p.get<std::string>("NeutrinoMCFlashModule", "NeutrinoMCFlash");
  _debug                   = p.get<bool>       ("DebugMode",             true);
  _opflash_producer_beam   = p.get<std::string>("BeamOpFlashProducer",   "simpleFlashBeam");
  _flash_trange_start      = p.get<double>     ("FlashVetoTimeStart",    3);
  _flash_trange_end        = p.get<double>     ("FlashVetoTimeEnd",      5);
  _use_genie_info          = p.get<bool>       ("UseGENIEInfo",          true); 
    
  _mgr.Configure(p.get<flashana::Config_t>("FlashMatchConfig"));

  if (_debug) {
    art::ServiceHandle<art::TFileService> fs;
    _tree1 = fs->make<TTree>("flashmatchtree","");
    _tree1->Branch("run",             &_run,                             "run/I");
    _tree1->Branch("subrun",          &_subrun,                          "subrun/I");
    _tree1->Branch("event",           &_event,                           "event/I");
    _tree1->Branch("beam_flash_spec", "std::vector<double>",             &_beam_flash_spec);
    _tree1->Branch("hypo_flash_spec", "std::vector<std::vector<double>>",&_hypo_flash_spec);
    _tree1->Branch("numc_flash_spec", "std::vector<double>",             &_numc_flash_spec);
    _tree1->Branch("score",           "std::vector<double>",             &_score);
    _tree1->Branch("t0",              "std::vector<double>",             &_t0);
    _tree1->Branch("qll_xmin",        "std::vector<double>",             &_qll_xmin);
    _tree1->Branch("tpc_xmin",        "std::vector<double>",             &_tpc_xmin);
    _tree1->Branch("xfixed_hypo_spec","std::vector<double>",             &_xfixed_hypo_spec);
    _tree1->Branch("fv",              &_fv,                              "fv/I");
    _tree1->Branch("ccnc",            &_ccnc,                            "ccnc/I");
    _tree1->Branch("nupdg",           &_nupdg,                           "nupdg/I");
  }


  produces< std::vector<ubana::FlashMatch>>();
  produces< art::Assns<ubana::FlashMatch,   recob::Track>>();
  produces< art::Assns<ubana::FlashMatch,   recob::PFParticle>>();
}

void NeutrinoFlashMatch::produce(art::Event & e)
{

  if(_debug) std::cout << "NeutrinoFlashMatch starts." << std::endl;

  // Instantiate the output
  std::unique_ptr< std::vector<ubana::FlashMatch>>                   flashMatchTrackVector      (new std::vector<ubana::FlashMatch>);
  std::unique_ptr< art::Assns<ubana::FlashMatch, recob::Track>>      assnOutFlashMatchTrack     (new art::Assns<ubana::FlashMatch, recob::Track>     );
  std::unique_ptr< art::Assns<ubana::FlashMatch, recob::PFParticle>> assnOutFlashMatchPFParticle(new art::Assns<ubana::FlashMatch, recob::PFParticle>);

  ::art::ServiceHandle<geo::Geometry> geo;

  _mgr.Reset();
  _result.clear();
  if(_debug) _mgr.PrintConfig();

  // Get Beam Flashes from the ART event
  ::art::Handle<std::vector<recob::OpFlash>> beamflash_h;
  e.getByLabel(_opflash_producer_beam,beamflash_h);
  if( !beamflash_h.isValid() || beamflash_h->empty() ) {
    std::cerr << "Don't have good flashes." << std::endl;
    e.put(std::move(flashMatchTrackVector));
    e.put(std::move(assnOutFlashMatchTrack));
    e.put(std::move(assnOutFlashMatchPFParticle));
    return;
  }
  int nBeamFlashes = 0;

  for (size_t n = 0; n < beamflash_h->size(); n++) {

    auto const& flash = (*beamflash_h)[n];

    if(flash.Time() < _flash_trange_start || _flash_trange_end < flash.Time()) {
      continue;
    }
    nBeamFlashes++;

    // Construct a Flash_t
    ::flashana::Flash_t f;
    f.x = f.x_err = 0;
    f.y = flash.YCenter();
    f.z = flash.ZCenter();
    f.y_err = flash.YWidth();
    f.z_err = flash.ZWidth();
    f.pe_v.resize(geo->NOpDets());
    f.pe_err_v.resize(geo->NOpDets());
    for (unsigned int i = 0; i < f.pe_v.size(); i++) {
      unsigned int opdet = geo->OpDetFromOpChannel(i);
      f.pe_v[opdet] = flash.PE(i);
      f.pe_err_v[opdet] = sqrt(flash.PE(i));
    }
    f.time = flash.Time();
    f.idx = nBeamFlashes-1;
    beam_flashes.resize(nBeamFlashes);
    beam_flashes[nBeamFlashes-1] = f;

    if (_debug && nBeamFlashes==1) {
      _beam_flash_spec.resize(f.pe_v.size());
      _beam_flash_spec = f.pe_v;
    } 
  } // flash loop

  // Don't waste other time if there are no flashes in the beam spill
  if (nBeamFlashes == 0) {
    std::cout << "Zero beam flashes in this event." << std::endl;
    e.put(std::move(flashMatchTrackVector));
    e.put(std::move(assnOutFlashMatchTrack));
    e.put(std::move(assnOutFlashMatchPFParticle));
    return;
  }

  // For now let's only consider cases where we have only 1 flash in the beam spill
  if (nBeamFlashes > 1) {
    std::cout << "More than one beam flashes in this event." << std::endl;
    e.put(std::move(flashMatchTrackVector));
    e.put(std::move(assnOutFlashMatchTrack));
    e.put(std::move(assnOutFlashMatchPFParticle));
    return;
  }


  // Emplace flash to Flash Matching Manager
  ::flashana::Flash_t f = beam_flashes[0];
  _mgr.Emplace(std::move(f));


  if(_debug && !e.isRealData()){
    // Also save the neutrino MC flash
    ::art::Handle<std::vector<recob::OpFlash> > nuMcflash_h;
    e.getByLabel(_nuMcFlash_producer,nuMcflash_h);
    if( !nuMcflash_h.isValid() || nuMcflash_h->empty() ) {
      std::cerr << "Don't have neutrino MC flashes." << std::endl;
      //e.put(std::move(flashMatchTrackVector));
      //e.put(std::move(assnOutFlashMatchTrack));
      //e.put(std::move(assnOutFlashMatchPFParticle));
      //return;
    } else {

      auto const& flash = (*nuMcflash_h)[0];
      _numc_flash_spec.resize(geo->NOpDets());
      for (unsigned int i = 0; i < geo->NOpDets(); i++) {
        unsigned int opdet = geo->OpDetFromOpChannel(i);
        _numc_flash_spec[opdet] = flash.PE(i);
      }
    }
  }



  // ********************
  // Construct TPC Objects
  // ********************

  std::vector<lar_pandora::TrackVector     > track_v_v; // the TPC obj as vector of tracks
  std::vector<lar_pandora::PFParticleVector> pfp_v_v;   // the TPC obj as vector of PFP
  lar_pandora::PFParticlesToSpacePoints pfp_to_spacept; // map from PFP to SpacePoint
  lar_pandora::SpacePointsToHits spacept_to_hits;       // map from SpacePoint to Hit

  std::vector<flashana::Flash_t> xfixed_hypo_v;
  std::vector<double> xfixed_chi2_v, xfixed_ll_v;

  //UBXSecHelper::GetTPCObjects(pfParticleList, pfParticleToTrackMap, pfParticleToVertexMap, pfp_v_v, track_v_v);
  UBXSecHelper::GetTPCObjects(e, _pfp_producer, _track_producer, pfp_v_v, track_v_v, pfp_to_spacept, spacept_to_hits);

  if(_debug) std::cout << "For this event we have " << track_v_v.size() << " pandora slices." << std::endl;

  xfixed_hypo_v.resize(track_v_v.size());
  xfixed_chi2_v.resize(track_v_v.size());
  xfixed_ll_v.resize(track_v_v.size());

  for (unsigned int tpcObj = 0; tpcObj < track_v_v.size(); tpcObj++) {

    lar_pandora::TrackVector tpcObjTrk_v      = track_v_v[tpcObj];
    lar_pandora::PFParticleVector tpcObjPfp_v = pfp_v_v[tpcObj];

    // Get QCluster for this TPC Object
    flashana::QCluster_t qcluster = this->GetQCluster(tpcObjPfp_v, pfp_to_spacept, spacept_to_hits);
    qcluster = this->GetQCluster(tpcObjTrk_v);

    qcluster.idx = tpcObj;

    // Emplace the QCluster to the FlashMatching Manager
    _mgr.Emplace(std::move(qcluster));

    double chi2, ll;
    xfixed_hypo_v[tpcObj] = this->Trial(tpcObjTrk_v, beam_flashes[0], chi2, ll);
    xfixed_chi2_v[tpcObj] = chi2;
    xfixed_ll_v[tpcObj] = ll;
  }

  if(_debug) std::cout << "Finished emplacing beam flash and tpc objects" << std::endl;

  // ********************
  // Run Flash Matching
  // ********************

  _result = _mgr.Match();


  // ********************
  // Save the results
  // ********************
  
  if(_debug) {
    std::cout << "Number of matches: " << _result.size() << std::endl;
    _hypo_flash_spec.resize(_result.size());
    _run    = e.id().run();
    _subrun = e.id().subRun();
    _event  = e.id().event();
  }

  _score.resize(_result.size());
  _t0.resize(_result.size()); 
  _qll_xmin.resize(_result.size());
  _tpc_xmin.resize(_result.size());

  for(_matchid=0; _matchid < (int)(_result.size()); ++_matchid) {

    auto const& match = _result[_matchid];

    _flashid         = match.flash_id;
    _score[_matchid] = match.score;

    auto const& flash = _mgr.FlashArray()[_flashid];
    _t0[_matchid] = flash.time;

    if(_debug) std::cout << "For this match, the score is " << match.score << std::endl;

    // Get the TPC obj 
    lar_pandora::TrackVector      track_v = track_v_v[match.tpc_id];
    lar_pandora::PFParticleVector pfp_v   = pfp_v_v[match.tpc_id];

    // Get hypo spec
    if(_debug){
      _hypo_flash_spec[_matchid].resize(geo->NOpDets());
      for(size_t pmt=0; pmt<_hypo_flash_spec[_matchid].size(); ++pmt) _hypo_flash_spec[_matchid][pmt] = match.hypothesis[pmt];
    }
    _xfixed_hypo_spec = xfixed_hypo_v[_matchid].pe_v; 
    _xfixed_chi2      = xfixed_chi2_v[_matchid];
    _xfixed_ll        = xfixed_ll_v[_matchid];

    // Save x position
    _qll_xmin[_matchid] = match.tpc_point.x;

    _tpc_xmin[_matchid] = 1.e4;
    for(auto const& pt : _mgr.QClusterArray()[match.tpc_id]) {
      if(pt.x < _tpc_xmin[_matchid]) _tpc_xmin[_matchid] = pt.x;
    }

    ubana::FlashMatch fm;
    fm.SetScore               ( _score[_matchid] );
    fm.SetTPCX                ( _tpc_xmin[_matchid] );
    fm.SetEstimatedX          ( _qll_xmin[_matchid] );
    fm.SetHypoFlashSpec       ( _hypo_flash_spec[_matchid] );
    fm.SetRecoFlashSpec       ( _beam_flash_spec );
    //fm.SetMCFlashSpec         ( _numc_flash_spec );
    fm.SetXFixedHypoFlashSpec ( _xfixed_hypo_spec );
    fm.SetXFixedChi2          ( _xfixed_chi2 );
    fm.SetXFixedLl            ( _xfixed_ll );

    flashMatchTrackVector->emplace_back(std::move(fm));
    util::CreateAssn(*this, e, *flashMatchTrackVector, track_v, *assnOutFlashMatchTrack);
    util::CreateAssn(*this, e, *flashMatchTrackVector, pfp_v,   *assnOutFlashMatchPFParticle);
  }


  e.put(std::move(flashMatchTrackVector));
  e.put(std::move(assnOutFlashMatchTrack));
  e.put(std::move(assnOutFlashMatchPFParticle));



  if (_debug && !e.isRealData() && _use_genie_info) {
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
  }

  if (_debug) _tree1->Fill();

  if (_debug) std::cout << "NeutrinoFlashMatch ends." << std::endl;


  return;
}




//______________________________________________________________________________________________________________________________________
flashana::QCluster_t NeutrinoFlashMatch::GetQCluster(std::vector<art::Ptr<recob::PFParticle>> pfp_v, lar_pandora::PFParticlesToSpacePoints pfp_to_spacept, lar_pandora::SpacePointsToHits spacept_to_hits) {

  flashana::QCluster_t summed_qcluster;
  summed_qcluster.clear();


  for (auto pfp : pfp_v) {

    // Get the spacepoints
    lar_pandora::SpacePointVector spacepoints;
    auto iter = pfp_to_spacept.find(pfp);
    if (iter != pfp_to_spacept.end()) {
      spacepoints = iter->second;
      std::cout << "[FlashMatching] Number of spacepoints ass to pfp with id " << pfp->Self() << " is: " << spacepoints.size() << std::endl;
    } else {
      std::cout << "[FlashMatching] Can't find ass spacepoints for pfp with id: " << pfp->Self() << "(pdg " << pfp->PdgCode() << ")" << std::endl;
      //return summed_qcluster;
    }

    // Get the hits per spacepoint
    for (auto sp_pt : spacepoints) {
      auto it = spacept_to_hits.find(sp_pt);
      if (it != spacept_to_hits.end()) {
        //std::cout << "[FlashMatching] Hit found. Hit plane is " << (it->second)->View() << std::endl;
      }
    }
  }



  return summed_qcluster;
}




//______________________________________________________________________________________________________________________________________
flashana::QCluster_t NeutrinoFlashMatch::GetQCluster(std::vector<art::Ptr<recob::Track>> track_v) {

  flashana::QCluster_t summed_qcluster;
  summed_qcluster.clear();

  for (unsigned int trk = 0; trk < track_v.size(); trk++) {

    art::Ptr<recob::Track> trk_ptr = track_v.at(trk);

    ::geoalgo::Trajectory track_geotrj;
    track_geotrj.resize(trk_ptr->NumberTrajectoryPoints(),::geoalgo::Vector(0.,0.,0.));

    for (size_t pt_idx=0; pt_idx < trk_ptr->NumberTrajectoryPoints(); ++pt_idx) {
      auto const& pt = trk_ptr->LocationAtPoint(pt_idx);
      track_geotrj[pt_idx][0] = pt[0];
      track_geotrj[pt_idx][1] = pt[1];
      track_geotrj[pt_idx][2] = pt[2];
    }

    auto qcluster = ((flashana::LightPath*)(_mgr.GetCustomAlgo("LightPath")))->FlashHypothesis(track_geotrj);
    summed_qcluster += qcluster;

  } // track loop

  return summed_qcluster;
}

//______________________________________________________________________________________________________________________________________
flashana::Flash_t NeutrinoFlashMatch::Trial(std::vector<art::Ptr<recob::Track>> track_v, flashana::Flash_t flashBeam, double & _chi2, double & _ll) {

  flashana::QCluster_t summed_qcluster;
  summed_qcluster.clear();

  double t0 = flashBeam.time;

  for (unsigned int trk = 0; trk < track_v.size(); trk++) {

    art::Ptr<recob::Track> trk_ptr = track_v.at(trk);

    ::geoalgo::Trajectory track_geotrj;
    track_geotrj.resize(trk_ptr->NumberTrajectoryPoints(),::geoalgo::Vector(0.,0.,0.));

    for (size_t pt_idx=0; pt_idx < trk_ptr->NumberTrajectoryPoints(); ++pt_idx) {
      auto const& pt = trk_ptr->LocationAtPoint(pt_idx);
      track_geotrj[pt_idx][0] = pt[0] - t0*0.1114359;
      track_geotrj[pt_idx][1] = pt[1];
      track_geotrj[pt_idx][2] = pt[2];
    }

    auto qcluster = ((flashana::LightPath*)(_mgr.GetCustomAlgo("LightPath")))->FlashHypothesis(track_geotrj);
    summed_qcluster += qcluster;

  } // track loop

  flashana::Flash_t flashHypo;
  flashHypo.pe_v.resize(32);
  ((flashana::PhotonLibHypothesis*)(_mgr.GetAlgo(flashana::kFlashHypothesis)))->FillEstimate(summed_qcluster,flashHypo);

  double O, H;
  _ll = 0;
  _chi2 = 0;

  //   Loop over PMTs and construct log-likelihood
  for (int pmt = 0; pmt < 32; pmt++){
    O = beam_flashes[0].pe_v[pmt];
    H = flashHypo.pe_v[pmt];

    //std::cout << "O: " << O << std::endl;
    //std::cout << "H: " << H << std::endl;

    if (H==0) continue;

    _chi2 += std::pow((O - H), 2) / (H);
    _ll -= std::log10(TMath::Poisson(O,H));

    //std::cout << "_chi2: " << _chi2 << std::endl;
    //std::cout << "_ll:   " << _ll << std::endl;
  }

  //std::cout << "------------ ------------ >>>>> TRIAL, -loglikelihood is: " << _ll << std::endl;

  return flashHypo;
}


DEFINE_ART_MODULE(NeutrinoFlashMatch)
