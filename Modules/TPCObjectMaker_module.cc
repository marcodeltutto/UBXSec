////////////////////////////////////////////////////////////////////////
// Class:       TPCObjectMaker
// Plugin Type: producer (art v2_05_00)
// File:        TPCObjectMaker_module.cc
//
// Generated at Mon May 15 10:56:05 2017 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

/**
 * \class TPCObjectMaker
 *
 * \ingroup UBXSec
 *
 * \brief Art producer module that creates TPCObjects
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
 * Created on: Mon May 15 10:56:05 2017
 *
 */

/** \addtogroup UBXSec

   @{*/
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

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larcore/Geometry/Geometry.h"

#include "uboone/LLBasicTool/GeoAlgo/GeoVector.h"
#include "uboone/LLBasicTool/GeoAlgo/GeoAABox.h"

#include "uboone/UBXSec/DataTypes/TPCObject.h"
#include "uboone/UBXSec/Algorithms/McPfpMatch.h"
#include "uboone/UBXSec/Algorithms/UBXSecHelper.h"
#include "uboone/UBXSec/Algorithms/TPCObjectFilter.h"

#include <memory>

const simb::Origin_t NEUTRINO_ORIGIN = simb::kBeamNeutrino;
const simb::Origin_t COSMIC_ORIGIN   = simb::kCosmicRay;

namespace ubana {
  class TPCObjectMaker;
}


class ubana::TPCObjectMaker : public art::EDProducer {
public:
  explicit TPCObjectMaker(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TPCObjectMaker(TPCObjectMaker const &) = delete;
  TPCObjectMaker(TPCObjectMaker &&) = delete;
  TPCObjectMaker & operator = (TPCObjectMaker const &) = delete;
  TPCObjectMaker & operator = (TPCObjectMaker &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  /**
   *  @brief Gets all the PFPs for a single Pandora slice, given the neutrino PFP as input
   *
   *  @param pfParticleList the list of PFP
   *  @param particle the neutrino PFP, input
   *  @param pfp_v output, a vector of PFP (the TPC object) */
  void CollectPFP(lar_pandora::PFParticleVector pfParticleList, art::Ptr<recob::PFParticle> particle, lar_pandora::PFParticleVector &pfp_v);

  /**
   *  @brief Gets all the tracks and PFP for a single Pandora slice
   *
   *  @param pfParticleToTrackMap map from PFP to tracks
   *  @param pfParticleToShowerMap map from PFP to showers
   *  @param pfp_v input, a vector of PFP (the TPC object)
   *  @param track_v output, a vector of tracks (the TPC object)   
   *  @param shower_v output, a vector of showers (the TPC object) */
  void CollectTracksAndShowers(lar_pandora::PFParticlesToTracks pfParticleToTrackMap, lar_pandora::PFParticlesToShowers pfParticleToShowerMap, lar_pandora::PFParticleVector pfp_v, lar_pandora::TrackVector &track_v, lar_pandora::ShowerVector &shower_v);

  /**
   *  @brief Gets the pfp, track and shower multiplicity for a neutrino PFP
   *
   *  @param pfParticleList the list of PFP
   *  @param particle the input neutrino PFP
   *  @param pfp_v input, a vector of PFP (the TPC object)
   *  @param p output, multiplicity in number of PFPs
   *  @param t output, multiplicity in number of tracks
   *  @param s output, multiplicity in number of showers */
  void GetMultiplicity(lar_pandora::PFParticleVector pfParticleList, lar_pandora::PFParticleVector pfp_v, art::Ptr<recob::PFParticle> particle, int & p, int & t, int & s);

  /**
   *  @brief Constructs TPC objects using Pandora PFP slices
   *
   *  @param pfParticleList the list of PFP
   *  @param pfParticleToTrackMap map from PFP to tracks
   *  @param pfParticleToShowerMap map from PFP to showers
   *  @param pfParticleToVertexMap map from PFP to vertices
   *  @param _pfp_producer the PFP producer module
   *  @param pfp_v_v output, a vector of vector of PFP (a vector of TPC objects)
   *  @param track_v_v output, a vector of vector of tracks (a vector of TPC objects)   
   *  @param shower_v_v output, a vector of vector of showers (a vector of TPC objects) 
   *  @param p_v output, multiplicity in number of PFPs
   *  @param t_v output, multiplicity in number of tracks
   *  @param s_v output, multiplicity in number of showers */
  void GetTPCObjects(lar_pandora::PFParticleVector pfParticleList, lar_pandora::PFParticlesToTracks pfParticleToTrackMap, lar_pandora::PFParticlesToShowers pfParticleToShowerMap, lar_pandora::PFParticlesToVertices  pfParticleToVertexMap, std::vector<lar_pandora::PFParticleVector> & pfp_v_v, std::vector<lar_pandora::TrackVector> & track_v_v, std::vector<lar_pandora::ShowerVector> & shower_v_v, std::vector<int> & p_v, std::vector<int> & t_v, std::vector<int> & s_v);


  /**
   *  @brief Returns the nu PFP from a TPC object
   *
   *  @param pfp_v the TPC object (vector of PFP) */
  art::Ptr<recob::PFParticle> GetNuPFP(lar_pandora::PFParticleVector pfp_v);

private:

  ubana::McPfpMatch mcpfpMatcher;
  ubana::TPCObjectFilter *_tpcobj_filter;

  bool _is_mc;

  std::string _pfp_producer;
  std::string _vertexLabel;
  std::string _trackLabel;
  std::string _showerLabel;
  std::string _hitfinderLabel;
  std::string _geantModuleLabel;
  std::string _spacepointLabel;
  std::string _mcpHitAssLabel;

  bool _use_premade_ass;
  bool _do_filter;
  bool _pandora_cosmic_mode;
  bool _debug;
};


ubana::TPCObjectMaker::TPCObjectMaker(fhicl::ParameterSet const & p)
{
  _pfp_producer        = p.get<std::string>("PFParticleProducer"); 
  _vertexLabel         = p.get<std::string>("VertexProducer");
  _trackLabel          = p.get<std::string>("TrackProducer");
  _showerLabel         = p.get<std::string>("ShowerProducer");
  _hitfinderLabel      = p.get<std::string>("HitProducer");
  _geantModuleLabel    = p.get<std::string>("GeantModule");
  _spacepointLabel     = p.get<std::string>("SpacePointProducer");
  _mcpHitAssLabel      = p.get<std::string>("MCPHitAssProducer", "pandoraCosmicHitRemoval");

  _use_premade_ass     = p.get<bool>       ("UsePremadeMCPHitAss");
  _do_filter           = p.get<bool>       ("FilterObjects");
  _pandora_cosmic_mode = p.get<bool>       ("PandoraCosmicMode", "false"); 
  _debug               = p.get<bool>       ("Debug");

  if (_do_filter) _tpcobj_filter = new ubana::TPCObjectFilter();

  produces< std::vector<ubana::TPCObject>>();
  produces< art::Assns<ubana::TPCObject, recob::Track>>();
  produces< art::Assns<ubana::TPCObject, recob::Shower>>();
  produces< art::Assns<ubana::TPCObject, recob::PFParticle>>();
  produces< art::Assns<ubana::TPCObject, recob::Vertex>>();
}

void ubana::TPCObjectMaker::produce(art::Event & e){

  if (_debug) std::cout << "[TPCObjectMaker] Starts" << std::endl;
 
  _is_mc = !e.isRealData();

  if (_is_mc && _use_premade_ass) 
    mcpfpMatcher.Configure(e, _pfp_producer, _spacepointLabel, _hitfinderLabel, _geantModuleLabel, _mcpHitAssLabel, lar_pandora::LArPandoraHelper::kAddDaughters);
  else if (_is_mc)
    mcpfpMatcher.Configure(e, _pfp_producer, _spacepointLabel, _hitfinderLabel, _geantModuleLabel);

  // Instantiate the output
  std::unique_ptr< std::vector< ubana::TPCObject > >                tpcObjectVector        (new std::vector<ubana::TPCObject>);
  std::unique_ptr< art::Assns<ubana::TPCObject, recob::Track>>      assnOutTPCObjectTrack  (new art::Assns<ubana::TPCObject, recob::Track>      );
  std::unique_ptr< art::Assns<ubana::TPCObject, recob::Shower>>     assnOutTPCObjectShower (new art::Assns<ubana::TPCObject, recob::Shower>     );
  std::unique_ptr< art::Assns<ubana::TPCObject, recob::PFParticle>> assnOutTPCObjectPFP    (new art::Assns<ubana::TPCObject, recob::PFParticle> );
  std::unique_ptr< art::Assns<ubana::TPCObject, recob::Vertex>>     assnOutTPCObjectVertex (new art::Assns<ubana::TPCObject, recob::Vertex> );

  // Get the needed services
  //::art::ServiceHandle<cheat::BackTracker> bt;
  ::art::ServiceHandle<geo::Geometry> geo;

  // Is nc?
  bool _is_nc = false;
  art::Handle< std::vector<simb::MCTruth> > mctruth_h;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if (e.getByLabel("generator",mctruth_h)) {
    art::fill_ptr_vector(mclist, mctruth_h);
    int iList = 0; // 1 nu int per spill
    _is_nc = (mclist[iList]->GetNeutrino().CCNC() == 1);
  }

  // Use LArPandoraHelper functions to collect Pandora information
  lar_pandora::PFParticleVector pfParticleList;              //vector of PFParticles
  lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, pfParticleList);

  // Collect vertices, tracks and shower
  lar_pandora::VertexVector           allPfParticleVertices;
  lar_pandora::PFParticlesToVertices  pfParticleToVertexMap;
  lar_pandora::LArPandoraHelper::CollectVertices(e, _vertexLabel, allPfParticleVertices, pfParticleToVertexMap);
  lar_pandora::TrackVector            allPfParticleTracks;
  lar_pandora::PFParticlesToTracks    pfParticleToTrackMap;
  lar_pandora::LArPandoraHelper::CollectTracks(e, _trackLabel, allPfParticleTracks, pfParticleToTrackMap);
  lar_pandora::ShowerVector           allPfParticleShowers;
  lar_pandora::PFParticlesToShowers   pfParticleToShowerMap;
  lar_pandora::LArPandoraHelper::CollectShowers(e, _showerLabel, allPfParticleShowers, pfParticleToShowerMap);

  std::vector<lar_pandora::TrackVector     > track_v_v;
  std::vector<lar_pandora::ShowerVector    > shower_v_v;
  std::vector<lar_pandora::PFParticleVector> pfp_v_v;
  std::vector<int> p_v, t_v, s_v;

  this->GetTPCObjects(pfParticleList, pfParticleToTrackMap, pfParticleToShowerMap, pfParticleToVertexMap, pfp_v_v, track_v_v, shower_v_v, p_v, t_v, s_v);


  // Do the MCParticle to PFParticle matching

  lar_pandora::PFParticlesToMCParticles matched_pfp_to_mcp_maps;    // This is a map: MCParticle to matched PFParticle
  if (_is_mc) {
    mcpfpMatcher.GetRecoToTrueMatches(matched_pfp_to_mcp_maps);
  }



  // Loop over true particle and find the pfp with cosmic and neutrino origin 

  std::vector<art::Ptr<recob::PFParticle>> neutrinoOriginPFP;
  std::vector<art::Ptr<recob::PFParticle>> cosmicOriginPFP;
  std::vector<art::Ptr<recob::PFParticle>> cosmicStoppingOriginPFP;
  std::vector<art::Ptr<recob::PFParticle>> neutrinoStoppingOriginPFP;
  std::vector<art::Ptr<recob::PFParticle>> protonNCOriginPFP;
  std::vector<art::Ptr<recob::PFParticle>> pionNCOriginPFP;
  cosmicOriginPFP.clear();
  cosmicStoppingOriginPFP.clear();
  neutrinoStoppingOriginPFP.clear();
  protonNCOriginPFP.clear();
  pionNCOriginPFP.clear();


  if (!_is_mc) goto constructobjects;

  for (lar_pandora::PFParticlesToMCParticles::const_iterator iter1 = matched_pfp_to_mcp_maps.begin(), iterEnd1 = matched_pfp_to_mcp_maps.end();
      iter1 != iterEnd1; ++iter1) {

    art::Ptr<recob::PFParticle> pf_par = iter1->first;    // The PFParticle 
    art::Ptr<simb::MCParticle>  mc_par = iter1->second;   // The matched MCParticle 

    const art::Ptr<simb::MCTruth> mc_truth = UBXSecHelper::TrackIDToMCTruth(e, "largeant", mc_par->TrackId()); //bt->TrackIDToMCTruth(mc_par->TrackId());

    if (!mc_truth) {
      std::cerr << "[TPCObjectMaker] Problem with MCTruth pointer." << std::endl;
      continue;
    }

    // Cosmic origin
    if (mc_truth->Origin() == COSMIC_ORIGIN) {
      if (_debug) std::cout << "[TPCObjectMaker] PFP " << pf_par->Self() << " has cosmic origin" << std::endl;
      cosmicOriginPFP.emplace_back(pf_par);

      // Check if this is a stopping muon in the TPC 
      ::geoalgo::Vector mcpar_end(mc_par->EndX(), mc_par->EndY(), mc_par->EndZ());

      ::geoalgo::AABox tpc_vol(0., (-1.)*(geo->DetHalfHeight()), 0.,
                               geo->DetHalfWidth()*2., geo->DetHalfHeight(), geo->DetLength());

      if((mc_par->PdgCode() == 13 || mc_par->PdgCode() == -13) && tpc_vol.Contain(mcpar_end)) {
        cosmicStoppingOriginPFP.emplace_back(pf_par); 
      }

    }

    // Neutrino origin
    if (mc_truth->Origin() == NEUTRINO_ORIGIN) {
      if (_debug) std::cout << "[TPCObjectMaker] PFP " << pf_par->Self() << " has neutrino origin" << std::endl;
      neutrinoOriginPFP.emplace_back(pf_par);

      // Check if this is a stopping muon in the TPC
      ::geoalgo::Vector mcpar_end(mc_par->EndX(), mc_par->EndY(), mc_par->EndZ());

      ::geoalgo::AABox tpc_vol(0., (-1.)*(geo->DetHalfHeight()), 0.,
                               geo->DetHalfWidth()*2., geo->DetHalfHeight(), geo->DetLength());
            
      if((mc_par->PdgCode() == 13 || mc_par->PdgCode() == -13) && tpc_vol.Contain(mcpar_end)) {
        neutrinoStoppingOriginPFP.emplace_back(pf_par);
      }

      // Check if it's a proton or pion from NC interaction
      if(_is_nc) {
        if (mc_par->PdgCode() == 2212)
          protonNCOriginPFP.emplace_back(pf_par);
        if (mc_par->PdgCode() == 211 || mc_par->PdgCode() == -211) 
          pionNCOriginPFP.emplace_back(pf_par);
      }
    }
  }




  // Construct TPCObjects
  constructobjects:

  for (size_t i = 0; i < pfp_v_v.size(); i++){

    ::ubana::TPCObject obj;

    // Set tracks
    std::vector<recob::Track> trk_v;
    trk_v.clear();
    for (auto t : track_v_v[i]) trk_v.emplace_back((*t));
    obj.SetTracks(trk_v);

    // Set PFPs
    std::vector<recob::PFParticle> pfp_v;
    pfp_v.clear();
    for (auto p : pfp_v_v[i]) pfp_v.emplace_back((*p));
    obj.SetPFPs(pfp_v);

    // Set vertex
    std::vector<art::Ptr<recob::Vertex>> vtx_v;
    art::Ptr<recob::PFParticle> pfp = this->GetNuPFP(pfp_v_v[i]);
    auto iter = pfParticleToVertexMap.find(pfp);
    if (iter != pfParticleToVertexMap.end()) {
      obj.SetVertex(*(iter->second[0]));
      vtx_v.resize(1);
      vtx_v.at(0) = iter->second[0];
    }

    // Set origin
    ::ubana::TPCObjectOrigin origin = ubana::kUnknown;
    if (_is_mc)
      origin = UBXSecHelper::GetSliceOrigin(neutrinoOriginPFP, cosmicOriginPFP, pfp_v_v[i]); 
    obj.SetOrigin(origin);

    // Set origin extra
    ::ubana::TPCObjectOriginExtra origin_extra = ubana::kNotSet;
    if (_is_mc) {
      if (origin == ubana::kCosmicRay) origin_extra = UBXSecHelper::GetSliceOriginExtra_Stopping(cosmicStoppingOriginPFP, pfp_v_v[i]);
      else {
        origin_extra = UBXSecHelper::GetSliceOriginExtra_NC(protonNCOriginPFP, pionNCOriginPFP, pfp_v_v[i]);
        if (origin_extra == ubana::kNotSet) {
          origin_extra = UBXSecHelper::GetSliceOriginExtra_Stopping(neutrinoStoppingOriginPFP, pfp_v_v[i]);
        }
      }
    }
    obj.SetOriginExtra(origin_extra);

    // Set Multiplicity
    obj.SetMultiplicity(p_v[i], t_v[i], s_v[i]);

    tpcObjectVector->emplace_back(obj);
    util::CreateAssn(*this, e, *tpcObjectVector, track_v_v[i],  *assnOutTPCObjectTrack);
    util::CreateAssn(*this, e, *tpcObjectVector, shower_v_v[i], *assnOutTPCObjectShower);
    util::CreateAssn(*this, e, *tpcObjectVector, pfp_v_v[i],    *assnOutTPCObjectPFP);
    util::CreateAssn(*this, e, *tpcObjectVector, vtx_v,         *assnOutTPCObjectVertex);
  }



  // Put TPCObjects into the Event
  e.put(std::move(tpcObjectVector)); 
  e.put(std::move(assnOutTPCObjectTrack));
  e.put(std::move(assnOutTPCObjectShower));
  e.put(std::move(assnOutTPCObjectPFP));
  e.put(std::move(assnOutTPCObjectVertex));


  if (_debug) std::cout << "[TPCObjectMaker] Ends" << std::endl;
}



//_____________________________________________________________________________________
art::Ptr<recob::PFParticle> ubana::TPCObjectMaker::GetNuPFP(lar_pandora::PFParticleVector pfp_v){

  for (unsigned int pfp = 0; pfp < pfp_v.size(); pfp++) {

    if(lar_pandora::LArPandoraHelper::IsNeutrino(pfp_v.at(pfp))) {
      return pfp_v.at(pfp);
    }
    if(_pandora_cosmic_mode && pfp_v.at(pfp)->IsPrimary()) {
      return pfp_v.at(pfp);
    }
  }

  std::cout << "[TPCObjectMaker] No neutrino PFP found." << std::endl;

  art::Ptr<recob::PFParticle> temp;
  return temp;

}



//___________________________________________________________________________________________________
void ubana::TPCObjectMaker::GetTPCObjects(lar_pandora::PFParticleVector pfParticleList,
                                          lar_pandora::PFParticlesToTracks pfParticleToTrackMap,
                                          lar_pandora::PFParticlesToShowers pfParticleToShowerMap, 
                                          lar_pandora::PFParticlesToVertices  pfParticleToVertexMap,
                                          std::vector<lar_pandora::PFParticleVector> & pfp_v_v,
                                          std::vector<lar_pandora::TrackVector> & track_v_v,
                                          std::vector<lar_pandora::ShowerVector> & shower_v_v,
                                          std::vector<int> & p_v, std::vector<int> & t_v, std::vector<int> & s_v) {

  track_v_v.clear();
  shower_v_v.clear();
  pfp_v_v.clear();
  p_v.clear();
  t_v.clear();
  s_v.clear();


  if (_debug) std::cout << "[TPCObjectMaker] Getting TPC Objects..." << std::endl;

  for (unsigned int n = 0; n < pfParticleList.size(); ++n) {
    const art::Ptr<recob::PFParticle> particle = pfParticleList.at(n);

    bool is_main_pfp = lar_pandora::LArPandoraHelper::IsNeutrino(particle);
    if (_pandora_cosmic_mode) is_main_pfp = particle->IsPrimary();

    if(is_main_pfp) {
      if (_debug) std::cout << "[TPCObjectMaker] \t >>> Creating TPC Object " << track_v_v.size() << " <<<"<< std::endl;

      if (_debug) std::cout << "[TPCObjectMaker] \t PFP " << particle->Self() << " is the " << (_pandora_cosmic_mode ? "primary" : "neutrino") << " PFP." << std::endl;
      if (_debug) std::cout << "[TPCObjectMaker] \t The main PFP " << (particle->IsPrimary() ? "is" : "is not") << " a primary" << std::endl;

      lar_pandora::TrackVector track_v;
      lar_pandora::ShowerVector shower_v;
      lar_pandora::PFParticleVector pfp_v;
      int p, t, s;

      if (_debug) std::cout << "[TPCObjectMaker] \t Before filtering..." << std::endl;

      // Collect PFPs for this TPC object
      this->CollectPFP(pfParticleList, particle, pfp_v);

      // Collect Tracks and Showers for this TPC object
      this->CollectTracksAndShowers(pfParticleToTrackMap, pfParticleToShowerMap, pfp_v, // input
                                    track_v, shower_v);                                 // output

      // If filtering is on, filter the PFP for this TPC object
      if (_debug) std::cout << "[TPCObjectMaker] \t After filtering..." << std::endl;
      lar_pandora::PFParticleVector filtered_pfp_v;
      if(_tpcobj_filter && _do_filter) {

        filtered_pfp_v = _tpcobj_filter->Filter(pfp_v, pfParticleToTrackMap, pfParticleToShowerMap, pfParticleToVertexMap);

        pfp_v = filtered_pfp_v;

        this->CollectTracksAndShowers(pfParticleToTrackMap, pfParticleToShowerMap, pfp_v, // input
                                      track_v, shower_v);                                 // output
      }

      // Calculate multiplicity for this TPC object
      this->GetMultiplicity(pfParticleList, pfp_v, particle, p, t, s);



      if (_debug) std::cout << "[TPCObjectMaker] \t Number of pfp for this TPC object: "    << pfp_v.size()   << std::endl;
      for (auto pfp : pfp_v) {
        if (_debug) std::cout << "[TPCObjectMaker] \t \t PFP " << pfp->Self() << " with pdg " << pfp->PdgCode();
        auto it = pfParticleToVertexMap.find(pfp);
        if (it == pfParticleToVertexMap.end()) {
           if (_debug) std::cout << " and vertex [vertex not available for this PFP]" << std::endl;
        } else {
          double xyz[3];
          (it->second)[0]->XYZ(xyz);
          if (_debug) std::cout << " and vertex " << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
        }
      }
      if (_debug) {
        std::cout << "[TPCObjectMaker] \t Number of tracks for this TPC object:  " << track_v.size()  << std::endl;
        std::cout << "[TPCObjectMaker] \t Number of showers for this TPC object: " << shower_v.size() << std::endl;
        std::cout << "[TPCObjectMaker] \t Multiplicity (PFP) for this TPC object:    " << p << std::endl;
        std::cout << "[TPCObjectMaker] \t Multiplicity (Track) for this TPC object:  " << t << std::endl;
        std::cout << "[TPCObjectMaker] \t Multiplicity (Shower) for this TPC object: " << s << std::endl;
        std::cout << "[TPCObjectMaker]" << std::endl;
      }

      pfp_v_v.emplace_back(pfp_v);
      track_v_v.emplace_back(track_v);
      shower_v_v.emplace_back(shower_v);
      p_v.emplace_back(p);
      t_v.emplace_back(t);
      s_v.emplace_back(s);

    } // end if neutrino
  } // end pfp loop
}


//______________________________________________________________________________________________________________________________________
void ubana::TPCObjectMaker::CollectPFP(lar_pandora::PFParticleVector pfParticleList,
                                       art::Ptr<recob::PFParticle> particle,
                                       lar_pandora::PFParticleVector &pfp_v) {

  pfp_v.emplace_back(particle);

  // And their daughters
  const std::vector<size_t> &daughterIDs = particle->Daughters();
  if(daughterIDs.size() == 0) return;
  else {
    for (unsigned int m = 0; m < daughterIDs.size(); ++m) {
      const art::Ptr<recob::PFParticle> daughter = pfParticleList.at(daughterIDs.at(m));
      // Recursive call
      this->CollectPFP(pfParticleList, daughter, pfp_v);
    }
  }

}


//______________________________________________________________________________________________________________________________________
void ubana::TPCObjectMaker::CollectTracksAndShowers(lar_pandora::PFParticlesToTracks pfParticleToTrackMap,
                                                    lar_pandora::PFParticlesToShowers pfParticleToShowerMap,
                                                    lar_pandora::PFParticleVector pfp_v,
                                                    lar_pandora::TrackVector &track_v,
                                                    lar_pandora::ShowerVector &shower_v) {

  // Cleaning
  track_v.clear();
  shower_v.clear();

  // Loop over the PFPs
  for (auto pfp : pfp_v) {

    if (_debug) std::cout << "[TPCObjectMaker] \t PFP " << pfp->Self() << " which " << (pfp->IsPrimary() ? "is" : "is not") << " a primary" << std::endl;

    auto iter1 = pfParticleToTrackMap.find(pfp);
    if (iter1 != pfParticleToTrackMap.end()) {
      lar_pandora::TrackVector tracks = iter1->second;
      if (_debug) std::cout << "[TPCObjectMaker] \t PFP " << pfp->Self() << " has " << tracks.size() << " tracks ass." << std::endl;
      for (unsigned int trk = 0; trk < tracks.size(); trk++) {
        track_v.emplace_back(tracks[trk]);
      }
    }

    auto iter2 = pfParticleToShowerMap.find(pfp);
    if (iter2 != pfParticleToShowerMap.end()) {
      lar_pandora::ShowerVector showers = iter2->second;
      if (_debug) std::cout << "[TPCObjectMaker] \t PFP " << pfp->Self() << " has " << showers.size() << " showers ass." << std::endl;
      for (unsigned int s = 0; s < showers.size(); s++) {
        shower_v.emplace_back(showers[s]);
      }
    }
  }

}



//______________________________________________________________________________________________________________________________________
void ubana::TPCObjectMaker::GetMultiplicity(lar_pandora::PFParticleVector pfParticleList, 
                                            lar_pandora::PFParticleVector pfp_v,
                                            art::Ptr<recob::PFParticle> particle,
                                            int & p,
                                            int & t,
                                            int & s) {

  // Initialize
  p = t = s = -1;

  // Input PFP has to be a neutrino
  if (!lar_pandora::LArPandoraHelper::IsNeutrino(particle)) {

    // Check if it's a primary at least, and in case return, otherwise raise exception
    if (particle->IsPrimary()) {
      if(_pandora_cosmic_mode) {
        return;
      } else {
        std::cout << "[TPCObjectMaker] Using ubana::TPCObjectMaker::GetMultiplicity with a primary but not a neutrino PFP as input. Will return now." << std::endl;
        return;
      }
    } else {
      std::cerr << "[TPCObjectMaker] Using ubana::TPCObjectMaker::GetMultiplicity with a non neutrino PFP as input. Exiting now." << std::endl;
      throw std::exception();
    }

  }

  // Reset
  p = t = s = 0;

  const std::vector<size_t> &daughterIDs = particle->Daughters();
  if(daughterIDs.size() == 0) {
    if (_debug) std::cout << "[TPCObjectMaker] No daughters for this neutrino PFP." << std::endl;
    return;
  }
  else {
    for (unsigned int m = 0; m < daughterIDs.size(); ++m) {

      const art::Ptr<recob::PFParticle> daughter = pfParticleList.at(daughterIDs.at(m));

      bool found_in_tpcobj = false;
      for (auto pfp : pfp_v) {
        if (daughter == pfp)
          found_in_tpcobj = true;
      }
      if (!found_in_tpcobj) continue;

      p++;
      if (lar_pandora::LArPandoraHelper::IsTrack(daughter))  t++;
      if (lar_pandora::LArPandoraHelper::IsShower(daughter)) s++;
    }
  }

}



DEFINE_ART_MODULE(ubana::TPCObjectMaker)

  /** @} */ // end of doxygen group
