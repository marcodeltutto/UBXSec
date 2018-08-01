/**
 * \class UBXSecEvent
 *
 * \ingroup UBXSec
 *
 * \brief Data product to store a UBXSec Event
 * 
 *
 * \author Marco Del Tutto <marco.deltutto@physics.ox.ac.uk>
 *
 * \version 1.0
 *
 * \date 2017/10/08
 *
 * Contact: marco.deltutto@physics.ox.ac.uk
 *
 * Created on: Sunday, October 08, 2017 at 14:53:36
 *
 */


#ifndef UBXSecEvent_h
#define UBXSecEvent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>

using namespace std;

class UBXSecEvent /*: public TObject*/{
  public :
  
  // Declaration of leaf types
  Int_t           run; ///< Run number
  Int_t           subrun; ///< Subrun number
  Int_t           event; ///< Event number
  Bool_t          muon_is_reco; ///< Is true if the muon from the neutrino interaction is reconstructed
  Double_t        muon_reco_pur; ///< If reco, stores the reco muon purity
  Double_t        muon_reco_eff; ///< If reco, stores the reco muon efficiency
  Double_t        true_muon_mom; ///< If exists, stores the true muon momentum
  Double_t        true_muon_mom_matched; ///< True momentum of the MCP matched to the muon PFP
  Int_t           n_pfp; ///< Number of PFP in the event
  Int_t           n_pfp_primary; ///< Number of primary PFP in the event (neutrino PFP)
  Int_t           n_primary_cosmic_pfp; ///< Number of primary PFP in the event from pandoraCosmic (primaries before the removal)
  Int_t           nPFPtagged; ///< Not used
  Int_t           muon_is_flash_tagged; ///< Not used
  Double_t        muon_tag_score; ///< Not used
  Double_t        fm_score; ///< Not used
  Int_t           fv; ///< Is 1 if the true neutrino vertex is in the fiducial volume
  Int_t           fv_sce; ///< Is 1 if the true neutrino vertex is in the fiducial volume (takes into account the sce correction)
  Int_t           ccnc; ///< Is 0 if CC, 1 if NC
  Int_t           mode; ///< Iteraction mode: 0=Quasi-elastic or Elastic, 1=Resonant (RES), 2=DIS, 3=Coherent production, 10=MEC
  Int_t           nupdg; ///< Neutrino flavour (pdg code)
  Bool_t          is_signal; ///< Is trues if the event is a true numu cc in FV
  Double_t        nu_e; ///< Stores the true neutrino energy
  Double_t        lep_costheta; ///< Lepton true cosThata angle at start
  Double_t        lep_phi; ///< Lepton true Phi angle at start
  Int_t           genie_mult; ///< Number of stable GENIE final state particles
  Int_t           genie_mult_ch; ///< Number of stable charged GENIE final state particles
  Double_t        bnb_weight; ///< BNB correction weight to correct nue flux
  Bool_t          is_selected; ///< True if event passed numu cc inclusive selection
  Int_t           selected_slice; ///< The index of the selected slice

  Double_t        sce_corr_x; ///< Space charge correction to be applied to the true nu vertex (to be summed on x)
  Double_t        sce_corr_y; ///< Space charge correction to be applied to the true nu vertex (to be summed on y)
  Double_t        sce_corr_z; ///< Space charge correction to be applied to the true nu vertex (to be summed on z)

  Int_t           mc_muon_contained; ///< Is 1 if the true mc muon is fully contained
  Int_t           is_swtriggered; ///< Is true if the event passed the software trigger
  Double_t        vtx_resolution; ///< Stores the vertex resolution
  Int_t           n_tpcobj_nu_origin; ///< Number of TPCObjects with neutrino origin in the event
  Int_t           n_tpcobj_cosmic_origin; ///< Number of TPCObjects with cosmic origin in the event
  Int_t           nslices; ///< Stores the number of TPCObjects in the event
  vector<double>   slc_flsmatch_score; ///< Flash matching score (-9999 means failed to match)
  vector<double>   slc_flsmatch_qllx; ///< Estimated X position given by flash matching
  vector<double>   slc_flsmatch_tpcx; ///< TPC X position given by the flash time
  vector<double>   slc_flsmatch_t0; ///< Time of the beam spill flash
  vector<double>   slc_flsmatch_hypoz; ///< Z center of the hypothesis flash
  vector<double>   slc_flsmatch_xfixed_chi2; ///< Not used
  vector<double>   slc_flsmatch_xfixed_ll; ///< Not used
  vector<double>   slc_flsmatch_cosmic_score; ///< Not used
  vector<double>   slc_flsmatch_cosmic_t0; ///< Not used
  vector<double>   slc_nuvtx_x; ///< Reconstructed neutrino vertex X (cm)
  vector<double>   slc_nuvtx_y; ///< Reconstructed neutrino vertex Y (cm)
  vector<double>   slc_nuvtx_z; ///< Reconstructed neutrino vertex Z (cm)
  vector<int>      slc_nuvtx_fv; ///< Is 1 if the recon neutrino vertex is in the fiducial volume
  vector<double>   slc_vtxcheck_angle; ///< Angle between longest tracks (if extist, -9999 otherwise)
  vector<int>      slc_origin; ///< Origin of the slice (0: neutrino, 1: cosmic, 2: mixed) (MC only)
  vector<int>      slc_origin_extra; ///< Origin extra of the slice (-1: not set, 0: stopping muon, 1: acpt, 2: ncpion, 3: ncproton) (MC only)
  vector<int>      slc_nhits_u; ///< Number of hits on the U plane
  vector<int>      slc_nhits_v; ///< Number of hits on the V plane
  vector<int>      slc_nhits_w; ///< Number of hits on the W plane
  vector<double>   slc_longesttrack_length; ///< Length of the longest track
  vector<double>   slc_longesttrack_phi; ///< Phi angle of the longest track
  vector<double>   slc_longesttrack_theta; ///< Cos(theta) of the longest track
  vector<bool>     slc_longesttrack_iscontained; ///< Is true if the longest track if fully contained
  vector<double>   slc_longestshower_length; ///< Length of the longest shower
  vector<double>   slc_longestshower_phi; ///< Phi angle of the longest shower
  vector<double>   slc_longestshower_theta; ///< Cos(theta) of the longest shower
  vector<double>   slc_longestshower_openangle; ///< Opening angle of the longest shower
  vector<double>   slc_longestshower_startx; ///< Start oiint in X of the longest shower
  vector<double>   slc_longestshower_starty; ///< Start oiint in Y of the longest shower
  vector<double>   slc_longestshower_startz; ///< Start oiint in Z of the longest shower
  vector<int>      slc_acpt_outoftime; ///< Not used
  vector<int>      slc_crosses_top_boundary; ///< Is 1 if the track crosses the top roof of the TPC
  vector<int>      slc_nuvtx_closetodeadregion_u; ///< Is 1 if the recon nu vertex is close to a dead region on the U plane
  vector<int>      slc_nuvtx_closetodeadregion_v; ///< Is 1 if the recon nu vertex is close to a dead region on the V plane
  vector<int>      slc_nuvtx_closetodeadregion_w; ///< Is 1 if the recon nu vertex is close to a dead region on the W plane
  vector<double>   slc_kalman_chi2; ///< Not used
  vector<int>      slc_kalman_ndof; ///< Not used
  vector<bool>     slc_passed_min_track_quality; ///< Is true if the TPCObject passed track minimum quality requirements
  vector<bool>     slc_passed_min_vertex_quality; ///< Is true if the TPCObject passed track minimum vertex requirements
  vector<double>   slc_n_intime_pe_closestpmt; ///< Not used
  vector<double>   slc_maxdistance_vtxtrack; ///< Not used
  vector<bool>     slc_geocosmictag; ///< Is true if the TPCObject is through-going as a whole
  vector<bool>     slc_consistency; ///< Is false if the TPCObject is not a consistent candidate
  vector<double>   slc_consistency_score; ///< Temporary
  vector<int>      slc_npfp; ///< Number of PFP in the TPCObject
  vector<int>      slc_ntrack; ///< Number of tracks in the TPCObject
  vector<int>      slc_nshower; ///< Number of showers in the TPCObject
  vector<bool>     slc_iscontained; ///< Is true if the slice is fully contained
  vector<int>      slc_mult_pfp; ///< PFP multiplicity
  vector<int>      slc_mult_track; ///< Track multiplicity
  vector<int>      slc_mult_shower; ///< Shower multiplicity
  vector<int>      slc_mult_track_tolerance; ///< Track multiplicity considering tracks n cm close from the reco vtx (default n=5)

  vector<bool>     slc_muoncandidate_exists; ///< Is true if we found a muon candidate for the TPCObject
  vector<double>   slc_muoncandidate_length; ///< Track length for the muon candidate in the TPCObject
  vector<double>   slc_muoncandidate_phi; ///< Phi angle for the muon candidate in the TPCObject
  vector<double>   slc_muoncandidate_theta; ///< Cos(theta) for the muon candidate in the TPCObject
  vector<double>   slc_muoncandidate_mom_range; ///< Momentum (by range) of the muon candidate in the TPCObject
  vector<double>   slc_muoncandidate_mom_mcs; ///< Momentum (by MCS) of the muon candidate in the TPCObject
  vector<double>   slc_muoncandidate_mom_mcs_pi; ///<  Momentum (by MCS) of the muon candidate in the TPCObject (using pion hypo)
  vector<double>   slc_muoncandidate_mcs_ll; ///< -LL of the MCS fit
  vector<bool>     slc_muoncandidate_contained; ///< Is true if the muon candidate in the TPCObject is fully contained
  vector<double>   slc_muoncandidate_dqdx_trunc; ///< dqdx truncated mean for the muon candidate, plane 2
  vector<double>   slc_muoncandidate_dqdx_u_trunc; ///< dqdx truncated mean for the muon candidate, plane 0
  vector<double>   slc_muoncandidate_dqdx_v_trunc; ///< dqdx truncated mean for the muon candidate, plane 1
  vector<vector<double> > slc_muoncandidate_dqdx_v; ///< dqdx for every hit for the muon candidate
  vector<bool>     slc_muoncandidate_mip_consistency; ///< true if the muon candidate pass mip consistency cut
  vector<bool>     slc_muoncandidate_mip_consistency2;
  vector<int>      slc_muoncandidate_truepdg; ///< True pdg code of the candated muon track
  vector<int>      slc_muoncandidate_trueorigin; ///< True origin of the candidate muon track
  vector<double>   slc_muoncandidate_mcs_delta_ll; ///< Delta LL from MCS fit
  vector<double>   slc_muoncandidate_residuals_mean; ///< Mean of the residuals (between hits and track)
  vector<double>   slc_muoncandidate_residuals_std; ///< Standard Deviation of the residuals (between hits and track)
  vector<int>      slc_muoncandidate_wiregap; ///< Biggest wire gap found in track path
  vector<int>      slc_muoncandidate_wiregap_dead; ///< Number of dead wires in the biggest wire gap found in track path
  vector<double>   slc_muoncandidate_linearity; ///< Linearity of the hit the track is made out of
  vector<double>   slc_muoncandidate_perc_used_hits_in_cluster; ///< Number of used hits in the cluster to make the track
  vector<double>   slc_muoncandidate_maxscatteringangle; ///< Maximum scattering angle along track

  Int_t            nbeamfls; ///< Number of beam flashes in the event
  vector<double>   beamfls_time; ///< Time of the beam flash
  vector<double>   beamfls_pe; ///< PE of the beam flash
  vector<double>   beamfls_z; ///< Z centre of the beam flash
  Int_t            candidate_flash_time; ///< Not used
  Double_t         candidate_flash_z; ///< Z position of flash in beam spill
  Bool_t           no_mcflash_but_op_activity; ///< Not used
  vector<vector<double> > beamfls_spec; ///< PE per PMT of the beam flashe
  vector<double>   numc_flash_spec; ///< PE per PMT of the neutrino MC flash
  vector<vector<double> > slc_flshypo_xfixed_spec; ///< Not used
  vector<vector<double> > slc_flshypo_spec; ///< PE per PMT of the hypothesis flash for the TPCObjetc
  Int_t           nsignal; ///< Not used

  vector<double>   tvtx_x; ///< True neutrino vertex X (cm)
  vector<double>   tvtx_y; ///< True neutrino vertex Y (cm)
  vector<double>   tvtx_z; ///< True neutrino vertex Z (cm)
  
  Double_t        pot; ///< Not used

  Int_t evtwgt_genie_pm1_nfunc; ///< Number of functions used for GENIE reweighting (pm1sigma)
  vector<std::string> evtwgt_genie_pm1_funcname; ///< Names of the functions used for GENIE reweighting (pm1sigma)
  vector<int> evtwgt_genie_pm1_nweight; ///< Number of weights per function name used for GENIE reweighting (pm1sigma)
  vector<vector<double>> evtwgt_genie_pm1_weight; ///< Weights per function name used for GENIE reweighting (pm1sigma)

  Int_t evtwgt_genie_multisim_nfunc; ///< Number of functions used for GENIE reweighting (multisim)
  vector<std::string> evtwgt_genie_multisim_funcname; ///< Names of the functions used for GENIE reweighting (multisim)
  vector<int> evtwgt_genie_multisim_nweight; ///< Number of weights per function name used for GENIE reweighting (multisim)
  vector<vector<double>> evtwgt_genie_multisim_weight; ///< Weights per function name used for GENIE reweighting (multisim)

  Int_t evtwgt_genie_models_multisim_nfunc; ///< Number of functions used for GENIE Models reweighting (multisim)
  vector<std::string> evtwgt_genie_models_multisim_funcname; ///< Names of the functions used for GENIE Models reweighting (multisim)
  vector<int> evtwgt_genie_models_multisim_nweight; ///< Number of weights per function name used for GENIE Models reweighting (multisim)
  vector<vector<double>> evtwgt_genie_models_multisim_weight; ///< Weights per function name used for GENIE Models reweighting (multisim)
 
  Int_t evtwgt_flux_multisim_nfunc; ///< Number of functions used for FLUX reweighting (multisim)
  vector<std::string> evtwgt_flux_multisim_funcname; ///< Names of the functions used for FLUX reweighting (multisim)
  vector<int> evtwgt_flux_multisim_nweight; ///< Number of weights per function name used for FLUX reweighting (multisim)
  vector<vector<double>> evtwgt_flux_multisim_weight; ///< Weights per function name used for FLUX reweighting (multisim)
 
  int _default_value = -9999; ///< Default value 

  UBXSecEvent();
  virtual ~UBXSecEvent();
  void Init();
  void ResizeVectors(int); 
  void ResizeGenieTruthVectors(int); 
  void ResetGenieEventWeightVectorsPM1();
  void ResetGenieEventWeightVectorsMultisim();
  void ResetGenieModelsEventWeightVectorsMultisim();
  void ResetFluxEventWeightVectorsMultisim();

};

#endif
