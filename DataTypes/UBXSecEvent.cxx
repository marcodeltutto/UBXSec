#ifndef UBXSecEvent_cxx
#define UBXSecEvent_cxx

#include "UBXSecEvent.h"
//#include <TH2.h>
//#include <TStyle.h>
//#include <TCanvas.h>



UBXSecEvent::UBXSecEvent()
{
  Init();
}

UBXSecEvent::~UBXSecEvent()
{
}



void UBXSecEvent::Init()
{
  
  run = _default_value;
  subrun = _default_value;
  event = _default_value;
  muon_is_reco = false;
  muon_reco_pur = _default_value;
  muon_reco_eff = _default_value;
  true_muon_mom = _default_value;
  true_muon_mom_matched = _default_value;
  n_pfp = _default_value;
  n_pfp_primary = _default_value;
  n_primary_cosmic_pfp = _default_value;
  nPFPtagged = _default_value;
  muon_is_flash_tagged = _default_value;
  muon_tag_score = _default_value;
  fm_score = _default_value;
  fv = _default_value;
  ccnc = _default_value;
  mode = _default_value;
  nupdg = _default_value;
  is_signal = false;
  nu_e = _default_value;
  lep_costheta = _default_value; 
  lep_phi = _default_value;
  genie_mult = _default_value;
  genie_mult_ch = _default_value;
  bnb_weight = _default_value;
  is_selected = false;

  mc_muon_contained = _default_value;
  is_swtriggered = _default_value;
  vtx_resolution = _default_value;

  nslices = _default_value;
  n_tpcobj_nu_origin = _default_value;
  n_tpcobj_cosmic_origin = _default_value;
  
  nsignal = _default_value;
  pot = _default_value;
 
  no_mcflash_but_op_activity = false;

  ResizeVectors(0);

}

void UBXSecEvent::ResizeVectors(int vsize) {

  slc_flsmatch_score.resize(vsize, _default_value);
  slc_flsmatch_qllx.resize(vsize, _default_value);
  slc_flsmatch_tpcx.resize(vsize, _default_value);
  slc_flsmatch_t0.resize(vsize, _default_value);
  slc_flsmatch_hypoz.resize(vsize, _default_value);
  slc_flsmatch_xfixed_chi2.resize(vsize, _default_value);
  slc_flsmatch_xfixed_ll.resize(vsize, _default_value);
  slc_nuvtx_x.resize(vsize);
  slc_nuvtx_y.resize(vsize);
  slc_nuvtx_z.resize(vsize);
  slc_nuvtx_fv.resize(vsize);
  slc_vtxcheck_angle.resize(vsize);
  slc_origin.resize(vsize);
  slc_origin_extra.resize(vsize);
  slc_flshypo_xfixed_spec.resize(vsize);
  slc_flshypo_spec.resize(vsize);
  slc_nhits_u.resize(vsize, _default_value);
  slc_nhits_v.resize(vsize, _default_value);
  slc_nhits_w.resize(vsize, _default_value);
  slc_flsmatch_cosmic_score.resize(vsize, _default_value);
  slc_flsmatch_cosmic_t0.resize(vsize, _default_value);
  slc_longesttrack_length.resize(vsize, _default_value);
  slc_longesttrack_phi.resize(vsize, _default_value);
  slc_longesttrack_theta.resize(vsize, _default_value);
  slc_longesttrack_iscontained.resize(vsize, _default_value);
  slc_longestshower_length.resize(vsize, _default_value);
  slc_longestshower_phi.resize(vsize, _default_value);
  slc_longestshower_theta.resize(vsize, _default_value);
  slc_longestshower_openangle.resize(vsize, _default_value);
  slc_longestshower_startx.resize(vsize, _default_value);
  slc_longestshower_starty.resize(vsize, _default_value);
  slc_longestshower_startz.resize(vsize, _default_value);
  slc_muoncandidate_exists.resize(vsize, _default_value);
  slc_muoncandidate_length.resize(vsize, _default_value);
  slc_muoncandidate_phi.resize(vsize, _default_value);
  slc_muoncandidate_theta.resize(vsize, _default_value);
  slc_muoncandidate_mom_range.resize(vsize, _default_value);
  slc_muoncandidate_mom_mcs.resize(vsize, _default_value);
  slc_muoncandidate_mom_mcs_pi.resize(vsize, _default_value); 
  slc_muoncandidate_mcs_ll.resize(vsize, _default_value);
  slc_muoncandidate_contained.resize(vsize, _default_value);
  slc_muoncandidate_dqdx_trunc.resize(vsize, _default_value);
  slc_muoncandidate_dqdx_u_trunc.resize(vsize);
  slc_muoncandidate_dqdx_v_trunc.resize(vsize);
  slc_muoncandidate_dqdx_v.resize(vsize);
  slc_muoncandidate_mip_consistency.resize(vsize, true);
  slc_muoncandidate_mip_consistency2.resize(vsize, true);
  slc_muoncandidate_truepdg.resize(vsize, _default_value);
  slc_muoncandidate_trueorigin.resize(vsize, _default_value);
  slc_muoncandidate_mcs_delta_ll.resize(vsize, _default_value);
  slc_muoncandidate_residuals_mean.resize(vsize, _default_value);
  slc_muoncandidate_residuals_std.resize(vsize, _default_value);
  slc_muoncandidate_wiregap.resize(vsize, _default_value);
  slc_muoncandidate_wiregap_dead.resize(vsize, _default_value);
  slc_muoncandidate_linearity.resize(vsize, _default_value);
  slc_muoncandidate_perc_used_hits_in_cluster.resize(vsize, _default_value);
  slc_muoncandidate_maxscatteringangle.resize(vsize, _default_value);
  slc_acpt_outoftime.resize(vsize, _default_value);
  slc_crosses_top_boundary.resize(vsize, _default_value);
  slc_nuvtx_closetodeadregion_u.resize(vsize, _default_value);
  slc_nuvtx_closetodeadregion_v.resize(vsize, _default_value);
  slc_nuvtx_closetodeadregion_w.resize(vsize, _default_value);
  slc_kalman_chi2.resize(vsize, _default_value);
  slc_kalman_ndof.resize(vsize, _default_value);
  slc_passed_min_track_quality.resize(vsize, _default_value);
  slc_passed_min_vertex_quality.resize(vsize, _default_value);
  slc_n_intime_pe_closestpmt.resize(vsize, _default_value);
  slc_maxdistance_vtxtrack.resize(vsize, _default_value);
  slc_npfp.resize(vsize, _default_value);
  slc_ntrack.resize(vsize, _default_value);
  slc_nshower.resize(vsize, _default_value);
  slc_iscontained.resize(vsize, _default_value);
  slc_mult_pfp.resize(vsize, _default_value);
  slc_mult_track.resize(vsize, _default_value);
  slc_mult_shower.resize(vsize, _default_value);
  slc_mult_track_tolerance.resize(vsize, _default_value);
  slc_geocosmictag.resize(vsize, false);
  slc_consistency.resize(vsize, true);
  slc_consistency_score.resize(vsize, 0.);

  slc_othershowers_longest_length.resize(vsize, _default_value); 
  slc_othershowers_longest_startx.resize(vsize, _default_value); 
  slc_othershowers_longest_starty.resize(vsize, _default_value); 
  slc_othershowers_longest_startz.resize(vsize, _default_value); 
  slc_othershowers_longest_phi.resize(vsize, _default_value); 
  slc_othershowers_longest_theta.resize(vsize, _default_value); 
  slc_othershowers_longest_openangle.resize(vsize, _default_value); 

  slc_othershowers_forward_length.resize(vsize, _default_value); 
  slc_othershowers_forward_startx.resize(vsize, _default_value); 
  slc_othershowers_forward_starty.resize(vsize, _default_value); 
  slc_othershowers_forward_startz.resize(vsize, _default_value); 
  slc_othershowers_forward_phi.resize(vsize, _default_value); 
  slc_othershowers_forward_theta.resize(vsize, _default_value); 
  slc_othershowers_forward_openangle.resize(vsize, _default_value); 
  
  slc_othershowers_flashmatch_length.resize(vsize, _default_value); 
  slc_othershowers_flashmatch_startx.resize(vsize, _default_value); 
  slc_othershowers_flashmatch_starty.resize(vsize, _default_value); 
  slc_othershowers_flashmatch_startz.resize(vsize, _default_value); 
  slc_othershowers_flashmatch_phi.resize(vsize, _default_value); 
  slc_othershowers_flashmatch_theta.resize(vsize, _default_value); 
  slc_othershowers_flashmatch_openangle.resize(vsize, _default_value); 

}


void UBXSecEvent::ResetGenieEventWeightVectorsPM1() {

  evtwgt_genie_pm1_funcname.clear();
  evtwgt_genie_pm1_weight.clear();
  evtwgt_genie_pm1_nweight.clear();

}

void UBXSecEvent::ResetGenieEventWeightVectorsMultisim() {

  evtwgt_genie_multisim_funcname.clear();
  evtwgt_genie_multisim_weight.clear();
  evtwgt_genie_multisim_nweight.clear();

}

void UBXSecEvent::ResetFluxEventWeightVectorsMultisim() {

  evtwgt_flux_multisim_funcname.clear();
  evtwgt_flux_multisim_weight.clear();
  evtwgt_flux_multisim_nweight.clear();

}

#endif

