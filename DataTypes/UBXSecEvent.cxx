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
  muon_is_reco = _default_value;
  muon_reco_pur = _default_value;
  muon_reco_eff = _default_value;
  true_muon_mom = _default_value;
  true_muon_mom_matched = _default_value;
  nPFPtagged = _default_value;
  muon_is_flash_tagged = _default_value;
  muon_tag_score = _default_value;
  fm_score = _default_value;
  fv = _default_value;
  ccnc = _default_value;
  nupdg = _default_value;
  is_signal = _default_value;
  nu_e = _default_value;
  
  
  
  mc_muon_contained = _default_value;
  is_swtriggered = _default_value;
  vtx_resolution = _default_value;
  nslices = _default_value;
  
  
  nsignal = _default_value;
  pot = _default_value;
  
  ResizeVectors(0);


  // Set object pointer
  /*
   slc_flsmatch_score = 0;
   slc_flsmatch_qllx = 0;
   slc_flsmatch_tpcx = 0;
   slc_flsmatch_t0 = 0;
   slc_flsmatch_hypoz = 0;
   slc_flsmatch_xfixed_chi2 = 0;
   slc_flsmatch_xfixed_ll = 0;
   slc_flsmatch_cosmic_score = 0;
   slc_flsmatch_cosmic_t0 = 0;
   slc_nuvtx_x = 0;
   slc_nuvtx_y = 0;
   slc_nuvtx_z = 0;
   slc_nuvtx_fv = 0;
   slc_vtxcheck_angle = 0;
   slc_origin = 0;
   slc_origin_extra = 0;
   slc_nhits_u = 0;
   slc_nhits_v = 0;
   slc_nhits_w = 0;
   slc_longesttrack_length = 0;
   slc_longesttrack_phi = 0;
   slc_longesttrack_theta = 0;
   slc_longesttrack_iscontained = 0;
   slc_acpt_outoftime = 0;
   slc_crosses_top_boundary = 0;
   slc_nuvtx_closetodeadregion_u = 0;
   slc_nuvtx_closetodeadregion_v = 0;
   slc_nuvtx_closetodeadregion_w = 0;
   slc_kalman_chi2 = 0;
   slc_kalman_ndof = 0;
   slc_passed_min_track_quality = 0;
   slc_passed_min_vertex_quality = 0;
   slc_n_intime_pe_closestpmt = 0;
   slc_maxdistance_vtxtrack = 0;
   slc_npfp = 0;
   slc_ntrack = 0;
   slc_nshower = 0;
   slc_iscontained = 0;
   slc_mult_pfp = 0;
   slc_mult_track = 0;
   slc_mult_shower = 0;
   slc_mult_track_tolerance = 0;
   beamfls_time = 0;
   beamfls_pe = 0;
   beamfls_z = 0;
   beamfls_spec = 0;
   numc_flash_spec = 0;
   slc_flshypo_xfixed_spec = 0;
   slc_flshypo_spec = 0;
   */
  /*
   mctrk_start_x = 0;
   mctrk_start_y = 0;
   mctrk_start_z = 0;
   trk_start_x = 0;
   trk_start_y = 0;
   trk_start_z = 0;
   vtx_x = 0;
   vtx_y = 0;
   vtx_z = 0;
   tvtx_x = 0;
   tvtx_y = 0;
   tvtx_z = 0;
   */
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
  slc_muoncandidate_exists.resize(vsize, _default_value);
  slc_muoncandidate_length.resize(vsize, _default_value);
  slc_muoncandidate_phi.resize(vsize, _default_value);
  slc_muoncandidate_theta.resize(vsize, _default_value);
  slc_muoncandidate_mom_range.resize(vsize, _default_value);
  slc_muoncandidate_mom_mcs.resize(vsize, _default_value);
  slc_muoncandidate_contained.resize(vsize, _default_value);
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

}
#endif

