/**
 * \class UBXSecEvent
 *
 * \ingroup UBXSec
 *
 * \brief Data product to store a UBXSec Event
 * 
 *
 * \author $Author: Marco Del Tutto<marco.deltutto@physics.ox.ac.uk> $
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2017/10/08 $
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
  Int_t           run;
  Int_t           subrun;
  Int_t           event;
  Bool_t          muon_is_reco;
  Double_t        muon_reco_pur;
  Double_t        muon_reco_eff;
  Double_t        true_muon_mom;
  Double_t        true_muon_mom_matched;
  Int_t           nPFPtagged;
  Int_t           muon_is_flash_tagged;
  Double_t        muon_tag_score;
  Double_t        fm_score;
  Int_t           fv;
  Int_t           ccnc;
  Int_t           nupdg;
  Bool_t          is_signal;
  Double_t        nu_e;
  /*
   Double_t        recon_muon_start_x;
   Double_t        recon_muon_start_y;
   Double_t        recon_muon_start_z;
   Double_t        recon_muon_end_x;
   Double_t        recon_muon_end_y;
   Double_t        recon_muon_end_z;
   Double_t        mc_muon_start_x;
   Double_t        mc_muon_start_y;
   Double_t        mc_muon_start_z;
   Double_t        mc_muon_end_x;
   Double_t        mc_muon_end_y;
   Double_t        mc_muon_end_z;
   */
  Int_t           mc_muon_contained;
  Int_t           is_swtriggered;
  Double_t        vtx_resolution;
  Int_t           nslices;
  vector<double>   slc_flsmatch_score;
  vector<double>   slc_flsmatch_qllx;
  vector<double>   slc_flsmatch_tpcx;
  vector<double>   slc_flsmatch_t0;
  vector<double>   slc_flsmatch_hypoz;
  vector<double>   slc_flsmatch_xfixed_chi2;
  vector<double>   slc_flsmatch_xfixed_ll;
  vector<double>   slc_flsmatch_cosmic_score;
  vector<double>   slc_flsmatch_cosmic_t0;
  vector<double>   slc_nuvtx_x;
  vector<double>   slc_nuvtx_y;
  vector<double>   slc_nuvtx_z;
  vector<int>      slc_nuvtx_fv;
  vector<double>   slc_vtxcheck_angle;
  vector<int>      slc_origin;
  vector<int>      slc_origin_extra;
  vector<int>      slc_nhits_u;
  vector<int>      slc_nhits_v;
  vector<int>      slc_nhits_w;
  vector<double>   slc_longesttrack_length;
  vector<double>   slc_longesttrack_phi;
  vector<double>   slc_longesttrack_theta;
  vector<bool>     slc_longesttrack_iscontained;
  vector<int>      slc_acpt_outoftime;
  vector<int>      slc_crosses_top_boundary;
  vector<int>      slc_nuvtx_closetodeadregion_u;
  vector<int>      slc_nuvtx_closetodeadregion_v;
  vector<int>      slc_nuvtx_closetodeadregion_w;
  vector<double>   slc_kalman_chi2;
  vector<int>      slc_kalman_ndof;
  vector<bool>     slc_passed_min_track_quality;
  vector<bool>     slc_passed_min_vertex_quality;
  vector<double>   slc_n_intime_pe_closestpmt;
  vector<double>   slc_maxdistance_vtxtrack;
  vector<bool>     slc_geocosmictag;
  vector<int>      slc_npfp;
  vector<int>      slc_ntrack;
  vector<int>      slc_nshower;
  vector<bool>     slc_iscontained;
  vector<int>      slc_mult_pfp;
  vector<int>      slc_mult_track;
  vector<int>      slc_mult_shower;
  vector<int>      slc_mult_track_tolerance;
  vector<bool>     slc_muoncandidate_exists;
  vector<double>   slc_muoncandidate_length;
  vector<double>   slc_muoncandidate_phi;
  vector<double>   slc_muoncandidate_theta;
  vector<double>   slc_muoncandidate_mom_range;
  vector<double>   slc_muoncandidate_mom_mcs;
  vector<bool>     slc_muoncandidate_contained;
  Int_t            nbeamfls;
  vector<double>   beamfls_time;
  vector<double>   beamfls_pe;
  vector<double>   beamfls_z;
  Int_t            candidate_flash_time;
  Bool_t           no_mcflash_but_op_activity;
  vector<vector<double> > beamfls_spec;
  vector<double>   numc_flash_spec;
  vector<vector<double> > slc_flshypo_xfixed_spec;
  vector<vector<double> > slc_flshypo_spec;
  Int_t           nsignal;
  /*
   vector<double>   mctrk_start_x;
   vector<double>   mctrk_start_y;
   vector<double>   mctrk_start_z;
   vector<double>   trk_start_x;
   vector<double>   trk_start_y;
   vector<double>   trk_start_z;
   vector<double>   vtx_x;
   vector<double>   vtx_y;
   vector<double>   vtx_z;
   */
   vector<double>   tvtx_x;
   vector<double>   tvtx_y;
   vector<double>   tvtx_z;
  
  Double_t        pot;
 
  int _default_value = -9999;

  UBXSecEvent();
  virtual ~UBXSecEvent();
  void Init();
  void ResizeVectors(int); 

  //ClassDef(UBXSecEvent,1)  //Event Header
};

#endif
