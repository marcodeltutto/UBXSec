#ifndef STOPPINGMUONTAGGERHELPER_CXX
#define STOPPINGMUONTAGGERHELPER_CXX

#include "StoppingMuonTaggerHelper.h"
#include <iostream>
#include <algorithm>
#include <numeric> 
#include "TVector3.h"

namespace ubana {

  //_________________________________________________________________________________ 
  StoppingMuonTaggerHelper::StoppingMuonTaggerHelper()
  {
    if (_debug) std::cout << "StoppingMuonTaggerHelper Instantiated." << std::endl;

    _csvfile.open ("stopping_muon_tagger_helper.csv", std::ofstream::out | std::ofstream::trunc);
    _csvfile << "n,i,dqdx,dqdx_slider,linearity" << std::endl;
  }

  //_________________________________________________________________________________ 
  void StoppingMuonTaggerHelper::Configure(fhicl::ParameterSet const& pset)
  {

    _w2cm = pset.get< double > ( "WireToCmConstant", 0.3 );
    _t2cm = pset.get< double > ( "TimeToCmConstant", 0.0557 ); 
    _dqds_calib = pset.get< double > ( "GainCalib", 198 );
    _slider_window = pset.get< size_t > ( "SliderWindow", 10 );
    _max_allowed_hit_distance = pset.get< double > ( "MaxAllowedHitDistance", 15 );
    _slope_threshold = pset.get< double > ( "SlopeThreshold", 0.15 );  
    _hits_to_remove = pset.get< int > ( "HitsToRemove", 3 );
    _pre_post_window = pset.get< int > ( "PrePostWindow", 5 ); 
    _perc_diff_cut = pset.get< double > ( "PercDiffCut", 20 );
    _local_linearity_threshold = pset.get< double > ( "LocalLinerityThreshold", 0.85 );
    _min_muon_hits = pset.get< int > ( "MinMuonHits", 20 );
    _max_muon_hits = pset.get< int > ( "MaxMuonHits", 2000 );
    _min_michel_hits = pset.get< int > ( "MinMichelHits", 2 );
    _max_michel_hits = pset.get< int > ( "MaxMichelHits", 50 );
    _max_end_hits = pset.get< int > ( "MaxEndHits", 120 );  
    _debug = pset.get< bool > ( "DebugMode", false );

  }

  //_________________________________________________________________________________ 
  void StoppingMuonTaggerHelper::PrintConfig() {

    std::cout << "--- StoppingMuonTaggerHelper configuration:" << std::endl;
    std::cout << "---   _min_muon_hits               = " << _min_muon_hits << std::endl;
    std::cout << "---   _max_muon_hits               = " << _max_muon_hits << std::endl;
    std::cout << "---   _min_michel_hits             = " << _min_michel_hits << std::endl;
    std::cout << "---   _max_michel_hits             = " << _max_michel_hits << std::endl;
    std::cout << "---   _hits_to_remove              = " << _hits_to_remove << std::endl;
    std::cout << "---   _pre_post_window             = " << _pre_post_window << std::endl;
    std::cout << "---   _perc_diff_cut               = " << _perc_diff_cut << std::endl;
    std::cout << "---   _local_linearity_threshold   = " << _local_linearity_threshold << std::endl;
    std::cout << "---   _slope_threshold             = " << _slope_threshold << std::endl;
    std::cout << "---   _max_allowed_hit_distance    = " << _max_allowed_hit_distance << std::endl;
    std::cout << "---   _slider_window               = " << _slider_window << std::endl;
    std::cout << "---   _dqds_calib                  = " << _dqds_calib << std::endl;
    std::cout << "---   _w2cm                        = " << _w2cm << std::endl;
    std::cout << "---   _t2cm                        = " << _t2cm << std::endl;
  }

  //_________________________________________________________________________________ 
  void StoppingMuonTaggerHelper::Clear(){
     _s_hit_v.clear();
     _dqds_v.clear();
     _ds_v.clear();
     _dqds_slider.clear();

     _start_index = -1;
     _hits_ordered = false;
  }


  //_________________________________________________________________________________ 
  size_t StoppingMuonTaggerHelper::FilterByPlane(int plane_no) {

    std::vector<ubana::SimpleHit> new_vector;
    new_vector.clear();

    for (size_t element = 0; element < _s_hit_v.size(); element++) {

      if (_s_hit_v.at(element).plane == plane_no)
        new_vector.push_back(_s_hit_v.at(element));

    }

    std::swap(new_vector, _s_hit_v);

    return _s_hit_v.size();

  }


  //_________________________________________________________________________________
  size_t StoppingMuonTaggerHelper::FilterOnSingleWire() {

    if (!_hits_ordered) {
      std::cout << __PRETTY_FUNCTION__ << ": need to order hits first." << std::endl;
      throw std::exception();
    }

    if (_s_hit_v.size() < 4) {
      std::cout << __PRETTY_FUNCTION__ << ": _s_hit_v size is less than 4." << std::endl;
      return 0;
    }

    std::vector<ubana::SimpleHit> new_vector;
    new_vector.clear(); 

    std::vector<double> new_vector_ds;
    new_vector_ds.clear();

    std::vector<double> mean_v;
    mean_v.clear();

    std::vector<double> wire_v;

    for(const auto& window : this->get_windows(_s_hit_v, 4)) {

      for (auto s_h : window) {
        wire_v.push_back(s_h.wire);
      }

      mean_v.push_back(this->mean(wire_v));

      wire_v.clear();
    }

    new_vector.push_back(_s_hit_v.at(0));
    new_vector.push_back(_s_hit_v.at(1));
    new_vector_ds.push_back(_ds_v.at(0));
    new_vector_ds.push_back(_ds_v.at(1));

    for (size_t i = 2; i < mean_v.size()-1; i++) {

      //if (_debug) std::cout << "i: " << i 
      //                      << " wire : " << _s_hit_v.at(i).wire 
      //                      << " time : " << _s_hit_v.at(i).time 
      //                      << " mean: " << mean_v.at(i) << std::endl;

      if (std::abs(mean_v.at(i-1) - mean_v.at(i)) < 1    && 
          _s_hit_v.at(i-1).wire !=  _s_hit_v.at(i).wire  &&
          std::abs(mean_v.at(i)   - mean_v.at(i+1)) < 1  &&
          _s_hit_v.at(i).wire !=  _s_hit_v.at(i+1).wire    ) {

        std::cout << ">>>" << std::endl;
        if (_s_hit_v.at(i).integral > _s_hit_v.at(i+1).integral) {
          new_vector.push_back(_s_hit_v.at(i));
          new_vector_ds.push_back(_ds_v.at(i));
        } else {
          new_vector.push_back(_s_hit_v.at(i+1));
          new_vector_ds.push_back(_ds_v.at(i+1)); 
        }

        i++;

      } else {

        new_vector.push_back(_s_hit_v.at(i));
        new_vector_ds.push_back(_ds_v.at(i));
      }

    }

    new_vector.push_back(_s_hit_v.at(_s_hit_v.size()-1));
    new_vector_ds.push_back(_ds_v.at(_ds_v.size()-1)); 

    std::swap(new_vector, _s_hit_v);
    std::swap(new_vector_ds, _ds_v); 

    if (_s_hit_v.size() != _ds_v.size()) {
      std::cout << "_s_hit_v and _ds_v size mismatch" << std::endl;
      throw std::exception();
    }
    return _s_hit_v.size();
  }

  //_________________________________________________________________________________
  void StoppingMuonTaggerHelper::SetStartHit(double time, 
                                             double wire_no, 
                                             int plane_no) {


    if (plane_no != 2) {
      std::cout << "Plane not supported." << std::endl;
      return;
    }

    std::cout << "wire_no: " << wire_no << ", time: " << time << std::endl;

    TVector3 pt1(time, wire_no, 0);
    
    double min_dist = 1e9;
    int best_hit_id = -1;

    if (_debug) std::cout << "Simple hit vector size " << _s_hit_v.size() << std::endl;

    for (size_t i = 0; i < _s_hit_v.size(); i++) {

      auto sh = _s_hit_v.at(i);
      TVector3 pt2(sh.time, sh.wire, 0);
      double dist = (pt1-pt2).Mag();

      if (dist< min_dist) {
        best_hit_id = i;
        min_dist = dist;
      }
    }

    if (best_hit_id == -1) {
      std::cout << "Could not find start hit." << std::endl;
      return;
    }
   
 
    // The best hit may not be the start one, but one very close,
    // So let's find a border hit

    auto almost_best_hit = _s_hit_v.at(best_hit_id);

    if (_debug) std::cout << "Almost best hit has wire " << almost_best_hit.wire
                          << ", and time " << almost_best_hit.time*4 << std::endl;

    TVector3 pt0(almost_best_hit.time, almost_best_hit.wire, 0);

    // First create a map from wire to hit
    std::map<int,ubana::SimpleHit> wire_to_hit;

    for (size_t i = 0; i < _s_hit_v.size(); i++) {

      auto sh = _s_hit_v.at(i); 

      // If we never encountered this wire, just save it
      auto iter = wire_to_hit.find(sh.wire);
      if (iter == wire_to_hit.end()) {
        wire_to_hit[sh.wire] = sh;
        continue;
      }

      // Otherwise pick the one that is closer to the almost_best_hit
      auto previous_hit = iter->second;
      TVector3 pt1(previous_hit.time, previous_hit.wire, 0);
      double dist_previous = (pt0-pt1).Mag();
      TVector3 pt2(sh.time, sh.wire, 0);
      double dist_current = (pt0-pt2).Mag();

      if (dist_current < dist_previous)
        wire_to_hit[sh.wire] = sh;
      else 
        continue;
    }

    // Then try going to the left first
    int n_step_left = 0;
    int best_index_left = -1;
    for (int w = almost_best_hit.wire - 1; w > 0; w--) {

      auto iter = wire_to_hit.find(w);
      if (iter == wire_to_hit.end()) {
        break;
      }

      //ubana::SimpleHit temp = iter->second;
      auto it = std::find(_s_hit_v.begin(), _s_hit_v.end(), iter->second);
      best_index_left = it - _s_hit_v.begin();//w;//std::distance(_s_hit_v.begin(), it);
      TVector3 pt1(iter->second.time, iter->second.wire, 0);
      double dist = (pt0-pt1).Mag();
  
      if (dist > 4.) // Check this number! 
        break;

      n_step_left ++;

      if (_debug) std::cout << "Found hit on the left, wire " << iter->second.wire
                            << ", time " << iter->second.time*4 << std::endl;
    }

    // Then go rigth
    int n_step_right = 0;
    int best_index_right = -1;
    for (int w = almost_best_hit.wire + 1; w < 3456; w++) {

      auto iter = wire_to_hit.find(w);
      if (iter == wire_to_hit.end()) {
        break;
      }

      auto it = std::find(_s_hit_v.begin(), _s_hit_v.end(), iter->second);
      best_index_right = it - _s_hit_v.begin();//w;//std::distance(_s_hit_v.begin(), it);
      TVector3 pt1(iter->second.time, iter->second.wire, 0);
      double dist = (pt0-pt1).Mag(); 

      if (dist > 4.) // Check this number!  
        break; 

      n_step_right ++;

      if (_debug) std::cout << "Found hit on the right, wire " << iter->second.wire 
                            << ", time " << iter->second.time*4 << std::endl; 
    }

    _start_index = 0;

    if (n_step_left == 0 || n_step_right == 0) {
      _start_index = best_hit_id;
    } else if (n_step_left > n_step_right) {
      _start_index = best_index_right;
    } else if (n_step_right >= n_step_left) {
      _start_index = best_index_left;
    }


    if (_debug) std::cout << "[StoppingMuonTaggerHelper] Start hit set to w: " 
                          << _s_hit_v.at(_start_index).wire << ", and t: " 
                          << _s_hit_v.at(_start_index).time*4 << std::endl;


    return;
  }

  

  //_________________________________________________________________________________ 
  size_t StoppingMuonTaggerHelper::OrderHits() {

    if (_start_index < 0) {
      std::cout << "Start hit not set." << std::endl;
      return 0;
    }

    if ((size_t)_start_index >= _s_hit_v.size()) {
      std::cout << "Start hit not compatible with current hit vector." 
       << "Start hit index is " << _start_index 
       << ", while hit vector size is " << _s_hit_v.size() << std::endl;
      return 0;
    }


    ubana::SimpleHitVector new_vector;
    new_vector.clear();
    new_vector.reserve(_s_hit_v.size());

    _ds_v.clear();
    _ds_v.reserve(_s_hit_v.size());

    new_vector.push_back(_s_hit_v.at(_start_index));

    //if (_debug) {
    //  for (auto h : _s_hit_v) {
    //    std::cout << "BEFORE: " << h.wire << ", " << h.time*4 << std::endl;
    //  }
    //}

    _s_hit_v.erase(_s_hit_v.begin() + _start_index);

    double min_dist = 1e9; 
    int min_index = -1;

    while (_s_hit_v.size() != 0) {

      min_dist = 1e9;
      min_index = -1; 

      for (size_t i = 0; i < _s_hit_v.size(); i++){

        TVector3 pt1(new_vector.back().time, new_vector.back().wire, 0);
        TVector3 pt2(_s_hit_v.at(i).time,    _s_hit_v.at(i).wire,    0);
        double dist = (pt1 - pt2).Mag();

        if (dist < min_dist) {
          min_index = i;
          min_dist = dist;
        }
      }

      if (min_index < 0) {
         std::cout << "Logic fail." << std::endl;
         throw std::exception();
      }

      // Emplace the next hit in the new vector...
      if (min_dist < _max_allowed_hit_distance) {
        if (_debug) std::cout << "min_dist: " << min_dist <<std::endl;
        new_vector.push_back(_s_hit_v.at(min_index));
        _ds_v.push_back(min_dist);
      } else if (new_vector.size() > 5){
 
        // Calculate previous slope
        auto iter = new_vector.end();
        auto sh_3 = _s_hit_v.at(min_index);
        auto sh_2 = *(--iter);
        auto sh_1 = *(iter-5); // go 5 hits back
        double slope = (sh_2.time - sh_1.time) / (sh_2.wire - sh_1.wire);

        // Calculate the new slope
        double slope_new = (sh_3.time - sh_2.time) / (sh_3.wire - sh_2.wire);

        std::cout << "sh_1.wire: "<<sh_1.wire<<", sh_1.time: "<< sh_1.time << std::endl;
        std::cout << "sh_2.wire: "<<sh_2.wire<<", sh_2.time: "<< sh_2.time << std::endl;
        std::cout << "sh_3.wire: "<<sh_3.wire<<", sh_3.time: "<< sh_3.time << std::endl;
        if (_debug) std::cout << "Current slope : " << slope 
                              << " New slope: " << slope_new 
                              << " Diff: " << slope_new - slope << std::endl;

        // Check the next hit will be in a consecutive wire
        bool progressive_order = false;

        if (sh_1.wire < sh_2.wire) {
          if (sh_3.wire > sh_2.wire) {
            progressive_order = true;
          }
        }
        if (sh_2.wire < sh_1.wire) {
          if (sh_3.wire < sh_2.wire) {
            progressive_order = true;
          }
        }

        std::cout << "Progressive order? " << (progressive_order ? "YES" : "NO") << std::endl;

        // If the two slopes are close, than there is 
        // probably a dead region between the point.
        // If so, increase the min distance by half a meter
        // and add the hit.
        if (std::abs(slope_new - slope) < _slope_threshold &&
            min_dist < _max_allowed_hit_distance + 50 &&
            progressive_order) {

          new_vector.push_back(_s_hit_v.at(min_index)); 
          _ds_v.push_back(min_dist);

        } 

      }

      // ...and delete it from the old vector
      _s_hit_v.erase(_s_hit_v.begin() + min_index);

    }

    // For the last hit, use the last values of min_dist
    _ds_v.push_back(min_dist);

    // Now that the vector is ordered, reassing to the original one
    _s_hit_v = new_vector;

    //if (_debug) {
      //for (auto h : _s_hit_v) {
        //std::cout << "Ordered hits: " << h.wire << ", " << h.time*4 << std::endl;
      //}
    //}

    _hits_ordered = true;

    return _s_hit_v.size();
  }

  void StoppingMuonTaggerHelper::CalculatedQds() {

    if (!_hits_ordered) {
      std::cout << "Call to " << __FUNCTION__ << " without having ordered hits." << std::endl;
      throw std::exception();
    }
 
    if (_ds_v.size() != _s_hit_v.size()) {
      std::cout << __FUNCTION__ << ": ds size is different than hit vector size" << std::endl;
      throw std::exception();
    }
 
    _dqds_v.clear();
    _dqds_v.reserve(_s_hit_v.size());

    double ds = 1.;

    for (size_t i = 0; i < _s_hit_v.size()-1; i++) {

      TVector3 this_point(_s_hit_v.at(i).wire * _w2cm, _s_hit_v.at(i).time*4*_t2cm, 0);
      TVector3 next_point(_s_hit_v.at(i+1).wire * _w2cm, _s_hit_v.at(i+1).time*4*_t2cm, 0);
      ds = (this_point - next_point).Mag();

      _dqds_v.emplace_back(_s_hit_v.at(i).integral/ds * _dqds_calib);
      if (_debug) std::cout << "_dqds_v.back() " << _dqds_v.back() << std::endl;

    }

    // Finish with the last point
    _dqds_v.emplace_back(_s_hit_v.at(_s_hit_v.size()-1).integral/ds * _dqds_calib);

    return;
  }



  bool StoppingMuonTaggerHelper::PerformdQdsSlider() {

    if (_dqds_v.size() != _s_hit_v.size()) {
      std::cout << __FUNCTION__ << ": dqds size is different than hit vector size" << std::endl;
      throw std::exception();
    }

    if (_dqds_v.size() < _slider_window * 2) {
      std::cout << __FUNCTION__ << ": not enough hits" << std::endl;
      return false;
    }

    //size_t window = _slider_window;

    if (_slider_window % 2 != 0) {
      std::cout << __FUNCTION__ << ": _slider_window has to be even." << std::endl; 
      throw std::exception();
    }

    _dqds_slider.clear();
    _dqds_slider.reserve(_dqds_v.size());


    for(const auto& window : this->get_windows(_dqds_v, _slider_window)) {

      double median_dqds = this->GetTruncMedian(window);
      _dqds_slider.push_back(median_dqds);
      //if (_debug) std::cout << "dqds_slider value " << median_dqds << std::endl;


    }

/*

    for (size_t i = window/2; i < _dqds_v.size() - window/2; i++) {

      std::vector<double> temp_v(_dqds_v.begin() + i - window/2, _dqds_v.begin() + i + window/2);
      double median_dqds = this->GetTruncMedian(temp_v);
      _dqds_slider.at(i) = median_dqds;
      if (_debug) std::cout << "dqds_slider value " << median_dqds << " at i " << i << std::endl;
    }

    // Now do the edges
    for (size_t i = 0; i < window/2; i++) {
      std::vector<double> temp_v(_dqds_v.begin(), _dqds_v.begin() + window); 
      double median_dqds = this->GetTruncMedian(temp_v);
      _dqds_slider.at(i) = median_dqds;
      if (_debug) std::cout << "dqds_slider value " << median_dqds << " at i " << i << std::endl;
    }
    for (size_t i = _dqds_v.size() - window/2; i < _dqds_v.size(); i++) {
      std::vector<double> temp_v(_dqds_v.begin() + _dqds_v.size() - window/2, _dqds_v.begin() + _dqds_v.size());
      double median_dqds = this->GetTruncMedian(temp_v); 
      _dqds_slider.at(i) = median_dqds;
      if (_debug) std::cout << "dqds_slider value " << median_dqds << " at i " << i << std::endl;
    }

*/
    return true;

  }

  double StoppingMuonTaggerHelper::GetTruncMedian(std::vector<double> v) {

    if (v.size() > 2) {
      // Find and erase max element
      auto it_max = std::max_element(v.begin(), v.end());
      v.erase(it_max);

      // Find and erase min element
      auto it_min = std::min_element(v.begin(), v.end());
      v.erase(it_min);
    }

    double median = -1;

    size_t size = v.size();
    std::sort(v.begin(), v.end());
    if (size % 2 == 0){
      median = (v[size/2 - 1] + v[size/2]) / 2;
    }
    else{
      median = v[size/2];
    }

    return median;
  }

  void StoppingMuonTaggerHelper::PrintOnFile(int index) {

    if (_s_hit_v.size() != _linearity_v.size() || _s_hit_v.size() != _dqds_slider.size()) {
      return;
    }

    for (size_t i = 0; i < _dqds_slider.size(); i++) {
      _csvfile << index << "," 
               << i << "," 
               << _dqds_v.at(i) << "," 
               << _dqds_slider.at(i) << ", "
               << _linearity_v.at(i)
               << std::endl;
    }

    if (_debug) {
      int counter = 0;
      for (auto h : _s_hit_v) {
        std::cout << "i " << counter
                  << ", wire: " << h.wire
                  << ", time: " << h.time*4
                  << ", dqdx_slide: " << _dqds_slider.at(counter)
                  << ", linearity: " << _linearity_v.at(counter) << std::endl;
        counter++;
      }
    }

  }


  void StoppingMuonTaggerHelper::CalculateLocalLinearity() {

    if (_s_hit_v.size() < 2) {
      std::cout << "Cannot calculate linearity if less than 2 hits" << std::endl;
      return;
    }

    std::vector<double> R;
    R.reserve(_s_hit_v.size());

    std::vector<double> X;
    std::vector<double> Y;
    X.reserve(_slider_window);
    Y.reserve(_slider_window);

    for(const auto& window : this->get_windows(_s_hit_v, 20/*_slider_window*/)) {

      for(const auto& s_hit : window) {
        X.push_back(s_hit.wire); 
        Y.push_back(s_hit.time);
      }

      auto c  = cov(X,Y);
      auto sX = stdev(X);
      auto sY = stdev(Y);
      auto r  = std::abs(c/(sX * sY));

      if(_debug) {
        //std::cout << "c: "  << c << std::endl
        //          << "sX: " << sX <<  std::endl
        //          << "sY: " << sY <<  std::endl
        //          << "r: "  << r <<  std::endl;
        //std::cout << "Local Linearity: " << r << std::endl;
      }

      if(std::isnan(r)) r = 0.0; 
      R.push_back(r);
      
      X.clear(); Y.clear();
    }   
 
    //first and last points will be nan. Lets set them equal to the points just above and below
    R.at(0)            = R.at(1);
    R.at(R.size() - 1) = R.at(R.size() - 2);
    
    _linearity_v = R;

    _linearity_is_set = true;

    return;
  } 

  double StoppingMuonTaggerHelper::cov (const std::vector<double>& data1,
                                        const std::vector<double>& data2) const
  {
    if(data1.size() == 0) std::cout << __FUNCTION__ << "You have me nill to cov" << std::endl;
    if(data2.size() == 0) std::cout << __FUNCTION__ << "You have me nill to cov" << std::endl;

    double result = 0.0;
    auto   mean1  = mean(data1);
    auto   mean2  = mean(data2);
    
    for(size_t i = 0; i < data1.size(); ++i)
      result += (data1[i] - mean1)*(data2[i] - mean2);
    
    return result/((double)data1.size());
      
  }

  double StoppingMuonTaggerHelper::stdev(const std::vector<double>& data) const
  {
    if(data.size() == 0) std::cout << __FUNCTION__ << "You have me nill to stdev" << std::endl;

    double result = 0.0;
    auto    avg   = mean(data);
    for(const auto& d: data)
      result += (d - avg)*(d - avg);
    
    return sqrt(result/((double)data.size()));
  }

  double StoppingMuonTaggerHelper::mean(const std::vector<double>& data) const
  {
    if(data.size() == 0) std::cout << __FUNCTION__ << "You have me nill to mean" << std::endl;
        
    double result = 0.0;

    for(const auto& d : data) 
      result += d;
        
    return (result / ((double)data.size()));
  }


  bool StoppingMuonTaggerHelper::MakeDecision(::ubana::StopMuAlgoType algo) {

   switch (algo) {

     case kAlgoMichel:
       return this->IsStopMuMichel();
       break;

     case kAlgoBragg:
       return this->IsStopMuBragg();
       break;

     case kAlgoSimpleMIP:
       return this->IsSimpleMIP();
       break;

     case kAlgoCurvature:
       return this->IsStopMuCurvature();
       break;

     default:
       std::cout << "Algo not found." << std::endl;
       throw std::exception();
    }

    return false;

  }


  bool StoppingMuonTaggerHelper::IsStopMuMichel() {

    if (_dqds_slider.size() < (unsigned int) (_hits_to_remove * 2 + _pre_post_window * 2 + 6)) {
      if (_debug) std::cout << "Can't make decision, number of simple hits is " << _dqds_slider.size() 
                            << ", which is less then " << _hits_to_remove * 2 + _pre_post_window * 2 + 6<< std::endl;
      return false;
    }


    // Vertex must not be in the FV
    //if (_fv.InFV(_vertex))
    //  return false;


    // Find the hits with the maximum dqds, that one will be the hit
    // where the Bragg peak is. If this is a big cluster, is likely that 
    // there'll be a big delta ray that may fake a Bragg peak. In this case, 
    // look only at the end of the cluster to find the Bragg peak.
    int offset = 1; // at least start from 1 as first hit is usually bad
    if (_dqds_slider.size() > (size_t)_max_end_hits) {
      if (_debug) std::cout << "[IsStopMuMichel] Many hits in this cluseter ("
                            << _dqds_slider.size() << "), finding Bragg in the last "
                            << _max_end_hits << " hits." << std::endl;
      offset = _dqds_slider.size() - _max_end_hits;
    }
    int bragg_index;
    auto it_max = std::max_element(_dqds_slider.begin()+offset, _dqds_slider.end());
    bragg_index = it_max - _dqds_slider.begin();
    double bragg_dqds = *it_max;
    for (bool flag = true; flag && it_max != _dqds_slider.end(); it_max++) {
      if (*it_max < bragg_dqds) {
        bragg_index = --it_max - _dqds_slider.begin();
        flag = false;
      }
    }

    if (_debug) std::cout << "[IsStopMuMichel] Bragg peak hit index is " << bragg_index << std::endl;

    // Check that the number of muon hits are below the maximum allowed
    int n_muon_hits = bragg_index + 1;
    if (n_muon_hits > _max_muon_hits) {
      if (_debug) std::cout << "[IsStopMuMichel] Number of muon hits is " << n_muon_hits
                            << " which is above maximum allowed (" << _max_muon_hits << ")" << std::endl;
      return false;
    }

    // Check that the number of muon hits are above the minimum allowed
    if (n_muon_hits < _min_muon_hits) {
      if (_debug) std::cout << "[IsStopMuMichel] Number of muon hits is " << n_muon_hits
                            << " which is below minimum allowed (" << _min_muon_hits << ")" << std::endl;
      return false;
    }

    // Check that the local linearity is less than threshold in the Bragg region
    double bragg_local_linearity = _linearity_v.at(bragg_index);

    if (bragg_local_linearity > _local_linearity_threshold) {
      if (_debug) std::cout << "[IsStopMuMichel] Local linearity is " << bragg_local_linearity 
                            << " which is above threshold (" << _local_linearity_threshold << ")" << std::endl;
      return false;
    }

    // Check that the local linearity in the muon region (before Bragg) 
    // is above threshold. Exclude first hits as things can get funny
    // at the beginning. Also exclude last hits as we expect the muon 
    // to curve in the last region
    /*
    for (size_t i = offset; i < (size_t) bragg_index - _pre_post_window; i++) {
      if (_linearity_v.at(i) < _local_linearity_threshold) {
        if (_debug) std::cout << "[IsStopMuMichel] Local linearity at hit " << i << " (before Bragg) is " << _linearity_v.at(i)
                              << " which is below threshold (" << _local_linearity_threshold << ")" << std::endl;
        return false;
      }
    }
    */

    // Check that the photon hits are below the maximum allowed
    int n_michel_hits = _dqds_slider.size() - bragg_index - 1;
    if (n_michel_hits > _max_michel_hits) {
      if (_debug) std::cout << "[IsStopMuMichel] Number of Michel hits is " << n_michel_hits
                            << " which is above maximum allowed (" << _max_michel_hits << ")" << std::endl;
      return false;
    }

    // Check that the photon hits are above the minimum allowed
    if (n_michel_hits < _min_michel_hits) {
      if (_debug) std::cout << "[IsStopMuMichel] Number of Michel hits is " << n_michel_hits
                            << " which is below the minimum allowed (" << _min_michel_hits << ")" << std::endl;
      return false;
    }

    // Get mean of first and last hits
    std::vector<double> temp = _dqds_slider;

    // Remove first "_hits_to_remove" and last "_hits_to_remove" hits
    /*if (offset > 0) {
      temp.erase(temp.begin(), temp.begin() + offset);
    } else {
      temp.erase(temp.begin(), temp.begin() + _hits_to_remove);
    }*/
    std::cout << "_dqds_slider vector has size " << _dqds_slider.size() << std::endl;
    std::cout << "temp vector has size " << temp.size() << std::endl;
    temp.erase(temp.begin(), temp.begin() + bragg_index - (_pre_post_window + 5));
    temp.erase(temp.end() - _hits_to_remove, temp.end());

    std::cout << "temp vector has size " << temp.size() << std::endl;

    bragg_index = (_pre_post_window + 5);

    if (temp.size() <= (size_t)bragg_index) {
      std::cout << "Not enough hits." << std::endl;
      return false;
    }
    std::cout << "bragg index is " << bragg_index << std::endl;
    std::cout << "at bragg index temp vector is " << temp.at(bragg_index) << std::endl;

    for (size_t i = 0; i < temp.size(); i++) std::cout << i << ": temp = " << temp.at(i) << std::endl;

    double start_mean = std::accumulate(temp.begin(), temp.begin() + _pre_post_window, 0);
    start_mean /= _pre_post_window;

    double end_mean = std::accumulate(temp.end() - _pre_post_window, temp.end(), 0);
    end_mean /= _pre_post_window;

    int edge = bragg_index + 5;

    if (temp.size() - edge < (size_t) _pre_post_window) {
      int vector_size = (int) temp.size();
      if (_debug) std::cout << "[IsStopMuMichel] Few Michel hits, calculating average only on "
                            << vector_size - edge << " hits." << std::endl;
      end_mean = std::accumulate(temp.end() - (vector_size - edge), temp.end(), 0);
      end_mean /= (double) vector_size - edge;
    }

    double perc_diff = (start_mean - end_mean) / start_mean * 100.;

    if (_debug) std::cout << "[IsStopMuMichel] Start mean: " << start_mean 
                          << ", end mean " << end_mean << ", Perc diff is " << perc_diff << std::endl;

    if (perc_diff > _perc_diff_cut) {
      return true;
    } 

    return false;   
  }




  bool StoppingMuonTaggerHelper::IsStopMuBragg() {

    if (_dqds_slider.size() < (unsigned int) (_pre_post_window * 2)) {
      if (_debug) std::cout << "Can't make decision, number of simple hits is " << _dqds_slider.size() 
                            << ", which is less then " << _pre_post_window * 2 << std::endl;
      return false;
    }

    // Vertex must not be in the FV
    if (_fv.InFV(_vertex))
      return false;

    // Find the hits with the maximum dqds, that one will be the hit
    // where the Bragg peak is
    size_t bragg_index;
    auto it_max = std::max_element(_dqds_slider.begin(), _dqds_slider.end());
    bragg_index = it_max - _dqds_slider.begin();
    double bragg_dqds = *it_max;
    for (bool flag = true; flag && it_max != _dqds_slider.end(); it_max++) {
      if (*it_max < bragg_dqds) {
        bragg_index = --it_max - _dqds_slider.begin();
        flag = false;
      }
    }

    if (_debug) std::cout << "[IsStopMuBragg] Bragg peak hit index is " << bragg_index << std::endl;

    // Check that the number of muon hits are below the maximum allowed
    int n_muon_hits = bragg_index + 1;
    if (n_muon_hits > _max_muon_hits) {
      if (_debug) std::cout << "[IsStopMuBragg] Number of muon hits is " << n_muon_hits
                            << " which is above maximum allowed (" << _max_muon_hits << ")" << std::endl;
      return false;
    }

    // Check that the number of muon hits are above the minimum allowed
    if (n_muon_hits < _min_muon_hits) {
      if (_debug) std::cout << "[IsStopMuBragg] Number of muon hits is " << n_muon_hits
                            << " which is below minimum allowed (" << _min_muon_hits << ")" << std::endl;
      return false;
    }

    // In this case we are looking for events that don't have a Michel, 
    // so we want to ensure that the local linearity is not below threshold
    // in the Bragg region
    /*
    double bragg_local_linearity = _linearity_v.at(bragg_index);

    if (bragg_local_linearity < _local_linearity_threshold) {
      if (_debug) std::cout << "[IsStopMuBragg] Local linearity is " << bragg_local_linearity
                            << " which is less than threshold (" << _local_linearity_threshold << ")" << std::endl;
      return false;
    }
    */

    // We actually want that there is no kink in this cluster,
    // as we just want the muon to stop. But exclude first hits 
    // as things can get funny at the beginning
    for (size_t l = _hits_to_remove; l < _linearity_v.size() - _hits_to_remove; l++) {
      if (_linearity_v.at(l) < _local_linearity_threshold) {
        if (_debug) std::cout << "[IsStopMuBragg] Local linearity at hit " << l << " is " << _linearity_v.at(l)
                              << " which is less than threshold (" << _local_linearity_threshold << ")" << std::endl;
        return false;
      }
    }


    // Take firsts hits, then lasts hits
    bool _use_mean = false; 
    double start_mean, end_mean;

    if (_use_mean) { 

      start_mean = std::accumulate(_dqds_slider.begin(), _dqds_slider.begin() + _pre_post_window, 0);
      start_mean /= _pre_post_window;
      end_mean = std::accumulate(_dqds_slider.end() - _pre_post_window, _dqds_slider.end(), 0);
      end_mean /= _pre_post_window;

    } else {

      std::vector<double>::const_iterator first = _dqds_slider.begin();
      std::vector<double>::const_iterator last  = _dqds_slider.begin() + _pre_post_window;
      std::vector<double> temp(first, last);
      std::sort(temp.begin(), temp.end());

      double _n_hits_remove = 7.;

      start_mean = std::accumulate(temp.begin(), temp.begin() + _n_hits_remove, 0);
      start_mean /= _n_hits_remove;

      first = _dqds_slider.end() - _pre_post_window;
      last  = _dqds_slider.end();
      std::vector<double> temp2(first, last);
      std::sort(temp2.begin(), temp2.end());

      end_mean = std::accumulate(temp2.begin(), temp2.begin() + _n_hits_remove, 0);
      end_mean /= _n_hits_remove;
    }

    double perc_diff = (start_mean - end_mean) / start_mean * 100.;

    if (_debug) std::cout << "[IsStopMuBragg] Start mean: " << start_mean 
                          << ", end mean " << end_mean << ", Perc diff is " << perc_diff << std::endl;

    if (perc_diff < -_perc_diff_cut) {
      return true;
    }

    return false;

  }

  bool StoppingMuonTaggerHelper::IsSimpleMIP() {

    // This algo excludes the first and last hit

    // Check that the local linearity in never below threshold
    for (size_t i = 2; i < _linearity_v.size() - 2; i++) {

      if (_linearity_v.at(i) < _local_linearity_threshold) {

        std::cout << "[IsSimpleMIP] Local linearity at hit " << i 
                  << " is " << _linearity_v.at(i) 
                  << " which is below threshold (" 
                  << _local_linearity_threshold << ")." << std::endl;
        return false;

      }
    }

    // Check compatibility
    if (_dqds_v.size() != _dqds_slider.size()) {
      std::cout << "[IsSimpleMIP] _dqds_v vector size is " << _dqds_v.size() 
                << " which is different to _dqds_slider vector size, which is "
                << _dqds_slider.size() << "." << std::endl;
      throw std::exception();
    }

    // Check that we have enough hits
    if (_dqds_v.size() < 6 + 1 + 6 + 1) {
      std::cout << "[IsSimpleMIP] Not enough hits." << std::endl;
      return false;
    }

    // Now verify the first and last hits are flat in dqds
    std::vector<double> start_v (_dqds_v.begin() + 1, _dqds_v.begin() + 6);
    std::vector<double> end_v (_dqds_v.end() - 6, _dqds_v.end() - 1);

    //double std_start = stdev(start_v);
    //double std_end = stdev(end_v);
    //std::cout << "std_start: " << std_start << ", std_end: " << std_end << std::endl;

    bool good_start = true;
    bool good_end = true;

    for (auto q : start_v) {
      if (q < 40000 || q > 75000) {
        good_start = false;
      }
    }

    for (auto q : end_v) { 
      if (q < 40000 || q > 75000) {
        good_end = false;
      }
    }

    std::cout << "[IsSimpleMIP] Start is " << (good_start ? "GOOD" : "BAD") << std::endl;
    std::cout << "[IsSimpleMIP] End is " << (good_end ? "GOOD" : "BAD") << std::endl;


    // Now use the smoothed dqds to evaluate the start and 
    // end dqds average value
    double start_mean = 0, end_mean = 0;

    start_mean = std::accumulate(_dqds_slider.begin() + 1, 
                                 _dqds_slider.begin() + 6, 0);
    start_mean /= 5.;

    end_mean = std::accumulate(_dqds_slider.end() - 6,
                               _dqds_slider.end() - 1, 0);
    end_mean /= 5.;

    double perc_diff = (start_mean - end_mean) / start_mean * 100.; 

    if (_debug) std::cout << "[IsSimpleMIP Start mean: " << start_mean
                          << ", end mean " << end_mean << ", Perc diff is " << perc_diff << std::endl; 

    if (good_start && good_end && std::abs(perc_diff) < 10) {
      return true;
    }

    return false;

  }

  /*
  double StoppingMuonTaggerHelper::GetSTD(std::vector<double> v) {

    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();

    double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / v.size() - mean * mean);
  }
  */

  bool StoppingMuonTaggerHelper::IsStopMuCurvature() {

    if (_s_hit_v.size() < (unsigned int) (_pre_post_window * 2)) {
      if (_debug) std::cout << "Can't make decision, number of simple hits is " << _s_hit_v.size()
                            << ", which is less then " << _pre_post_window * 2 << std::endl;
      return false;
    }

    // Find the hits with the maximum dqds, that one will be the hit
    // where the Bragg peak is
    size_t bragg_index;
    auto it_max = std::max_element(_dqds_slider.begin(), _dqds_slider.end());
    bragg_index = it_max - _dqds_slider.begin();
    double bragg_dqds = *it_max;
    for (bool flag = true; flag && it_max != _dqds_slider.end(); it_max++) {
      if (*it_max < bragg_dqds) {
        bragg_index = --it_max - _dqds_slider.begin();
        flag = false;
      }
    }

    if (_debug) std::cout << "[IsStopMuCurvature] Bragg peak hit index is " << bragg_index << std::endl;

    if (bragg_index < (size_t) _pre_post_window * 2) {
      std::cout << "[IsStopMuCurvature] Not enough points"  << std::endl;
      return false;
    }

    // In this case we are looking for events that don't have a Michel, 
    // so we want to ensure that the local linearity is not below threshold
    // in the Bragg region
    double bragg_local_linearity = _linearity_v.at(bragg_index);

    if (bragg_local_linearity > _local_linearity_threshold) {
      if (_debug) std::cout << "[IsStopMuBragg] Local linearity is " << bragg_local_linearity
                            << " which is above threshold (" << _local_linearity_threshold << ")" << std::endl;
      return false;
    }

    /* Calulate curvature first up to the bragg hit
    // http://paulbourke.net/geometry/circlesphere/

    TVector3 pt0(0., 0., 0.);
    TVector3 pt1(0., 0., 0.);
    TVector3 pt2(0., 0., 0.);

    for (size_t i = 1; i < bragg_index - _pre_post_window - 1; i++) {
    
      pt0.SetX(_s_hit_v.at(i-1).wire);
      pt0.SetY(_s_hit_v.at(i-1).time);
 
      pt1.SetX(_s_hit_v.at(i).wire);
      pt1.SetY(_s_hit_v.at(i).time);

      pt2.SetX(_s_hit_v.at(i+1).wire);
      pt2.SetY(_s_hit_v.at(i+1).time);
     
      double curvature = this->GetCurvature(pt0, pt1, pt2);

      std::cout << "Curvature is " << curvature << std::endl;
    }
    */

    // Look at linearity before the bragg peak

    std::vector<double> lin;
    lin.clear();
    lin.insert(lin.end(), _linearity_v.begin() + bragg_index - 3 * _pre_post_window, 
                          _linearity_v.begin() + bragg_index - 3);

    double lin_mean = std::accumulate(lin.begin(), lin.end(), 0);
    lin_mean /= (double) lin.size();

    if (_debug) std::cout << "[IsStopMuBragg] Average of local linearity before the Bragg peak is " 
                          << lin_mean << std::endl;

    if (lin_mean < 0.75/*_linearity_curvature_threshold*/)
      return true;


    return false;
  }


  double StoppingMuonTaggerHelper::GetCurvature(TVector3 pt1, TVector3 pt2, TVector3 pt3) {

    if (pt1.Z() != 0 || pt2.Z() != 0 || pt3.Z() != 0) {
      std::cout << __PRETTY_FUNCTION__ << "Third component of input points has to be zero." << std::endl;
      throw std::exception();
    }

    double ma = ( pt2.Y() - pt1.Y() ) / ( pt2.X() - pt1.X() );
    double mb = ( pt3.Y() - pt2.Y() ) / ( pt3.X() - pt2.X() );

    double x = ( ma * mb * (pt1.Y() - pt3.Y()) + mb * (pt1.X() + pt2.X()) - ma * (pt2.X() + pt3.X()) ) / ( 2 * (mb - ma) );

    double y = (-1/ma) * ( x - (pt1.X() + pt2.X())/2 ) + (pt1.Y() + pt2.Y())/2;

    TVector3 center(x, y, 0.);

    double radius = (pt1 - center).Mag();

    if (radius == 0)
      return 1e9;

    return 1./radius;

  } 


  template<typename T>
  std::vector<std::vector<T>> StoppingMuonTaggerHelper::get_windows(const std::vector<T>& the_thing,
                                                                    const size_t window_size) const
  {

    // given a vector of values return a vector of the same length
    // with each element being a vector of the values of the local neighbors
    // of the element at position i in the original vector
    // input  : [0,1,2,3,4,5,6,...,...,N-3,N-2,N-1] (input vector of size N)
    // output  (assuming a value of 'w' below == 3):
    // 0th element: [0]
    // 1st element: [0,1,2]
    // 2nd element: [0,1,2,3,4]
    // jth element: [j-w,j-w+1,..,j+w-2,j+w-1]
    
    std::vector<std::vector<T>> data;
    
    auto w = window_size + 2;
    w = (unsigned int)((w - 1)/2);
    auto num = the_thing.size();
    
    data.reserve(num);
    
    for(size_t i = 1; i <= num; ++i) {
      std::vector<T> inner;
      inner.reserve(20);
      // if we are at the beginning of the vector (and risk accessing -1 elements)
      if(i < w)
        {
          for(size_t j = 0; j < 2 * (i%w) - 1; ++j)
            inner.push_back(the_thing[j]);
        }
      // if we are at the end of the vector (and risk going past it)
      else if (i > num - w + 1)
        {
          for(size_t j = num - 2*((num - i)%w)-1 ; j < num; ++j)
            inner.push_back(the_thing[j]);
        }
      // if we are in the middle of the waveform
      else
        {
          for(size_t j = i - w; j < i + w - 1; ++j)
            inner.push_back(the_thing[j]);
        }
      data.emplace_back(inner);
    }

    return data;
  
  }

}


#endif
