/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-24 08:44:27
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-07 22:25:28
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef MCDATA_H
#define MCDATA_H

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <stdexcept>
#include <Eigen/Core>

namespace mcdata {

using data_t = Eigen::Array<double,Eigen::Dynamic,1>;
using scalardata_t = Eigen::Array<double,1,1>;

/*----------------------DataBin class------------------*/
class DataBin 
{
public:
  DataBin(const unsigned& size=1);
  ~DataBin() {}
  void clear(void);
  void resize(const unsigned& size);
  bool add_sample(const data_t& sample);
  bool has_samples(void) const { return (num_samples_ > 0); }
  bool has_carry_over(void) const { return !waiting_sample_exist_; }
  const unsigned& num_samples(void) const { return num_samples_; }
  const data_t& carry(void) const { return carry_; } 
  bool have_new_samples(void) const { return num_samples_!=num_samples_last_; }
  void finalize(void) const;
  const data_t& mean(void) const { finalize(); return mean_; }
  const data_t& stddev(void) const { finalize(); return stddev_; }
private:
  unsigned size_;
  unsigned num_samples_{0};
  mutable unsigned num_samples_last_{0};
  data_t ssum_;
  data_t sumsq_;
  data_t carry_;
  data_t waiting_sample_;
  bool waiting_sample_exist_{false};
  data_t MinusOne_;
  mutable data_t mean_;
  mutable data_t stddev_;
};

/*----------------------mcdata class------------------*/
class MC_Data : private std::vector<DataBin>
{
public:
  MC_Data() {}
  MC_Data(const std::string& name, const unsigned& size=1) { init(name,size); }
  ~MC_Data() {}
  virtual void init(const std::string& name, const unsigned& size=1);
  virtual void resize(const unsigned& size);
  void clear(void);
  void add_sample(const data_t& sample);
  void add_sample(const double& sample);
  void operator<<(const data_t& sample);
  void operator<<(const double& sample);
  const unsigned& num_samples(void) const { return top_bin->num_samples(); }
  void finalize(void) const;
  const std::string& name(void) const { return name_; }
  unsigned size(void) const { return mean_.size(); }
  const data_t& mean_data(void) const; 
  double mean(void) const; 
  const double& mean(const int& n) const; 
  const data_t& stddev_data(void) const; 
  double stddev(void) const; // { return top_bin->stddev(); } 
  //const data_t& tau(void) const { finalize(); return tau_; } 
  const double& stddev(const int& n) const;
  const double& tau(void) const;
  std::string result_str(const int& n=0) const; 
  std::string conv_str(const int& n=0) const; 
  const MC_Data& with_statistic(void) const { show_statistic_=true; return *this; }
  void show_statistic(std::ostream& os=std::cout) const;
  friend std::ostream& operator<<(std::ostream& os, const MC_Data& obs);
private:
  std::string name_;
  static const unsigned max_binlevel_default_ = 20;
  static const unsigned good_sample_size_ = 30;
  std::vector<DataBin>::iterator top_bin;
  std::vector<DataBin>::iterator end_bin;
  mutable unsigned dcorr_level_;
  mutable data_t mean_;
  mutable data_t stddev_;
  mutable double tau_;
  mutable bool show_statistic_;
  mutable std::string error_converged_;
  mutable std::string convergence_str_;

  void find_conv_and_tau(const unsigned& n=0) const;
  void check_convergence(const std::vector<double>& xv, const std::vector<double>& yv) const;
};



} // end namespace mcdata

#endif
