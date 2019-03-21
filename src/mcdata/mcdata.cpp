/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-24 08:54:44
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-23 00:06:47
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./mcdata.h"

namespace mcdata {

DataBin::DataBin(const unsigned& size) 
{
  resize(size);
  clear();
}

void DataBin::clear(void) 
{
  num_samples_ = 0;
  num_samples_last_ = 0;
  ssum_.setZero();
  sumsq_.setZero();
  carry_.setZero();
  waiting_sample_.setZero();
  waiting_sample_exist_ = false;
  MinusOne_.setConstant(-1.0);
  mean_.setZero();
  stddev_.setZero();
}

void DataBin::resize(const unsigned& size) 
{
  size_ = size;
  ssum_.resize(size);
  sumsq_.resize(size);
  carry_.resize(size);
  waiting_sample_.resize(size);
  MinusOne_.resize(size);
  mean_.resize(size);
  stddev_.resize(size);
  this->clear();
}

bool DataBin::add_sample(const data_t& new_sample) 
{
  num_samples_++;
  ssum_ += new_sample;
  sumsq_ += new_sample * new_sample;
  if (waiting_sample_exist_) {
    carry_ = (waiting_sample_ + new_sample)*0.5;
    waiting_sample_exist_ = false;
  }
  else {
    waiting_sample_ = new_sample;
    waiting_sample_exist_ = true;
  }
  return !waiting_sample_exist_;
}

void DataBin::finalize(void) const
{
  if (num_samples_last_ != num_samples_) {
    num_samples_last_ = num_samples_;
    if (num_samples_ > 1) {
      mean_ = ssum_/num_samples_;
      data_t variance_ = sumsq_/num_samples_ - mean_ * mean_;
      stddev_ = variance_/(num_samples_-1);
      stddev_ = stddev_.sqrt();
    }
    else {
      mean_ = ssum_;
      stddev_ = MinusOne_;  // to mean infinity
    }
  }
}

/*----------------------MC_Data class------------------*/
void MC_Data::init(const std::string& name, const unsigned& size) 
{
  std::vector<DataBin>::clear();
  for (unsigned i=0; i<max_binlevel_default_; ++i) 
    this->push_back(DataBin(size));
  top_bin = this->begin();
  end_bin = this->end();
  name_ = name;
  this->clear();
}

void MC_Data::resize(const unsigned& size) 
{
  for (auto& bin : *this) bin.resize(size);
  this->clear();
}

void MC_Data::clear(void) 
{
  for (auto& bin : *this) bin.clear();
  dcorr_level_ = 0;
  mean_ = top_bin->mean();
  stddev_ = top_bin->stddev();
  tau_ = -1.0;
  show_statistic_ = false;
  error_converged_ = "NOT_CONVD";
  convergence_str_ = "NULL";
}

void MC_Data::add_sample(const data_t& sample)
{
  auto this_bin = top_bin;  
  data_t new_sample(sample);
  while (this_bin->add_sample(new_sample)) {
    new_sample = this_bin->carry();
    //if (this_bin++ == end_bin) break; // wrong logic?
    if (++this_bin == end_bin) break;
  }
}

void MC_Data::add_sample(const double& sample)
{
  data_t new_sample(1);
  new_sample << sample;
  add_sample(new_sample);
}

void MC_Data::operator<<(const data_t& sample) {
  add_sample(sample);
}

void MC_Data::operator<<(const double& sample)
{
  data_t new_sample(1);
  new_sample << sample;
  add_sample(new_sample);
}

const data_t& MC_Data::mean_data(void) const 
{
  this->finalize(); return mean_;
}

double MC_Data::mean(void) const 
{ 
  this->finalize(); return mean_.sum();
}

const double& MC_Data::mean(const int& n) const 
{ 
  this->finalize(); return mean_(n);
}

const data_t& MC_Data::stddev_data(void) const 
{
  this->finalize(); return stddev_;
}

double MC_Data::stddev(void) const 
{ 
  this->finalize(); return stddev_.sum();
} 

const double& MC_Data::stddev(const int& n) const 
{ 
  this->finalize(); return stddev_(n);
}

const double& MC_Data::tau(void) const 
{
  this->find_conv_and_tau(0);
  return tau_;
}

void MC_Data::finalize(void)  const
{ 
  if (top_bin->have_new_samples()) {
    // mean
    mean_ = top_bin->mean();
    // stddev is the 'bin stddev' of the lowest depth bin with a 'good sample size'
    stddev_ = top_bin->stddev();
    dcorr_level_ = 0;
    for (auto bin=this->begin(); bin!=this->end(); ++bin) {
      if (bin->num_samples() < good_sample_size_) break;
      stddev_ = bin->stddev();
      dcorr_level_++;
    }
  }
}

void MC_Data::find_conv_and_tau(const unsigned& n) const 
{ 
  this->finalize();
  std::vector<double> xv, yv;
  unsigned level_n = dcorr_level_;
  if (dcorr_level_ == 2) level_n = dcorr_level_-2;
  else if (dcorr_level_ >= 3) level_n = dcorr_level_-3;
  for (unsigned i=level_n; i<=dcorr_level_; ++i) {
    xv.push_back(static_cast<double>(i));
    data_t stddev = this->operator[](i).stddev();
    yv.push_back(static_cast<double>(stddev(n)));
  }
  // assuming convergence
  double stddev_0 = static_cast<double>(this->operator[](0).stddev()(n));
  double stddev_d = static_cast<double>(this->operator[](dcorr_level_).stddev()(n));
  if (stddev_0 < 1.0E-12) tau_ = -1.0;
  else {
    double r = stddev_d/stddev_0;
    tau_ = 0.5*( r*r - 1.0);
    tau_ = std::abs(tau_);
  }
  // 'tau' gets reset in case of non-convergence
  this->check_convergence(xv, yv);
  if (error_converged_ == "NOT_CONVD") tau_ = -1.0;
}

void MC_Data::check_convergence(const std::vector<double>& xv, const std::vector<double>& yv) const 
{
  error_converged_ = "NOT_CONVD";
  if (xv.size() >= 3) {
    // slope of a least square fit straight line through these points
    double a1 = static_cast<double>(xv.size());
    double b1 = 0.0;
    for (const auto& x : xv) b1 += x;
    double c1 = 0.0;
    for (const auto& y : yv) c1 += y;
    double a2 = b1;
    double b2 = 0.0;
    for (const auto& x : xv) b2 += x*x;
    double c2 = 0.0;
    for (unsigned i=0; i<xv.size(); ++i) c2 += xv[i]*yv[i];
    double slope = (c2 - c1*a2/a1)/(b2 - b1*a2/a1);
    //if (std::abs(slope) < 1.0E-4) {
    if (std::abs(slope) < 0.1 * yv.back()) error_converged_ = "CONVERGED";
  }
} 

std::string MC_Data::result_str(const int& n) const 
{ 
  std::ostringstream os;
  double mean, stddev;
  if (n<0) {
    mean = this->mean();
    stddev = this->stddev();
  }
  else {
    mean = this->mean(n);
    stddev = this->stddev(n);
  }
  os << std::scientific << std::uppercase << std::setprecision(6) << std::right;
  os << std::setw(15) << mean << std::fixed << std::setw(10) << stddev;
  return os.str();
}

std::string MC_Data::conv_str(const int& n) const 
{ 
  this->find_conv_and_tau(n);
  std::ostringstream os;
  os << std::scientific << std::uppercase << std::setprecision(6) << std::right;
  os << std::setw(10) << this->num_samples();
  os << std::setw(11) << error_converged_;
  os << std::setw(7) << std::setprecision(2) << std::right << std::fixed << tau_;
  os << std::resetiosflags(std::ios_base::floatfield) << std::nouppercase; 
  return os.str();
}

void MC_Data::show_statistic(std::ostream& os) const
{
  auto bin = top_bin;
  os << name_ << " Statistic:\n";
  unsigned i = 0;
  while (bin->num_samples()>1) {
    os <<"bin-"<<std::setw(2)<<i++<<": "<< std::right<<std::setw(8)<<bin->num_samples(); 
    data_t mean = bin->mean();
    data_t stddev = bin->stddev();
    for (int n=0; n<mean.rows(); ++n) {
      os << std::scientific << std::uppercase << std::setw(14) << mean(n); 
      os << std::setw(6) << " (+/-) " << std::fixed << std::setw(10) << stddev(n); 
    }
    os << "\n";
    bin++;
  }
  os << std::right << std::resetiosflags(std::ios_base::floatfield) << std::nouppercase; 
}

std::ostream& operator<<(std::ostream& os, const MC_Data& obs)
{
  using namespace std;
  streamsize dp = cout.precision(); 
  os << scientific << uppercase << setprecision(6);
  os << "#" << std::string(36, '-') << "\n";
  os << "# " << obs.name() << ": " << "(samples = " << obs.num_samples() << ")\n";
  os << left << setw(8) << "# i" << setw(16) << "mean" << setw(12) << "err" << "\n";
  os << "#" << std::string(36, '-') << "\n";
  auto mean = obs.mean_data();
  auto stddev = obs.stddev_data();
  for (int i=0; i<mean.rows(); ++i) {
    os << "  ";
    os << scientific << setw(6) << i << setw(16) << mean(i); // << setw(6) << "(+/-)";
    os << fixed << setw(12) << stddev(i) << "\n";
  }
  if (obs.show_statistic_) {
    os << "\n";
    obs.show_statistic(os);
    obs.show_statistic_ = false;
  }
  os << right << resetiosflags(ios_base::floatfield) << nouppercase << setprecision(dp);
  return os;
}


} // end namespace mcdata
