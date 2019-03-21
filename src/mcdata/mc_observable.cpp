/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-05-06 11:31:00
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-16 09:32:43
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./mc_observable.h"
#include <locale>
//#include <boost/algorithm/string.hpp>

namespace mcdata {

MC_Observable::MC_Observable() 
  : MC_Data() 
{
  num_dataset_=0;
  avg_stddev_.setZero();
  avg_tau_ = 0.0;
}

MC_Observable::MC_Observable(const std::string& name, const unsigned& size,
  const bool& replace_mode) 
{
  this->init(name,size);
  replace_mode_ = replace_mode;
} 

void MC_Observable::init(const std::string& name, const unsigned& size)
{
  MC_Data::init(name, size);
  name_ = name;
  num_dataset_=0;
  avg_mcdata_.init(name, size);
  avg_stddev_.setZero();
  avg_tau_ = 0.0;
  //elem_names_ = {name};
  elem_names_.resize(size);
  // file name
  fname_ = name_;
  // lowercase
  std::locale loc;
  for (int i=0; i<fname_.size(); ++i) fname_[i]=std::tolower(fname_[i],loc);
  auto pos = fname_.find('^');
  if (pos != std::string::npos) fname_.erase(pos,1);
  fname_ = "res_"+fname_+".txt";
}

void MC_Observable::resize(const unsigned& size)
{
  MC_Data::resize(size);
  avg_mcdata_.resize(size);
  elem_names_.resize(size);
  for (unsigned i=0; i<size; ++i) 
    elem_names_[i] = "elem" + std::to_string(i);
}

void MC_Observable::resize(const unsigned& size, const std::vector<std::string>& elem_names)
{
  MC_Data::resize(size);
  avg_mcdata_.resize(size);
  elem_names_ = elem_names;
}

/*void MC_Observable::check_on(const input::Parameters& inputs, const bool& replace_mode) 
{
  int no_warn;
  is_on_ = inputs.set_value(name(), false, no_warn);
  replace_mode_ = replace_mode;
}*/

void MC_Observable::print_heading(const std::string& header,
  const std::vector<std::string>& xvars) 
{
  if (!is_on()) return;
  if (heading_printed_) return;
  if (!replace_mode_) return;
  if (!fs_.is_open()) open_file();
  fs_ << header;
  fs_ << "# Results: " << name() << "\n";
  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << "# ";
  fs_ << std::left;
  //fs_ << std::setw(14)<<xvar_name;
  for (const auto& p : xvars) fs_ << std::setw(14)<<p.substr(0,14);
  // total value
  if (MC_Data::size()>1 && have_total_)
    fs_ << std::setw(14)<<"Total"<<std::setw(11)<<"err";
  for (const auto& name : elem_names_) 
    fs_ << std::setw(14)<<name<<std::setw(11)<<"err";
  //fs_ << std::setw(9)<<"samples";
  fs_ << std::setw(9)<<"samples"<<std::setw(12)<<"converged"<<std::setw(6)<<"tau";
  fs_ << std::endl;
  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << std::flush;
  heading_printed_ = true;
  close_file();
}

void MC_Observable::print_result(const std::vector<double>& xpvals) 
{
  if (!is_on()) return;
  if (!fs_.is_open()) open_file();
  fs_ << std::right;
  fs_ << std::scientific << std::uppercase << std::setprecision(6);
  for (const auto& p : xpvals) 
    fs_ << std::setw(14) << p;
  // total value
  //if (MC_Data::size()>1)
  if (MC_Data::size()>1 && have_total_)
    fs_ << MC_Data::result_str(-1); 
  for (unsigned i=0; i<MC_Data::size(); ++i) 
    fs_ << MC_Data::result_str(i); 
  fs_ << MC_Data::conv_str(0); //.substr(0,10); 
  fs_ << std::endl;
  fs_ << std::flush;
  close_file();
} 

void MC_Observable::open_file(void) 
{
  if (fs_.is_open()) return;
  if (replace_mode_) {
    fs_.open(fname_);
    replace_mode_ = false;
  }
  else fs_.open(fname_, std::ios::app);
  if (!fs_.is_open()) 
    throw std::runtime_error("Observable::open_file: file open failed");
}

void MC_Observable::close_file(void) 
{
  fs_.close();
}

void MC_Observable::save_result(void)
{
  avg_mcdata_ << MC_Data::mean_data();
  avg_stddev_ += MC_Data::stddev_data();
  avg_tau_ += MC_Data::tau();
  ++num_dataset_;
}


} // end namespace mcdata


