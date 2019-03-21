/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-03-20 11:50:30
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-03-21 11:38:00
*----------------------------------------------------------------------------*/
// File: sysconfig.cpp
#include <iomanip>
#include "sysconfig.h"

void SysConfig::init(const lattice_id& lid, const lattice_size& size, const wf_id& wid)
{
  lattice_.construct(lid,size);
  // one body part of the wavefunction
  num_sites_ = lattice_.num_sites();
  basis_state_.init(num_sites_);
  hole_doping_ = 0.0;
  wf_.init(wid, lattice_, hole_doping_);
  num_upspins_ = wf_.num_upspins();
  num_dnspins_ = wf_.num_dnspins();
  basis_state_.init_spins(num_upspins_,num_dnspins_);

  // variational parameters
  num_wf_params_ = wf_.num_vparams();
  num_total_vparams_ = num_wf_params_;
  vparams_.resize(num_total_vparams_);

  // work arrays
  psi_mat_.resize(num_upspins_,num_dnspins_);
  psi_inv_.resize(num_upspins_,num_dnspins_);
  psi_row_.resize(num_dnspins_);
  psi_col_.resize(num_upspins_);
  inv_row_.resize(num_upspins_);
}

int SysConfig::build(const RealVector& vparams)
{
  wf_.compute(lattice_, vparams, 0);
  return 0;
}

int SysConfig::init_state(void)
{
  // try for a well condictioned amplitude matrix
  basis_state_.set_random();
  int num_attempt = 0;
  while (true) {
    wf_.get_amplitudes(psi_mat_,basis_state_.upspin_sites(), basis_state_.dnspin_sites());
    // reciprocal conditioning number
    Eigen::JacobiSVD<ComplexMatrix> svd(psi_mat_);
    // reciprocal cond. num = smallest eigenval/largest eigen val
    double rcond = svd.singularValues()(svd.singularValues().size()-1)/svd.singularValues()(0);
    if (std::isnan(rcond)) rcond = 0.0; 
    if (rcond>1.0E-15) break;
    //std::cout << "rcondition number = "<< rcond << "\n";
    // try new basis state
    basis_state_.set_random();
    if (++num_attempt > 1000) {
      throw std::underflow_error("*SysConfig::init: configuration wave function ill conditioned.");
    }
  }
  psi_inv_ = psi_mat_.inverse();
  // reset run parameters
  num_updates_ = 0;
  refresh_cycle_ = 100;
  num_proposed_moves_ = 0;
  num_accepted_moves_ = 0;
  return 0;
}

int SysConfig::update_state(void)
{
  for (int n=0; n<num_upspins_; ++n) do_upspin_hop();
  for (int n=0; n<num_dnspins_; ++n) do_dnspin_hop();
  //for (int n=0; n<num_exchange_moves_; ++n) do_spin_exchange();
  num_updates_++;
  if (num_updates_ % refresh_cycle_ == 0) {
    psi_inv_ = psi_mat_.inverse();
  }
  //std::cout << basis_state_ << "\n"; getchar();
  return 0;
}

int SysConfig::do_upspin_hop(void)
{
  if (basis_state_.gen_upspin_hop()) {
    int upspin = basis_state_.which_upspin();
    int to_site = basis_state_.which_site();
    wf_.get_amplitudes(psi_row_, to_site, basis_state_.dnspin_sites());
    amplitude_t det_ratio = psi_row_.cwiseProduct(psi_inv_.col(upspin)).sum();
    if (std::abs(det_ratio) < 1.0E-12) {
      // for safety
      basis_state_.undo_last_move();
      return 0; 
    } 
    auto weight_ratio = det_ratio;
    double transition_proby = std::norm(weight_ratio);
    num_proposed_moves_++;
    if (basis_state_.rng().random_real()<transition_proby) {
      num_accepted_moves_++;
      // upddate state
      basis_state_.commit_last_move();
      // update amplitudes
      inv_update_upspin(upspin,psi_row_,det_ratio);
    }
    else {
      basis_state_.undo_last_move();
    }
  } 
  return 0;
}

int SysConfig::do_dnspin_hop(void)
{
  if (basis_state_.gen_dnspin_hop()) {
    int dnspin = basis_state_.which_dnspin();
    int to_site = basis_state_.which_site();
    wf_.get_amplitudes(psi_col_, basis_state_.upspin_sites(), to_site);
    amplitude_t det_ratio = psi_col_.cwiseProduct(psi_inv_.row(dnspin)).sum();
    if (std::abs(det_ratio) < 1.0E-12) { // for safety
      basis_state_.undo_last_move();
      return 0; 
    } 
    amplitude_t weight_ratio = det_ratio;
    double transition_proby = std::norm(weight_ratio);
    num_proposed_moves_++;
    if (basis_state_.rng().random_real()<transition_proby) {
      num_accepted_moves_++;
      // upddate state
      basis_state_.commit_last_move();
      // update amplitudes
      inv_update_dnspin(dnspin,psi_col_,det_ratio);
    }
    else {
      basis_state_.undo_last_move();
    }
  } 
  return 0;
}

int SysConfig::inv_update_upspin(const int& upspin, const ColVector& psi_row, 
  const amplitude_t& det_ratio)
{
  psi_mat_.row(upspin) = psi_row;
  amplitude_t ratio_inv = amplitude_t(1.0)/det_ratio;
  for (int i=0; i<upspin; ++i) {
    amplitude_t beta = ratio_inv*psi_row.cwiseProduct(psi_inv_.col(i)).sum();
    psi_inv_.col(i) -= beta * psi_inv_.col(upspin);
  }
  for (int i=upspin+1; i<num_upspins_; ++i) {
    amplitude_t beta = ratio_inv*psi_row.cwiseProduct(psi_inv_.col(i)).sum();
    psi_inv_.col(i) -= beta * psi_inv_.col(upspin);
  }
  psi_inv_.col(upspin) *= ratio_inv;
  return 0;
}

int SysConfig::inv_update_dnspin(const int& dnspin, const RowVector& psi_col, 
  const amplitude_t& det_ratio)
{
  psi_mat_.col(dnspin) = psi_col;
  amplitude_t ratio_inv = amplitude_t(1.0)/det_ratio;
  for (int i=0; i<dnspin; ++i) {
    amplitude_t beta = ratio_inv*psi_col_.cwiseProduct(psi_inv_.row(i)).sum();
    psi_inv_.row(i) -= beta * psi_inv_.row(dnspin);
  }
  for (int i=dnspin+1; i<num_dnspins_; ++i) {
    amplitude_t beta = ratio_inv*psi_col_.cwiseProduct(psi_inv_.row(i)).sum();
    psi_inv_.row(i) -= beta * psi_inv_.row(dnspin);
  }
  psi_inv_.row(dnspin) *= ratio_inv;
  return 0;
}

void SysConfig::print_stats(std::ostream& os) const
{
  std::streamsize dp = std::cout.precision(); 
  double accept_ratio = 100.0*double(num_accepted_moves_)/(num_proposed_moves_);
  os << "--------------------------------------\n";
  os << " total mcsteps = " << num_updates_ <<"\n";
  os << std::fixed << std::showpoint << std::setprecision(1);
  os << " acceptance ratio = " << accept_ratio << " %\n";
  os << "--------------------------------------\n";
  // restore defaults
  os << std::resetiosflags(std::ios_base::floatfield) << std::setprecision(dp);
}

double SysConfig::get_energy(void) const
{
  // hopping energy
  //double bond_sum = 0.0;
  //for (int i=0; i<lattice.num_bonds(); ++i) {
  //  int src = lattice.bond(i).src();
  //  int tgt = lattice.bond(i).tgt();
  //}

  return 1.0;
}

