/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-03-20 11:50:30
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-03-20 12:46:54
*----------------------------------------------------------------------------*/
// File: sysconfig.cpp
#include "sysconfig.h"

void SysConfig::init(const Lattice& lattice)
{
	// one body part of the wavefunction
	num_sites_ = lattice.num_sites();
	basis_state_.init(num_sites_);
	hole_doping_ = 0.0;
	wf_.init(wf_id::BCS, lattice, hole_doping_);
	num_upspins_ = wf_.num_upspins();
	num_dnspins_ = wf_.num_dnspins();
	psi_mat_.resize(num_upspins_,num_dnspins_);
	psi_inv_.resize(num_upspins_,num_dnspins_);
	// initial configuration
  basis_state_.init_spins(num_upspins_,num_dnspins_);
  // try for a well condictioned amplitude matrix
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

  // work arrays
  psi_row_.resize(num_dnspins_);
  psi_col_.resize(num_upspins_);
  inv_row_.resize(num_upspins_);

  // update parameters
  num_updates_ = 0;
  refresh_cycle_ = 100;
}

void SysConfig::update_state(void)
{

}


