/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-03-20 11:50:30
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-03-20 11:51:26
*----------------------------------------------------------------------------*/
// File: sysconfig.h
#ifndef SYSCONFIG_H
#define SYSCONFIG_H

#include "lattice.h"
#include "wavefunction.h"
#include "basis.h"

using amplitude_t = std::complex<double>;

class SysConfig
{
public:
	SysConfig() {}
	SysConfig(const Lattice& lattice) { init(lattice); }
	~SysConfig() {}
	void init(const Lattice& lattice);
	int build(const Lattice& lattice, const RealVector& vparams);
	int init_state(void);
	int update_state(void);
	const int& num_vparams(void) const { return num_total_vparams_; }
  void print_stats(std::ostream& os=std::cout) const;
private:
	int num_sites_;
	int num_upspins_;
	int num_dnspins_;
	double hole_doping_;
	FockBasis basis_state_;
	Wavefunction wf_;
	ComplexMatrix psi_mat_;
	ComplexMatrix psi_inv_;
	// variational parameters
	int num_total_vparams_;
	int num_wf_params_;
	RealVector vparams_;

	// work arrays
  mutable ColVector psi_row_;
  mutable RowVector psi_col_;
  mutable RowVector inv_row_;

	// update parameters_
  int num_updates_;
  int refresh_cycle_;
  int num_proposed_moves_;
  int num_accepted_moves_;

  int do_upspin_hop(void);
  int do_dnspin_hop(void);
  int inv_update_upspin(const int& upspin, const ColVector& psi_row, 
    const std::complex<double>& det_ratio);
  int inv_update_dnspin(const int& dnspin, const RowVector& psi_col, 
    const std::complex<double>& det_ratio);
};


#endif