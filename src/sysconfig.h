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

class SysConfig
{
public:
	SysConfig() {}
	SysConfig(const Lattice& lattice) { init(lattice); }
	~SysConfig() {}
	void init(const Lattice& lattice);
	void update_state(void);
private:
	int num_sites_;
	int num_upspins_;
	int num_dnspins_;
	double hole_doping_;
	FockBasis basis_state_;
	Wavefunction wf_;
	ComplexMatrix psi_mat_;
	ComplexMatrix psi_inv_;
	// work arrays
  mutable ColVector psi_row_;
  mutable RowVector psi_col_;
  mutable RowVector inv_row_;

	// update parameters_
  int num_updates_;
  int refresh_cycle_;
};


#endif