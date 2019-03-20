/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-03-19 14:22:06
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-03-19 14:22:53
*----------------------------------------------------------------------------*/
// File: wavefunction.h
#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <complex>
#include <Eigen/Eigenvalues>
#include "./constants.h"
#include "./matrix.h"
#include "./lattice.h"

enum class wf_id {FEARMISEA, BCS};

class Wavefunction 
{
public:
	Wavefunction() {}
  Wavefunction(const wf_id& id, const Lattice& lattice, const double& hole_doping=0.0)
  	{ init(id, lattice, hole_doping); }
  ~Wavefunction() {}
  void init(const wf_id& id, const Lattice& lattice, const double& hole_doping=0.0);
  void compute(const Lattice& lattice, const RealVector& vparams, 
    const int& start_pos, const bool& psi_gradient=false);
  const int& num_upspins(void) const { return num_upspins_; }
  const int& num_dnspins(void) const { return num_dnspins_; }
  const int& num_vparams(void) const { return num_vparams_; }
  const double& hole_doping(void) const { return hole_doping_; }
  void get_amplitudes(ComplexMatrix& ampl_mat, const std::vector<int>& row,  
    const std::vector<int>& col) const;
  void get_amplitudes(ColVector& ampl_vec, const int& irow,  
    const std::vector<int>& col) const;
  void get_amplitudes(RowVector& ampl_vec, const std::vector<int>& row,
    const int& icol) const;
  void get_amplitudes(std::complex<double>& elem, const int& irow, const int& jcol) const;
  //void get_gradients(Matrix& psi_grad, const int& n, 
  //  const std::vector<int>& row, const std::vector<int>& col) const;
private:
	wf_id id_;
  int num_sites_;
  int num_spins_;
  int num_upspins_;
  int num_dnspins_;
  int num_vparams_;
  double hole_doping_;
  double band_filling_;
  double ch_potential_;
  RealVector vparams_;
  ComplexMatrix psi_;
  std::vector<RealMatrix> psi_gradient_;
  //bool have_gradient_{false};
  // matrices & solvers
	void set_particle_num(const double& hole_doping);
  void compute_BCS(const Lattice& lattice, const RealVector& vparams, 
    const int& start_pos, const bool& psi_gradient=false);
};


#endif