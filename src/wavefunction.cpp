/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-03-19 14:22:06
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-03-19 23:34:38
*----------------------------------------------------------------------------*/
// File: wavefunction.cpp
#include "wavefunction.h"

void Wavefunction::init(const wf_id& id, const Lattice& lattice, const double& hole_doping)
{
	id_ = id;
	num_sites_ = lattice.num_sites();
	switch (id_) {
		case wf_id::BCS: 
      num_vparams_ = 1;
      vparams_.resize(num_vparams_);
			break;
		default: 
			throw std::range_error("This wavefunction not implemented\n");
			break;
	}
  set_particle_num(hole_doping);
}

void Wavefunction::set_particle_num(const double& hole_doping) 
{
  hole_doping_ = hole_doping;
  band_filling_ = 1.0-hole_doping_;
  if (id_==wf_id::BCS) {
    int n = static_cast<int>(std::round(0.5*band_filling_*num_sites_));
    if (n<0 || n>num_sites_) throw std::range_error("Wavefunction:: hole doping out-of-range");
    num_upspins_ = n;
    num_dnspins_ = num_upspins_;
    num_spins_ = num_upspins_ + num_dnspins_;
    band_filling_ = static_cast<double>(2*n)/num_sites_;
  }
  else{
    int n = static_cast<int>(std::round(band_filling_*num_sites_));
    if (n<0 || n>2*num_sites_) throw std::range_error("Wavefunction:: hole doping out-of-range");
    num_spins_ = n;
    num_dnspins_ = num_spins_/2;
    num_upspins_ = num_spins_ - num_dnspins_;
    band_filling_ = static_cast<double>(n)/num_sites_;
  }
  hole_doping_ = 1.0 - band_filling_;
}

void Wavefunction::compute(const Lattice& lattice, const RealVector& vparams, 
    const int& start_pos, const bool& psi_gradient)
{
  switch (id_) {
    case wf_id::BCS: 
      compute_BCS(lattice, vparams, start_pos, psi_gradient);
      break;
    default: 
      throw std::range_error("This wavefunction not implemented\n");
      break;
  }
}

void Wavefunction::compute_BCS(const Lattice& lattice, const RealVector& vparams, 
    const int& start_pos, const bool& psi_gradient)
{
  if (lattice.id() != lattice_id::SQUARE) {
    throw std::range_error("BCS wavefunction is not implemented for this lattice\n");
  }

  double Delta = vparams(start_pos);
  // BCS pair amplitudes for one band system 
  /*
  for (int k=0; k<lattice.num_kpoints(); ++k) {
    Vector3d kvec = lattice.kpoint(k);
    //----------------------------------
    //std::cout << "** Hack at BCS_State\n";
    //ek = -4.0*(std::cos(kvec[0])+std::cos(kvec[1]));
    //delta_k = 0.5 * (std::cos(kvec[0])-std::cos(kvec[1])); 
    //----------------------------------
    double deltak_sq = std::norm(delta_k);
    double ek_plus_Ek = ek + std::sqrt(ek*ek + 4.0*deltak_sq);
    if (std::sqrt(deltak_sq)<1.0E-12 && ek<0.0) {
      phi_k[k](0,0) = large_number_ * std::exp(ii()*std::arg(delta_k));
    }
    else {
      phi_k[k](0,0) = 2.0*delta_k/ek_plus_Ek;
    }
  }
  */
}











