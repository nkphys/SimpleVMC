/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-03-19 14:22:06
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-03-24 11:48:35
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
  psi_.resize(num_sites_,num_sites_);
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
  if (lattice.id() == lattice_id::SQUARE) {

    // Chemical potential (non-interacting system) 
    ch_potential_ = 0.0;
    std::vector<double> ek;
    for (int k=0; k<lattice.num_kpoints(); ++k) {
      Vector3d kvec = lattice.kpoint(k);
      double cos_kx = std::cos(kvec[0]);
      double cos_ky = std::cos(kvec[1]);
      double e = -2.0*(cos_kx+cos_ky);
      ek.push_back(e);
    }
    std::sort(ek.begin(),ek.end());
    if (num_upspins_ < num_sites_) {
      ch_potential_ = 0.5*(ek[num_upspins_-1]+ek[num_upspins_]);
    }
    else {
      ch_potential_ = ek[num_upspins_-1];
    }

    //------------BCS wavefunction for SQUARE lattice-----------
    double delta_sc = vparams(start_pos);
    double large_number = 1.0E+4;
    // k-space pair amplitudes 'phi_k' 
    RealVector phi_k(lattice.num_kpoints());
    for (int k=0; k<lattice.num_kpoints(); ++k) {
      Vector3d kvec = lattice.kpoint(k);
      double cos_kx = std::cos(kvec[0]);
      double cos_ky = std::cos(kvec[1]);
      double ek = -2.0*(cos_kx+cos_ky) - ch_potential_;
      if (std::abs(delta_sc) < 1.0E-12) {
        // wavefunction reduces to FEARMISEA
        if (ek < 0.0) phi_k[k] = 1.0;
        else phi_k[k] = 0.0;
      }
      else {
        double deltak = delta_sc*(cos_kx-cos_ky); 
        double deltak_sq = deltak*deltak;
        if (std::sqrt(deltak_sq)<1.0E-12 && ek<0.0) {
          // singular k-points
          phi_k[k] = large_number; 
        }
        else {
          double eps_k = std::sqrt(ek*ek + deltak_sq);
          phi_k[k] = deltak/(ek+eps_k);
        }
      }
      //std::cout << "phi_k["<<k<<"] = "<<phi_k[k]<<"\n"; getchar();
    }
    // pair amplitudes in lattice space
    for (int i=0; i<num_sites_; ++i) {
      Vector3d R_i = lattice.site(i).cell_coord();
      for (int j=0; j<num_sites_; ++j) {
        Vector3d R_j = lattice.site(j).cell_coord();
        Vector3d R_ij = R_i-R_j;
        std::complex<double> ksum=0.0;
        for (int k=0; k<lattice.num_kpoints(); ++k) {
          Vector3d kvec = lattice.kpoint(k);
          std::complex<double> exp_kdor=std::exp(II*kvec.dot(R_ij));
          ksum += phi_k[k] * exp_kdor;
        }
        psi_(i,j) = ksum/double(lattice.num_kpoints());
        //std::cout << "phi["<<i<<","<<j<<"] = "<<psi_(i,j)<<"\n"; getchar();
      }
    }
  }
  else {
    throw std::range_error("BCS wavefunction is not implemented for this lattice\n");
  }
}

void Wavefunction::get_amplitudes(ComplexMatrix& ampl_mat, const std::vector<int>& row, 
  const std::vector<int>& col) const
{
  for (int i=0; i<row.size(); ++i) {
    for (int j=0; j<col.size(); ++j) {
      ampl_mat(i,j) = psi_(row[i],col[j]);
    }
  }
}

void Wavefunction::get_amplitudes(ColVector& ampl_vec, const int& irow,  
    const std::vector<int>& col) const
{
  for (int j=0; j<col.size(); ++j)
    ampl_vec[j] = psi_(irow,col[j]);
}

void Wavefunction::get_amplitudes(RowVector& ampl_vec, const std::vector<int>& row,
    const int& icol) const
{
  for (int j=0; j<row.size(); ++j)
    ampl_vec[j] = psi_(row[j],icol);
}

void Wavefunction::get_amplitudes(std::complex<double>& elem, const int& irow, 
  const int& jcol) const
{
  elem = psi_(irow,jcol);
}









