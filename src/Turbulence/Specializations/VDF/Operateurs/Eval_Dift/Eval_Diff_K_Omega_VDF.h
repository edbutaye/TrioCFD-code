/****************************************************************************
* Copyright (c) 2023, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
//////////////////////////////////////////////////////////////////////////////
//
// File:        Eval_Diff_K_Omega_VDF.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Operateurs/Eval_Dift
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Eval_Diff_K_Omega_VDF_included
#define Eval_Diff_K_Omega_VDF_included

#include <Champ_Fonc_base.h>
#include <Champ_Uniforme.h>
#include <Champ_base.h>
#include <TRUST_Ref.h>
#include <K_Omega_constants.h>

class Eval_Diff_K_Omega_VDF
{
public:
  virtual ~Eval_Diff_K_Omega_VDF() { }

  Eval_Diff_K_Omega_VDF(double Sigma_K = SIGMA_K, double Sigma_Omega = SIGMA_OMEGA ) : Sigma_K_(Sigma_K) , Sigma_Omega_(Sigma_Omega), db_diffusivite_(-123.)
  {
    Sigma_[0]=Sigma_K;
    Sigma_[1]=Sigma_Omega;
  }

  inline void associer_diff_turb(const Champ_Fonc_base& diffu) { diffusivite_turbulente_ = diffu; }

  inline void associer_mvolumique(const Champ_base& mvol)
  {
    masse_volumique_ = mvol;
    dv_mvol_.ref(mvol.valeurs());
  }

  inline void associer_Sigma_K_Omega(double Sigma_K,double Sigma_Omega)
  {
    Sigma_K_ = Sigma_K;
    Sigma_Omega_ = Sigma_Omega;
    Sigma_[0] = Sigma_K;
    Sigma_[1] = Sigma_Omega;
  }
  inline void associer_Sigma_K_Omega_SST(double Sigma_K1, double Sigma_K2, double Sigma_Omega1, double Sigma_Omega2)
  {
    Sigma_K_SST_[0] = Sigma_K1;
    Sigma_K_SST_[1] = Sigma_K2;
    Sigma_Omega_SST_[0] = Sigma_Omega1;
    Sigma_Omega_SST_[1] = Sigma_Omega2;
  }
  inline void associer(const Champ_base& diffu)
  {
    diffusivite_ =  diffu;
    if (sub_type(Champ_Uniforme, diffu))  db_diffusivite_ = diffu.valeurs()(0,0);
  }

  inline void associer_tab_F1(const DoubleTab& tabF1)
  {
    tab_F1_.ref(tabF1);
  }
  inline void raise_is_SST_or_BSL() { is_SST_or_BSL_=1; }
  inline void initialize_tab_F1(const int nb_elem_tot)
  {
    tab_F1_.resize(nb_elem_tot);
    tab_F1_=0;
  }

  inline virtual void mettre_a_jour()
  {
    dv_diffusivite_turbulente_.ref(diffusivite_turbulente_->valeurs());
    if (sub_type(Champ_Uniforme, diffusivite_.valeur())) db_diffusivite_ = diffusivite_->valeurs()(0,0);
  }

  inline const Champ_Fonc_base& diffusivite_turbulente() const { return diffusivite_turbulente_.valeur(); }
  inline const Champ_base& diffusivite() const { return diffusivite_.valeur(); }

  // Pour CRTP !
  inline int get_ind_Fluctu_Term() const { throw; }
  inline double get_equivalent_distance(int boundary_index,int local_face) const { throw; }
  inline double get_dv_mvol(const int i) const { return dv_mvol_[i]; }

protected:
  double Sigma_K_, Sigma_Omega_, db_diffusivite_, Sigma_[2];
  double Sigma_K_SST_[2], Sigma_Omega_SST_[2];
  DoubleVect dv_diffusivite_turbulente_, dv_mvol_;
  DoubleTab tab_F1_;
  int is_SST_or_BSL_=0;
  OBS_PTR(Champ_Fonc_base) diffusivite_turbulente_;
  OBS_PTR(Champ_base) masse_volumique_, diffusivite_;
};

#endif /* Eval_Diff_K_Omega_VDF_included */
