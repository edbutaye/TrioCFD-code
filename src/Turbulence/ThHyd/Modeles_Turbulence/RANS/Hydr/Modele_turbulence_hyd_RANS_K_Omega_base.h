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
// File:        Modele_turbulence_hyd_RANS_K_Omega_base.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Modele_turbulence_hyd_RANS_K_Omega_base_included
#define Modele_turbulence_hyd_RANS_K_Omega_base_included

#include <Modele_turbulence_hyd_2_eq_base.h>
#include <Modele_turbulence_hyd_RANS_Gen.h>
#include <K_Omega_constants.h>

class Transport_K_Omega_base;

/*! @brief Classe Modele_turbulence_hyd_RANS_K_Omega_base Classe de base des modeles de type RANS_komega
 *
 * @sa Modele_turbulence_hyd_base
 */
class Modele_turbulence_hyd_RANS_K_Omega_base: public Modele_turbulence_hyd_2_eq_base
{
  Declare_base_sans_constructeur(Modele_turbulence_hyd_RANS_K_Omega_base);
public:

  Modele_turbulence_hyd_RANS_K_Omega_base()
  {
    Prandtl_K_ = PRANDTL_K;
    Prandtl_Omega_ = PRANDTL_OMEGA;
    model_variant_ = "SST";
    Sigma_K_ = SIGMA_K;
    Sigma_Omega_ = SIGMA_OMEGA;
    Sigma_K1_ = SIGMA_K1;
    Sigma_K2_ = SIGMA_K2;
    Sigma_OMEGA1_ = SIGMA_OMEGA1;
    Sigma_OMEGA2_ = SIGMA_OMEGA2;
  }

  void set_param(Param& param) override;
  void completer() override;
  int sauvegarder(Sortie& os) const override;
  int reprendre(Entree& is) override;
  void verifie_loi_paroi() override { };
  std::vector<YAML_data> data_a_sauvegarder() const override;

  const Champ_base& get_champ(const Motcle& nom) const override;
  bool has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const override;
  bool has_champ(const Motcle& nom) const override;
  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override;

  inline const Motcle& get_model_variant() const { return model_variant_; }

  const double& get_Sigma_K() const { return Sigma_K_; }
  const double& get_Sigma_Omega() const { return Sigma_Omega_; }
  const double& get_Sigma_K1() const { return Sigma_K1_; }
  const double& get_Sigma_K2() const { return Sigma_K2_; }
  const double& get_Sigma_Omega1() const { return Sigma_OMEGA1_; }
  const double& get_Sigma_Omega2() const { return Sigma_OMEGA2_; }



protected:
  Motcle model_variant_; // default model will be k-omega STD
  static constexpr double CST_A1 = 0.31;
  double Sigma_K_, Sigma_Omega_;
  double Sigma_K1_, Sigma_K2_, Sigma_OMEGA1_, Sigma_OMEGA2_;

};

#endif /* Modele_turbulence_hyd_RANS_K_Omega_base_included */
