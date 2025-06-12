/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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
// File:        Modele_turbulence_hyd_K_Eps_Bicephale.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_turbulence_hyd_K_Eps_Bicephale_included
#define Modele_turbulence_hyd_K_Eps_Bicephale_included

#include <Modele_Fonc_Bas_Reynolds_Base.h>
#include <Transport_K_ou_Eps.h>

/*! @brief Classe Modele_turbulence_hyd_K_Eps_Bicephale Cette classe represente le modele de turbulence (k,eps) pour les
 *
 *     equations de Navier-Stokes ou les 2 equations de k et eps sont gerees separement du point de vue informatique.
 *
 * @sa Modele_turbulence_hyd_base Modele_turbulence_hyd_LES_base
 */
class Modele_turbulence_hyd_K_Eps_Bicephale: public Modele_turbulence_hyd_RANS_Bicephale_base, public Modele_turbulence_hyd_RANS_Gen<Modele_turbulence_hyd_K_Eps_Bicephale>
{
  Declare_instanciable(Modele_turbulence_hyd_K_Eps_Bicephale);
public:
  void set_param(Param& param) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  int preparer_calcul() override;
  void verifie_loi_paroi() override;
  bool initTimeStep(double dt) override;
  void mettre_a_jour(double) override;

  inline OWN_PTR(Modele_Fonc_Bas_Reynolds_Base)& associe_modele_fonction() { return mon_modele_fonc_; }
  inline const OWN_PTR(Modele_Fonc_Bas_Reynolds_Base)& associe_modele_fonction() const { return mon_modele_fonc_; }

  const Champ_base& get_champ(const Motcle& nom) const override;
  bool has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const override;
  bool has_champ(const Motcle& nom) const override;
  void get_noms_champs_postraitables(Noms& nom, Option opt = NONE) const override;
  virtual Champ_Fonc_base& calculer_viscosite_turbulente(double temps);

protected:
  OWN_PTR(Modele_Fonc_Bas_Reynolds_Base) mon_modele_fonc_;
  void fill_turbulent_viscosity_tab(const int , const DoubleTab&, const DoubleTab& , const DoubleTab& ,const DoubleTab& , const DoubleTab&  , DoubleTab& );
};


#endif /* Modele_turbulence_hyd_K_Eps_Bicephale_included */
