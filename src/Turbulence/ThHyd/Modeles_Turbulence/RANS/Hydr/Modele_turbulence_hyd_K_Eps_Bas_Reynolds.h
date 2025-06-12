/****************************************************************************
* Copyright (c) 2019, CEA
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

#ifndef Modele_turbulence_hyd_K_Eps_Bas_Reynolds_included
#define Modele_turbulence_hyd_K_Eps_Bas_Reynolds_included

#include <Transport_K_Eps_Bas_Reynolds.h>
#include <Modele_turbulence_hyd_K_Eps.h>



/*! @brief class Modele_turbulence_hyd_K_Eps_Bas_Reynolds
 *
 *  Decrire ici la classe Modele_turbulence_hyd_K_Eps_Bas_Reynolds
 *  Cette classe represente le modele de turbulence (k,eps) pour de faible nombre de Reynolds.
 *
 */
class Modele_turbulence_hyd_K_Eps_Bas_Reynolds: public Modele_turbulence_hyd_RANS_K_Eps_base, public Modele_turbulence_hyd_RANS_Gen<Modele_turbulence_hyd_K_Eps_Bas_Reynolds>
{
  Declare_instanciable(Modele_turbulence_hyd_K_Eps_Bas_Reynolds);
public:
  void set_param(Param& param) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  int preparer_calcul() override;
  bool initTimeStep(double dt) override;
  void mettre_a_jour(double) override;
  void completer() override;
  void controler() override;
  const Champ_base& get_champ(const Motcle& nom) const override;
  bool has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const override;
  bool has_champ(const Motcle& nom) const override;
  void get_noms_champs_postraitables(Noms& nom, Option opt = NONE) const override;
  inline Transport_K_Eps_Bas_Reynolds& get_set_eq_transport() override;
  inline const Transport_K_Eps_Bas_Reynolds& get_eq_transport() const override;

  void imprimer(Sortie&) const override { /* Don nothing */ }
  virtual Champ_Fonc_base& calculer_viscosite_turbulente(double temps);

private:
  void fill_turbulent_viscosity_tab(const int , const DoubleTab&, const DoubleTab& , DoubleTab& );
};

inline Transport_K_Eps_Bas_Reynolds& Modele_turbulence_hyd_K_Eps_Bas_Reynolds::get_set_eq_transport()
{
  Transport_K_Eps_Bas_Reynolds& eq_transport = ref_cast(Transport_K_Eps_Bas_Reynolds,ptr_eq_transport_.valeur());
  return eq_transport;
}

inline const Transport_K_Eps_Bas_Reynolds& Modele_turbulence_hyd_K_Eps_Bas_Reynolds::get_eq_transport() const
{
  const Transport_K_Eps_Bas_Reynolds& eq_transport = ref_cast(Transport_K_Eps_Bas_Reynolds,ptr_eq_transport_.valeur());
  return eq_transport;
}

#endif /* Modele_turbulence_hyd_K_Eps_Bas_Reynolds_included */
