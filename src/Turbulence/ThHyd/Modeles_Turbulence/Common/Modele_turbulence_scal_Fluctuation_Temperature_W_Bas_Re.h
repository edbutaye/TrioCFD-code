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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Modele_turbulence_scal_Fluctuation_Temperature_W_Bas_Re.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/Common/Scal
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_turbulence_scal_Fluctuation_Temperature_W_Bas_Re_included
#define Modele_turbulence_scal_Fluctuation_Temperature_W_Bas_Re_included

#include <Modele_turbulence_scal_Fluctuation_Temperature_W.h>
#include <Transport_Fluctuation_Temperature_W_Bas_Re.h>
#include <Modele_Fonc_Bas_Reynolds_Thermique_Base.h>

class Modele_turbulence_scal_Fluctuation_Temperature_W_Bas_Re :  public Modele_turbulence_scal_Fluctuation_Temperature_W
{
  Declare_instanciable(Modele_turbulence_scal_Fluctuation_Temperature_W_Bas_Re);

public:
  void completer() override;
  int preparer_calcul() override;
  void mettre_a_jour(double ) override;
  //void associer_eqn(const Equation_base&);

  inline Champ_Inc_base& Fluctu_Temperature() override;
  inline const Champ_Inc_base& Fluctu_Temperature() const override;

  inline Transport_Fluctuation_Temperature_W& equation_Fluctu() override;
  inline const Transport_Fluctuation_Temperature_W& equation_Fluctu() const override;

  inline Modele_Fonc_Bas_Reynolds_Thermique_Base& associe_modele_fonction();
  inline const Modele_Fonc_Bas_Reynolds_Thermique_Base& associe_modele_fonction() const;
  Champ_Fonc_base& calculer_diffusivite_turbulente() override;
  void set_param(Param&) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  //////////////////////////////////////////////////////
  //Methode creer_champ pas codee a surcharger si necessaire
  const Champ_base& get_champ(const Motcle& nom) const override;
  bool has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const override;
  bool has_champ(const Motcle& nom) const override;
  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override;
  /////////////////////////////////////////////////////

private :

//Entree& lire(const Motcle&, Entree&);
  OBS_PTR(Transport_Fluctuation_Temperature_W_Bas_Re) eqn_transport_Fluctu_Temp;
  OWN_PTR(Modele_Fonc_Bas_Reynolds_Thermique_Base) mon_modele_fonc;


protected :
  // nous n'avons plus alpha_turb = visco_turb/Prdt_turb
};

inline Transport_Fluctuation_Temperature_W& Modele_turbulence_scal_Fluctuation_Temperature_W_Bas_Re::equation_Fluctu()
{
  return eqn_transport_Fluctu_Temp.valeur();
}

inline const Transport_Fluctuation_Temperature_W& Modele_turbulence_scal_Fluctuation_Temperature_W_Bas_Re::equation_Fluctu() const
{
  return eqn_transport_Fluctu_Temp.valeur();
}

inline const Champ_Inc_base& Modele_turbulence_scal_Fluctuation_Temperature_W_Bas_Re::Fluctu_Temperature() const
{
  return eqn_transport_Fluctu_Temp->inconnue();
}

inline Champ_Inc_base& Modele_turbulence_scal_Fluctuation_Temperature_W_Bas_Re::Fluctu_Temperature()
{
  return eqn_transport_Fluctu_Temp->inconnue();
}

inline Modele_Fonc_Bas_Reynolds_Thermique_Base& Modele_turbulence_scal_Fluctuation_Temperature_W_Bas_Re::associe_modele_fonction()
{
  return mon_modele_fonc.valeur();
}

inline const Modele_Fonc_Bas_Reynolds_Thermique_Base& Modele_turbulence_scal_Fluctuation_Temperature_W_Bas_Re::associe_modele_fonction() const
{
  return mon_modele_fonc.valeur();
}

#endif
