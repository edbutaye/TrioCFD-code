/****************************************************************************
* Copyright (c) 2018, CEA
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
// File:        Modele_turbulence_hyd_RANS_Bicephale_base.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_RANS_Bicephale_base.h>
#include <Transport_K_ou_Eps_base.h>
#include <Param.h>

Implemente_base_sans_constructeur(Modele_turbulence_hyd_RANS_Bicephale_base, "Modele_turbulence_hyd_RANS_Bicephale_base", Modele_turbulence_hyd_2_eq_base);
// XD mod_turb_hyd_rans_bicephale modele_turbulence_hyd_deriv mod_turb_hyd_rans_bicephale -1 Class for RANS turbulence model for Navier-Stokes equations.

Sortie& Modele_turbulence_hyd_RANS_Bicephale_base::printOn(Sortie& is) const { return Modele_turbulence_hyd_base::printOn(is); }

Entree& Modele_turbulence_hyd_RANS_Bicephale_base::readOn(Entree& is) { return Modele_turbulence_hyd_base::readOn(is); }

void Modele_turbulence_hyd_RANS_Bicephale_base::set_param(Param& param)
{
  Modele_turbulence_hyd_2_eq_base::set_param(param);
  param.ajouter("eps_min", &EPS_MIN_); // XD_ADD_P double Lower limitation of epsilon (default value 1.e-10).
  param.ajouter("eps_max", &EPS_MAX_); // XD_ADD_P double Upper limitation of epsilon (default value 1.e+10).
  param.ajouter("PRANDTL_K", &Prandtl_K_); // XD_ADD_P double Keyword to change the Prk value (default 1.0).
  param.ajouter("PRANDTL_EPS", &Prandtl_Eps_); // XD_ADD_P double Keyword to change the Pre value (default 1.3)
}

void Modele_turbulence_hyd_RANS_Bicephale_base::completer()
{
  eqn_transp_K().completer();
  eqn_transp_Eps().completer();
  verifie_loi_paroi();
}

const Champ_base& Modele_turbulence_hyd_RANS_Bicephale_base::get_champ(const Motcle& nom) const
{
  try
    {
      return Modele_turbulence_hyd_base::get_champ(nom);
    }
  catch (Champs_compris_erreur&)
    {
    }

  for (int i = 0; i < nombre_d_equations(); i++)
    {
      try
        {
          return equation_k_eps(i).get_champ(nom);
        }
      catch (Champs_compris_erreur&)
        {
        }
    }
  throw Champs_compris_erreur();
}

void Modele_turbulence_hyd_RANS_Bicephale_base::get_noms_champs_postraitables(Noms& nom, Option opt) const
{
  Modele_turbulence_hyd_base::get_noms_champs_postraitables(nom, opt);

  for (int i = 0; i < nombre_d_equations(); i++)
    equation_k_eps(i).get_noms_champs_postraitables(nom, opt);
}

/*! @brief Sauvegarde le modele de turbulence sur un flot de sortie.
 *
 * (en vue d'une reprise)
 *     Sauvegarde le type de l'objet et
 *     les equations de transport K-epsilon associees.
 *
 * @param (Sortie& os) un flot de sortie
 * @return (int) code de retour propage de: Transport_K_ou_Eps::sauvegarder(Sortie&)
 */
int Modele_turbulence_hyd_RANS_Bicephale_base::sauvegarder(Sortie& os) const
{

  Modele_turbulence_hyd_base::sauvegarder(os);
  return (eqn_transp_K().sauvegarder(os) + eqn_transp_Eps().sauvegarder(os));
}

/*! @brief Reprise du modele a partir d'un flot d'entree.
 *
 * Si l'equation portee par l'objet est non nulle
 *     on effectue une reprise "bidon".
 *
 * @param (Entree& is) un flot d'entree
 * @return (int) code de retour propage de: Transport_K_ou_Eps::reprendre(Entree& is) ou 1 si la reprise est bidon.
 */
int Modele_turbulence_hyd_RANS_Bicephale_base::reprendre(Entree& is)
{
  Modele_turbulence_hyd_base::reprendre(is);
  if (mon_equation_.non_nul())
    return (eqn_transp_K().reprendre(is) * eqn_transp_Eps().reprendre(is));
  else
    return reprendre_generique(is);
}

/*! @brief Associe la seconde equation en parametre au modele de turbulence.
 *
 * @param (Equation_base& eqn) la seconde equation a laquelle l'objet s'associe
 */
void Modele_turbulence_hyd_RANS_Bicephale_base::associer_seconde_eqn(const Equation_base& eqn)
{
  ma_seconde_equation_ = eqn;
}
