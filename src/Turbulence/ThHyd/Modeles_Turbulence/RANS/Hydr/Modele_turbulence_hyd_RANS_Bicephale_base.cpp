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
// XD mod_turb_hyd_rans_bicephale mod_turb_hyd_rans mod_turb_hyd_rans_bicephale -1 Class for RANS turbulence model for Navier-Stokes equations.

Sortie& Modele_turbulence_hyd_RANS_Bicephale_base::printOn(Sortie& is) const { return Modele_turbulence_hyd_base::printOn(is); }

Entree& Modele_turbulence_hyd_RANS_Bicephale_base::readOn(Entree& is) { return Modele_turbulence_hyd_base::readOn(is); }

void Modele_turbulence_hyd_RANS_Bicephale_base::set_param(Param& param)
{
  Modele_turbulence_hyd_2_eq_base::set_param(param);
  param.ajouter_non_std("List_transport_equations", (this), Param::REQUIRED);
  param.ajouter("eps_min", &EPS_MIN_); // XD_ADD_P double Lower limitation of epsilon (default value 1.e-20).
  param.ajouter("eps_max", &EPS_MAX_); // XD_ADD_P double Upper limitation of epsilon (default value 1.e+10).
  param.ajouter("PRANDTL_K|Sigma_K", &Sigma_K_); // XD_ADD_P double Keyword to change the Prk|Sigma_K value (default 1.0). See https://turbmodels.larc.nasa.gov/ke-chien.html
  param.ajouter("PRANDTL_EPS|Sigma_Eps", &Sigma_Eps_); // XD_ADD_P double Keyword to change the Pre|Sigma_Eps value (default 1.3).
}

int Modele_turbulence_hyd_RANS_Bicephale_base::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot == "list_transport_equations")
    {
      Motcle dump;
      is >> dump; // read "{" keyword
      Motcle word,which_eq,name_transport_eq;

      is>>word;

      while (word != "}")
        {
          which_eq=word;
          is >> word;
          name_transport_eq=word;
          if (which_eq=="K_EQUATION")
            {
              ptr_eq_transport_K_.typer(name_transport_eq);
              Transport_K_ou_Eps_base& eq_transport_K=ref_cast(Transport_K_ou_Eps_base, ptr_eq_transport_K_.valeur());
              eq_transport_K.set_K_equation(); // to do BEFORE is >> ptr_eq_transport_K_.valeur();
              ptr_eq_transport_K_.valeur().associer_modele_turbulence(*this);
              is >> ptr_eq_transport_K_.valeur();
            }
          else if (which_eq=="EPS_EQUATION")
            {
              ptr_eq_transport_Eps_.typer(name_transport_eq);
              ptr_eq_transport_Eps_.valeur().associer_modele_turbulence(*this);
              is >> ptr_eq_transport_Eps_.valeur();
            }
          else
            {
              Cerr << "Error in Modele_turbulence_hyd_RANS_Bicephale_base::lire_motcle_non_standard"
                   << finl;
              Cerr << "Unknown equation type: " << which_eq << finl;
              Process::exit();
            }
          is>>word;
        }
    }
  else
    return Modele_turbulence_hyd_2_eq_base::lire_motcle_non_standard(mot, is);
  return 1;
}

void Modele_turbulence_hyd_RANS_Bicephale_base::completer()
{
  get_set_eq_transp_K().completer();
  get_set_eq_transp_Eps().completer();
  verifie_loi_paroi();
}

bool Modele_turbulence_hyd_RANS_Bicephale_base::has_champ(const Motcle& nom, OBS_PTR(Champ_base)& ref_champ) const
{
  if (Modele_turbulence_hyd_base::has_champ(nom, ref_champ))
    return true;

  for (int i = 0; i < nombre_d_equations(); i++)
    if (get_equation_k_eps(i).has_champ(nom, ref_champ))
      return true;

  return false; /* rien trouve */
}

bool Modele_turbulence_hyd_RANS_Bicephale_base::has_champ(const Motcle& nom) const
{
  if (Modele_turbulence_hyd_base::has_champ(nom))
    return true;

  for (int i = 0; i < nombre_d_equations(); i++)
    if (get_equation_k_eps(i).has_champ(nom))
      return true;

  return false; /* rien trouve */
}

const Champ_base& Modele_turbulence_hyd_RANS_Bicephale_base::get_champ(const Motcle& nom) const
{
  OBS_PTR(Champ_base) ref_champ;

  if (Modele_turbulence_hyd_base::has_champ(nom, ref_champ))
    return ref_champ;

  for (int i = 0; i < nombre_d_equations(); i++)
    if (get_equation_k_eps(i).has_champ(nom, ref_champ))
      return ref_champ;

  throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));
}

void Modele_turbulence_hyd_RANS_Bicephale_base::get_noms_champs_postraitables(Noms& nom, Option opt) const
{
  Modele_turbulence_hyd_base::get_noms_champs_postraitables(nom, opt);

  for (int i = 0; i < nombre_d_equations(); i++)
    get_equation_k_eps(i).get_noms_champs_postraitables(nom, opt);
}

/*! @brief for PDI IO: retrieve name, type and dimensions of the fields to save/restore
 *
 */
std::vector<YAML_data> Modele_turbulence_hyd_RANS_Bicephale_base::data_a_sauvegarder() const
{
  std::vector<YAML_data> data = Modele_turbulence_hyd_base::data_a_sauvegarder();
  std::vector<YAML_data> eqn_transp = get_eq_transp_K().data_a_sauvegarder();
  data.insert(data.end(), eqn_transp.begin(), eqn_transp.end());
  return data;
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
  return (get_eq_transp_K().sauvegarder(os) + get_eq_transp_Eps().sauvegarder(os));
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
    return (get_set_eq_transp_K().reprendre(is) * get_set_eq_transp_Eps().reprendre(is));
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

void Modele_turbulence_hyd_RANS_Bicephale_base::controler()
{
  get_set_eq_transp_K().controler_variable();
  get_set_eq_transp_Eps().controler_variable();
}

/*! @brief Renvoie le champ inconnue K du modele de turbulence Cette inconnue est portee
 *
 *     par l equation de transport K porte par le modele.
 *     (version const)
 *
 * @return (Champ_Inc_base&) le champ inconnue (K)
 */
const Champ_Inc_base& Modele_turbulence_hyd_RANS_Bicephale_base::get_K() const
{
  return ptr_eq_transport_K_.valeur().inconnue();
}

/*! @brief Renvoie le champ inconnue K du modele de turbulence Cette inconnue est portee
 *
 *     par l equation de transport K porte par le modele.
 *
 * @return (Champ_Inc_base&) le champ inconnue (K)
 */
Champ_Inc_base& Modele_turbulence_hyd_RANS_Bicephale_base::get_set_K()
{
  return ptr_eq_transport_K_.valeur().inconnue();
}

/*! @brief Renvoie le champ inconnue epsilon du modele de turbulence Cette inconnue est portee
 *
 *     par l equation de transport epsilon porte par le modele.
 *     (version const)
 *
 * @return (Champ_Inc_base&) le champ inconnue (epsilon)
 */
const Champ_Inc_base& Modele_turbulence_hyd_RANS_Bicephale_base::get_Eps() const
{
  return ptr_eq_transport_Eps_.valeur().inconnue();
}

/*! @brief Renvoie le champ inconnue epsilon du modele de turbulence Cette inconnue est portee
 *
 *     par l equation de transport epsilon porte par le modele.
 *
 * @return (Champ_Inc_base&) le champ inconnue (epsilon)
 */
Champ_Inc_base& Modele_turbulence_hyd_RANS_Bicephale_base::get_set_Eps()
{
  return ptr_eq_transport_Eps_.valeur().inconnue();
}

/*! @brief Renvoie l equation d evolution de K du modele de turbulence
 *
 * @return (Transport_K_ou_Eps&) equation (K)
 */
Transport_K_ou_Eps_base& Modele_turbulence_hyd_RANS_Bicephale_base::get_set_eq_transp_K()
{
  Transport_K_ou_Eps_base& eq_transport_K=ref_cast(Transport_K_ou_Eps_base, ptr_eq_transport_K_.valeur());
  return eq_transport_K;
}

/*! @brief Renvoie l equation d evolution de K du modele de turbulence (version const)
 *
 * @return (Transport_K_ou_Eps&) equation (K)
 */
const Transport_K_ou_Eps_base& Modele_turbulence_hyd_RANS_Bicephale_base::get_eq_transp_K() const
{
  const Transport_K_ou_Eps_base& eq_transport_K=ref_cast(Transport_K_ou_Eps_base, ptr_eq_transport_K_.valeur());
  return eq_transport_K;
}

/*! @brief Renvoie l equation d evolution de epsilon du modele de turbulence
 *
 * @return (Transport_K_ou_Eps&) equation (epsilon)
 */
Transport_K_ou_Eps_base& Modele_turbulence_hyd_RANS_Bicephale_base::get_set_eq_transp_Eps()
{
  Transport_K_ou_Eps_base& eq_transport_Eps=ref_cast(Transport_K_ou_Eps_base, ptr_eq_transport_Eps_.valeur());
  return eq_transport_Eps;
}

/*! @brief Renvoie l equation d evolution de epsilon du modele de turbulence (version const)
 *
 * @return (Transport_K_ou_Eps&) equation (epsilon)
 */
const Transport_K_ou_Eps_base& Modele_turbulence_hyd_RANS_Bicephale_base::get_eq_transp_Eps() const
{
  const Transport_K_ou_Eps_base& eq_transport_Eps=ref_cast(Transport_K_ou_Eps_base, ptr_eq_transport_Eps_.valeur());
  return eq_transport_Eps;
}
