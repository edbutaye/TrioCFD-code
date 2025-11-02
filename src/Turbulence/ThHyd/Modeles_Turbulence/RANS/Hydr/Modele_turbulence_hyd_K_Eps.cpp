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
// File:        Modele_turbulence_hyd_K_Eps.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_scal_base.h>
#include <Modele_turbulence_hyd_K_Eps.h>
#include <Champ_Inc_P0_base.h>
#include <Schema_Temps_base.h>
#include <communications.h>
#include <Champ_Uniforme.h>
#include <Probleme_base.h>
#include <Perf_counters.h>
#include <Fluide_base.h>
#include <TRUSTTrav.h>
#include <Param.h>
#include <Debog.h>

Implemente_instanciable(Modele_turbulence_hyd_K_Eps, "Modele_turbulence_hyd_K_Epsilon", Modele_turbulence_hyd_RANS_K_Eps_base);
// XD k_epsilon mod_turb_hyd_rans_keps k_epsilon -1 Turbulence model (k-eps).

Sortie& Modele_turbulence_hyd_K_Eps::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}

Entree& Modele_turbulence_hyd_K_Eps::readOn(Entree& s)
{
  return Modele_turbulence_hyd_RANS_K_Eps_base::readOn(s);
}

void Modele_turbulence_hyd_K_Eps::set_param(Param& param)
{
  Modele_turbulence_hyd_RANS_K_Eps_base::set_param(param);
  param.ajouter_non_std("Modele_Fonc_Bas_Reynolds", (this)); // XD_ADD_P modele_fonction_bas_reynolds_base This keyword is used to set the bas Reynolds model used.
  param.ajouter("CMU", &LeCmu_); // XD_ADD_P double Keyword to modify the Cmu constant of k-eps model : Nut=Cmu*k*k/eps Default value is 0.09
  param.ajouter("PRANDTL_K|Sigma_K", &Sigma_K_); // XD_ADD_P double Keyword to change the Prk|Sigma_K value (default 1.0). See https://turbmodels.larc.nasa.gov/ke-chien.html
  param.ajouter("PRANDTL_EPS|Sigma_Eps", &Sigma_Eps_); // XD_ADD_P double Keyword to change the Pre|Sigma_Eps value (default 1.3).
}

int Modele_turbulence_hyd_K_Eps::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot == "Modele_Fonc_Bas_Reynolds")
    {
      Cerr << "Lecture du modele bas reynolds associe " << finl;
      Modele_Fonc_Bas_Reynolds_Base::typer_lire_Modele_Fonc_Bas_Reynolds(mon_modele_fonc_, get_eq_transport(), is);
      Cerr << "mon_modele_fonc.que_suis_je() avant discretisation " << mon_modele_fonc_.que_suis_je() << finl;
      mon_modele_fonc_->discretiser();
      Cerr << "mon_modele_fonc.que_suis_je() " << mon_modele_fonc_->que_suis_je() << finl;
      mon_modele_fonc_->lire_distance_paroi();
      return 1;
    }
  else
    return Modele_turbulence_hyd_RANS_K_Eps_base::lire_motcle_non_standard(mot, is);
}

/*! @brief Calcule la viscosite turbulente au temps demande.
 *
 * @param (double temps) le temps auquel il faut calculer la viscosite
 * @return (Champ_Fonc_base&) la viscosite turbulente au temps demande
 * @throws erreur de taille de visco_turb_K_eps
 */
Champ_Fonc_base& Modele_turbulence_hyd_K_Eps::calculer_viscosite_turbulente(double temps)
{
  const Champ_base& chK_Eps = get_eq_transport().inconnue();
  Nom type = chK_Eps.que_suis_je();
  const Domaine_Cl_dis_base& le_dom_Cl_dis = get_eq_transport().domaine_Cl_dis();
  const DoubleTab& tab_K_Eps = chK_Eps.valeurs();
  DoubleTab& visco_turb = la_viscosite_turbulente_->valeurs();

  // K_Eps(i,0) = K au noeud i
  // K_Eps(i,1) = Epsilon au noeud i
  // int n = tab_K_Eps.dimension(0);
  int n = tab_K_Eps.dimension(0);
  if (n < 0)
    {
      if (sub_type(Champ_Inc_P0_base, chK_Eps))
        n = get_eq_transport().domaine_dis().domaine().nb_elem();
      else
        {
          Cerr << "Unsupported K_Eps field in Modele_turbulence_hyd_K_Eps::calculer_viscosite_turbulente" << finl;
          Process::exit(-1);
        }
    }

  DoubleTrav Fmu, D(tab_K_Eps.dimension(0)), Cmu(tab_K_Eps.dimension(0));
  if (mon_modele_fonc_.non_nul())
    {
      // pour avoir nu en incompressible et mu en QC
      // et non comme on a divise K et eps par rho (si on est en QC)
      // on veut toujours nu
      const OWN_PTR(Champ_Don_base) ch_visco = ref_cast(Fluide_base,get_eq_transport().milieu()).viscosite_cinematique();
      const Champ_Don_base& ch_visco_cin = ref_cast(Fluide_base,get_eq_transport().milieu()).viscosite_cinematique();

      const DoubleTab& tab_visco = ch_visco_cin.valeurs();
      Fmu.resize(tab_K_Eps.dimension(0));
      const Domaine_dis_base& le_dom_dis = get_eq_transport().domaine_dis();

      mon_modele_fonc_->Calcul_Fmu(Fmu, le_dom_dis, le_dom_Cl_dis, tab_K_Eps, ch_visco);
      int is_Cmu_constant = mon_modele_fonc_->Calcul_is_Cmu_constant();
      if (is_Cmu_constant == 0)
        {
          const DoubleTab& vitesse = mon_equation_->inconnue().valeurs();
          mon_modele_fonc_->Calcul_Cmu(Cmu, le_dom_dis, le_dom_Cl_dis, vitesse, tab_K_Eps, EPS_MIN_);

          /*Paroi*/
          Nom lp = get_eq_transport().modele_turbulence().loi_paroi().que_suis_je();
          if (lp != "negligeable_VEF")
            {
              DoubleTab visco_tab(visco_turb.dimension_tot(0));
              assert(sub_type(Champ_Uniforme,ch_visco_cin));
              visco_tab = tab_visco(0, 0);
              const int idt = mon_equation_->schema_temps().nb_pas_dt();
              const DoubleTab& tab_paroi = loi_paroi().Cisaillement_paroi();
              mon_modele_fonc_->Calcul_Cmu_Paroi(Cmu, le_dom_dis, le_dom_Cl_dis, visco_tab, visco_turb, tab_paroi, idt, vitesse, tab_K_Eps, EPS_MIN_);
            }
        }
    }

  // dans le cas d'un domaine nul on doit effectuer le dimensionnement
  double non_prepare = 1;
  Debog::verifier("Modele_turbulence_hyd_K_Eps::calculer_viscosite_turbulente la_viscosite_turbulente before", la_viscosite_turbulente_->valeurs());
  Debog::verifier("Modele_turbulence_hyd_K_Eps::calculer_viscosite_turbulente tab_K_Eps", tab_K_Eps);
  if (visco_turb.size() == n)
    non_prepare = 0.;
  non_prepare = mp_max(non_prepare);

  if (non_prepare == 1)
    {
      OWN_PTR(Champ_Inc_base) visco_turb_au_format_K_eps;
      visco_turb_au_format_K_eps.typer(type);
      DoubleTab& visco_turb_K_eps = complete_viscosity_field(n, get_eq_transport().domaine_dis(), visco_turb_au_format_K_eps);

      if (visco_turb_K_eps.size() != n)
        {
          Cerr << "visco_turb_K_eps size is " << visco_turb_K_eps.size() << " instead of " << n << finl;
          Process::exit();
        }
      // A la fin de cette boucle, le tableau visco_turb_K_eps contient les valeurs de la viscosite turbulente
      // au centre des faces du maillage.
      fill_turbulent_viscosity_tab(n, tab_K_Eps, Cmu, Fmu, D, visco_turb_K_eps);

      // On connait donc la viscosite turbulente au centre des faces de chaque element
      // On cherche maintenant a interpoler cette viscosite turbulente au centre des elements.
      la_viscosite_turbulente_->affecter(visco_turb_au_format_K_eps.valeur());
      Debog::verifier("Modele_turbulence_hyd_K_Eps::calculer_viscosite_turbulente visco_turb_au_format_K_eps", visco_turb_au_format_K_eps.valeur());
    }
  else
    fill_turbulent_viscosity_tab(n, tab_K_Eps, Cmu, Fmu, D, visco_turb);

  la_viscosite_turbulente_->changer_temps(temps);
  Debog::verifier("Modele_turbulence_hyd_K_Eps::calculer_viscosite_turbulente la_viscosite_turbulente after", la_viscosite_turbulente_->valeurs());
  return la_viscosite_turbulente_;
}

void Modele_turbulence_hyd_K_Eps::fill_turbulent_viscosity_tab(const int n, const DoubleTab& tab_K_Eps, const DoubleTab& tab_Cmu, const DoubleTab& tab_Fmu, const DoubleTab& tab_D, DoubleTab& tab_turbulent_viscosity)
{
  const double EPS_MIN = EPS_MIN_;
  const double le_cmu = LeCmu_;
  const bool has_low_re_model = mon_modele_fonc_.non_nul();
  const int is_cmu_constant = has_low_re_model ? mon_modele_fonc_->Calcul_is_Cmu_constant() : 0;

  CDoubleTabView K_Eps = tab_K_Eps.view_ro();
  CDoubleArrView Cmu = static_cast<const ArrOfDouble&>(tab_Cmu).view_ro();
  CDoubleArrView Fmu = static_cast<const ArrOfDouble&>(tab_Fmu).view_ro();
  CDoubleArrView D = static_cast<const ArrOfDouble&>(tab_D).view_ro();
  DoubleArrView turbulent_viscosity = static_cast<ArrOfDouble&>(tab_turbulent_viscosity).view_wo();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), n, KOKKOS_LAMBDA(const int i)
  {
    double eps = K_Eps(i, 1);
    double k = K_Eps(i, 0);
    double value = eps <= EPS_MIN ? 0 : k * k / (eps + D(i));
    turbulent_viscosity[i] = has_low_re_model ? (Fmu(i) * (is_cmu_constant ? le_cmu : Cmu(i)) * value) : le_cmu * value;
  });
  end_gpu_timer(__KERNEL_NAME__);
}

int Modele_turbulence_hyd_K_Eps::preparer_calcul()
{
  get_set_eq_transport().preparer_calcul();
  Modele_turbulence_hyd_RANS_K_Eps_base::preparer_calcul();
  // GF pour initialiser la loi de paroi thermique en TBLE
  if (equation().probleme().nombre_d_equations() > 1)
    {
      const RefObjU& modele_turbulence = equation().probleme().equation(1).get_modele(TURBULENCE);
      if (sub_type(Modele_turbulence_scal_base, modele_turbulence.valeur()))
        {
          Turbulence_paroi_scal_base& loi_paroi_T = ref_cast_non_const(Modele_turbulence_scal_base,modele_turbulence.valeur()).loi_paroi();
          loi_paroi_T.init_lois_paroi();
        }
    }

  calculate_limit_viscosity<MODELE_TYPE::K_EPS>(get_set_unknown(), LeCmu_);
  Debog::verifier("Modele_turbulence_hyd_K_Eps::preparer_calcul la_viscosite_turbulente", la_viscosite_turbulente_->valeurs());
  return 1;
}

bool Modele_turbulence_hyd_K_Eps::initTimeStep(double dt)
{
  return get_set_eq_transport().initTimeStep(dt);
}

/*! @brief Effectue une mise a jour en temps du modele de turbulence.
 *
 * Met a jour l'equation de transport K-epsilon,
 *     calcule la loi de paroi et la viscosite turbulente
 *     au nouveau temps.
 *
 * @param (double temps) le temps de mise a jour
 */
void Modele_turbulence_hyd_K_Eps::mettre_a_jour(double temps)
{
  Schema_Temps_base& sch = get_set_eq_transport().schema_temps();
  get_set_eq_transport().domaine_Cl_dis().mettre_a_jour(temps);
  if (!get_eq_transport().equation_non_resolue())
    {
      statistics().end_count(STD_COUNTERS::update_variables);
      sch.faire_un_pas_de_temps_eqn_base(get_set_eq_transport());
      statistics().begin_count(STD_COUNTERS::update_variables, statistics().get_last_opened_counter_level()+1);
    }
  get_set_eq_transport().mettre_a_jour(temps);

  statistics().begin_count(STD_COUNTERS::turbulent_viscosity, statistics().get_last_opened_counter_level()+1);
  Debog::verifier("Modele_turbulence_hyd_K_Eps::mettre_a_jour la_viscosite_turbulente before", la_viscosite_turbulente_->valeurs());
  calculate_limit_viscosity<MODELE_TYPE::K_EPS>(get_set_unknown(), LeCmu_);
  Debog::verifier("Modele_turbulence_hyd_K_Eps::mettre_a_jour la_viscosite_turbulente after", la_viscosite_turbulente_->valeurs());
  statistics().end_count(STD_COUNTERS::turbulent_viscosity);
}

bool Modele_turbulence_hyd_K_Eps::has_champ(const Motcle& nom, OBS_PTR(Champ_base)& ref_champ) const
{
  if (Modele_turbulence_hyd_RANS_K_Eps_base::has_champ(nom, ref_champ))
    return true;

  if (mon_modele_fonc_.non_nul())
    if (mon_modele_fonc_->has_champ(nom, ref_champ))
      return true;

  return false; /* rien trouve */
}

bool Modele_turbulence_hyd_K_Eps::has_champ(const Motcle& nom) const
{
  if (Modele_turbulence_hyd_RANS_K_Eps_base::has_champ(nom))
    return true;

  if (mon_modele_fonc_.non_nul())
    if (mon_modele_fonc_->has_champ(nom))
      return true;

  return false; /* rien trouve */
}

const Champ_base& Modele_turbulence_hyd_K_Eps::get_champ(const Motcle& nom) const
{
  OBS_PTR(Champ_base) ref_champ;

  if (Modele_turbulence_hyd_RANS_K_Eps_base::has_champ(nom, ref_champ))
    return ref_champ;

  if (mon_modele_fonc_.non_nul())
    if (mon_modele_fonc_->has_champ(nom, ref_champ))
      return ref_champ;

  throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));
}

void Modele_turbulence_hyd_K_Eps::get_noms_champs_postraitables(Noms& nom, Option opt) const
{
  Modele_turbulence_hyd_RANS_K_Eps_base::get_noms_champs_postraitables(nom, opt);
  if (mon_modele_fonc_.non_nul())
    mon_modele_fonc_->get_noms_champs_postraitables(nom, opt);

}
void Modele_turbulence_hyd_K_Eps::verifie_loi_paroi()
{
  Nom lp = loipar_->que_suis_je();
  if (lp == "negligeable_VEF" || lp == "negligeable_VDF")
    if (!associe_modele_fonction().non_nul())
      {
        Cerr << "The turbulence model of type " << que_suis_je() << finl;
        Cerr << "must not be used with a wall law of type negligeable or with a modele_function." << finl;
        Cerr << "Another wall law must be selected with this kind of turbulence model." << finl;
      }
}


Transport_K_Eps& Modele_turbulence_hyd_K_Eps::get_set_eq_transport()
{
  Transport_K_Eps& eq_transport = ref_cast(Transport_K_Eps, ptr_eq_transport_.valeur());
  return eq_transport;
}

const Transport_K_Eps& Modele_turbulence_hyd_K_Eps::get_eq_transport() const
{
  const Transport_K_Eps& eq_transport = ref_cast(Transport_K_Eps, ptr_eq_transport_.valeur());
  return eq_transport;
}

void Modele_turbulence_hyd_K_Eps::controler()
{
  get_set_eq_transport().controler_K_Eps();
}
