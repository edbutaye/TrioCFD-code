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
// File:        Source_Transport_K_Omega_VEF_Face.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Sources
//
//////////////////////////////////////////////////////////////////////////////

#include <Discretisation_tools.h>
#include <Operateur_Grad.h>
#include <TRUSTTabs_forward.h>
#include <Source_Transport_K_Omega_VEF_Face.h>
#include <Modele_turbulence_hyd_K_Omega.h>
#include <Transport_K_Omega.h>
#include <Navier_Stokes_std.h>
#include <Pb_Hydraulique_Turbulent.h>
#include <Milieu_base.h>
#include <VEF_discretisation.h>
#include <Domaine_VEF.h>
#include <Domaine_Cl_VEF.h>
#include <K_Omega_Eps_constants.h>
#include <Fluide_base.h>

Implemente_instanciable_sans_constructeur(Source_Transport_K_Omega_VEF_Face,
                                          "Source_Transport_K_Omega_VEF_P1NC",
                                          Source_Transport_K_Omega_VEF_Face_base);

Sortie& Source_Transport_K_Omega_VEF_Face::printOn(Sortie& s) const
{
  return s << que_suis_je();
}

Entree& Source_Transport_K_Omega_VEF_Face::readOn(Entree& is)
{
  Source_Transport_K_Omega_VEF_Face_base::verifier_pb_komega(mon_equation->probleme(),
                                                             que_suis_je());
  Source_Transport_K_Omega_VEF_Face_base::readOn(is);
  return (is);
}

void Source_Transport_K_Omega_VEF_Face::get_noms_champs_postraitables(Noms& nom,
                                                                      Option opt) const
{
  Source_Transport_K_Omega_VEF_Face_base::get_noms_champs_postraitables(nom, opt);
  Noms noms_compris = champs_compris_.liste_noms_compris();
  noms_compris.add("grad_k_grad_omega");

  if (opt == DESCRIPTION)
    Cerr << " Source_Transport_K_Omega_VEF_Face : " << noms_compris << finl;
  else
    nom.add(noms_compris);
}

void Source_Transport_K_Omega_VEF_Face::creer_champ(const Motcle& nom)
{
  Source_Transport_K_Omega_VEF_Face_base::creer_champ(nom);
  const VEF_discretisation& disc = ref_cast(VEF_discretisation, equation().discretisation());
  Noms noms(1), unites(1);

  if (grad_k_omega_.est_nul())
    {

      noms[0] = "grad_k_grad_omega";
      disc.discretiser_champ("champ_elem", equation().domaine_dis(), scalaire,
                             noms , unites, 1, 0, grad_k_omega_);
      champs_compris_.ajoute_champ(grad_k_omega_);
    }
  if (grad_k_elem_.est_nul())
    {
      noms[0] = "grad_k";
      disc.discretiser_champ("champ_elem", equation().domaine_dis(), scalaire,
                             noms , unites, dimension, 0, grad_k_elem_);
      champs_compris_.ajoute_champ(grad_k_elem_);
    }
  if (grad_omega_elem_.est_nul())
    {
      noms[0] = "grad_omega";
      disc.discretiser_champ("champ_elem", equation().domaine_dis(), scalaire,
                             noms , unites, dimension, 0, grad_omega_elem_);
      champs_compris_.ajoute_champ(grad_omega_elem_);
    }

  if (production_k_face_.est_nul())
    {
      noms[0] = "production_k";
      disc.discretiser_champ("champ_face", equation().domaine_dis(), scalaire,
                             noms , unites, 1, 0, production_k_face_);
      champs_compris_.ajoute_champ(production_k_face_);
    }
  if (production_omega_face_.est_nul())
    {
      noms[0] = "production_omega";
      disc.discretiser_champ("champ_face", equation().domaine_dis(), scalaire,
                             noms, unites, 1, 0, production_omega_face_);
      champs_compris_.ajoute_champ(production_omega_face_);
    }
  if (dissipation_k_face_.est_nul())
    {
      noms[0] = "dissipation_k";
      disc.discretiser_champ("champ_face", equation().domaine_dis(), scalaire,
                             noms, unites, 1, 0, dissipation_k_face_);
      champs_compris_.ajoute_champ(dissipation_k_face_);
    }
  if (dissipation_omega_face_.est_nul())
    {
      noms[0] = "dissipation_omega";
      disc.discretiser_champ("champ_face", equation().domaine_dis(), scalaire,
                             noms, unites, 1, 0, dissipation_omega_face_);
      champs_compris_.ajoute_champ(dissipation_omega_face_);
    }
  if (cross_diffusion_k_omega_face_.est_nul())
    {
      noms[0] = "cross_diffusion_k_omega";
      disc.discretiser_champ("champ_face", equation().domaine_dis(), scalaire,
                             noms, unites, 1, 0, cross_diffusion_k_omega_face_);
      champs_compris_.ajoute_champ(cross_diffusion_k_omega_face_);
    }

}

void Source_Transport_K_Omega_VEF_Face::associer_pb(const Probleme_base& pb)
{
  Source_Transport_K_Omega_VEF_Face_base::associer_pb(pb);
  eqn_K_Omega = ref_cast(Transport_K_Omega, equation());
  turbulence_model = ref_cast(Modele_turbulence_hyd_K_Omega, eqn_K_Omega->modele_turbulence());
}

const DoubleTab& Source_Transport_K_Omega_VEF_Face::get_visc_turb() const
{
  return eqn_K_Omega->modele_turbulence().viscosite_turbulente().valeurs();
}

const DoubleTab& Source_Transport_K_Omega_VEF_Face::get_cisaillement_paroi() const
{
  return turbulence_model->loi_paroi().Cisaillement_paroi();
}

const DoubleTab& Source_Transport_K_Omega_VEF_Face::get_K_pour_production() const
{
  return eqn_K_Omega->inconnue().valeurs();
}

const Nom Source_Transport_K_Omega_VEF_Face::get_type_paroi() const
{
  return turbulence_model->loi_paroi().que_suis_je();
}

void Source_Transport_K_Omega_VEF_Face::compute_blending_F1(DoubleTab& tab_gradKgradOmega) const
{
  const Navier_Stokes_std& eq_ns = ref_cast(Navier_Stokes_std, eqn_K_Omega->modele_turbulence().equation());
  const double visc_cst = eq_ns.fluide().viscosite_cinematique().valeurs()(0,0);

  DoubleTab& tab_F1 = ref_cast_non_const(DoubleTab, turbulence_model->get_tabF1());
  DoubleTab& tab_F2 = ref_cast_non_const(DoubleTab, turbulence_model->get_tabF2());
  DoubleTrav tab_visc_face(le_dom_VEF->nb_faces_tot()); // dimension_tot(0)
  const double& cst_min_cd_komega = turbulence_model->get_expert_mode().get_cst_cd_komega();
  const bool is_dilatable = eq_ns.fluide().is_dilatable();

  if (is_dilatable)
    elem_to_face(le_dom_VEF.valeur(), eq_ns.fluide().viscosite_cinematique().valeurs(), tab_visc_face);

  CDoubleArrView visc_face;
  if (is_dilatable) visc_face = static_cast<ArrOfDouble&>(tab_visc_face).view_ro();
  CDoubleTabView K_Omega = eqn_K_Omega->inconnue().valeurs().view_ro();
  CDoubleArrView distmin = static_cast<const ArrOfDouble&>(le_dom_VEF->y_faces()).view_ro();
  CDoubleArrView gradKgradOmega = static_cast<const ArrOfDouble&>(tab_gradKgradOmega).view_ro();
  DoubleArrView F1 = static_cast<ArrOfDouble&>(tab_F1).view_wo();
  DoubleArrView F2 = static_cast<ArrOfDouble&>(tab_F2).view_wo();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), le_dom_VEF->nb_faces(), KOKKOS_LAMBDA(const int face)
  {
    double const dmin = Kokkos::fmax(distmin(face), 1e-12);
    double const enerK = K_Omega(face, 0);
    double const omega = K_Omega(face, 1);
    double const kinematic_viscosity = is_dilatable ? visc_face(face) : visc_cst;
    double const tmp1 = Kokkos::sqrt(enerK)/(BETA_K*omega*dmin);
    double const tmp2 = 500.0*kinematic_viscosity/(omega*dmin*dmin);
    double const maxval = Kokkos::fmax(2*SIGMA_OMEGA2*gradKgradOmega(face)/omega, cst_min_cd_komega);
    double const tmp3 = 4.0*SIGMA_OMEGA2*enerK/(maxval*dmin*dmin);

    double const arg1 = Kokkos::fmin(Kokkos::fmax(tmp1, tmp2), tmp3); // Common name of the variable
    F1(face) = Kokkos::tanh(arg1*arg1*arg1*arg1);

    double const arg2 = Kokkos::fmax(2.*tmp1, tmp2);
    F2(face) = Kokkos::tanh(arg2*arg2);
  });
  end_gpu_timer(__KERNEL_NAME__);

  tab_F1.echange_espace_virtuel();
  tab_F2.echange_espace_virtuel();
}

double Source_Transport_K_Omega_VEF_Face::blender(double const val1, double const val2,
                                                  int const face) const
{
  const DoubleTab& F1 = turbulence_model->get_tabF1();
  return F1(face)*val1 + (1 - F1(face))*val2;
}


void Source_Transport_K_Omega_VEF_Face::compute_cross_diffusion(DoubleTab& tab_gradKgradOmega) const
{
  // Same structure than Source_WC_Chaleur_VEF.cpp in TRUST

  // = First, store K and Omega in two separated Tab for the gradient computation, not ideal
  const DoubleTab& tab_K_Omega = eqn_K_Omega->inconnue().valeurs();
  int size = tab_K_Omega.dimension_tot(0);
  DoubleTrav tab_enerK(size); // field on faces
  DoubleTrav tab_omega(size); // field on faces
  CDoubleTabView K_Omega = tab_K_Omega.view_ro();
  DoubleArrView enerK = static_cast<ArrOfDouble&>(tab_enerK).view_wo();
  DoubleArrView omega = static_cast<ArrOfDouble&>(tab_omega).view_wo();
  DoubleArrView gradKgradOmega = static_cast<ArrOfDouble&>(tab_gradKgradOmega).view_rw();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), le_dom_VEF->nb_faces_tot(), KOKKOS_LAMBDA(const int num_face)
  {
    enerK(num_face) = K_Omega(num_face, 0);
    omega(num_face) = K_Omega(num_face, 1);
    gradKgradOmega(num_face) = 0;
  });
  end_gpu_timer(__KERNEL_NAME__);

  // = Definition of the gradients
  // Get number of componants, for gradient resize
  // We use the velocity field to resize the gradient as the velocity is on faces
  const Navier_Stokes_std& eqHyd = ref_cast(Navier_Stokes_std, mon_equation->probleme().equation(0));
  const DoubleTab& velocity_field_face = eqHyd.vitesse().valeurs(); // Velocity on faces
  const int nbr_velocity_components = velocity_field_face.dimension(1); // dimension_tot ?

  DoubleTab& tab_gradK_elem = ref_cast_non_const(DoubleTab, grad_k_elem_->valeurs());
  DoubleTab& tab_gradOmega_elem = ref_cast_non_const(DoubleTab, grad_omega_elem_->valeurs());

  // Compute the two gradients
  const Operateur_Grad& Op_Grad_komega = eqn_K_Omega->gradient_operator_komega();
  Op_Grad_komega.calculer(tab_enerK, tab_gradK_elem);
  Op_Grad_komega.calculer(tab_omega, tab_gradOmega_elem);

  CDoubleTabView gradK_elem = tab_gradK_elem.view_ro();
  CDoubleTabView gradOmega_elem = tab_gradOmega_elem.view_ro();
  DoubleArrView chmp_post = static_cast<ArrOfDouble&>(ref_cast_non_const(DoubleTab, grad_k_omega_->valeurs())).view_rw();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), le_dom_VEF->nb_elem_tot(), KOKKOS_LAMBDA(const int num_elem)
  {
    for (int ncompo = 0; ncompo < nbr_velocity_components; ++ncompo)
      chmp_post(num_elem) += gradK_elem(num_elem, ncompo) * gradOmega_elem(num_elem, ncompo);
  });
  end_gpu_timer(__KERNEL_NAME__);

  // Interpolate from elem to face
  const Domaine_dis_base& domaine_dis = mon_equation->inconnue().domaine_dis_base();
  const Domaine_VF& domaine = ref_cast(Domaine_VF, domaine_dis);

  DoubleTrav tab_gradK_face(velocity_field_face.dimension_tot(0), nbr_velocity_components);
  Discretisation_tools::cells_to_faces(domaine, tab_gradK_elem, tab_gradK_face);

  DoubleTrav tab_gradOmega_face(velocity_field_face.dimension_tot(0), nbr_velocity_components);
  Discretisation_tools::cells_to_faces(domaine, tab_gradOmega_elem, tab_gradOmega_face);

  // Dot Product gradKgradOmega
  CDoubleTabView gradK_face = tab_gradK_face.view_ro();
  CDoubleTabView gradOmega_face = tab_gradOmega_face.view_ro();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), le_dom_VEF->nb_faces_tot(), KOKKOS_LAMBDA(const int num_face)
  {
    for (int ncompo = 0; ncompo < nbr_velocity_components; ++ncompo)
      gradKgradOmega(num_face) +=
        gradK_face(num_face, ncompo) * gradOmega_face(num_face, ncompo);
  });
  end_gpu_timer(__KERNEL_NAME__);
}

void Source_Transport_K_Omega_VEF_Face::fill_resu_k_omega(const DoubleVect& tab_volumes_entrelaces,
                                                          const DoubleTrav& tab_ProdK,
                                                          const DoubleTab& tab_gradKgradOmega,
                                                          DoubleTab& tab_resu) const
{
  const double LeK_MIN = eqn_K_Omega->modele_turbulence().get_K_MIN();
  const double& gamma1 = turbulence_model->get_expert_mode().get_gamma1();
  const double& gamma2 = turbulence_model->get_expert_mode().get_gamma2();
  const bool is_SST_or_BSL = turbulence_model->is_SST_or_BSL();

  CDoubleTabView K_Omega = eqn_K_Omega->inconnue().valeurs().view_ro();
  CDoubleArrView volumes_entrelaces = tab_volumes_entrelaces.view_ro();
  CDoubleArrView ProdK = static_cast<const ArrOfDouble&>(tab_ProdK).view_ro();
  CDoubleArrView gradKgradOmega = static_cast<const ArrOfDouble&>(tab_gradKgradOmega).view_ro();
  CDoubleArrView F1 = static_cast<const ArrOfDouble&>(turbulence_model->get_tabF1()).view_ro();
  DoubleArrView production_omega_face = static_cast<ArrOfDouble&>(ref_cast_non_const(DoubleTab,production_omega_face_->valeurs())).view_wo();
  DoubleArrView production_k_face = static_cast<ArrOfDouble&>(ref_cast_non_const(DoubleTab,production_k_face_->valeurs())).view_wo();
  DoubleArrView dissipation_k_face = static_cast<ArrOfDouble&>(ref_cast_non_const(DoubleTab,dissipation_k_face_->valeurs())).view_wo();
  DoubleArrView dissipation_omega_face = static_cast<ArrOfDouble&>(ref_cast_non_const(DoubleTab,dissipation_omega_face_->valeurs())).view_wo();
  DoubleArrView cross_diffusion_k_omega_face = static_cast<ArrOfDouble&>(ref_cast_non_const(DoubleTab,cross_diffusion_k_omega_face_->valeurs())).view_wo();
  DoubleTabView resu = tab_resu.view_rw();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), le_dom_VEF->nb_faces(), KOKKOS_LAMBDA(const int face)
  {
    const double tke = K_Omega(face, 0);
    const double omega = K_Omega(face, 1);
    production_k_face(face) = ProdK(face);
    dissipation_k_face(face) = BETA_K*tke*omega;
    resu(face, 0) += (ProdK(face) - dissipation_k_face(face))*volumes_entrelaces(face);

    if (tke >= LeK_MIN)
      {
        const double cALPHA = is_SST_or_BSL
                              ? F1(face)*gamma1 + (1 - F1(face))*gamma2
                              : ALPHA_OMEGA;
        const double cBETA = is_SST_or_BSL
                             ? F1(face)*BETA1 + (1 - F1(face))*BETA2
                             : BETA_OMEGA;
        double cSIGMA;
        if (is_SST_or_BSL)
          cSIGMA = 2*(1 - F1(face)*SIGMA_OMEGA2);
        else
          cSIGMA = (gradKgradOmega(face) > 0) ? 0.125 : 0.;

        double contrib {0};
        production_omega_face(face) = cALPHA * ProdK(face) * omega / tke; // production
        dissipation_omega_face(face) = cBETA * omega * omega; // dissipation
        cross_diffusion_k_omega_face(face) =  cSIGMA / omega * gradKgradOmega(face);

        contrib += production_omega_face(face); // production
        contrib -= dissipation_omega_face(face); // dissipation
        contrib += cross_diffusion_k_omega_face(face); // cross diffusion
        contrib *= volumes_entrelaces(face);
        resu(face, 1) += contrib;
      }
  });
  end_gpu_timer(__KERNEL_NAME__);
}

DoubleTab& Source_Transport_K_Omega_VEF_Face::ajouter(DoubleTab& resu) const
{
  return Source_Transport_K_Omega_VEF_Face_base::ajouter_komega(resu);
}

void Source_Transport_K_Omega_VEF_Face::contribuer_a_avec(const DoubleTab& a,
                                                          Matrice_Morse& tab_matrice) const
{
  const double LeK_MIN = eqn_K_Omega->modele_turbulence().get_K_MIN();
  const DoubleTab& tab_K_Omega = equation().inconnue().valeurs();

  DoubleTrav tab_production_TKE {le_dom_VEF->nb_faces()};

  const Domaine_Cl_VEF& domaine_Cl_VEF = ref_cast(Domaine_Cl_VEF,
                                                  eq_hydraulique->domaine_Cl_dis());
  const DoubleTab& visco_turb = get_visc_turb(); // voir les classes filles
  const DoubleTab& velocity = eq_hydraulique->inconnue().valeurs();
  const DoubleTab& TKE = get_K_pour_production(); // voir les classes filles
  const bool& deactivate_production_limiter =  turbulence_model->get_expert_mode().get_deactivate_production_limiter();
  calculer_terme_production_K(le_dom_VEF.valeur(), domaine_Cl_VEF, tab_production_TKE,
                              TKE, velocity, visco_turb,
                              _interpolation_viscosite_turbulente,
                              _coefficient_limiteur,
                              deactivate_production_limiter);

  DoubleTrav tab_gradKgradOmega(le_dom_VEF->nb_faces_tot());
  compute_cross_diffusion(tab_gradKgradOmega);
  if (turbulence_model->is_SST_or_BSL())
    compute_blending_F1(tab_gradKgradOmega);
  const double& gamma1 = turbulence_model->get_expert_mode().get_gamma1();
  const double& gamma2 = turbulence_model->get_expert_mode().get_gamma2();
  const bool is_SST_or_BSL = turbulence_model->is_SST_or_BSL();

  CDoubleTabView K_Omega = tab_K_Omega.view_ro();
  CDoubleArrView porosite_face = eqn_K_Omega->milieu().porosite_face().view_ro();
  CDoubleArrView volumes_entrelaces = le_dom_VEF->volumes_entrelaces().view_ro();
  CDoubleArrView production_TKE = static_cast<const ArrOfDouble&>(tab_production_TKE).view_ro();
  CDoubleArrView gradKgradOmega = static_cast<const ArrOfDouble&>(tab_gradKgradOmega).view_ro();
  CDoubleArrView F1 = is_SST_or_BSL ? static_cast<const ArrOfDouble&>(turbulence_model->get_tabF1()).view_ro() : CDoubleArrView();
  Matrice_Morse_View matrice;
  matrice.set(tab_matrice);
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), tab_K_Omega.dimension(0), KOKKOS_LAMBDA(const int face)
  {
    if (K_Omega(face, 0) >= LeK_MIN)
      {
        const double tke = K_Omega(face, 0);
        const double omega = K_Omega(face, 1);

        const double volporo = porosite_face(face) * volumes_entrelaces(face);

        const double coef_k = BETA_K*omega*volporo;
        matrice.add(face*2, face*2, coef_k);

        const double cALPHA = is_SST_or_BSL
                              ? F1(face)*gamma1 + (1 - F1(face))*gamma2
                              : ALPHA_OMEGA;

        const double cBETA = is_SST_or_BSL
                             ? F1(face)*BETA1 + (1 - F1(face))*BETA2
                             : BETA_OMEGA;

        const double cSIGMA = is_SST_or_BSL
                              ? 2*(1 - F1(face)*SIGMA_OMEGA2)
                              : (gradKgradOmega(face) > 0)*0.125;

        const double coef_omega = (-cALPHA*production_TKE(face)/tke
                                   + cBETA*omega
                                   - cSIGMA/(omega*omega)*gradKgradOmega(face) ) * volporo;
        matrice.add(face*2 + 1, face*2 + 1, coef_omega);
      }
  });
  end_gpu_timer(__KERNEL_NAME__);
}
