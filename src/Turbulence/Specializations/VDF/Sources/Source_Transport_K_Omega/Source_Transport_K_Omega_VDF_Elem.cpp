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
// File:        Source_Transport_K_Omega_VDF_Elem.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Sources
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_K_Omega_VDF_Elem.h>
#include <Modele_turbulence_hyd_K_Omega.h>
#include <Pb_Hydraulique_Turbulent.h>
#include <Discretisation_tools.h>
#include <Navier_Stokes_std.h>
#include <Milieu_base.h>
#include <TRUSTTrav.h>
#include <VDF_discretisation.h>

Implemente_instanciable_sans_constructeur(Source_Transport_K_Omega_VDF_Elem,
                                          "Source_Transport_K_Omega_VDF_P0_VDF",
                                          Source_Transport_K_Omega_VDF_Elem_base);

Sortie& Source_Transport_K_Omega_VDF_Elem::printOn(Sortie& s) const { return s << que_suis_je() ; }

Entree& Source_Transport_K_Omega_VDF_Elem::readOn(Entree& is)
{
  Source_Transport_K_Omega_VDF_Elem_base::verifier_pb_komega(mon_equation->probleme(), que_suis_je());
  return Source_Transport_K_Omega_VDF_Elem_base::readOn(is);
}

void Source_Transport_K_Omega_VDF_Elem::get_noms_champs_postraitables(Noms& nom,
                                                                      Option opt) const
{
  Source_Transport_K_Omega_VDF_Elem_base::get_noms_champs_postraitables(nom, opt);
  Noms noms_compris = champs_compris_.liste_noms_compris();
  noms_compris.add("grad_k_grad_omega");

  if (opt == DESCRIPTION)
    Cerr << " Source_Transport_K_Omega_VEF_Face : " << noms_compris << finl;
  else
    nom.add(noms_compris);
}

void Source_Transport_K_Omega_VDF_Elem::creer_champ(const Motcle& nom)
{
  Source_Transport_K_Omega_VDF_Elem_base::creer_champ(nom);
  const VDF_discretisation& disc = ref_cast(VDF_discretisation, equation().discretisation());
  Noms noms(1), unites(1);

  if (grad_k_omega_elem_.est_nul())
    {
      noms[0] = "grad_k_grad_omega";
      disc.discretiser_champ("champ_elem", equation().domaine_dis(), scalaire,
                             noms , unites, 1, 0, grad_k_omega_elem_);
      champs_compris_.ajoute_champ(grad_k_omega_elem_);
    }
  if (grad_k_face_.est_nul())
    {
      noms[0] = "grad_k";
      disc.discretiser_champ("vitesse", equation().domaine_dis(), scalaire,
                             noms , unites, dimension, 0, grad_k_face_);
      champs_compris_.ajoute_champ(grad_k_face_);
    }
  if (grad_omega_face_.est_nul())
    {
      noms[0] = "grad_omega";
      disc.discretiser_champ("vitesse", equation().domaine_dis(), scalaire,
                             noms , unites, dimension, 0, grad_omega_face_);
      champs_compris_.ajoute_champ(grad_omega_face_);
    }
  if (grad_k_elem_.est_nul())
    {
      noms[0] = "grad_k_elem";
      disc.discretiser_champ("champ_elem", equation().domaine_dis(), scalaire,
                             noms , unites, dimension, 0, grad_k_elem_);
      champs_compris_.ajoute_champ(grad_k_elem_);
    }
  if (grad_omega_elem_.est_nul())
    {
      noms[0] = "grad_omega_elem";
      disc.discretiser_champ("champ_elem", equation().domaine_dis(), scalaire,
                             noms , unites, dimension, 0, grad_omega_elem_);
      champs_compris_.ajoute_champ(grad_omega_elem_);
    }
  if (production_k_elem_.est_nul())
    {
      noms[0] = "production_k";
      disc.discretiser_champ("champ_elem", equation().domaine_dis(), scalaire,
                             noms , unites, 1, 0, production_k_elem_);
      champs_compris_.ajoute_champ(production_k_elem_);
    }
  if (production_omega_elem_.est_nul())
    {
      noms[0] = "production_omega";
      disc.discretiser_champ("champ_elem", equation().domaine_dis(), scalaire,
                             noms, unites, 1, 0, production_omega_elem_);
      champs_compris_.ajoute_champ(production_omega_elem_);
    }
  if (dissipation_k_elem_.est_nul())
    {
      noms[0] = "dissipation_k";
      disc.discretiser_champ("champ_elem", equation().domaine_dis(), scalaire,
                             noms, unites, 1, 0, dissipation_k_elem_);
      champs_compris_.ajoute_champ(dissipation_k_elem_);
    }
  if (dissipation_omega_elem_.est_nul())
    {
      noms[0] = "dissipation_omega";
      disc.discretiser_champ("champ_elem", equation().domaine_dis(), scalaire,
                             noms, unites, 1, 0, dissipation_omega_elem_);
      champs_compris_.ajoute_champ(dissipation_omega_elem_);
    }
  if (cross_diffusion_k_omega_elem_.est_nul())
    {
      noms[0] = "cross_diffusion_k_omega";
      disc.discretiser_champ("champ_elem", equation().domaine_dis(), scalaire,
                             noms, unites, 1, 0, cross_diffusion_k_omega_elem_);
      champs_compris_.ajoute_champ(cross_diffusion_k_omega_elem_);
    }
}


void Source_Transport_K_Omega_VDF_Elem::associer_pb(const Probleme_base& pb)
{
  Source_Transport_K_Omega_VDF_Elem_base::associer_pb(pb);
  eqn_K_Omega = ref_cast(Transport_K_Omega, equation());
  turbulence_model = ref_cast(Modele_turbulence_hyd_K_Omega, eqn_K_Omega->modele_turbulence());
}

const DoubleTab& Source_Transport_K_Omega_VDF_Elem::get_visc_turb() const
{
  return eqn_K_Omega->modele_turbulence().viscosite_turbulente().valeurs();
}

void Source_Transport_K_Omega_VDF_Elem::calculer_terme_production(const Champ_Face_VDF& vitesse, const DoubleTab& visco_turb, const DoubleTab& vit, DoubleVect& P) const
{
  const DoubleTab& K_Omega = eqn_K_Omega->inconnue().valeurs();
  if (axi) calculer_terme_production_K_Axi(le_dom_VDF.valeur(), vitesse, P, K_Omega, visco_turb);
  else calculer_terme_production_K_for_komega(le_dom_VDF.valeur(), le_dom_Cl_VDF.valeur(), P, K_Omega, vit, vitesse, visco_turb);
}

void Source_Transport_K_Omega_VDF_Elem::fill_resu(const DoubleVect& P, DoubleTab& resu) const
{
  const DoubleVect& volumes = le_dom_VDF->volumes();
  const DoubleVect& porosite_vol = le_dom_Cl_VDF->equation().milieu().porosite_elem();
  const DoubleTab& K_Omega = eqn_K_Omega->inconnue().valeurs();
  const double LeK_MIN = eqn_K_Omega->modele_turbulence().get_K_MIN();
  const DoubleTab& gradKgradOmega_elem = grad_k_omega_elem_->valeurs();
  DoubleTab& production_omega_elem = ref_cast_non_const(DoubleTab,production_omega_elem_->valeurs());
  DoubleTab& dissipation_k_elem = ref_cast_non_const(DoubleTab,dissipation_k_elem_->valeurs());
  DoubleTab& dissipation_omega_elem = ref_cast_non_const(DoubleTab,dissipation_omega_elem_->valeurs());
  DoubleTab& cross_diffusion_k_omega_elem = ref_cast_non_const(DoubleTab,cross_diffusion_k_omega_elem_->valeurs());

  for (int elem = 0; elem < le_dom_VDF->nb_elem(); ++elem)
    {
      const double tke = K_Omega(elem, 0);
      const double omega = K_Omega(elem, 1);
      double volporo = volumes(elem)*porosite_vol(elem);
      dissipation_k_elem(elem) = BETA_K*tke*omega;
      resu(elem, 0) += (P(elem) - dissipation_k_elem(elem))*volporo;

      if (tke >= LeK_MIN)
        {
          const double cALPHA = turbulence_model->is_SST_or_BSL()
                                ? blender(GAMMA1, GAMMA2, elem)
                                : ALPHA_OMEGA;
          const double cBETA = turbulence_model->is_SST_or_BSL()
                               ? blender(BETA1, BETA2, elem)
                               : BETA_OMEGA;

          double cSIGMA;
          if (turbulence_model->is_SST_or_BSL())
            cSIGMA = 2*(1 - turbulence_model->get_tabF1()(elem)*SIGMA_OMEGA2);
          else
            cSIGMA = (gradKgradOmega_elem(elem) > 0) ? 0.125 : 0.;

          double contrib { 0 };
          production_omega_elem(elem) = cALPHA * P(elem) * omega / tke; // production
          dissipation_omega_elem(elem) = cBETA * omega * omega; // dissipation
          cross_diffusion_k_omega_elem(elem) =  cSIGMA / omega * gradKgradOmega_elem(elem);

          contrib += production_omega_elem(elem); //production
          contrib -= dissipation_omega_elem(elem); // dissipation
          contrib += cross_diffusion_k_omega_elem(elem); // cross diffusion
          contrib *= volporo;
          resu(elem, 1) += contrib;
        }
    }
}

void Source_Transport_K_Omega_VDF_Elem::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  Source_Transport_K_Omega_VDF_Elem_base::ajouter_komega(secmem);

  const std::string& nom_inco = equation().inconnue().le_nom().getString();
  Matrice_Morse* mat = matrices.count(nom_inco) ? matrices.at(nom_inco) : nullptr;
  if(!mat) return;

  const DoubleTab& K_Omega = equation().inconnue().valeurs();
  const DoubleVect& porosite = le_dom_Cl_VDF->equation().milieu().porosite_elem(), &volumes = le_dom_VDF->volumes();
  const int size = K_Omega.dimension(0);
  // on implicite le -eps et le -eps^2/k
  // cAlan : impliciter omega ?
  const DoubleTab& gradKgradOmega_elem = grad_k_omega_elem_->valeurs();
  const DoubleVect& production_TKE = production_k_elem_->valeurs();
  for (int c = 0; c < size; ++c)
    {
      // -eps*vol  donne +vol dans la bonne case
      if (K_Omega(c, 0) > DMINFLOAT)
        {
          const double tke = K_Omega(c, 0);
          const double omega = K_Omega(c, 1);

          // cAlan : a adapter
          double coef_k = porosite(c)*volumes(c)*BETA_K*omega;
          (*mat)(c*2, c*2) += coef_k;
          const double cALPHA = turbulence_model->is_SST_or_BSL()
                                ? blender(GAMMA1, GAMMA2, c)
                                : ALPHA_OMEGA;

          const double cBETA = turbulence_model->is_SST_or_BSL()
                               ? blender(BETA1, BETA2, c)
                               : BETA_OMEGA;

          double cSIGMA;
          if (turbulence_model->is_SST_or_BSL())
            cSIGMA = 2*(1 - turbulence_model->get_tabF1()(c)*SIGMA_OMEGA2);
          else
            cSIGMA = (gradKgradOmega_elem(c) > 0) ? 0.125 : 0.;

          const double coef_omega = (-cALPHA*production_TKE(c)/tke
                                     + cBETA*omega
                                     - cSIGMA/(omega*omega)*gradKgradOmega_elem(c) ) * porosite(c)*volumes(c);
          (*mat)(c*2 + 1, c*2 + 1) += coef_omega;
        }
    }
}

void Source_Transport_K_Omega_VDF_Elem::compute_cross_diffusion() const
{
  const int nb_elem_tot = le_dom_VDF->nb_elem_tot();
  const int nb_faces_tot = le_dom_VDF->nb_faces_tot();
  const int nb_faces = le_dom_VDF->nb_faces();

  const int nb_elem = le_dom_VDF->nb_elem();
  DoubleTab& gradKgradOmega_elem = ref_cast_non_const(DoubleTab, grad_k_omega_elem_->valeurs());
  const DoubleTab& K_Omega = eqn_K_Omega->inconnue().valeurs();
  DoubleTab enerK; // field on elem
  DoubleTab omega; // field on elem

  enerK.resize(nb_elem_tot,1);
  omega.resize(nb_elem_tot,1);
  for (int elem = 0; elem < nb_elem_tot; ++elem)
    {
      enerK(elem,0) = K_Omega(elem, 0);
      omega(elem,0) = K_Omega(elem, 1);
      gradKgradOmega_elem(elem) = 0;
    }

  DoubleTab gradK_face_tmp; // field on faces
  DoubleTab  gradOmega_face_tmp; // field on faces
  gradK_face_tmp.resize(nb_faces_tot,1);
  gradOmega_face_tmp.resize(nb_faces_tot,1);

  // Compute the two gradients
  const Operateur_Grad& Op_Grad_komega = eqn_K_Omega->gradient_operator_komega();
  Op_Grad_komega.calculer(enerK, gradK_face_tmp);
  Op_Grad_komega.calculer(omega, gradOmega_face_tmp);

  DoubleTab& tab_gradK_face = ref_cast_non_const(DoubleTab, grad_k_face_->valeurs());
  DoubleTab& tab_gradOmega_face = ref_cast_non_const(DoubleTab, grad_omega_face_->valeurs());
  for (int face = 0; face < nb_faces; ++face)
    {
      tab_gradK_face(face) = gradK_face_tmp(face, 0);
      tab_gradOmega_face(face) = gradOmega_face_tmp(face, 0);
    }
  tab_gradK_face.echange_espace_virtuel();
  tab_gradOmega_face.echange_espace_virtuel();

  // Interpolate from faces to elem
  const Domaine_dis_base& domaine_dis = mon_equation->inconnue().domaine_dis_base();
  /*DoubleTab gradK_elem, gradOmega_elem;
  gradK_elem.resize(nb_elem_tot,dimension);
  gradOmega_elem.resize(nb_elem_tot,dimension);
  */
  DoubleTab& tab_gradK_elem = ref_cast_non_const(DoubleTab, grad_k_elem_->valeurs());
  DoubleTab& tab_gradOmega_elem = ref_cast_non_const(DoubleTab, grad_omega_elem_->valeurs());

  grad_k_face_->valeur_aux_centres_de_gravite(domaine_dis.domaine(), tab_gradK_elem);
  grad_omega_face_->valeur_aux_centres_de_gravite(domaine_dis.domaine(), tab_gradOmega_elem);

  // Dot Product gradKgradOmega
  // EB :  no need to update virtual cells
  for (int elem = 0; elem < nb_elem; ++elem)
    for (int dim=0; dim<dimension; ++dim)
      gradKgradOmega_elem(elem) += tab_gradK_elem(elem,dim) * tab_gradOmega_elem(elem,dim);
}

void Source_Transport_K_Omega_VDF_Elem::compute_blending_F1() const
{
  const DoubleTab& K_Omega = eqn_K_Omega->inconnue().valeurs();
  const DoubleTab& kinematic_viscosity = get_visc_turb();
  const DoubleTab& distmin = le_dom_VDF->y_elem(); // Minimum distance to the edge
  DoubleTab& gradKgradOmega_elem = ref_cast_non_const(DoubleTab, grad_k_omega_elem_->valeurs());

  DoubleTab& tabF1 = ref_cast_non_const(DoubleTab, turbulence_model->get_tabF1());
  DoubleTab& tabF2 = ref_cast_non_const(DoubleTab, turbulence_model->get_tabF2());

  for (int elem = 0; elem < le_dom_VDF->nb_elem(); ++elem)
    {
      double const dmin = std::max(distmin(elem), 1e-12);
      double const enerK = K_Omega(elem, 0);
      double const omega = K_Omega(elem, 1);

      double const tmp1 = sqrt(enerK)/(BETA_K*omega*dmin);
      double const tmp2 = 500.0*kinematic_viscosity(elem)/(omega*dmin*dmin);
      double const maxval = std::max(2*SIGMA_OMEGA2*gradKgradOmega_elem(elem)/omega, 1e-20);
      double const tmp3 = 4.0*SIGMA_OMEGA2*enerK/(maxval*dmin*dmin);

      double const arg1 = std::min(std::max(tmp1, tmp2), tmp3); // Common name of the variable
      tabF1(elem) = std::tanh(arg1*arg1*arg1*arg1);
      double const arg2 = std::max(2.*tmp1, tmp2);
      tabF2(elem) = std::tanh(arg2*arg2);
    }
  tabF1.echange_espace_virtuel();
  tabF2.echange_espace_virtuel();
}

double Source_Transport_K_Omega_VDF_Elem::blender(double const val1, double const val2,
                                                  int const elem) const
{
  const DoubleTab& F1 = turbulence_model->get_tabF1();
  return F1(elem)*val1 + (1 - F1(elem))*val2;
}

