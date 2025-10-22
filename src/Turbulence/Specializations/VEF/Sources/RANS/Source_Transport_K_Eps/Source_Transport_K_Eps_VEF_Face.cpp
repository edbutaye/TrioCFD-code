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
// File:        Source_Transport_K_Eps_VEF_Face.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Sources
//
//////////////////////////////////////////////////////////////////////////////

#include <TRUSTTabs_forward.h>
#include <Source_Transport_K_Eps_VEF_Face.h>
#include <Modele_turbulence_hyd_K_Eps.h>
#include <VEF_discretisation.h>
#include <Transport_K_Eps.h>
#include <Milieu_base.h>
#include <Domaine_VEF.h>

Implemente_instanciable_sans_constructeur(Source_Transport_K_Eps_VEF_Face,
                                          "Source_Transport_K_Eps_VEF_P1NC",
                                          Source_Transport_VEF_Face_base);

Sortie& Source_Transport_K_Eps_VEF_Face::printOn(Sortie& s) const
{
  return s << que_suis_je();
}

Entree& Source_Transport_K_Eps_VEF_Face::readOn(Entree& is)
{
  Source_Transport_VEF_Face_base::verifier_pb_keps(mon_equation->probleme(), que_suis_je());
  return Source_Transport_VEF_Face_base::readOn(is);
}

void Source_Transport_K_Eps_VEF_Face::associer_pb(const Probleme_base& pb)
{
  Source_Transport_VEF_Face_base::associer_pb(pb);
  mon_eq_transport_K_Eps = ref_cast(Transport_K_Eps, equation());
}

void Source_Transport_K_Eps_VEF_Face::creer_champ(const Motcle& nom)
{
  Source_Transport_VEF_Face_base::creer_champ(nom);
  const VEF_discretisation& disc = ref_cast(VEF_discretisation, equation().discretisation());
  Noms noms(1), unites(1);

  if (production_k_face_.est_nul())
    {
      noms[0] = "production_k";
      disc.discretiser_champ("champ_face", equation().domaine_dis(), scalaire,
                             noms , unites, 1, 0, production_k_face_);
      champs_compris_.ajoute_champ(production_k_face_);
    }
}

void Source_Transport_K_Eps_VEF_Face::get_noms_champs_postraitables(Noms& nom,
                                                                    Option opt) const
{
  Source_Transport_VEF_Face_base::get_noms_champs_postraitables(nom, opt);
  Noms noms_compris = champs_compris_.liste_noms_compris();
  noms_compris.add("production_k");

  if (opt == DESCRIPTION)
    Cerr << " Source_Transport_K_Omega_VEF_Face : " << noms_compris << finl;
  else
    nom.add(noms_compris);
}

const DoubleTab& Source_Transport_K_Eps_VEF_Face::get_visc_turb() const
{
  return mon_eq_transport_K_Eps->modele_turbulence().viscosite_turbulente().valeurs();
}

const DoubleTab& Source_Transport_K_Eps_VEF_Face::get_cisaillement_paroi() const
{
  const Modele_turbulence_hyd_K_Eps& mod  = ref_cast(Modele_turbulence_hyd_K_Eps, mon_eq_transport_K_Eps->modele_turbulence());
  return mod.loi_paroi().Cisaillement_paroi();
}

const DoubleTab& Source_Transport_K_Eps_VEF_Face::get_K_pour_production() const
{
  return mon_eq_transport_K_Eps->inconnue().valeurs();
}

const OWN_PTR(Modele_Fonc_Bas_Reynolds_Base)& Source_Transport_K_Eps_VEF_Face::get_modele_fonc_bas_reyn() const
{
  return ref_cast(Modele_turbulence_hyd_K_Eps,
                  mon_eq_transport_K_Eps->modele_turbulence()).associe_modele_fonction();
}

void Source_Transport_K_Eps_VEF_Face::calcul_tabs_bas_reyn(const DoubleTrav& P, const DoubleTab& vit, const DoubleTab& visco_turb, const Champ_Don_base& ch_visco_cin,
                                                           const Champ_base& ch_visco_cin_ou_dyn, DoubleTab& D, DoubleTab& E, DoubleTab& F1, DoubleTab& F2) const
{
  const DoubleTab& K_eps = mon_eq_transport_K_Eps->inconnue().valeurs();
  get_modele_fonc_bas_reyn()->Calcul_D(D, mon_eq_transport_K_Eps->domaine_dis(),
                                       mon_eq_transport_K_Eps->domaine_Cl_dis(),
                                       vit, K_eps, ch_visco_cin);
  D.echange_espace_virtuel();

  get_modele_fonc_bas_reyn()->Calcul_E(E, mon_eq_transport_K_Eps->domaine_dis(),
                                       mon_eq_transport_K_Eps->domaine_Cl_dis(),
                                       vit, K_eps, ch_visco_cin, visco_turb);
  E.echange_espace_virtuel();

  get_modele_fonc_bas_reyn()->Calcul_F1(F1, mon_eq_transport_K_Eps->domaine_dis(),
                                        mon_eq_transport_K_Eps->domaine_Cl_dis(),
                                        P, K_eps, ch_visco_cin_ou_dyn);
  get_modele_fonc_bas_reyn()->Calcul_F2(F2, D, mon_eq_transport_K_Eps->domaine_dis(),
                                        K_eps, ch_visco_cin_ou_dyn);
}

const Nom Source_Transport_K_Eps_VEF_Face::get_type_paroi() const
{
  const Modele_turbulence_hyd_K_Eps& mod  = ref_cast(Modele_turbulence_hyd_K_Eps,mon_eq_transport_K_Eps->modele_turbulence());
  return mod.loi_paroi().que_suis_je();
}

void Source_Transport_K_Eps_VEF_Face::calcul_tenseur_reyn(const DoubleTab& visco_turb,
                                                          const DoubleTab& gradient_elem,
                                                          DoubleTab& Re) const
{
  get_modele_fonc_bas_reyn()->calcul_tenseur_Re(visco_turb, gradient_elem, Re);
}

void Source_Transport_K_Eps_VEF_Face::fill_resu_bas_rey(const DoubleVect& tab_volumes_entrelaces,
                                                        const DoubleTrav& tab_prod, const DoubleTab& tab_D,
                                                        const DoubleTab& tab_E, const DoubleTab& tab_F1,
                                                        const DoubleTab& tab_F2, DoubleTab& tab_resu) const
{
  const double LeK_MIN = mon_eq_transport_K_Eps->modele_turbulence().get_K_MIN();
  const double _C1_ = C1;
  const double _C2_ = C2;
  CDoubleTabView K_eps = mon_eq_transport_K_Eps->inconnue().valeurs().view_ro();
  CDoubleArrView volumes_entrelaces = tab_volumes_entrelaces.view_ro();
  CDoubleArrView prod = static_cast<const ArrOfDouble&>(tab_prod).view_ro();
  CDoubleArrView D = static_cast<const ArrOfDouble&>(tab_D).view_ro();
  CDoubleArrView E = static_cast<const ArrOfDouble&>(tab_E).view_ro();
  CDoubleArrView F1 = static_cast<const ArrOfDouble&>(tab_F1).view_ro();
  CDoubleArrView F2 = static_cast<const ArrOfDouble&>(tab_F2).view_ro();
  DoubleTabView resu = tab_resu.view_rw();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), le_dom_VEF->nb_faces(), KOKKOS_LAMBDA(const int fac)
  {
    const double tke = K_eps(fac, 0);
    const double eps = K_eps(fac, 1);
    resu(fac, 0) += (prod(fac) - eps - D(fac))*volumes_entrelaces(fac);
    if (tke >= LeK_MIN)
      resu(fac, 1) += ((_C1_*prod(fac)*F1(fac) - _C2_*eps*F2(fac))*eps/tke + E(fac))*volumes_entrelaces(fac);
  });
  end_gpu_timer(__KERNEL_NAME__);
}

void Source_Transport_K_Eps_VEF_Face::fill_resu(const DoubleVect& tab_volumes_entrelaces,
                                                const DoubleTrav& tab_prod, DoubleTab& tab_resu) const
{
  const double LeK_MIN = mon_eq_transport_K_Eps->modele_turbulence().get_K_MIN();
  const double _C1_ = C1;
  const double _C2_ = C2;
  CDoubleTabView K_eps = mon_eq_transport_K_Eps->inconnue().valeurs().view_ro();
  CDoubleArrView volumes_entrelaces = tab_volumes_entrelaces.view_ro();
  CDoubleArrView prod = static_cast<const ArrOfDouble&>(tab_prod).view_ro();
  DoubleArrView production_k_face = static_cast<ArrOfDouble&>(ref_cast_non_const(DoubleTab,production_k_face_->valeurs())).view_wo();
  DoubleTabView resu = tab_resu.view_rw();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), le_dom_VEF->nb_faces(), KOKKOS_LAMBDA(const int fac)
  {
    const double tke = K_eps(fac, 0);
    const double eps = K_eps(fac, 1);
    production_k_face(fac)=prod(fac);
    resu(fac, 0) += (prod(fac) - eps)*volumes_entrelaces(fac);
    if (K_eps(fac, 0) >= LeK_MIN)
      resu(fac, 1) += (_C1_*prod(fac) - _C2_*eps)*eps/tke*volumes_entrelaces(fac);
  });
  end_gpu_timer(__KERNEL_NAME__);
}

DoubleTab& Source_Transport_K_Eps_VEF_Face::ajouter(DoubleTab& resu) const
{
  return Source_Transport_VEF_Face_base::ajouter_keps(resu);
}

void Source_Transport_K_Eps_VEF_Face::contribuer_a_avec(const DoubleTab& a,
                                                        Matrice_Morse& matrice_morse) const
{
  const double LeK_MIN = mon_eq_transport_K_Eps->modele_turbulence().get_K_MIN();
  const double _C2_ = C2;
  CDoubleArrView porosite_face = mon_eq_transport_K_Eps->milieu().porosite_face().view_ro();
  CDoubleArrView volumes_entrelaces = le_dom_VEF->volumes_entrelaces().view_ro();
  CDoubleTabView K_eps = mon_eq_transport_K_Eps->inconnue().valeurs().view_ro();
  Matrice_Morse_View matrice;
  matrice.set(matrice_morse);
  // on implicite le -eps et le -eps^2/k
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), K_eps.extent(0), KOKKOS_LAMBDA(const int face)
  {
    if (K_eps(face, 0) >= LeK_MIN) // -eps*vol  donne +vol dans la bonne case
      {
        const double tke = K_eps(face, 0);
        const double eps = K_eps(face, 1);
        const double volporo = porosite_face(face) * volumes_entrelaces(face);

        const double coef_k = eps / tke * volporo;
        matrice.add(face * 2, face * 2, coef_k);

        const double coef_eps = _C2_ * eps / tke * volporo;
        matrice.add(face * 2 + 1, face * 2 + 1, coef_eps);
      }
  });
  end_gpu_timer(__KERNEL_NAME__);
}
