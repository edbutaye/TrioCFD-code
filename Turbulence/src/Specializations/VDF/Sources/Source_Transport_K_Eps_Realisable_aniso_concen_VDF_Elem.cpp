/****************************************************************************
* Copyright (c) 2022, CEA
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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Source_Transport_K_Eps_Realisable_aniso_concen_VDF_Elem.cpp
// Directory : $TURBULENCE_ROOT/src/Specializations/VDF/Sources
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_K_Eps_Realisable_aniso_concen_VDF_Elem.h>
#include <Convection_Diffusion_Temperature.h>
#include <Convection_Diffusion_Concentration.h>
#include <Modele_turbulence_scal_base.h>
#include <Fluide_base.h>
#include <Probleme_base.h>
#include <Champ_Uniforme.h>
#include <Champ_Face.h>
#include <Zone_VDF.h>
#include <Zone_Cl_VDF.h>
#include <Modele_turbulence_hyd_K_Eps_Realisable.h>
#include <DoubleTrav.h>
#include <Fluide_Quasi_Compressible.h>
#include <Convection_Diffusion_Temperature_Turbulent.h>
#include <Ref_Transport_K_Eps_Realisable.h>
#include <Constituant.h>


Implemente_instanciable_sans_constructeur(Source_Transport_K_Eps_Realisable_aniso_concen_VDF_Elem,"Source_Transport_K_Eps_Realisable_aniso_concen_VDF_P0_VDF",Source_Transport_K_Eps_Realisable_VDF_Elem);

Sortie& Source_Transport_K_Eps_Realisable_aniso_concen_VDF_Elem::printOn(Sortie& s) const
{
  return s << que_suis_je() ;
}


Entree& Source_Transport_K_Eps_Realisable_aniso_concen_VDF_Elem::readOn(Entree& is)
{
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  Motcle motlu;

  is >> motlu;
  if (motlu != accolade_ouverte)
    {
      Cerr << "On attendait { pour commencer a lire les constantes de Source_Transport_K_Eps_Realisable_aniso_concen (peut-estre choisi automatiquement en raison de la nature du probleme)" << finl;
      exit();
    }
  Cerr << "Lecture des constantes de Source_Transport_K_Eps_Realisable_aniso_concen (peut-estre choisi automatiquement en raison de la nature du probleme)" << finl;
  Motcles les_mots(2);
  {
    les_mots[0] = "C2_eps";
    les_mots[1] = "C3_eps";
  }
  is >> motlu;
  while (motlu != accolade_fermee)
    {
      int rang=les_mots.search(motlu);
      switch(rang)
        {
        case 0 :
          {
            is >> C2_;
            break;
          }
        case 1 :
          {
            is >> C3_;
            break;
          }
        default :
          {
            Cerr << "On ne comprend pas le mot cle : " << motlu << "dans Source_Transport_K_Eps_Realisable_aniso_concen (peut-estre choisi automatiquement en raison de la nature du probleme)" << finl;
            exit();
          }
        }

      is >> motlu;
    }
  return is;
}




void Source_Transport_K_Eps_Realisable_aniso_concen_VDF_Elem::associer_pb(const Probleme_base& pb)
{
  if (pb.nombre_d_equations()<2)
    {
      Cerr<<"The K_Eps_Realisable source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"for a "<<pb.que_suis_je()<<" problem."<<finl;
    }

  const Equation_base& eqn = pb.equation(1);
  const Milieu_base& milieu = pb.equation(0).milieu();
  const Fluide_base& fluide = ref_cast(Fluide_base,milieu);

  if (sub_type(Fluide_Quasi_Compressible,fluide))
    {
      Cerr<<"The K_Eps_Realisable source term "<<que_suis_je()<<" cannot be activated"<<finl;
      Cerr<<"with a "<<milieu.que_suis_je()<<" medium."<<finl;
      exit();
    }
  Source_Transport_K_Eps_Realisable_VDF_Elem::associer_pb(pb);

  const Convection_Diffusion_Concentration& eqn_c =
    ref_cast(Convection_Diffusion_Concentration,eqn);
  eq_concentration = eqn_c;
  if (!fluide.beta_c().non_nul())
    {
      Cerr << "You forgot to define beta_co field in the fluid." << finl;
      Cerr << "It is mandatory when using the K-Eps_Realisable model (buoyancy effects)." << finl;
      Cerr << "If you don't want buoyancy effects, then specify: beta_co champ_uniforme 1 0." << finl;
      exit();
    }
  beta_c = fluide.beta_c();
  gravite = fluide.gravite();
}


DoubleTab& Source_Transport_K_Eps_Realisable_aniso_concen_VDF_Elem::ajouter(DoubleTab& resu) const
{
  Source_Transport_K_Eps_Realisable_VDF_Elem::ajouter(resu);
  //
  //// Plutost que de calculer P, on appelle Source_Transport_K_Eps_Realisable_VDF_Elem::ajouter(resu)
  //// et on ajoute directement G
  ////
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const Zone_Cl_VDF& zcl_VDF_co = ref_cast(Zone_Cl_VDF,eq_concentration->zone_Cl_dis().valeur());
  const DoubleTab& K_eps = eqn_keps_Rea->inconnue().valeurs();
  const DoubleTab& concen = eq_concentration->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire =
    ref_cast(Modele_turbulence_scal_base,eq_concentration->get_modele(TURBULENCE).valeur());
  const DoubleTab& diffu_turb  = le_modele_scalaire.diffusivite_turbulente().valeurs();
  const Champ_Uniforme& ch_beta_concen = ref_cast(Champ_Uniforme, beta_c->valeur());
  const DoubleVect& g = gravite->valeurs();
  const DoubleVect& volumes = zone_VDF.volumes();
  const DoubleVect& porosite_vol = zone_VDF.porosite_elem();
  int nb_elem = zone_VDF.nb_elem();
  int nb_consti = eq_concentration->constituant().nb_constituants();

  // Ajout d'un espace virtuel au tableau G
  DoubleVect G;
  zone_VDF.zone().creer_tableau_elements(G);

  if (nb_consti == 1)
    {
      double d_beta_c = ch_beta_concen(0,0);
      calculer_terme_destruction_K(zone_VDF,zcl_VDF_co,G,
                                   concen,diffu_turb,d_beta_c,g);
    }
  else
    {
      const DoubleVect& d_beta_c = ch_beta_concen.valeurs();
      calculer_terme_destruction_K(zone_VDF,zcl_VDF_co,G,
                                   concen,diffu_turb,d_beta_c,g,
                                   nb_consti);
    }

  double C1_loc=1.44; // C1 value is not a constant in Realizable K-Epsilon model but here, we take the default value of C1 used in standard K-Epsilon, as proposed by litterature
  double C3_loc;
  double LeK_MIN = eqn_keps_Rea->modele_turbulence().get_LeK_MIN();
  //const Mod_turb_hyd_RANS& mod_turb_RANS = ref_cast(Mod_turb_hyd_RANS,eq_hydraulique->modele_turbulence().valeur());
  //double LeK_MIN = mod_turb_RANS.get_LeK_MIN() ;
  for (int elem=0; elem<nb_elem; elem++)
    {

      resu(elem,0) += G(elem)*volumes(elem)*porosite_vol(elem);

      if (K_eps(elem,0) >= LeK_MIN)
        {
          C3_loc = C3_ ;
          if ( G(elem) > 0. ) C3_loc = 0. ;
          resu(elem,1) += C1_loc*(1.-C3_loc)*G(elem) *volumes(elem)*porosite_vol(elem)
                          * K_eps(elem,1)/K_eps(elem,0);
        }
    }
  return resu;
}

DoubleTab& Source_Transport_K_Eps_Realisable_aniso_concen_VDF_Elem::calculer(DoubleTab& resu) const
{

  resu = 0.;
  return ajouter(resu);

}









