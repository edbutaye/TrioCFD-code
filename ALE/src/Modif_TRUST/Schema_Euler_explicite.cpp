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
// File:        Schema_Euler_explicite.cpp
// Directory:   $TRUST_ROOT/src/Kernel/Schemas_Temps
// Version:     /main/29
//
//////////////////////////////////////////////////////////////////////////////

#include <Schema_Euler_explicite.h>
#include <Equation.h>

#include <Debog.h>
#include <Probleme_base.h>
#include <Domaine.h>
#include <Domaine_ALE.h>

Implemente_instanciable(Schema_Euler_explicite,"Schema_euler_explicite|Scheme_euler_explicit",Schema_Temps_base);


// Description:
//    Simple appel a: Schema_Temps_base::printOn(Sortie& )
//    Ecrit le schema en temps sur un flot de sortie.
// Precondition:
// Parametre: Sortie& s
//    Signification: un flot de sortie
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Sortie&
//    Signification: le flot de sortie modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
Sortie& Schema_Euler_explicite::printOn(Sortie& s) const
{
  return  Schema_Temps_base::printOn(s);
}


// Description:
//    Lit le schema en temps a partir d'un flot d'entree.
//    Simple appel a: Schema_Temps_base::readOn(Entree& )
// Precondition:
// Parametre: Entree& s
//    Signification: un flot d'entree
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Entree&
//    Signification: le flot d'entree modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
Entree& Schema_Euler_explicite::readOn(Entree& s)
{
  return Schema_Temps_base::readOn(s) ;
}


////////////////////////////////
//                            //
// Caracteristiques du schema //
//                            //
////////////////////////////////


// Description:
//    Renvoie le nombre de valeurs temporelles a conserver.
//    Ici : n et n+1, donc 2.
int Schema_Euler_explicite::nb_valeurs_temporelles() const
{
  return 2 ;
}

// Description:
//    Renvoie le nombre de valeurs temporelles futures.
//    Ici : n+1, donc 1.
int Schema_Euler_explicite::nb_valeurs_futures() const
{
  return 1 ;
}

// Description:
//    Renvoie le le temps a la i-eme valeur future.
//    Ici : t(n+1)
double Schema_Euler_explicite::temps_futur(int i) const
{
  assert(i==1);
  return temps_courant()+pas_de_temps();
}

// Description:
//    Renvoie le le temps le temps que doivent rendre les champs a
//    l'appel de valeurs()
//    Ici : t(n+1)
double Schema_Euler_explicite::temps_defaut() const
{
  return temps_courant()+pas_de_temps();
}

/////////////////////////////////////////
//                                     //
// Fin des caracteristiques du schema  //
//                                     //
/////////////////////////////////////////

// Description:
//    Effectue un pas de temps d'Euler explicite
//    sur l'equation passee en parametre.
int Schema_Euler_explicite::faire_un_pas_de_temps_eqn_base(Equation_base& eqn)
{
  DoubleTab& present = eqn.inconnue().valeurs(); // Un
  DoubleTab& futur   = eqn.inconnue().futur();   // Un+1
  DoubleTab dudt(futur);
  Debog::verifier("Schema_Euler_explicite::faire_un_pas_de_temps_eqn_base -futur avant", futur);

  // Boundary conditions applied on Un+1:
  eqn.zone_Cl_dis()->imposer_cond_lim(eqn.inconnue(),temps_courant()+pas_de_temps());

  // On tourne la roue pour que les operateurs utilisent les champs au temps futur
  eqn.inconnue().avancer();
  eqn.derivee_en_temps_inco(dudt);
  eqn.inconnue().reculer();
  dudt.echange_espace_virtuel();
  Debog::verifier("Schema_Euler_explicite::faire_un_pas_de_temps_eqn_base -dudt", dudt);

  // Un+1=Un+dt_*dU/dt
  futur=dudt;
  futur*=dt_;
  futur.echange_espace_virtuel();
  Debog::verifier("Schema_Euler_explicite::faire_un_pas_de_temps_eqn_base -futur apres", futur);

  //Adding ALE Jacobians. Jacobians are renewed in Navier_Stokes_std::corriger_derivee_impl().
  // In ALE Un+1=(Jn/Jn+1)*Un+dt_*dU/dt
  Probleme_base& problem=pb_base();
  if (problem.domaine().que_suis_je()=="Domaine_ALE")
    {
      Domaine_ALE& domaineALE=ref_cast(Domaine_ALE, problem.domaine());
      DoubleTab ALEjacobian_Old=domaineALE.getOldJacobian();
      DoubleTab ALEjacobian_New=domaineALE.getNewJacobian();

      for (int num_face=0; num_face<(futur.size()/dimension); num_face++)
        {
          for (int dim=0; dim<dimension; dim++)
            {
              futur(num_face,dim)+=present(num_face,dim)*(ALEjacobian_Old(num_face,dim)/ALEjacobian_New(num_face,dim));
            }
        }
      futur.echange_espace_virtuel();
      Debog::verifier("Schema_Euler_explicite::faire_un_pas_de_temps_eqn_base -futur apres ALE", futur);
    }
  else  //No ALE
    {
      futur+=present;
    }

  eqn.zone_Cl_dis()->imposer_cond_lim(eqn.inconnue(),temps_courant()+pas_de_temps());
  update_critere_statio(dudt, eqn);

  return 1;
}
