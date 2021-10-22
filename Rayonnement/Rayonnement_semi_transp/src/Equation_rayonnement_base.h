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
// File:        Equation_rayonnement_base.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Equation_rayonnement_base_included
#define Equation_rayonnement_base_included

#include <Operateur_Diff.h>
#include <Matrice_Morse.h>
#include <Ref_Fluide_base.h>
#include <Ref_Modele_rayo_semi_transp.h>
#include <Operateur_Grad.h>
class Motcle;
class Milieu_base;
class Fluide_base;

class Equation_rayonnement_base: public Equation_base
{
  Declare_base_sans_constructeur(Equation_rayonnement_base);
public:

  Equation_rayonnement_base();
  void set_param(Param& titi);
  int lire_motcle_non_standard(const Motcle&, Entree&);
  virtual bool initTimeStep(double dt);
  virtual bool solve();
  void associer_milieu_base(const Milieu_base&);
  inline void associer_fluide(const Fluide_base&);
  void associer_modele_rayonnement(const Modele_rayo_semi_transp&);
  Milieu_base& milieu();
  const Milieu_base& milieu() const;
  inline const Modele_rayo_semi_transp&  Modele() const;
  inline Modele_rayo_semi_transp& Modele();
  int nombre_d_operateurs() const;
  const Operateur& operateur(int) const;
  Operateur& operateur(int);
  virtual const Champ_Inc& inconnue() const;
  virtual Champ_Inc& inconnue();
  void discretiser();
  inline Fluide_base& fluide();
  inline const Fluide_base& fluide() const;

  // pas de flux calcule correctement par les operateurs...
  inline int impr(Sortie& os) const
  {
    return 1;
  };

  virtual int nb_colonnes_tot()=0;
  virtual int nb_colonnes()=0;
  void Mat_Morse_to_Mat_Bloc(Matrice& matrice_tmp);
  void dimensionner_Mat_Bloc_Morse_Sym(Matrice& matrice_tmp);


  const Discretisation_base& discretisation() const;

  Operateur_Grad& operateur_gradient();
  const Operateur_Grad& operateur_gradient() const;
  virtual void typer_op_grad()=0;

  void associer_pb_base(const Probleme_base& pb);

  virtual void resoudre(double temps)=0;
  virtual void assembler_matrice()=0;


  /////////////////////////////////////////////////////

  virtual void modifier_matrice()=0;
  virtual void evaluer_cl_rayonnement(double temps)=0;

  void completer();

protected:

  REF(Fluide_base) le_fluide;
  REF(Modele_rayo_semi_transp) le_modele;

  Operateur_Diff terme_diffusif;

  Champ_Inc irradiance_;
  SolveurSys solveur;
  Matrice_Morse la_matrice;

  Operateur_Grad gradient;


};


// Description:
//    Renvoie le modele de rayonnement semi transparent
// Precondition :
// Parametre :
//    Signification :
//    Contraintes :
//    Acces :
// Retour
//    Signification :
//    Contraintes :
// Exception :
// Effets de bord :
// Postcondition :
// Postcondition :
inline Modele_rayo_semi_transp& Equation_rayonnement_base::Modele()
{
  return le_modele.valeur();
}


// Description:
//    Renvoie le modele de rayonnement semi transparent
// Precondition :
// Parametre :
//    Signification :
//    Contraintes :
//    Acces :
// Retour
//    Signification :
//    Contraintes : reference constante
// Exception :
// Effets de bord :
// Postcondition :
// Postcondition :
inline const Modele_rayo_semi_transp& Equation_rayonnement_base::Modele() const
{
  return le_modele.valeur();
}


// Description:
//    Associe un fluide incompressible semi transparent a l'equation.
// Precondition:
// Parametre: Fluide_base& un_fluide
//    Signification: le fluide incompressible semi transparent a associer
//    Valeurs par defaut:
//    Contraintes: reference constante
//    Acces:
// Retour:
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
// Postcondition: l'equation a un fluide associe
inline void Equation_rayonnement_base::associer_fluide(const Fluide_base& un_fluide)
{
  le_fluide = un_fluide;
}


// Description:
//    renvoie le fluide semi transparent associe a l'equation de rayonnement
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: le fluide associe a l'equation de rayonnement
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
inline const Fluide_base& Equation_rayonnement_base::fluide() const
{
  return le_fluide.valeur();
}


// Description:
//    renvoie le fluide semi transparent associe a l'equation de rayonnement
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: le fluide associe a l'equation de rayonnement
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
inline Fluide_base& Equation_rayonnement_base::fluide()
{
  return le_fluide.valeur();
}


#endif
