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
// File:        Echange_externe_impose_rayo_semi_transp.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src/VDF
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////


#include <Echange_externe_impose_rayo_semi_transp.h>
#include <Motcle.h>

Implemente_instanciable(Echange_externe_impose_rayo_semi_transp,"Paroi_Echange_externe_impose_rayo_semi_transp",Echange_externe_impose);



/*! @brief
 *
 * @param (Sortie& os) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Echange_externe_impose_rayo_semi_transp::printOn(Sortie& os) const
{
  return os;
}


/*! @brief Lecture des parametres de la condition Echange_impose Lecture de l'emissivite de la paroi
 *
 *     Lecture du coefficient A
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Echange_externe_impose_rayo_semi_transp::readOn(Entree& is)
{
  Motcle motlu;
  Motcles les_motcles(2);
  {
    les_motcles[0] = "h_imp";
    les_motcles[1] = "T_ext";
  }

  int ind = 0;
  while (ind < 2)
    {
      is >> motlu;
      int rang = les_motcles.search(motlu);

      switch(rang)
        {
        case 0:
          {
            is >> h_imp_;
            break;
          }
        case 1:
          {
            is >> T_ext();
            break;
          }
        default:
          {
            Cerr << "Erreur a la lecture de la condition aux limites de type "<<finl;
            Cerr << "Echange_externe_impose_rayo_semi_transp " << finl;
            Cerr << "On attendait " << les_motcles << "a la place de " <<  motlu << finl;
            exit();
          }
        }
      ind++;
    }

  if(local_min_vect(h_imp_->valeurs()) < 1.e9)
    {
      Cerr<<"Erreur sur l'utilisation de la condition a la limite"<<finl;
      Cerr<<"Echange_externe_impose_rayo_semi_transp. Celle ci ne peut"<<finl;
      Cerr<<"etre utilisee pour un probleme de rayonnement que pour "<<finl;
      Cerr<<"imposer une temperature sur une paroi"<<finl;
      exit();
    }

  return is;
}



void Echange_externe_impose_rayo_semi_transp::completer()
{
  Echange_externe_impose::completer();
}



const Cond_lim_base& Echange_externe_impose_rayo_semi_transp::la_cl() const
{
  return (*this);
}


Champ_front_base& Echange_externe_impose_rayo_semi_transp::temperature_bord()
{
  return T_ext();
}

void Echange_externe_impose_rayo_semi_transp::calculer_temperature_bord(double temps)
{
  // La temperature de paroi etant directement donnee par le champ_front
  // T_ext, il n'y a rien a calculer ici
  ;
}
