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
// File:        Frontiere_ouverte_temperature_imposee_rayo_semi_transp.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src
// Version:     /main/10
//
//////////////////////////////////////////////////////////////////////////////

#include <Frontiere_ouverte_temperature_imposee_rayo_semi_transp.h>

Implemente_instanciable(Frontiere_ouverte_temperature_imposee_rayo_semi_transp,"Frontiere_ouverte_temperature_imposee_rayo_semi_transp",Entree_fluide_temperature_imposee);

/*! @brief
 *
 * @param (Sortie& os) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Frontiere_ouverte_temperature_imposee_rayo_semi_transp::printOn(Sortie& os) const
{
  return os;
}

/*! @brief Lecture des parametres de la condition Neumann_paroi
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Frontiere_ouverte_temperature_imposee_rayo_semi_transp::readOn(Entree& is)
{
  return Entree_fluide_temperature_imposee::readOn(is);
}

const Cond_lim_base& Frontiere_ouverte_temperature_imposee_rayo_semi_transp::la_cl() const
{
  return (*this);
}

void Frontiere_ouverte_temperature_imposee_rayo_semi_transp::calculer_temperature_bord(double temps)
{
  // Cette methode ne fait rien car la temperature de paroi est
  // directement donnee par T_ext()
  ;
}

void Frontiere_ouverte_temperature_imposee_rayo_semi_transp::completer()
{
  Entree_fluide_temperature_imposee::completer();
}
