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
// File:        Equation_rayonnement.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src
// Version:     /main/7
//
//////////////////////////////////////////////////////////////////////////////
//
//
//
//

#include <Equation_rayonnement.h>
#include <Discretisation_base.h>

Implemente_deriv(Equation_rayonnement_base);
Implemente_instanciable(Equation_rayonnement,"Equation_rayonnement",DERIV(Equation_rayonnement_base));


Sortie& Equation_rayonnement::printOn(Sortie& os) const
{
  return DERIV(Equation_rayonnement_base)::printOn(os);
}

Entree& Equation_rayonnement::readOn(Entree& is)
{
  return DERIV(Equation_rayonnement_base)::readOn(is);
}

void Equation_rayonnement::typer(Nom& type, const Equation_base& eqn)
{
  type = "Eq_rayo_semi_transp_";
  Nom disc = eqn.discretisation().que_suis_je();
  // les operateurs des sources sont communs aux discretisations VEF et VEFP1B
  if(disc=="VEFPreP1B")
    disc="VEF";

  type+=disc;

  Cerr << type << finl;
  DERIV(Equation_rayonnement_base)::typer(type);
}
