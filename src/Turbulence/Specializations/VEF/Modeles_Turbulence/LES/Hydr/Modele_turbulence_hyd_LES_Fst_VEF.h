/****************************************************************************
* Copyright (c) 2024, CEA
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
// File:        Modele_turbulence_hyd_LES_Fst_VEF.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Modeles_Turbulence/LES/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_turbulence_hyd_LES_Fst_VEF_included
#define Modele_turbulence_hyd_LES_Fst_VEF_included

#include <Modele_turbulence_hyd_LES_VEF_base.h>

/*! @brief classe Modele_turbulence_hyd_LES_Fst_VEF Cette classe correspond a la mise en oeuvre du modele sous
 *
 *  maille de fonction de structure en VEF, en utilisant la formulation
 *  nu_t = 0.02587*delta_x^2*sqrt(2*Sij*Sij+wi*wi)
 *
 *  .SECTION  voir aussi
 *  Modele_turbulence_hyd_LES_base
 *
 */
class Modele_turbulence_hyd_LES_Fst_VEF: public Modele_turbulence_hyd_LES_VEF_base
{
  Declare_instanciable_sans_constructeur(Modele_turbulence_hyd_LES_Fst_VEF);
public:
  Modele_turbulence_hyd_LES_Fst_VEF();

protected:
  Champ_Fonc& calculer_viscosite_turbulente() override;
  virtual void calculer_racine();

  double C1_;
  DoubleVect Racine_;
};

#endif /* Modele_turbulence_hyd_LES_Fst_VEF_included */