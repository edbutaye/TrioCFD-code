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
/*! @brief class Op_Diff_K_Eps_VEF_Face Cette classe represente l'operateur de diffusion turbulente
 *
 *    La discretisation est VEF
 *   Les methodes pour l'implicite sont codees.
 *
 */
#ifndef Op_Diff_K_Eps_VEF_base_included
#define Op_Diff_K_Eps_VEF_base_included

#include <Op_Dift_VEF_base.h>
#include <Op_VEF_Face.h>
#include <TRUST_Ref.h>
#include <K_Omega_Eps_constants.h>

//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Op_Diff_K_Eps_VEF_base
//
//////////////////////////////////////////////////////////////////////////////

class Op_Diff_K_Eps_VEF_base : public Op_Dift_VEF_base
{

  Declare_instanciable_sans_constructeur(Op_Diff_K_Eps_VEF_base);

  Op_Diff_K_Eps_VEF_base(double Sigma_K = SIGMA_K_KEPS, double Sigma_Eps = SIGMA_EPS_KEPS ) : Sigma_K_(Sigma_K) , Sigma_Eps_(Sigma_Eps)
  {
    Sigma_[0]=Sigma_K;
    Sigma_[1]=Sigma_Eps;
  }

public:

  void completer() override;

protected:
  double Sigma_K_;
  double Sigma_Eps_;
  double Sigma_[2];

};

#endif
