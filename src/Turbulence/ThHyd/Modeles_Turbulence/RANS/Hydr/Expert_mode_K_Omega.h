/****************************************************************************
* Copyright (c) 2025, CEA
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
// File:        Expert_mode_K_Omega.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Expert_mode_K_Omega_included
#define Expert_mode_K_Omega_included
#include <Param.h>

enum class Menter_version_ { ORIGINAL_1994, MODIFIED_2003 };

class Expert_mode_K_Omega: public Objet_U
{
  Declare_instanciable_sans_constructeur(Expert_mode_K_Omega);
public:
  Expert_mode_K_Omega()
  {
    deactivate_production_limiter_=false;
    menter_version_=Menter_version_::ORIGINAL_1994;
  }

  void set_param(Param& param);
  int lire_motcle_non_standard(const Motcle& mot, Entree& is ) override;
  const Menter_version_& get_menter_version() const { return menter_version_; }
  const bool& get_deactivate_production_limiter() const { return deactivate_production_limiter_; }
  const double get_gamma1() const;
  const double get_gamma2() const;

protected:
  bool deactivate_production_limiter_;
  Menter_version_ menter_version_;
};
#endif
