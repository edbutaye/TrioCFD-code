/****************************************************************************
* Copyright (c) 2023, CEA
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
// File:        Expert_mode_K_Omega.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Expert_mode_K_Omega.h>
#include <Param.h>
#include <K_Omega_Eps_constants.h>

Implemente_instanciable_sans_constructeur(Expert_mode_K_Omega,"Expert_mode_K_Omega",Objet_U);

Sortie& Expert_mode_K_Omega::printOn(Sortie& os) const
{
  Cerr << "Error : ::printOn is not implemented." << finl;
  return os;
}

Entree& Expert_mode_K_Omega::readOn(Entree& is)
{
  Param p(que_suis_je());
  set_param(p);
  p.lire_avec_accolades_depuis(is);
  return is;
}

void Expert_mode_K_Omega::set_param(Param& param)
{
  param.ajouter_flag("deactivate_production_limiter",	&deactivate_production_limiter_); // XD_ADD_P bool Deactivate the production limiter in the k equation (default value false).
  param.ajouter_non_std("menter_version", (this)); // XD_ADD_P motcle Version of Menter SST model: "ORIGINAL_1994" or "MODIFIED_2003" (default value "ORIGINAL_1994").
}

int Expert_mode_K_Omega::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="menter_version")
    {
      Motcle menter_version_str;
      is >> menter_version_str;
      if (menter_version_str == "ORIGINAL_1994")
        menter_version_ = Menter_version::ORIGINAL_1994;
      else if (menter_version_str == "MODIFIED_2003")
        menter_version_ = Menter_version::MODIFIED_2003;
      else
        {
          Cerr << "Error: unknown Menter version: " << menter_version_str << finl;
          Cerr
              << "Available versions are: ORIGINAL_1994 and MODIFIED_2003."
              << finl;
          Process::exit();
        }
      return 1;
    }
  return 1;
}

const double& Expert_mode_K_Omega::get_gamma1() const
{
  if (menter_version_ == Menter_version::ORIGINAL_1994)
    return GAMMA1_1994;
  else
    return GAMMA1_2003;
}
const double& Expert_mode_K_Omega::get_gamma2() const
{
  if (menter_version_ == Menter_version::ORIGINAL_1994)
    return GAMMA2_1994;
  else
    return GAMMA2_2003;
}
const double& Expert_mode_K_Omega::get_cst_production_limiter() const
{
  if (menter_version_ == Menter_version::ORIGINAL_1994)
    return CST_PRODUCTION_LIMITER_1994;
  else
    return CST_PRODUCTION_LIMITER_2003;
}
const double& Expert_mode_K_Omega::get_cst_cd_komega() const
{
  if (menter_version_ == Menter_version::ORIGINAL_1994)
    return CST_MIN_CD_KOMEGA_1994;
  else
    return CST_MIN_CD_KOMEGA_2003;
}
