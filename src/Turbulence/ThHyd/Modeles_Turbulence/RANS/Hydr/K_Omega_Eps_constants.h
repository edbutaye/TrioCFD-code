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

#ifndef K_Omega_Eps_Constants
#define K_Omega_Eps_Constants

#include <cmath>

// Constants for the classic k-omega model Wilcox 1988
// See https://turbmodels.larc.nasa.gov/wilcox.html
static constexpr double BETA_K = 0.09; // Cmu or BETA_STAR, but clearer with _K
static constexpr double BETA_OMEGA = 3./40.; // BETA
static constexpr double PRANDTL_K = 1./2.; // SIGMA_STAR 1/Prandtl_K_ = sigma_k
static constexpr double PRANDTL_OMEGA = 1./2.; // SIGMA
static constexpr double SIGMA_K = PRANDTL_K;
static constexpr double SIGMA_OMEGA = PRANDTL_OMEGA;
static constexpr double ALPHA_OMEGA = 5./9.; // ALPHA

// Constants for the k-omega SST model
// See https://turbmodels.larc.nasa.gov/sst.html
static constexpr double SIGMA_K1 = 0.85;
static constexpr double SIGMA_K2 = 1.0;
static constexpr double SIGMA_OMEGA1 = 0.5;
static constexpr double SIGMA_OMEGA2 = 0.856;
static constexpr double BETA1 = 0.075;
static constexpr double BETA2 = 0.0828;
static constexpr double KAPPA = 0.41;
static const double GAMMA1_1994 = BETA1/BETA_K - SIGMA_OMEGA1*KAPPA*KAPPA/sqrt(BETA_K);
static const double GAMMA2_1994 = BETA2/BETA_K - SIGMA_OMEGA2*KAPPA*KAPPA/sqrt(BETA_K);
static constexpr double GAMMA1_2003 = 5./9.;// 5./9. -> Menter_2003. old : BETA1/BETA_K - SIGMA_OMEGA1*KAPPA*KAPPA/sqrt(BETA_K);
static constexpr double GAMMA2_2003 = 0.44;// 0.44 -> Menter_2003. old :BETA2/BETA_K - SIGMA_OMEGA2*KAPPA*KAPPA/sqrt(BETA_K);
static constexpr double CST_PRODUCTION_LIMITER_1994=20.;
static constexpr double CST_PRODUCTION_LIMITER_2003=10.;
static constexpr double CST_MIN_CD_KOMEGA_1994=1.e-20;
static constexpr double CST_MIN_CD_KOMEGA_2003=1.e-10;

// Constants for the k-epsilon model
static constexpr double SIGMA_K_KEPS = 1.0;
static constexpr double SIGMA_EPS_KEPS = 1.3;

#endif
