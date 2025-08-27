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


#ifndef Flux_parietal_bilineaire_included
#define Flux_parietal_bilineaire_included


#include <Flux_parietal_base.h>
#include <Saturation_base.h>
/*! @brief : class Flux_parietal_bilineaire
 *
 * <Description of class Flux_parietal_bilineaire>
 */

class Flux_parietal_bilineaire : public Flux_parietal_base
{

  Declare_instanciable( Flux_parietal_bilineaire ) ;
public:
  void qp(const input_t& input, output_t& output) const override;
  int T_at_wall() const override { return 0; } // Use T cell center
  double fraction_flux_vap(const input_t& in, double phi, double dTf_phi, double dp_phi, double dTp_phi, double& dTf_eps, double& dp_eps, double& dTp_eps) const;

private:
  double h_mono_ = -123., h_diph_ = -123.; // Coefficient d'échange thermique dans la zone monophasique /zone diphasique /
  double coeff_ = -123.; // Paramètre d'osv
  double b_ = 10000.0; // abscisse à l'origine de la fonction dans la zone diphasique, 10kW/m2 de base pour l'héritage du stage
  int n_l = -1, n_g = -1; // indice des phases d'injection du flux
  const Saturation_base *sat = nullptr; //saturation entre liquide et gas (si elle existe)

};

#endif /* Flux_parietal_bilineaire_included */
