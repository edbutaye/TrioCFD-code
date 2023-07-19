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
/////////////////////////////////////////////////////////////////////////////
//
// File      : IJK_Thermal.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Thermal.h>

Implemente_instanciable( IJK_Thermal, "IJK_Thermal", DERIV(IJK_Thermal_base) ) ;

Sortie& IJK_Thermal::printOn( Sortie& os ) const
{
  DERIV(IJK_Thermal_base)::printOn( os );
  return os;
}

Entree& IJK_Thermal::readOn( Entree& is )
{
  Cerr<<"Read and Cast IJK_Thermal :"<<finl;
  Motcle word;
  is>>word;
  Nom type = "IJK_Thermal_";
  Motcles thermal_words(3);
  {
    thermal_words[0] = "subresolution";
    thermal_words[1] = "multiplesubresolutions";
    thermal_words[2] = "onefluid";
  }
  int thermal_rank = thermal_words.search(word);
  switch(thermal_rank)
    {
    case 0 :
      {
        type += "Subresolution";
        break;
      }
    case 1 :
      {
        type += "Multiple_Subresolutions";
        break;
      }
    case 2 :
      {
        type += "Onefluid";
        break;
      }
    default :
      {
        Cerr<<"ERROR : Thermal problems already implemented are:"<<finl;
        Cerr<<thermal_words<<finl;
        abort();
      }
    }
  typer(type);
  // valeur().readOn(is);
  is >> valeur();
  return is;
  //  thermal_reader_->readOn(is);
  //  Cerr << thermal_reader_->set_type_thermal_problem() << finl;
  //  typer_problem();
  //  DERIV(IJK_Thermal_base)::readOn( is );
}
