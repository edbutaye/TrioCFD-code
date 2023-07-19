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
// File      : IJK_Thermal.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_Thermal_included
#define IJK_Thermal_included

#include <IJK_Thermal_base.h>
#include <IJK_Thermal_Reader.h>
#include <TRUST_Deriv.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_Thermal
//
// <Description of class IJK_Thermal>
//
/////////////////////////////////////////////////////////////////////////////

class IJK_Thermal : public DERIV(IJK_Thermal_base)
{
  Declare_instanciable( IJK_Thermal );

public :
//  IJK_Thermal_Reader* thermal_reader_;
//  inline void set_thermal_reader(IJK_Thermal_Reader& thermal_reader);
//  inline void typer_problem();
//  inline void typer_problem(std::string string_type);
  inline int initialize(const IJK_Splitting& splitting, const int idx);
  inline void update_thermal_properties();
  inline void euler_time_step(const double timestep);
  inline void associer(const IJK_FT_double& ijk_ft);
  inline void sauvegarder_temperature(Nom& lata_name, int idx);
  inline double compute_timestep(const double timestep,
                                 const double dxmin) const;

};

//inline void IJK_Thermal::set_thermal_reader(IJK_Thermal_Reader& thermal_reader)
//{
//  thermal_reader_ = &thermal_reader;
//}
//
//inline void IJK_Thermal::typer_problem(std::string string_type)
//{
//  typer(string_type.c_str());
//}
//
//inline void IJK_Thermal::typer_problem()
//{
////  Cout << "Test" << valeur().get_type_thermal_problem() << finl;
////  Cerr << "Test" << valeur().get_type_thermal_problem() << finl;
//  typer(thermal_reader_->set_type_thermal_problem().c_str());
//}

inline int IJK_Thermal::initialize(const IJK_Splitting& splitting, const int idx)
{
  return valeur().initialize(splitting, idx);
}

inline void IJK_Thermal::update_thermal_properties()
{
  return valeur().update_thermal_properties();
}

inline void IJK_Thermal::euler_time_step(const double timestep)
{
  valeur().euler_time_step(timestep);
}

inline void IJK_Thermal::associer(const IJK_FT_double& ijk_ft)
{
  valeur().associer(ijk_ft);
}

inline void IJK_Thermal::sauvegarder_temperature(Nom& lata_name, int idx)
{
  valeur().sauvegarder_temperature(lata_name, idx);
}

inline double IJK_Thermal::compute_timestep(const double timestep, const double dxmin) const
{
  return valeur().compute_timestep(timestep, dxmin);
}

#endif /* IJK_Thermal_included */
