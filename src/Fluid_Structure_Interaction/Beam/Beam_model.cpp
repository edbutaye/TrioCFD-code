/****************************************************************************
* Copyright (c) 2021, CEA
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
// File      : Beam_model.cpp
// Directory : $BEAM_MODEL_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Beam_model.h>
#include <Domaine_ALE.h>
#include <TRUSTVects.h>
#include <TRUSTTabs.h>
#include <Probleme_base.h>
#include <Domaine_VEF.h>
#include <Domaine_Cl_VEF.h>
#include <Frontiere.h>
#include <MD_Vector_tools.h>
#include <MD_Vector_std.h>
#include <Faces.h>
#include <fstream>
#include <iostream>
#include <math.h>       /* sin */
#include <TRUST_Ref.h>
#include <Domaine.h>

#define PI 3.14159265

using namespace std;


Implemente_instanciable_sans_constructeur_ni_destructeur(Beam_model, "Beam_model", Interprete_geometrique_base ) ;

// Syntaxe:
//  Beam_model NOMDOMAINE {
//    nb_beam
//    Name
//     nb_modes number of modes
//     longitudinal_axis x|y|z
//     bendingDirection x|y|z
//     nb_planes 1|2
//     NewmarkTimeScheme MA|FD|HHT
//     Mass_and_stiffness_file_name
//     Absc_file_name
//     Modal_deformation_file_name nb_modes files
//     [ Young_Module ]
//     { Rho_beam ]
//     [ BaseCenterCoordinates position  ]
//     [ CI_file_name ]
//     [ Output_position_1D nb_points  position ]
//     [ Output_position_3D nb_points  position ]
//     [ Restart_file_name file ]
//
//  }
// XD Beam_model interprete Beam_model 0 Reduced mechanical model: a beam model. Resolution based on a modal analysis. Temporal discretization: Newmark or Hilber-Hughes-Taylor (HHT)
// XD  attr dom ref_domaine dom 0 domain name
// XD  attr bloc bloc_lecture_beam_model bloc 0 not_set

// XD bloc_lecture_beam_model objet_lecture nul 0 bloc
// XD  attr aco chaine(into=["{"]) aco 0 Opening curly bracket.
// XD  attr nb_beam chaine(into=["nb_beam"]) nb_beam 0 Keyword to specify the number of beams
// XD  attr nb_beam_val entier nb_beam_val 0 Number of beams
// XD  attr Name chaine(into=["name"]) BeamName 0 keyword to specify the Name of the beam (the name must match with the name of the edge in the fluid domain)
// XD  attr Name_of_beam chaine Name_of_beam 0 keyword to specify the Name of the beam (the name must match with the name of the edge in the fluid domain)
// XD  attr bloc bloc_poutre bloc 0 not_set
// XD  attr Name2 chaine(into=["name"]) BeamName2 1 keyword to specify the Name of the beam (the name must match with the name of the edge in the fluid domain)
// XD  attr Name_of_beam2 chaine Name_of_beam2 1 keyword to specify the Name of the beam (the name must match with the name of the edge in the fluid domain)
// XD  attr bloc2 bloc_poutre bloc2 1 not_set
// XD  attr acof chaine(into=["}"]) acof 0 Closing curly bracket.

// XD NewmarkTimeScheme_deriv objet_lecture NewmarkTimeScheme_deriv -1 Solve the beam dynamics. Selection of time integration scheme.
// XD NewmarkTimeScheme_HHT NewmarkTimeScheme_deriv HHT 0 HHT alpha (Hilber-Hughes-Taylor, alpha usually -0.1 ) time integration scheme.
// XD  attr alpha floattant alpha 1 usually, alpha is set to -0.1
// XD NewmarkTimeScheme_MA NewmarkTimeScheme_deriv MA 0  MA (Newmark mean acceleration) time integration scheme.
// XD NewmarkTimeScheme_FD NewmarkTimeScheme_deriv FD 0 FD (Newmark finite differences) time integration scheme. Warning: this scheme is conditionally stable. The time step should satisfy the corresponding stability constraint, but this implementation does not automatically enforce it.The Newmark finite difference scheme is retained primarily for advanced users and benchmarking purposes.

// XD bloc_poutre objet_lecture nul 1 Read poutre bloc
// XD  attr nb_modes entier n 0 Number of modes. If there are two planes, indicate the total for both planes.
// XD  attr longitudinal_axis chaine dir 0 x, y, z. Axis along the length of the beam
// XD  attr bendingDirection chaine dir_bending 0 x, y, z . Direction of  bending. If two planes, give the first one, the second direction is determine automatically
// XD  attr nb_planes entier nplanes 0 Number of planes used in the beam dynamic model
// XD  attr NewmarkTimeScheme NewmarkTimeScheme_deriv NewmarkTimeScheme 0 Solve the beam dynamics. Time integration scheme: choice between MA (Newmark mean acceleration),  FD (Newmark finite differences), and HHT alpha (Hilber-Hughes-Taylor, alpha usually -0.1 )
// XD  attr Mass_and_stiffness_file_name chaine  Mass_and_stiffness_file_name 0 Name of the file containing the diagonal modal mass, stiffness, and damping matrices. The file must contain nb_modes lines. When several bending planes are defined, the modes must be ordered plane-by-plane. The first modes_per_plane lines correspond to plane 0, the next modes_per_plane lines correspond to plane 1.
// XD  attr Absc_file_name chaine Absc_file_name 0 Name of the file containing the coordinates of the Beam
// XD  attr Modal_deformation_file_name listchaine  Modal_deformation_file_name 0 Name of the file containing the modal deformation of the Beam (Two-column modal deformation file: W(X) modal displacement and the rotation dW/dX modal rotation).
// XD  attr Young_Module floattant young 1 Young Module
// XD  attr Rho_beam floattant rho 1 Beam density
// XD  attr BaseCenterCoordinates listf pos_center 1 position of the base center coordinates on the Beam
// XD  attr CI_file_name chaine CI_file_name 1 Name of the file containing the initial condition of the Beam. The initial condition file must contain nb_modes displacement values, one per line. When several bending planes are defined, the modes must be ordered plane-by-plane. The first modes_per_plane lines correspond to plane 0, the next modes_per_plane lines correspond to plane 1.
// XD  attr Absc_file_name chaine Absc_file_name 0 Name of the file containing the coordinates of the Beam
// XD  attr Restart_file_name chaine Restart_file_name 1 SaveBeamForRestart.txt file to restart the calculation. The file must contain nb_modes displacement values, one per line. When several bending planes are defined, the modes must be ordered plane-by-plane. The first modes_per_plane lines correspond to plane 0, the next modes_per_plane lines correspond to plane 1.
// XD  attr Absc_file_name chaine Absc_file_name 0 Name of the file containing the coordinates of the Beam
// XD  attr Output_position_1D list pt1d 1 nb_points  position Post-traitement of specific points on the Beam
// XD  attr Output_position_3D listpoints pt3d 1 nb_points  position Post-traitement of specific points on the 3d FSI boundary

Beam_model::Beam_model()
{

  nbModes_=0;
  nb_planes_=1;
  longitudinal_axis_=-1;
  bending_dir_.resize(1);
  bending_dir_=-1;
  young_=200.e+9;
  rho_ = 8100.;
  alpha_= 0.;
  beta_= 0.25;
  gamma_=0.5;
  temps_ =0.;
  dt_stab_ = 1.e+8;
  output_position_1D_.resize(0);
  output_position_1D_=0.;
  output_position_3D_.resize(0);
  output_position_3D_=0.;
  fluidForceOnBeam_.resize(0);
  fluidForceOnBeam_=0.;
  tempsComputeForceOnBeam_=0.;
  x0_=0.;
  y0_=0.;
  z0_=0.;
  mass_.reset();
  stiffness_.reset();
  damping_.reset();
  abscissa_.reset();
  qSpeed_.reset();
  qAcceleration_.reset();
  qDisplacement_.reset();
  output_position_1D_.reset();
  output_position_3D_.reset();
  beamName_= Nom();
}
Beam_model::~Beam_model()
{

}

Sortie& Beam_model::printOn( Sortie& os ) const
{
  return Interprete::printOn(os);
}

Entree& Beam_model::readOn( Entree& is )
{
  return Interprete::readOn(is);
}

Entree&  Beam_model::interpreter_(Entree& is)
{
  associer_domaine(is);
  Domaine_ALE& dom=ref_cast(Domaine_ALE, domaine());
  dom.reading_beam_model(is);
  return is;
}

void Beam_model::readInputMassStiffnessFiles(Nom& masse_and_stiffness_file_name)
{
  Cerr << "Reading beam model coefficients from file: " << masse_and_stiffness_file_name << finl;


  if (nbModes_ % nb_planes_ != 0)
    {
      Cerr << "ERROR: invalid configuration. The number of modes (" << nbModes_
           << ") must be divisible by the number of planes (" << nb_planes_ << ")." << finl;
      Cerr << "Each plane must have the same number of modes." << finl;
      Process::exit();
    }

  int modes_per_plane= nbModes_/nb_planes_;

  mass_.resize(modes_per_plane, nb_planes_);
  stiffness_.resize(modes_per_plane, nb_planes_);
  damping_.resize(modes_per_plane, nb_planes_);

  const std::string filename(masse_and_stiffness_file_name);
  std::ifstream monFlux(filename.c_str());

  if (!monFlux)
    {
      Cerr << "ERROR: Unable to open file '" << filename << "'." << finl;
      Process::exit();
    }

  std::string line;
  int index = 0;
  int lineNumber = 0;

  while (std::getline(monFlux, line))
    {
      ++lineNumber;

      // Skip empty or whitespace-only lines
      bool empty = true;
      for (char c : line)
        {
          if (!std::isspace(static_cast<unsigned char>(c))) { empty = false; break; }
        }
      if (empty) continue;

      // Skip comment lines starting with '#' or header lines (first line non-numeric)
      std::istringstream issCheck(line);
      double dummy;
      if (!(issCheck >> dummy)) continue;

      if (index >= nbModes_)
        {
          Cerr << "ERROR: File '" << filename << "' contains more data lines than nb_modes = "
               << nbModes_ << " (at line " << lineNumber << ")." << finl;
          Process::exit();
        }

      // Parse exactly three values
      std::istringstream iss(line);
      double mass, stiffness, damping;

      if (!(iss >> mass >> stiffness >> damping))
        {
          Cerr << "ERROR: In file '" << filename << "', line " << lineNumber
               << " does not contain exactly three numerical values." << finl;
          Cerr << "Line content: [" << line << "]" << finl;
          Process::exit();
        }

      // Check if a fourth value exists
      double extra;
      if (iss >> extra)
        {
          Cerr << "ERROR: In file '" << filename << "', line " << lineNumber
               << " contains more than three values." << finl;
          Cerr << "Line content: [" << line << "]" << finl;
          Process::exit();
        }

      if (std::abs(mass) <= 1.e-16 || std::abs(stiffness) <= 1.e-16)
        {
          Cerr << "ERROR: Invalid mass or stiffness in file '" << filename
               << "' at line " << lineNumber << "." << finl;
          Cerr << "mass = " << mass << ", stiffness = " << stiffness
               << ", damping = " << damping << finl;
          Process::exit();
        }
      // Store values
      int p = index / modes_per_plane;      // plane index
      int m = index % modes_per_plane;      // mode index
      mass_(m, p) = mass;
      stiffness_(m, p) = stiffness;
      damping_(m, p) = damping;
      ++index;
    }

  if (index < nbModes_)
    {
      Cerr << "ERROR: File '" << filename << "' contains too few data lines. "
           << "Expected nb_modes = " << nbModes_ << ", but only " << index << " valid lines were read." << finl;
      Process::exit();
    }

  Cerr << "Beam coefficients successfully read from '" << filename
       << "' (" << index << " modes)." << finl;
}

void Beam_model::readInputCIFile(Nom& CI_file_name)
{
  Cerr << "Reading initial beam conditions from file: " << CI_file_name << finl;

  if (nbModes_ % nb_planes_ != 0)
    {
      Cerr << "ERROR: invalid configuration. nb_modes (" << nbModes_
           << ") must be divisible by nb_planes (" << nb_planes_ << ")." << finl;
      Process::exit();
    }

  int modes_per_plane = nbModes_ / nb_planes_;

  qSpeed_.resize(modes_per_plane, nb_planes_);
  qAcceleration_.resize(modes_per_plane, nb_planes_);
  qDisplacement_.resize(modes_per_plane, nb_planes_);
  fluidForceOnBeam_.resize(modes_per_plane, nb_planes_);

  qSpeed_ = 0.;
  qAcceleration_ = 0.;
  qDisplacement_ = 0.;
  fluidForceOnBeam_ = 0.;

  const std::string filename(CI_file_name);
  std::ifstream monFlux(filename.c_str());

  if (!monFlux)
    {
      Cerr << "ERROR: Unable to open file '" << filename << "'." << finl;
      Process::exit();
    }

  std::string line;
  int index = 0;
  int lineNumber = 0;

  while (std::getline(monFlux, line))
    {
      ++lineNumber;

      // Skip empty or whitespace-only
      bool empty = true;
      for (char c : line)
        if (!std::isspace(static_cast<unsigned char>(c))) { empty = false; break; }
      if (empty) continue;

      // Skip comments
      if (line[0] == '#') continue;

      // Skip non-numeric header lines
      std::istringstream issCheck(line);
      double dummy;
      if (!(issCheck >> dummy)) continue;

      if (index >= nbModes_)
        {
          Cerr << "ERROR: File '" << filename
               << "' contains more data than expected (nb_modes = " << nbModes_
               << "), at line " << lineNumber << "." << finl;
          Process::exit();
        }

      std::istringstream iss(line);
      double displacement;

      if (!(iss >> displacement))
        {
          Cerr << "ERROR: Invalid displacement value at line "
               << lineNumber << " in file '" << filename << "'." << finl;
          Cerr << "Line content: [" << line << "]" << finl;
          Process::exit();
        }

      // Ensure only one value per line
      double extra;
      if (iss >> extra)
        {
          Cerr << "ERROR: Too many values at line " << lineNumber
               << " in file '" << filename << "'." << finl;
          Cerr << "Line content: [" << line << "]" << finl;
          Process::exit();
        }

      int p = index / modes_per_plane;  // plane index
      int m = index % modes_per_plane;  // mode index

      qDisplacement_(m, p) = displacement;

      // Compute initial acceleration
      qAcceleration_(m, p) =
        fluidForceOnBeam_(m, p) - (stiffness_(m, p) / mass_(m, p)) * displacement;

      ++index;
    }

  // Check file length
  if (index < nbModes_)
    {
      Cerr << "ERROR: File '" << filename << "' contains too few data lines. "
           << "Expected " << nbModes_ << ", but found " << index << "." << finl;
      Process::exit();
    }

  Cerr << "Initial beam conditions successfully read from '" << filename
       << "' (" << index << " values)." << finl;
}


void Beam_model::readRestartFile(Nom& Restart_file_name)
{
  Cerr << "Reading restart file: " << Restart_file_name << finl;

  int modes_per_plane = nbModes_ / nb_planes_;

  // Resize 2D arrays
  qDisplacement_.resize(modes_per_plane, nb_planes_);
  qSpeed_.resize(modes_per_plane, nb_planes_);
  qAcceleration_.resize(modes_per_plane, nb_planes_);

  const std::string filename(Restart_file_name);
  std::ifstream monFlux(filename.c_str());

  if (!monFlux)
    {
      Cerr << "ERROR: Unable to open restart file '" << filename << "'." << finl;
      Process::exit();
    }

  std::string line;
  int index = 0;
  int lineNumber = 0;

  while (std::getline(monFlux, line))
    {
      ++lineNumber;


      // Read exactly 4 values
      std::istringstream iss(line);
      double temps, displacement, speed, acceleration;

      if (!(iss >> temps >> displacement >> speed >> acceleration))
        {
          Cerr << "ERROR: Line " << lineNumber << " in file '" << filename
               << "' does not contain exactly 4 numerical values." << finl;
          Cerr << "Line content: [" << line << "]" << finl;
          Process::exit();
        }

      // Check for extra values
      double extra;
      if (iss >> extra)
        {
          Cerr << "ERROR: Line " << lineNumber << " in file '" << filename
               << "' contains more than 4 values." << finl;
          Cerr << "Line content: [" << line << "]" << finl;
          Process::exit();
        }

      int plane = index / modes_per_plane;
      int mode = index % modes_per_plane;

      temps_ = temps;
      qDisplacement_(mode, plane) = displacement;
      qSpeed_(mode, plane) = speed;
      qAcceleration_(mode, plane) = acceleration;

      ++index;
    }

  // Check that all modes were read
  if (index < nbModes_)
    {
      Cerr << "ERROR: File '" << filename << "' contains too few valid data lines. "
           << "Expected nb_modes = " << nbModes_ << ", but only " << index << " were read." << finl;
      Process::exit();
    }

  tempsComputeForceOnBeam_ = temps_;
  Cerr << "Restart file successfully read from '" << filename
       << "' (" << index << " modes, " << nb_planes_ << " planes)." << finl;
}

void Beam_model::readInputAbscFiles(Nom& absc_file_name)
{
  Cerr << "Reading beam coordinates from file: " << absc_file_name << finl;

  const std::string filename(absc_file_name);
  std::ifstream monFlux(filename.c_str());
  if (!monFlux)
    {
      Cerr << "ERROR: Unable to open file '" << filename << "'." << finl;
      Process::exit();
    }

  // ===== First pass: validate lines and count valid numeric values =====
  std::string line;
  int lineNumber = 0;
  int nValues = 0;

  while (std::getline(monFlux, line))
    {
      ++lineNumber;

      // Skip empty lines
      if (line.find_first_not_of(" \t\n\r") == std::string::npos) continue;

      // Skip comment lines
      if (line[0] == '#') continue;

      // Skip non-numeric headers
      std::istringstream iss(line);
      double absc;
      if (!(iss >> absc)) continue;

      // Check that there’s exactly one value
      double extra;
      if (iss >> extra)
        {
          Cerr << "ERROR: In file '" << filename << "', line " << lineNumber
               << " contains more than one value." << finl;
          Cerr << "Line content: [" << line << "]" << finl;
          Process::exit();
        }

      ++nValues; // count valid numeric lines
    }

  if (nValues == 0)
    {
      Cerr << "ERROR: No valid numeric values found in file '" << filename << "'." << finl;
      Process::exit();
    }

  // ===== Resize the DoubleVect to hold all values =====
  abscissa_.resize(nValues);

  // ===== Second pass: read and store values =====
  monFlux.clear();
  monFlux.seekg(0);

  lineNumber = 0;
  int i = 0;
  while (std::getline(monFlux, line))
    {
      ++lineNumber;

      // Skip empty lines
      if (line.find_first_not_of(" \t\n\r") == std::string::npos) continue;

      // Skip comment lines
      if (line[0] == '#') continue;

      // Skip non-numeric headers
      std::istringstream issCheck(line);
      double absc;
      if (!(issCheck >> absc)) continue;

      // Read the value
      std::istringstream iss(line);
      if (!(iss >> absc))
        {
          Cerr << "ERROR: Unable to read numeric value in file '" << filename
               << "' at line " << lineNumber << "." << finl;
          Cerr << "Line content: [" << line << "]" << finl;
          Process::exit();
        }

      abscissa_[i] = absc;
      ++i;
    }

  Cerr << "Beam coordinates successfully read from '" << filename
       << "' (" << abscissa_.size() << " points)." << finl;
}

void Beam_model::readInputModalDeformation(Noms& modal_deformation_file_name)
{
  if (nbModes_ == 0)
    {
      Cerr << "Error: no deformation modes defined. Use 'nb_modes'.\n";
      Process::exit();
    }

  if (nbModes_ != modal_deformation_file_name.size())
    {
      Cerr << "Error: mismatch: nbModes = " << nbModes_
           << " but " << modal_deformation_file_name.size()
           << " files provided.\n";
      Process::exit();
    }

  Cerr << "Reading beam modal deformation from: "
       << modal_deformation_file_name << finl;

  const int nPoints         = abscissa_.size();
  const int modes_per_plane = nbModes_ / nb_planes_;

  std::vector<DoubleVect> tmp_u(nbModes_);
  std::vector<DoubleVect> tmp_R(nbModes_);


  for (int global_mode = 0; global_mode < nbModes_; ++global_mode)
    {
      const std::string filename(modal_deformation_file_name[global_mode]);
      std::ifstream monFlux(filename.c_str());
      if (!monFlux)
        {
          Cerr << "ERROR: Cannot open '" << filename << "'" << finl;
          Process::exit();
        }

      DoubleVect u_vec(nPoints);
      DoubleVect R_vec(nPoints);

      std::string lineContent;
      int lineNumber = 0;
      int k = 0;

      while (std::getline(monFlux, lineContent))
        {
          ++lineNumber;

          if (lineContent.find_first_not_of(" \t\r\n") == std::string::npos) continue;
          if (lineContent[0] == '#') continue;

          std::istringstream issCheck(lineContent);
          double dummy;
          if (!(issCheck >> dummy)) continue;

          std::istringstream iss(lineContent);
          double disp, rot;
          if (!(iss >> disp >> rot))
            {
              Cerr << "ERROR: '" << filename << "', line "
                   << lineNumber << " : could not read two numeric values.\n";
              Process::exit();
            }

          // Check for extra tokens
          double extra;
          if (iss >> extra)
            {
              Cerr << "ERROR: '" << filename << "', line "
                   << lineNumber << " : line contains more than two values.\n";
              Cerr << "Line: [" << lineContent << "]\n";
              Process::exit();
            }

          if (k >= nPoints)
            {
              Cerr << "ERROR: file '" << filename
                   << "' has too many data lines (expected " << nPoints << ").\n";
              Process::exit();
            }

          u_vec(k) = disp;
          R_vec(k) = rot;
          ++k;
        }

      if (k != nPoints)
        {
          Cerr << "ERROR: file '" << filename
               << "' contains " << k << " valid data rows, expected "
               << nPoints << ".\n";
          Process::exit();
        }

      tmp_u[global_mode] = u_vec;
      tmp_R[global_mode] = R_vec;

      monFlux.close();
    }

  u_.vide();
  R_.vide();


  for (int mode = 0; mode < modes_per_plane; ++mode)
    {
      DoubleTab Ut(nPoints, nb_planes_);
      DoubleTab Rt(nPoints, nb_planes_);

      for (int plane = 0; plane < nb_planes_; ++plane)
        {
          int global = plane * modes_per_plane + mode;
          for (int p = 0; p < nPoints; ++p)
            {
              Ut(p, plane) = tmp_u[global](p);
              Rt(p, plane) = tmp_R[global](p);
            }
        }
      u_.add(Ut);
      R_.add(Rt);
    }

  Cerr << "Modal deformation successfully loaded for "
       << nbModes_ << " modes over " << nb_planes_ << " planes." << finl;
}



void Beam_model::initialization(double displacement)
{
  int modes_per_plane = nbModes_ / nb_planes_;

  qSpeed_.resize(modes_per_plane, nb_planes_);
  qAcceleration_.resize(modes_per_plane, nb_planes_);
  qDisplacement_.resize(modes_per_plane, nb_planes_);
  fluidForceOnBeam_.resize(modes_per_plane, nb_planes_);

  qSpeed_ = 0.;
  qAcceleration_ = 0.;
  qDisplacement_ = displacement;
  fluidForceOnBeam_ = 0.;
}
void Beam_model::initialization()
{
  int modes_per_plane = nbModes_ / nb_planes_;

  qSpeed_.resize(modes_per_plane, nb_planes_);
  qAcceleration_.resize(modes_per_plane, nb_planes_);
  qDisplacement_.resize(modes_per_plane, nb_planes_);
  fluidForceOnBeam_.resize(modes_per_plane, nb_planes_);

  qSpeed_ = 0.;
  qAcceleration_ = 0.;
  qDisplacement_ = 0.;
  fluidForceOnBeam_ = 0.;
}

void Beam_model::interpolationWeights(const double& x, const double& y, const double& z, int& i, int& j, double& alpha, double& betha) const
{

  double h = abscissa_[1] -abscissa_[0]; //1d mesh pitch
  int abscissa_size = abscissa_.size();
  double s = 0.0;

  switch (longitudinal_axis_)
    {
    case 0:
      s = x;
      break;
    case 1:
      s = y;
      break;
    default:
      s = z;
      break;
    }

  i = int(s / h);
  j = (i + 1 < abscissa_size) ? i + 1 : i;

  //linear interpolation between points i and j
  if (i==j)
    {
      alpha=1.;
      betha=0.;
    }
  else if(abs(abscissa_[i] - s)< 1.e-4)
    {
      alpha=1.;
      betha=0.;
    }
  else if (abs(abscissa_[j] - s)< 1.e-4)
    {
      alpha=0.;
      betha=1.;
    }
  else
    {
      alpha = (abscissa_[j] - s)/h;
      betha = (s - abscissa_[i])/h;
      if(alpha <0.)
        {
          alpha=0.;
          betha=1.;
        }
      else if (betha < 0.)
        {
          alpha=1.;
          betha=0.;
        }

    }
}


double Beam_model::interpolationPhiOnThe3DSurface(const double& x, const double& y, const double& z, const int& comp, const DoubleTab& u) const
{
  double alpha=0.0, betha=0.0 ;
  int i=0, j=0;
  interpolationWeights(x,y,z, i, j, alpha, betha);
  double phi=alpha*u(i, comp) + betha*u(j, comp);

  return phi;
}

DoubleTab Beam_model::interpolationOnThe3DSurface(const double& x, const double& y, const double& z, const DoubleTab& u, const DoubleTab& R) const
{

  double alpha=0.0, betha=0.0 ;
  int i=0, j=0;
  interpolationWeights(x,y,z, i, j, alpha, betha);

  DoubleTab u_interp(nb_planes_), R_interp(nb_planes_);
  for (int plane = 0; plane < nb_planes_; ++plane)
    {
      u_interp[plane]=alpha*u(i, plane) + betha*u(j, plane);
      R_interp[plane]=alpha*R(i, plane) + betha*R(j, plane);
    }

  DoubleTab phi(3, nb_planes_);   // 3 components: axial + 2 transverse
  phi = 0.;

  DoubleVect pos(3);
  pos(0) = x - x0_;  // relative X position from beam axis
  pos(1) = y - y0_;  // relative Y position from beam axis
  pos(2) = z - z0_;  // relative Z position from beam axis

  for (int plane = 0; plane < nb_planes_; ++plane)
    {
      int e_long  = longitudinal_axis_;      // main axis of the beam (0=x,1=y,2=z)
      int e_trans = bending_dir_[plane];     // bending direction for this plane

      double W   = u_interp[plane];          // transverse displacement
      double Rot = R_interp[plane];          // derivative of displacement (rotation)

      // indices of the remaining axes for computing axial displacement
      int idx1 = (e_long + 1) % 3;
      int idx2 = (e_long + 2) % 3;

      if (e_trans == idx1)
        {
          phi(e_long, plane) =  Rot * pos(idx2);
          phi(idx2, plane)   = 0.0; // transverse component perpendicular to bending
          phi(idx1, plane)   = W;   // transverse along bending
        }
      else if (e_trans == idx2)
        {
          phi(e_long, plane) = -Rot * pos(idx1);
          phi(idx1, plane)   = 0.0; // transverse component perpendicular to bending
          phi(idx2, plane)   = W;   // transverse along bending
        }
      else
        {
          //bending direction equals longitudinal axis
          phi(e_long, plane) = 0.0;
          phi(idx1, plane)   = 0.0;
          phi(idx2, plane)   = 0.0;
        }

    }
  return phi;
}

void Beam_model::saveBeamForRestart() const
{
  if (!je_suis_maitre()) return;

  const Nom filename = beamName_ + "SaveBeamForRestart.txt";
  std::ofstream ofs_sauve(filename, std::ofstream::out | std::ofstream::trunc);
  ofs_sauve.precision(32);

  if (!ofs_sauve)
    {
      Cerr << "ERROR: Unable to open file '" << filename << "' for writing restart data." << finl;
      Process::exit();
    }

  int modes_per_plane = nbModes_ / nb_planes_;

  for (int plane = 0; plane < nb_planes_; ++plane)
    {
      for (int mode = 0; mode < modes_per_plane; ++mode)
        {
          ofs_sauve << temps_ << "  "
                    << qDisplacement_(mode, plane) << " "
                    << qSpeed_(mode, plane) << " "
                    << qAcceleration_(mode, plane) << std::endl;
        }
    }

  ofs_sauve.close();
  Cerr << "Beam state saved for restart in '" << filename << "' ("
       << nbModes_ << " modes, " << nb_planes_ << " planes)." << finl;
}


void Beam_model::printOutputBeam1D(bool first_writing) const
{
  if (!je_suis_maitre()) return;

  int nb_output_points = output_position_1D_.size();
  int modes_per_plane = nbModes_ / nb_planes_;

  for (int plane = 0; plane < nb_planes_; ++plane)
    {
      //Compute 1D fields
      DoubleVect displacement(nb_output_points);
      DoubleVect velocity(nb_output_points);
      DoubleVect acceleration(nb_output_points);

      displacement = 0.;
      velocity     = 0.;
      acceleration = 0.;

      for (int mode = 0; mode < modes_per_plane; ++mode)
        {
          const DoubleTab& u = u_(mode);

          for (int k = 0; k < nb_output_points; ++k)
            {
              int idx = int(output_position_1D_[k]);
              displacement[k] += qDisplacement_(mode, plane) * u(idx, plane);
              velocity[k]     += qSpeed_(mode, plane)        * u(idx, plane);
              acceleration[k] += qAcceleration_(mode, plane)  *u(idx, plane);

            }
        }

      // ======== DISPLACEMENT FILE ========
      {
        Nom filename = beamName_;
        filename += "_Displacement1D_plane_";
        filename += Nom(plane);
        filename += ".out";

        SFichier out;
        out.ouvrir(filename, first_writing ? ios::out : ios::app);
        out.setf(ios::scientific);

        if (first_writing)
          {
            out << "# Plane " << plane
                << " Beam 1D displacement: time x y z at points ";
            for (int k = 0; k < nb_output_points; ++k)
              out << output_position_1D_[k] << " ";
            out << finl;
          }

        out << temps_ << " ";
        for (int k = 0; k < nb_output_points; ++k)
          out << displacement[k] << " ";
        out << finl;
      }

      // ======== VELOCITY FILE ========
      {
        Nom filename = beamName_;
        filename += "_Velocity1D_plane_";
        filename += Nom(plane);
        filename += ".out";

        SFichier out;
        out.ouvrir(filename, first_writing ? ios::out : ios::app);
        out.setf(ios::scientific);

        if (first_writing)
          {
            out << "# Plane " << plane
                << " Beam 1D velocity: time x y z at points ";
            for (int k = 0; k < nb_output_points; ++k)
              out << output_position_1D_[k] << " ";
            out << finl;
          }

        out << temps_ << " ";
        for (int k = 0; k < nb_output_points; ++k)
          out << velocity[k] << " ";
        out << finl;
      }

      // ======== ACCELERATION FILE ========
      {
        Nom filename = beamName_;
        filename += "_Acceleration1D_plane_";
        filename += Nom(plane);
        filename += ".out";

        SFichier out;
        out.ouvrir(filename, first_writing ? ios::out : ios::app);
        out.setf(ios::scientific);

        if (first_writing)
          {
            out << "# Plane " << plane
                << " Beam 1D acceleration: time x y z at points ";
            for (int k = 0; k < nb_output_points; ++k)
              out << output_position_1D_[k] << " ";
            out << finl;
          }

        out << temps_ << " ";
        for (int k = 0; k < nb_output_points; ++k)
          out << acceleration[k] << " ";
        out << finl;
      }
    }
}



void Beam_model::printOutputBeam3D(bool first_writing) const
{
  if (!je_suis_maitre()) return;
  const int nb_output_points = output_position_3D_.dimension(0);
  DoubleTab displacement(nb_output_points, 3);
  DoubleTab velocity(nb_output_points, 3);
  DoubleTab acceleration(nb_output_points, 3);

  displacement = 0.;
  velocity = 0.;
  acceleration = 0.;

  const int modes_per_plane = qDisplacement_.dimension(0);
  const int nb_planes = qDisplacement_.dimension(1);

  DoubleTab phi3D(3, nb_planes);

  for (int p = 0; p < nb_planes; p++)
    {
      for (int m = 0; m < modes_per_plane; m++)
        {
          const DoubleTab& u = u_(m);
          const DoubleTab& R = R_(m);

          for (int k = 0; k < nb_output_points; k++)
            {

              phi3D = interpolationOnThe3DSurface(output_position_3D_(k,0),
                                                  output_position_3D_(k,1),
                                                  output_position_3D_(k,2),
                                                  u, R);

              for (int i = 0; i < 3; i++)
                {
                  displacement(k,i) += qDisplacement_(m,p) * phi3D(i,p);
                  velocity(k,i)     += qSpeed_(m,p)        * phi3D(i,p);
                  acceleration(k,i) += qAcceleration_(m,p) * phi3D(i,p);
                }
            }
        }
    }



  Nom filename_disp(beamName_ + "_Displacement3D.out");
  Nom filename_speed(beamName_ + "_Velocity3D.out");
  Nom filename_acc(beamName_ + "_Acceleration3D.out");

  if (!displacement_out_3d_.is_open())
    {
      displacement_out_3d_.ouvrir(filename_disp, (first_writing ? ios::out : ios::app));
      displacement_out_3d_.setf(ios::scientific);
    }
  if (!speed_out_3d_.is_open())
    {
      speed_out_3d_.ouvrir(filename_speed, (first_writing ? ios::out : ios::app));
      speed_out_3d_.setf(ios::scientific);
    }
  if (!acceleration_out_3d_.is_open())
    {
      acceleration_out_3d_.ouvrir(filename_acc, (first_writing ? ios::out : ios::app));
      acceleration_out_3d_.setf(ios::scientific);
    }

  // --- HEADER ---
  if (first_writing)
    {
      displacement_out_3d_ << "# 3D beam displacement: time, (x,y,z) at each output point ";
      speed_out_3d_        << "# 3D beam velocity: time, (x,y,z) at each output point ";
      acceleration_out_3d_ << "# 3D beam acceleration: time, (x,y,z) at each output point ";

      for (int k = 0; k < nb_output_points; k++)
        displacement_out_3d_ << "( " << output_position_3D_(k,0) << " " << output_position_3D_(k,1) << " " << output_position_3D_(k,2) << " ) ";

      displacement_out_3d_ << finl;

      for (int k = 0; k < nb_output_points; k++)
        speed_out_3d_ << "( " << output_position_3D_(k,0) << " " << output_position_3D_(k,1) << " " << output_position_3D_(k,2) << " ) ";

      speed_out_3d_ << finl;

      for (int k = 0; k < nb_output_points; k++)
        acceleration_out_3d_ << "( " << output_position_3D_(k,0) << " " << output_position_3D_(k,1) << " " << output_position_3D_(k,2) << " ) ";

      acceleration_out_3d_ << finl;
    }

  // --- WRITE VALUES ---
  displacement_out_3d_ << temps_ << " ";
  speed_out_3d_        << temps_ << " ";
  acceleration_out_3d_ << temps_ << " ";

  for (int k = 0; k < nb_output_points; k++)
    {
      for (int i = 0; i < 3; i++)
        {
          displacement_out_3d_  << displacement(k,i) << " ";
          speed_out_3d_         << velocity(k,i)     << " ";
          acceleration_out_3d_  << acceleration(k,i) << " ";
        }
    }

  displacement_out_3d_  << finl;
  speed_out_3d_         << finl;
  acceleration_out_3d_  << finl;

}

void Beam_model::printOutputFluidForceOnBeam(bool first_writing) const
{
  if (!je_suis_maitre())
    return;

  int modes_per_plane = nbModes_ / nb_planes_;

  for (int plane = 0; plane < nb_planes_; ++plane)
    {
      // Build the filename for this plane
      Nom filename = beamName_;
      filename += "_ModalForceFluide1D_plane_";
      filename += Nom(plane);
      filename += ".out";

      // Open the file: overwrite if first writing, append otherwise
      SFichier out;
      out.ouvrir(filename, (first_writing ? ios::out : ios::app));
      out.setf(ios::scientific);

      // ===== Header =====
      if (first_writing)
        {
          out << "# Modal fluid force for plane " << plane << ": time ";
          for (int m = 0; m < modes_per_plane; ++m)
            out << "(mode " << m + 1 << ") ";
          out << finl;
        }
      else
        {
          // ===== Values =====
          out << temps_ << " ";
          for (int m = 0; m < modes_per_plane; ++m)
            out << fluidForceOnBeam_(m, plane) << " ";
          out << finl;
        }
      out.close();
    }
}



void Beam_model::setCenterCoordinates(const double& x0,const double& y0, const double& z0)
{
  x0_=x0;
  y0_=y0;
  z0_=z0;
}

DoubleTab& Beam_model::getVelocity(const double& tps, const double& dt)
{
  if(dt == 0.)
    {
      return qSpeed_;
    }
  else if(temps_!=tps) //   update qSpeed_ only once per time step!
    {
      temps_=tps;
      return NewmarkScheme(dt);

    }
  else
    return qSpeed_;
}

//Solve the beam dynamics. Time integration scheme
DoubleTab& Beam_model::NewmarkScheme(const double& dt)
{
  double squareDt = dt * dt;
  int modes_per_plane = nbModes_ / nb_planes_;

  // Loop over planes
  for (int plane = 0; plane < nb_planes_; ++plane)
    {
      for (int mode = 0; mode < modes_per_plane; ++mode)
        {
          // Predict displacement and speed
          double qDispl_pred = qDisplacement_(mode, plane) + dt * qSpeed_(mode, plane)
                               + squareDt * (0.5 - beta_) * qAcceleration_(mode, plane);
          double qSpeed_pred = qSpeed_(mode, plane) + dt * (1. - gamma_) * qAcceleration_(mode, plane);

          // Newmark coefficients
          double coeff1 = mass_(mode, plane)
                          + dt * gamma_ * (1. + alpha_) * damping_(mode, plane)
                          + squareDt * beta_ * (1. + alpha_) * stiffness_(mode, plane);
          double coeff2 = (1. + alpha_) * damping_(mode, plane) * qSpeed_pred;
          double coeff3 = (1. + alpha_) * stiffness_(mode, plane) * qDispl_pred;
          double coeff4 = alpha_ * stiffness_(mode, plane) * qDisplacement_(mode, plane);
          double coeff5 = alpha_ * damping_(mode, plane) * qSpeed_(mode, plane);

          // Update acceleration, displacement, and speed
          qAcceleration_(mode, plane) = (fluidForceOnBeam_(mode, plane) - coeff2 - coeff3 + coeff4 + coeff5) / coeff1;
          qDisplacement_(mode, plane) = qDispl_pred + squareDt * beta_ * qAcceleration_(mode, plane);
          qSpeed_(mode, plane) = qSpeed_pred + dt * gamma_ * qAcceleration_(mode, plane);
        }
    }

  saveBeamForRestart();
  if (output_position_1D_.size() > 0) printOutputBeam1D();
  if (output_position_3D_.size() > 0) printOutputBeam3D();

  return qSpeed_;
}



void Beam_model::read_beam(Entree& is)
{
  Motcle open_brace("{");
  Motcle close_brace("}");
  Motcle motlu;
  Nom nomlu;
  Nom masse_and_stiffness_file_name;
  Noms phi_file_name;
  Nom absc_file_name;
  Nom CI_file_name="none";
  Nom Restart_file_name="none";
  int var_int=0;
  int nb_modes=0;
  int nb_output_points_1D=0;
  DoubleVect output_position_1D(nb_output_points_1D);
  int nb_output_points_3D=0;
  DoubleTab output_position_3D(nb_output_points_3D,0);
  double var_double;
  is >> motlu;
  if (motlu != open_brace)
    {
      Cerr << "Error while reading Beam\n";
      Cerr << "Expected a { but found: "<< motlu;
      Process::exit();
    }
  while(1)
    {
      // reading a boundary name or }
      is >> nomlu;
      motlu=nomlu;
      if(motlu=="nb_modes")
        {
          is >> nb_modes;
          nbModes_=nb_modes;
          Cerr << "Number of modes : " <<  nbModes_ << finl;
        }
      if(motlu=="nb_planes")
        {
          is >> var_int;

          // Validate that number of planes is either 1 or 2
          if (var_int != 1 && var_int != 2)
            {
              Cerr << "ERROR: invalid number of beam planes: " << var_int << finl;
              Cerr << "Valid values are: 1 or 2." << finl;
              Process::exit();
            }
          nb_planes_=var_int;
          Cerr << "Number of planes : " <<  nb_planes_ << finl;
        }
      if(motlu=="longitudinal_axis")
        {
          is >> nomlu;
          if (nomlu == "x" || nomlu == "X")
            var_int = 0;
          else if (nomlu == "y" || nomlu == "Y")
            var_int = 1;
          else if (nomlu == "z" || nomlu == "Z")
            var_int = 2;
          else
            {
              Cerr << "ERROR: invalid main axis: " << nomlu
                   << "'. Valid options are: x, y, or z." << finl;
              Process::exit();
            }
          longitudinal_axis_=var_int;
          Cerr << "Longitudinal axis : " <<  longitudinal_axis_ << finl;
        }
      if(motlu=="bendingDirection")
        {
          is >> nomlu;
          if (nomlu == "x" || nomlu == "X")
            var_int = 0;
          else if (nomlu == "y" || nomlu == "Y")
            var_int = 1;
          else if (nomlu == "z" || nomlu == "Z")
            var_int = 2;
          else
            {
              Cerr << "ERROR: invalid direction of bending: '" << nomlu
                   << "'. Valid options are: x, y, or z." << finl;
              Process::exit();
            }
          bending_dir_[0]=var_int;
          Cerr << "Bending direction : " <<  bending_dir_ << finl;
        }
      if(motlu=="BaseCenterCoordinates")
        {
          is >> var_double;
          double x=var_double;
          is >> var_double;
          double y=var_double;
          is >> var_double;
          double z=var_double;
          setCenterCoordinates(x,y,z);

        }
      if(motlu=="Young_Module")
        {
          is >> var_double;
          young_=var_double;
          Cerr << "Young module : " <<  young_ << finl;
        }
      if(motlu=="Rho_beam")
        {
          is >> var_double;
          rho_=var_double;
          Cerr << "Rho beam : " <<  rho_ << finl;
        }
      if(motlu=="Mass_and_stiffness_file_name")
        {
          is >> nomlu;
          masse_and_stiffness_file_name=nomlu;

        }
      if(motlu=="Absc_file_name")
        {
          is >> nomlu;
          absc_file_name=nomlu;

        }
      if(motlu=="CI_file_name")
        {
          is >> nomlu;
          CI_file_name=nomlu;

        }

      if(motlu=="Modal_deformation_file_name")
        {
          is >> var_int;
          for (int i=0; i< var_int; i ++)
            {
              is >>nomlu;
              phi_file_name.add(nomlu);
            }
        }
      if(motlu=="NewmarkTimeScheme")
        {
          is >> nomlu;
          double alpha=-0.1; //default value
          if(nomlu =="HHT")
            {
              is>>alpha;
            }
          setTimeScheme(nomlu, alpha);

        }
      if(motlu=="Output_position_1D")
        {
          Cerr << "Beam Output_position_1D "<< finl;
          is>>nb_output_points_1D;
          output_position_1D.resize(nb_output_points_1D);
          double poz;
          for (int i=0; i< nb_output_points_1D; i ++)
            {
              is >>poz;
              output_position_1D[i]=poz;
            }

          setOutputPosition1D(output_position_1D);
        }
      if(motlu=="Output_position_3D")
        {
          Cerr << "Beam Output_position_3D "<< finl;
          is>>nb_output_points_3D;
          output_position_3D.resize(nb_output_points_3D, 3);
          double poz=0.;
          for (int i=0; i< nb_output_points_3D; i ++)
            {
              for (int j=0; j< 3; j ++)
                {
                  is >>poz;
                  output_position_3D(i,j)=poz;
                }
            }
          setOutputPosition3D(output_position_3D);
        }
      if (motlu=="Restart_file_name")
        {
          is >> nomlu;
          Restart_file_name=nomlu;
        }

      if (motlu == close_brace)
        break;
    }

  // Warning: Do NOT change the order of these function calls. The correct execution of the code depends on this sequence.
  if (longitudinal_axis_==-1)
    {
      Cerr << "ERROR: you must specify the axis along the length of the beam."<<finl;
      Process::exit();
    }
  if (bending_dir_[0] == -1)
    {
      Cerr << "ERROR: you must specify the bending direction."<<finl;
      Process::exit();
    }

  // Check that the first bending direction is not the longitudinal axis
  if (bending_dir_[0] == longitudinal_axis_)
    {
      Cerr << "ERROR: invalid bendingDirection: it must be different from the longitudinal_axis."<<finl;
      Process::exit();
    }

  // If there are two bending planes, determine the second direction automatically
  if (nb_planes_ == 2)
    {
      bending_dir_.resize(2);

      // The second bending direction is the remaining axis among {0,1,2}
      for (int axis = 0; axis < 3; ++axis)
        {
          if (axis != longitudinal_axis_ && axis != bending_dir_[0])
            {
              bending_dir_[1] = axis;
              break;
            }
        }
    }

  // Warning: Do NOT change the order of these function calls. The correct execution of the code depends on this sequence.
  readInputAbscFiles(absc_file_name);
  readInputMassStiffnessFiles(masse_and_stiffness_file_name);
  readInputModalDeformation(phi_file_name);
  if(CI_file_name!="none")
    {
      readInputCIFile(CI_file_name);
    }
  else
    {
      initialization();
    }
  if(Restart_file_name!="none")
    {
      readRestartFile(Restart_file_name);
    }
  else
    {
      if(je_suis_maitre())
        {
          bool first_writing=true;
          printOutputFluidForceOnBeam(first_writing);
          if (nb_output_points_1D>0)
            printOutputBeam1D(first_writing);
          if (nb_output_points_3D>0)
            printOutputBeam3D(first_writing);
        }
    }
}

