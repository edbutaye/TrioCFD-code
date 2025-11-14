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
// XD  attr nb_modes entier n 0 Number of modes
// XD  attr longitudinal_axis chaine dir 0 x, y, z. Axis along the length of the beam
// XD  attr bendingDirection chaine dir_bending 0 x, y, z . Direction of  bending
// XD  attr nb_planes entier nplanes 0 Number of planes used in the beam dynamic model
// XD  attr NewmarkTimeScheme NewmarkTimeScheme_deriv NewmarkTimeScheme 0 Solve the beam dynamics. Time integration scheme: choice between MA (Newmark mean acceleration),  FD (Newmark finite differences), and HHT alpha (Hilber-Hughes-Taylor, alpha usually -0.1 )
// XD  attr Mass_and_stiffness_file_name chaine  Mass_and_stiffness_file_name 0 Name of the file containing the diagonal modal mass, stiffness, and damping matrices.
// XD  attr Absc_file_name chaine Absc_file_name 0 Name of the file containing the coordinates of the Beam
// XD  attr  Modal_deformation_file_name listchaine  Modal_deformation_file_name 0 Name of the file containing the modal deformation of the Beam (mandatory if different from 0. 0. 0.)
// XD  attr Young_Module floattant young 1 Young Module
// XD  attr Rho_beam floattant rho 1 Beam density
// XD  attr BaseCenterCoordinates listf pos_center 1 position of the base center coordinates on the Beam
// XD  attr CI_file_name chaine CI_file_name 1 Name of the file containing the initial condition of the Beam
// XD  attr Restart_file_name chaine Restart_file_name 1 SaveBeamForRestart.txt file to restart the calculation
// XD  attr Output_position_1D list pt1d 1 nb_points  position Post-traitement of specific points on the Beam
// XD  attr Output_position_3D listpoints pt3d 1 nb_points  position Post-traitement of specific points on the 3d FSI boundary

Beam_model::Beam_model()
{

  nbModes_=0;;
  longitudinal_axis_=0;
  bending_dir_=1;
  nb_planes_=1;
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

    mass_.resize(nbModes_);
    stiffness_.resize(nbModes_);
    damping_.resize(nbModes_);

    const std::string filename(masse_and_stiffness_file_name);
    std::ifstream monFlux(filename.c_str());

    if (!monFlux)
    {
        Cerr << "ERROR: Unable to open file '" << filename << "'." << finl;
        Process::exit();
    }

    std::string line;
    int i = 0;
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

        if (i >= nbModes_)
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
        mass_[i] = mass;
        stiffness_[i] = stiffness;
        damping_[i] = damping;
        ++i;
    }

    if (i < nbModes_)
    {
        Cerr << "ERROR: File '" << filename << "' contains too few data lines. "
             << "Expected nb_modes = " << nbModes_ << ", but only " << i << " valid lines were read." << finl;
        Process::exit();
    }

    Cerr << "Beam coefficients successfully read from '" << filename
         << "' (" << i << " modes)." << finl;
}


void Beam_model::readInputCIFile(Nom& CI_file_name)
{
    Cerr << "Reading initial condition beam from file: " << CI_file_name << finl;

    qSpeed_.resize(nbModes_);
    qAcceleration_.resize(nbModes_);
    qDisplacement_.resize(nbModes_);
    fluidForceOnBeam_.resize(nbModes_);

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
    int i = 0;
    int lineNumber = 0;

    while (std::getline(monFlux, line))
    {
        ++lineNumber;

        // Skip empty lines
        bool empty = true;
        for (char c : line) { if (!std::isspace(static_cast<unsigned char>(c))) { empty = false; break; } }
        if (empty) continue;

        // Skip comment lines starting with '#'
        if (line[0] == '#') continue;

        // Skip non-numeric header lines
        std::istringstream issCheck(line);
        double dummy;
        if (!(issCheck >> dummy)) continue;

        if (i >= nbModes_)
        {
            Cerr << "ERROR: File '" << filename << "' contains more data lines than nb_modes = "
                 << nbModes_ << " (at line " << lineNumber << ")." << finl;
            Process::exit();
        }

        std::istringstream iss(line);
        double displacement;
        if (!(iss >> displacement))
        {
            Cerr << "ERROR: Unable to read displacement value in file '" << filename
                 << "' at line " << lineNumber << "." << finl;
            Cerr << "Line content: [" << line << "]" << finl;
            Process::exit();
        }

        // Check if a second value exists (should be only one per line)
        double extra;
        if (iss >> extra)
        {
            Cerr << "ERROR: In file '" << filename << "', line " << lineNumber
                 << " contains more than one value." << finl;
            Cerr << "Line content: [" << line << "]" << finl;
            Process::exit();
        }

        qDisplacement_[i] = displacement;
        // Initialize acceleration from fluid force and stiffness
        qAcceleration_[i] = fluidForceOnBeam_[i] - (stiffness_[i] / mass_[i]) * qDisplacement_[i];

        ++i;
    }

    // Check that the file contains exactly nbModes_ values
    if (i < nbModes_)
    {
        Cerr << "ERROR: File '" << filename << "' contains too few valid data lines. "
             << "Expected nb_modes = " << nbModes_ << ", but only " << i << " were read." << finl;
        Process::exit();
    }

    Cerr << "Initial condition successfully read from '" << filename
         << "' (" << i << " modes)." << finl;
}

void Beam_model::readRestartFile(Nom& Restart_file_name)
{
    Cerr << "Reading restart file: " << Restart_file_name << finl;
    const std::string filename(Restart_file_name);
    std::ifstream monFlux(filename.c_str());

    if (!monFlux)
    {
        Cerr << "Beam_model::readRestartFile" << finl;
        Cerr << "ERROR: Unable to open the restart file '" << filename << "'." << finl;
        Process::exit();
    }

    std::string line;
    int i = 0;
    int lineNumber = 0;

    while (std::getline(monFlux, line))
    {
        ++lineNumber;

        // Skip empty lines
        bool empty = true;
        for (char c : line) { if (!std::isspace(static_cast<unsigned char>(c))) { empty = false; break; } }
        if (empty) continue;

        // Skip comment lines starting with '#'
        if (line[0] == '#') continue;

        // Skip non-numeric header lines
        std::istringstream issCheck(line);
        double dummy;
        if (!(issCheck >> dummy)) continue;

        if (i >= nbModes_)
        {
            Cerr << "ERROR: File '" << filename << "' contains more data lines than nb_modes = "
                 << nbModes_ << " (at line " << lineNumber << ")." << finl;
            Process::exit();
        }

        std::istringstream iss(line);
        double temps, displacement, speed, acceleration;

        if (!(iss >> temps >> displacement >> speed >> acceleration))
        {
            Cerr << "ERROR: In file '" << filename << "', line " << lineNumber
                 << " does not contain exactly 4 numerical values." << finl;
            Cerr << "Line content: [" << line << "]" << finl;
            Process::exit();
        }

        // Check for extra values
        double extra;
        if (iss >> extra)
        {
            Cerr << "ERROR: In file '" << filename << "', line " << lineNumber
                 << " contains more than 4 values." << finl;
            Cerr << "Line content: [" << line << "]" << finl;
            Process::exit();
        }

        // Store values
        temps_ = temps;
        qDisplacement_[i] = displacement;
        qSpeed_[i] = speed;
        qAcceleration_[i] = acceleration;

        ++i;
    }

    // Check count
    if (i < nbModes_)
    {
        Cerr << "ERROR: File '" << filename << "' contains too few valid data lines. "
             << "Expected nb_modes = " << nbModes_ << ", but only " << i << " were read." << finl;
        Process::exit();
    }

    tempsComputeForceOnBeam_ = temps_;
    Cerr << "Restart file successfully read from '" << filename
         << "' (" << i << " modes)." << finl;
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
        Cerr << "Error: no deformation modes defined. At least one deformation mode is required.\n";
        Cerr << "Use the 'nb_modes' keyword to specify the number of modes.\n";
        Process::exit();
    }

    if (nbModes_ != modal_deformation_file_name.size())
    {
        Cerr << "Error: mismatch between the number of modal deformation files and nb_modes.\n"
             << "Expected " << nbModes_ << " files, but found " << modal_deformation_file_name.size() << ".\n"
             << "Adjust the inputs and restart the program.\n";
        Process::exit();
    }

    Cerr << "Reading beam modal deformation coefficients from files: " << modal_deformation_file_name << finl;

    for (int count = 0; count < nbModes_; ++count)
    {
        const std::string filename(modal_deformation_file_name[count]);
        std::ifstream monFlux(filename.c_str());
        if (!monFlux)
        {
            Cerr << "ERROR: Unable to open the file '" << filename << "'." << finl;
            Process::exit();
        }

        int nPoints = abscissa_.size();
        DoubleTab u(nPoints, 3);
        DoubleTab R(nPoints, 3);

        std::string lineContent;
        int lineNumber = 0;
        int line = 0;

        while (std::getline(monFlux, lineContent))
        {
            ++lineNumber;

            // Skip empty lines
            if (lineContent.find_first_not_of(" \t\n\r") == std::string::npos) continue;

            // Skip comment lines
            if (lineContent[0] == '#') continue;

            // Skip non-numeric headers
            std::istringstream issCheck(lineContent);
            double dummy;
            if (!(issCheck >> dummy)) continue;

            // Parse exactly 6 values
            std::istringstream iss(lineContent);
            double ux, uy, uz, rx, ry, rz;
            if (!(iss >> ux >> uy >> uz >> rx >> ry >> rz))
            {
                Cerr << "ERROR: Unable to read 6 numeric values in file '" << filename
                     << "' at line " << lineNumber << "." << finl;
                Cerr << "Line content: [" << lineContent << "]" << finl;
                Process::exit();
            }

            // Check for extra values
            double extra;
            if (iss >> extra)
            {
                Cerr << "ERROR: Line " << lineNumber << " in file '" << filename
                     << "' contains more than 6 values." << finl;
                Cerr << "Line content: [" << lineContent << "]" << finl;
                Process::exit();
            }

            // Store data
            u(line, 0) = ux; u(line, 1) = uy; u(line, 2) = uz;
            R(line, 0) = rx; R(line, 1) = ry; R(line, 2) = rz;

            ++line;
        }

        monFlux.close();

        //  Check that we read exactly nPoints lines
        if (line != nPoints)
        {
            Cerr << "ERROR: File '" << filename << "' contains " << line
                 << " valid data lines, but " << nPoints << " lines were expected." << finl;
            Process::exit();
        }

        // Add to modal arrays
        u_.add(u);
        R_.add(R);
    }

    Cerr << "All modal deformation files successfully read." << finl;
}


void Beam_model::initialization(double displacement)
{
  qSpeed_.resize(nbModes_);
  qAcceleration_.resize(nbModes_);
  qDisplacement_.resize(nbModes_);
  fluidForceOnBeam_.resize(nbModes_);
  qSpeed_=0.;
  qAcceleration_=0.;
  qDisplacement_=displacement;
  fluidForceOnBeam_=0.;

}
void Beam_model::initialization()
{
  qSpeed_.resize(nbModes_);
  qAcceleration_.resize(nbModes_);
  qDisplacement_.resize(nbModes_);
  fluidForceOnBeam_.resize(nbModes_);
  qSpeed_=0.;
  qAcceleration_=0.;
  qDisplacement_=0.;
  fluidForceOnBeam_=0.;
}

DoubleVect Beam_model::interpolationOnThe3DSurface(const double& x, const double& y, const double& z, const DoubleTab& u, const DoubleTab& R) const
{
  DoubleVect phi(3);
  phi=0.;
  double h = abscissa_[1] -abscissa_[0]; //1d mesh pitch
  int abscissa_size = abscissa_.size();
  double s=0.;
  double xs=x;
  double ys=y;
  double zs=z;
  if (longitudinal_axis_== 0)
    {
      s = xs;
      xs=0.;
    }
  else if (longitudinal_axis_== 1)
    {
      s = ys;
      ys=0.;
    }

  else
    {
      s = zs;
      zs=0.;
    }

  int i, j ;
  i = int(s/h);
  if((i+1) < abscissa_size)
    {
      j= i+1;
    }
  else
    {
      j=i;
    }
  double ux, uy, uz, Rx, Ry, Rz;

  //linear interpolation between points i and j
  double alpha, betha ;
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

  ux=alpha*u(i, 0) + betha*u(j, 0);
  uy=alpha*u(i, 1) + betha*u(j, 1);
  uz=alpha*u(i, 2) + betha*u(j, 2);
  Rx=alpha*R(i, 0) + betha*R(j, 0);
  Ry=alpha*R(i, 1) + betha*R(j, 1);
  Rz=alpha*R(i, 2) + betha*R(j, 2);

  phi[0] =ux + Ry*(zs-z0_) -Rz*(ys-y0_);
  phi[1] =uy + Rz*(xs-x0_) -Rx*(zs-z0_);
  phi[2] =uz + Rx*(ys-y0_) -Ry*(xs-x0_);



  return phi;
}

void Beam_model::saveBeamForRestart() const
{

  if (je_suis_maitre())
    {
      std::ofstream ofs_sauve;
      ofs_sauve.open (beamName_+"SaveBeamForRestart.txt", std::ofstream::out | std::ofstream::trunc);
      ofs_sauve.precision(32);
      for(int j=0; j < nbModes_; j++)
        {
          ofs_sauve<<temps_<<"  "<<qDisplacement_[j]<<" "<<qSpeed_[j]<<" "<<qAcceleration_[j]<<" "<<endl;
        }
      ofs_sauve.close();
    }

}

void Beam_model::printOutputBeam1D(bool first_writing) const
{

  if (je_suis_maitre())
    {
      int nb_output_points= output_position_1D_.size();
      DoubleTab displacement(nb_output_points,3);
      DoubleTab velocity(nb_output_points,3);
      DoubleTab acceleration(nb_output_points,3);
      displacement=0.;
      velocity=0.;
      acceleration=0.;
      for(int j=0; j < nbModes_; j++)
        {
          const DoubleTab& u=u_(j);
          for(int k=0; k<nb_output_points; k++)
            {
              for(int i=0; i<3; i++)
                {
                  displacement(k, i) += qDisplacement_[j]*u(int(output_position_1D_[k]),i);
                  velocity(k, i) += qSpeed_[j]*u(int(output_position_1D_[k]),i);
                  acceleration(k, i) += qAcceleration_[j]*u(int(output_position_1D_[k]),i);
                }
            }
        }
      Nom filename_disp(beamName_);
      filename_disp+="_Displacement1D.out";
      Nom filename_speed(beamName_);
      filename_speed+="_Velocity1D.out";
      Nom filename_acc(beamName_);
      filename_acc+="_Acceleration1D.out";
      if (!displacement_out_1d_.is_open())
        {
          displacement_out_1d_.ouvrir(filename_disp, (first_writing?ios::out:ios::app));
          displacement_out_1d_.setf(ios::scientific);
        }
      if (!speed_out_1d_.is_open())
        {
          speed_out_1d_.ouvrir(filename_speed, (first_writing?ios::out:ios::app));
          speed_out_1d_.setf(ios::scientific);
        }
      if (!acceleration_out_1d_.is_open())
        {
          acceleration_out_1d_.ouvrir(filename_acc, (first_writing?ios::out:ios::app));
          acceleration_out_1d_.setf(ios::scientific);
        }
      // comments are added to the file header
      if (first_writing)
        {
          displacement_out_1d_<< "# Printing Beam 1D displacement: time  values of x y z -component at points ";
          speed_out_1d_<< "# Printing Beam 1D velocity: time  values of x y z -component at points ";
          acceleration_out_1d_<< "# Printing Beam 1D acceleration: time  values of x y z -component at points ";
          for(int k=0; k<nb_output_points; k++)
            {
              displacement_out_1d_<<output_position_1D_[k]<<" ";
              speed_out_1d_<<output_position_1D_[k]<<" ";
              acceleration_out_1d_<<output_position_1D_[k]<<" ";
            }
          displacement_out_1d_<<finl;
          speed_out_1d_<<finl;
          acceleration_out_1d_<<finl;
        }
      displacement_out_1d_<< temps_<< " ";
      speed_out_1d_<< temps_<< " ";
      acceleration_out_1d_<< temps_<< " ";
      for(int k=0; k<nb_output_points; k++)
        {
          for(int i=0; i<3; i++)
            {
              displacement_out_1d_<<displacement(k, i)<<" ";
              speed_out_1d_<<velocity(k, i)<<" ";
              acceleration_out_1d_<<acceleration(k, i)<<" ";
            }
        }
      displacement_out_1d_<<finl;
      speed_out_1d_<<finl;
      acceleration_out_1d_<<finl;
    }
}
void Beam_model::printOutputBeam3D(bool first_writing) const
{

  if (je_suis_maitre())
    {
      int nb_output_points= output_position_3D_.dimension(0);
      DoubleTab displacement(nb_output_points,3);
      DoubleTab velocity(nb_output_points,3);
      DoubleTab acceleration(nb_output_points,3);
      displacement=0.;
      velocity=0.;
      acceleration=0.;
      DoubleVect phi3D(3);
      for(int j=0; j < nbModes_; j++)
        {
          const DoubleTab& u=u_(j);
          const DoubleTab& R=R_(j);
          for(int k=0; k<nb_output_points; k++)
            {
              phi3D=interpolationOnThe3DSurface(output_position_3D_(k,0),output_position_3D_(k,1),output_position_3D_(k,2), u, R);
              for(int i=0; i<3; i++)
                {
                  displacement(k, i) += qDisplacement_[j]*phi3D[i];
                  velocity(k, i) += qSpeed_[j]*phi3D[i];
                  acceleration(k, i) += qAcceleration_[j]*phi3D[i];
                }
            }
        }

      Nom filename_disp(beamName_);
      filename_disp+="_Displacement3D.out";
      Nom filename_speed(beamName_);
      filename_speed+="_Velocity3D.out";
      Nom filename_acc(beamName_);
      filename_acc+="_Acceleration3D.out";
      if (!displacement_out_3d_.is_open())
        {
          displacement_out_3d_.ouvrir(filename_disp, (first_writing?ios::out:ios::app));
          displacement_out_3d_.setf(ios::scientific);
        }
      if (!speed_out_3d_.is_open())
        {
          speed_out_3d_.ouvrir(filename_speed, (first_writing?ios::out:ios::app));
          speed_out_3d_.setf(ios::scientific);
        }
      if (!acceleration_out_3d_.is_open())
        {
          acceleration_out_3d_.ouvrir(filename_acc, (first_writing?ios::out:ios::app));
          acceleration_out_3d_.setf(ios::scientific);
        }
      // comments are added to the file header
      if (first_writing)
        {
          displacement_out_3d_<< "# Printing Beam 3D displacement: time  values of x y z -component at points ";
          speed_out_3d_<< "# Printing Beam 3D velocity: time  values of x y z -component at points ";
          acceleration_out_3d_<< "# Printing Beam 3D acceleration: time  values of x y z -component at points ";
          for(int k=0; k<nb_output_points; k++)
            {
              displacement_out_3d_<<"( ";
              speed_out_3d_<<"( ";
              acceleration_out_3d_<<"( ";
              for(int i=0; i<3; i++)

                {
                  displacement_out_3d_<<output_position_3D_(k, i)<<" ";
                  speed_out_3d_<<output_position_3D_(k, i)<<" ";
                  acceleration_out_3d_<<output_position_3D_(k, i)<<" ";
                }
              displacement_out_3d_<<")";
              speed_out_3d_<<")";
              acceleration_out_3d_<<")";
            }
          displacement_out_3d_<<finl;
          speed_out_3d_<<finl;
          acceleration_out_3d_<<finl;
        }
      displacement_out_3d_<< temps_<< " ";
      speed_out_3d_<< temps_<< " ";
      acceleration_out_3d_<< temps_<< " ";
      for(int k=0; k<nb_output_points; k++)
        {
          for(int i=0; i<3; i++)
            {
              displacement_out_3d_<<displacement(k, i)<<" ";
              speed_out_3d_<<velocity(k, i)<<" ";
              acceleration_out_3d_<<acceleration(k, i)<<" ";
            }
        }
      displacement_out_3d_<<finl;
      speed_out_3d_<<finl;
      acceleration_out_3d_<<finl;
    }
}

void Beam_model::printOutputFluidForceOnBeam(bool first_writing) const
{
  if (je_suis_maitre()) // Write the result in the ModalForceFluide1D.txt file
    {
      Nom filename(beamName_);
      filename+="_ModalForceFluide1D.out";
      if (!fluidForceOnBeam_out_.is_open())
        {
          fluidForceOnBeam_out_.ouvrir(filename, (first_writing?ios::out:ios::app));
          fluidForceOnBeam_out_.setf(ios::scientific);
        }
      if (first_writing)
        {
          fluidForceOnBeam_out_<< "# Printing modal 1D fluid force: time mode ";
          for(int nbmodes=0; nbmodes<nbModes_; nbmodes++)
            fluidForceOnBeam_out_<<nbmodes+1<<" ";
          fluidForceOnBeam_out_<<finl;
        }
      fluidForceOnBeam_out_<< temps_<< " ";
      for(int nbmodes=0; nbmodes<nbModes_; nbmodes++)
        fluidForceOnBeam_out_<<fluidForceOnBeam_[nbmodes]<<" ";
      fluidForceOnBeam_out_<<endl;
    }
}

void Beam_model::setCenterCoordinates(const double& x0,const double& y0, const double& z0)
{

  x0_=x0;
  y0_=y0;
  z0_=z0;
}

DoubleVect& Beam_model::getVelocity(const double& tps, const double& dt)
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
DoubleVect& Beam_model::NewmarkScheme (const double& dt)
{

  double squareDt=dt*dt;
  for(int j=0; j < nbModes_; j++)
    {
      double qDispl_pred= qDisplacement_[j] + dt*qSpeed_[j] + squareDt*(0.5-beta_)*qAcceleration_[j];
      double qSpeed_pred= qSpeed_[j] + dt*(1-gamma_)*qAcceleration_[j];

      double coeff1 = mass_[j] + dt*gamma_*(1.+ alpha_)*damping_[j] + squareDt*beta_*(1.+ alpha_)*stiffness_[j];
      double coeff2 = (1. + alpha_)*damping_[j]*qSpeed_pred;
      double coeff3 = (1. + alpha_)*stiffness_[j]*qDispl_pred;
      double coeff4 = alpha_*stiffness_[j]*qDisplacement_[j];
      double coeff5 = alpha_*damping_[j]*qSpeed_[j];


      qAcceleration_[j]=(fluidForceOnBeam_[j] - coeff2 - coeff3 + coeff4 + coeff5)/coeff1;
      qDisplacement_[j] = qDispl_pred + squareDt*beta_*qAcceleration_[j];
      qSpeed_[j] = qSpeed_pred + dt*gamma_*qAcceleration_[j];
    }



  saveBeamForRestart();
  if(output_position_1D_.size()>0) printOutputBeam1D();
  if(output_position_3D_.size()>0) printOutputBeam3D();

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
    	            bending_dir_=var_int;
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

