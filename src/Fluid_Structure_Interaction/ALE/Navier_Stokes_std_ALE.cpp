/****************************************************************************
* Copyright (c) 2022, CEA
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
// File:        Navier_Stokes_std_ALE.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/ALE/src/New
//
//////////////////////////////////////////////////////////////////////////////

#include <Navier_Stokes_std_ALE.h>
#include <Probleme_base.h>
#include <Domaine_ALE.h>
#include <Schema_Temps_base.h>
#include <Op_Conv_ALE_VEF.h>
#include <Op_Grad_VEF_P1B_Face.h>
#include <Discretisation_base.h>
#include <TRUSTTrav.h>
#include <Debog.h>
#include <Discret_Thyd.h>
#include <EcritureLectureSpecial.h>
#include <Avanc.h>
#include <TRUST_2_PDI.h>
#include <Discret_Thyd.h>
#include <Champ_Inc_base.h>


Implemente_instanciable(Navier_Stokes_std_ALE,"Navier_Stokes_standard_ALE",Navier_Stokes_std);
// XD Navier_Stokes_std_ALE navier_stokes_standard Navier_Stokes_std_ALE -1 Resolution of hydraulic Navier-Stokes eq. on mobile domain (ALE)

Sortie& Navier_Stokes_std_ALE::printOn(Sortie& os ) const
{
  return Navier_Stokes_std::printOn(os);
}

Entree& Navier_Stokes_std_ALE::readOn(Entree& is )
{
  return Navier_Stokes_std::readOn(is);
}


/*! @brief for PDI IO: retrieve name, type and dimensions of the fields to save/restore
 *
 */
std::vector<YAML_data> Navier_Stokes_std_ALE::data_a_sauvegarder() const
{
  std::vector<YAML_data> data = Navier_Stokes_std::data_a_sauvegarder();
  int nb_dim = vitesse().valeurs().nb_dim(); // Initialize with same discretization

  // Helper lambda to create YAML_data with uppercase name
  auto make_yaml = [&](const std::string& suffix) -> YAML_data
  {
    std::string name = probleme().le_nom().getString() + suffix;
    std::transform(name.begin(), name.end(), name.begin(), ::toupper);
    return YAML_data(name, "double", nb_dim);
  };

  // Add base data
  data.push_back(make_yaml("_JacobianOld"));
  data.push_back(make_yaml("_JacobianNew"));
  data.push_back(make_yaml("_MeshCoords"));

  const Domaine_ALE& dom_ale = ref_cast(Domaine_ALE, probleme().domaine());
  if (dom_ale.getMeshMotionModel() == 1)
    {
      // Add ALE-specific data
      data.push_back(make_yaml("_MeshDisplacement"));
      data.push_back(make_yaml("_MeshVelocity"));
      data.push_back(make_yaml("_MeshAcceleration"));
      data.push_back(make_yaml("_MeshPosition"));
      data.push_back(make_yaml("_MeshReferenceConfiguration"));
      data.push_back(make_yaml("_MeshTransformationGradient"));
      data.push_back(make_yaml("_MeshStress"));
    }

  return data;
}

int Navier_Stokes_std_ALE::sauvegarder(Sortie& os) const
{
  int bytes = 0, a_faire, special;
  bytes += Navier_Stokes_std::sauvegarder(os);
  EcritureLectureSpecial::is_ecriture_special(special, a_faire);

  const Domaine_ALE& dom_ale = ref_cast(Domaine_ALE, probleme().domaine());
  const Discret_Thyd& dis = ref_cast(Discret_Thyd, discretisation());
  double temps = schema_temps().temps_courant();
  Noms noms(1), unit(1);
  unit[0] = "none";

  auto create_champ = [&](const std::string& name, int nbComp, bool elem=false, const DoubleTab& val = DoubleTab()) -> OWN_PTR(Champ_Inc_base)
  {
    OWN_PTR(Champ_Inc_base) champ;
    noms[0] = name;
    dis.discretiser_champ(elem ? "champ_elem" : "champ_sommets",
                          domaine_dis(),
                          vectoriel,
                          noms,
                          unit,
                          nbComp,
                          1,
                          temps,
                          champ);
    champ->associer_eqn(*this);
    champ->valeurs() = val;
    return champ;
  };

  if (a_faire)
    {
      // Jacobians and mesh coordinates
      auto JacobianOld = create_champ("JacobianOld", vitesse().valeurs().nb_dim(), false, dom_ale.getOldJacobian());
      auto JacobianNew = create_champ("JacobianNew", vitesse().valeurs().nb_dim(), false, dom_ale.getNewJacobian());
      auto meshCoords  = create_champ("meshCoords", dimension, false, domaine_dis().domaine().les_sommets());

      if (special && Process::nproc() > 1)
        Cerr << "ATTENTION: For a parallel calculation, the field Jacobian is not saved in xyz format ..." << finl;
      else
        {
          bytes += JacobianOld->sauvegarder(os);
          bytes += JacobianNew->sauvegarder(os);
          bytes += meshCoords->sauvegarder(os);
        }

      if (dom_ale.getMeshMotionModel() == 1)
        {
          int nbn = (dimension == 2) ? 3 : 4;
          int nSymSize = (dimension == 2) ? 5 : 9;
          int symSize = (dimension == 2) ? 4 : 6;

          auto meshDisplacement           = create_champ("meshDisplacement", dimension, false, dom_ale.getMeshDisplacement());
          auto meshVelocity               = create_champ("meshVelocity", dimension, false, dom_ale.getMeshVelocity());
          auto meshAcceleration           = create_champ("meshAcceleration", dimension, false, dom_ale.getMeshAcceleration());
          auto meshPosition               = create_champ("meshPosition", dimension, false, dom_ale.getMeshPosition());
          auto meshReferenceConfiguration = create_champ("meshReferenceConfiguration", dimension*nbn, true, dom_ale.getMeshReferenceConfiguration());
          auto meshTransformationGradient = create_champ("meshTransformationGradient", nSymSize, true, dom_ale.getMeshTransformationGradient());
          auto meshStress                 = create_champ("meshStress", symSize, true, dom_ale.getMeshStress());

          auto save_field = [&](OWN_PTR(Champ_Inc_base)& champ)
          {
            if (!(special && Process::nproc() > 1))
              bytes += champ->sauvegarder(os);
          };

          save_field(meshDisplacement);
          save_field(meshVelocity);
          save_field(meshAcceleration);
          save_field(meshPosition);
          save_field(meshReferenceConfiguration);
          save_field(meshTransformationGradient);
          save_field(meshStress);
        }
    }
  else if (TRUST_2_PDI::is_PDI_checkpoint())
    {
      auto write_tab = [](const DoubleTab& d, const Nom& name) -> int
      {
        TRUST_2_PDI pdi_interface;
        pdi_interface.share_TRUSTTab_dimensions(d, name, 1);
        if (d.dimension_tot(0))
          pdi_interface.TRUST_start_sharing(name.getString(), d.addr());
        else
          {
            ArrOfDouble garbage(d.nb_dim());
            pdi_interface.TRUST_start_sharing(name.getString(), garbage.addr());
          }
        return 8 * d.size_array();
      };

      bytes += write_tab(dom_ale.getOldJacobian(), (probleme().le_nom() + "_JacobianOld").majuscule());
      bytes += write_tab(dom_ale.getNewJacobian(), (probleme().le_nom() + "_JacobianNew").majuscule());
      bytes += write_tab(domaine_dis().domaine().les_sommets(), (probleme().le_nom() + "_meshCoords").majuscule());

      if (dom_ale.getMeshMotionModel() == 1)
        {
          bytes += write_tab(dom_ale.getMeshDisplacement(), (probleme().le_nom() + "_meshDisplacement").majuscule());
          bytes += write_tab(dom_ale.getMeshVelocity(), (probleme().le_nom() + "_meshVelocity").majuscule());
          bytes += write_tab(dom_ale.getMeshAcceleration(), (probleme().le_nom() + "_meshAcceleration").majuscule());
          bytes += write_tab(dom_ale.getMeshPosition(), (probleme().le_nom() + "_meshPosition").majuscule());
          bytes += write_tab(dom_ale.getMeshReferenceConfiguration(), (probleme().le_nom() + "_meshReferenceConfiguration").majuscule());
          bytes += write_tab(dom_ale.getMeshTransformationGradient(), (probleme().le_nom() + "_meshTransformationGradient").majuscule());
          bytes += write_tab(dom_ale.getMeshStress(), (probleme().le_nom() + "_meshStress").majuscule());
        }
    }

  return bytes;
}



int Navier_Stokes_std_ALE::reprendre(Entree& is)
{
  // Base class resumption
  Navier_Stokes_std::reprendre(is);

  Domaine_ALE& dom_ale = ref_cast(Domaine_ALE, probleme().domaine());
  const Discret_Thyd& dis = ref_cast(Discret_Thyd, discretisation());
  double temps = schema_temps().temps_courant();
  Noms noms(1), unit(1);
  unit[0] = "none";

  auto create_field = [&](const std::string& name, int nbComp, bool elem=false) -> OWN_PTR(Champ_Inc_base)
  {
    OWN_PTR(Champ_Inc_base) champ;
    noms[0] = name;
    dis.discretiser_champ(elem ? "champ_elem" : "champ_sommets",
                          domaine_dis(),
                          vectoriel,
                          noms,
                          unit,
                          nbComp,
                          1,
                          temps,
                          champ);
    champ->associer_eqn(*this);

    // Create the tag for resumption
    Nom tag(champ->le_nom());
    tag += champ->que_suis_je();
    tag += probleme().domaine().le_nom();
    tag += Nom(temps, probleme().reprise_format_temps());
    return champ;
  };

  // Jacobians
  auto JacobianOld = create_field("JacobianOld", vitesse().valeurs().nb_dim());
  auto JacobianNew = create_field("JacobianNew", vitesse().valeurs().nb_dim());

  // Mesh coordinates
  auto meshCoords = create_field("meshCoords", dimension);

  if (EcritureLectureSpecial::is_lecture_special() && Process::nproc() > 1)
    {
      Cerr << "Error in Navier_Stokes_std_ALE::reprendre !" << finl;
      Cerr << "Use the sauv file to resume a parallel Navier_Stokes_std_ALE calculation (Jacobian is required) ... " << finl;
      Process::exit();
    }
  else
    {
      if(!TRUST_2_PDI::is_PDI_restart()) avancer_fichier(is, Nom(JacobianOld->le_nom() + JacobianOld->que_suis_je() + probleme().domaine().le_nom() + Nom(temps, probleme().reprise_format_temps())));
      JacobianOld->reprendre(is);

      if(!TRUST_2_PDI::is_PDI_restart()) avancer_fichier(is, Nom(JacobianNew->le_nom() + JacobianNew->que_suis_je() + probleme().domaine().le_nom() + Nom(temps, probleme().reprise_format_temps())));
      JacobianNew->reprendre(is);

      if(!TRUST_2_PDI::is_PDI_restart()) avancer_fichier(is, Nom(meshCoords->le_nom() + meshCoords->que_suis_je() + probleme().domaine().le_nom() + Nom(temps, probleme().reprise_format_temps())));
      meshCoords->reprendre(is);
    }

  // Set resumption field
  dom_ale.resumptionJacobianCoords(JacobianOld->valeurs(), JacobianNew->valeurs(), meshCoords->valeurs());
  dom_ale.updateMetrics(domaine_dis(), probleme());

  // Mesh motion fields
  if(dom_ale.getMeshMotionModel() == 1)
    {
      int nbn = (dimension == 2) ? 3 : 4;
      int nSymSize = (dimension == 2) ? 5 : 9;
      int symSize = (dimension == 2) ? 4 : 6;

      auto meshDisplacement          = create_field("meshDisplacement", dimension);
      auto meshVelocity              = create_field("meshVelocity", dimension);
      auto meshAcceleration          = create_field("meshAcceleration", dimension);
      auto meshPosition              = create_field("meshPosition", dimension);
      auto meshReferenceConfiguration= create_field("meshReferenceConfiguration", dimension*nbn, true);
      auto meshTransformationGradient= create_field("meshTransformationGradient", nSymSize, true);
      auto meshStress                = create_field("meshStress", symSize, true);

      auto resume_field = [&](OWN_PTR(Champ_Inc_base)& champ)
      {
        if(!TRUST_2_PDI::is_PDI_restart())
          {
            Nom tag(champ->le_nom());
            tag += champ->que_suis_je();
            tag += probleme().domaine().le_nom();
            tag += Nom(temps, probleme().reprise_format_temps());
            avancer_fichier(is, tag);
          }
        champ->reprendre(is);
      };

      resume_field(meshDisplacement);
      resume_field(meshVelocity);
      resume_field(meshAcceleration);
      resume_field(meshPosition);
      resume_field(meshReferenceConfiguration);
      resume_field(meshTransformationGradient);
      resume_field(meshStress);

      dom_ale.resumptionStructuralDynamicsMesh(
        temps,
        meshDisplacement->valeurs(),
        meshVelocity->valeurs(),
        meshAcceleration->valeurs(),
        meshPosition->valeurs(),
        meshReferenceConfiguration->valeurs(),
        meshTransformationGradient->valeurs(),
        meshStress->valeurs()
      );
    }

  return 1;
}



void Navier_Stokes_std_ALE::renewing_jacobians( DoubleTab& derivee )
{
  //Renewing ALE Jacobians
  int TimeStepNr=probleme().schema_temps().nb_pas_dt();
  Domaine_ALE& dom_ale=ref_cast(Domaine_ALE, probleme().domaine());

  DoubleTab New_ALEjacobian_Old=dom_ale.getNewJacobian(); //New  value for ALEjacobian_old
  DoubleTab New_ALEjacobian_New(New_ALEjacobian_Old);

  Op_Conv_ALE_VEF& opALEforJacob=ref_cast(Op_Conv_ALE_VEF, terme_convectif.valeur());
  opALEforJacob.calculateALEjacobian(New_ALEjacobian_New); //New  value for ALEjacobian_new
  dom_ale.update_ALEjacobians(New_ALEjacobian_Old, New_ALEjacobian_New, TimeStepNr); // Update new values of ALEjacobian_old and ALEjacobian_new saved in Domaine_ALE

  //End of renewing ALE Jacobians

  Nom discr=discretisation().que_suis_je();
  if (discr != "VEFPreP1B")
    {
      Cerr<<"volume_entrelace_Cl used in the mass matrix is wrong for the ALE treatment on the boundaries"<<finl;
      Cerr<<"(vol_entrelace_Cl=0 for some of the boundaries indeed a correction is necessary) :"<<finl;
      Cerr<<"the VEFPreP1B discretization must be used to avoid this problem. "<<finl;
      Process::exit();
    }
  Cerr << "Adding ALE contribution..." << finl;
  Op_Conv_ALE& opale=ref_cast(Op_Conv_ALE, terme_convectif.valeur());
  DoubleTrav ALE(derivee); // copie de la structure, initialise a zero
  opale.ajouterALE(la_vitesse->valeurs(), ALE);
  ALE.echange_espace_virtuel();
  solveur_masse->appliquer(ALE);
  ALE.echange_espace_virtuel();
  derivee+=ALE; // M-1(F + ALEconvectiveTerm - BtP(n))=derivee_withALEconvectiveTerm
  derivee.echange_espace_virtuel();
  //Cerr << "ALE => norme(derivee) = " << mp_norme_vect(derivee) << finl;
  Debog::verifier("derivee_pression Navier_Stokes_std::corriger_derivee_impl",derivee);
}

void Navier_Stokes_std_ALE::div_ale_derivative( DoubleTrav& deriveeALE, double timestep, DoubleTab& derivee, DoubleTrav& secmemP )
{
  Domaine_ALE& dom_ale=ref_cast(Domaine_ALE, probleme().domaine());
  DoubleTab ALEjacobian_Old=dom_ale.getOldJacobian();
  DoubleTab ALEjacobian_New=dom_ale.getNewJacobian();
  DoubleTab& vitesse_faces_ALE= dom_ale.vitesse_faces();
  DoubleTab term_Jacobian_ratio_U_n(la_vitesse->valeurs());               // For (Jacobian n/Jacobian n+1) * Un / timestep , initialized every iteration with with new la_vitesse values.

  for (int num_face=0; num_face<(vitesse_faces_ALE.size()/dimension); num_face++)
    {
      for (int dim=0; dim<dimension; dim++)
        {
          term_Jacobian_ratio_U_n(num_face,dim)*=(ALEjacobian_Old(num_face,dim)/(ALEjacobian_New(num_face,dim)*timestep));
          deriveeALE(num_face,dim)=term_Jacobian_ratio_U_n(num_face,dim)+derivee(num_face,dim); //(J_{n}/J_{n+1})*(U_{n}/timestep)+derivee_out
        }
    }
  deriveeALE.echange_espace_virtuel();

  divergence.calculer(deriveeALE, secmemP); //Div((J_{n}/J_{n+1})*(U_{n}/timestep)+derivee_out)
  secmemP *= -1; // car div =-B
  secmemP.echange_espace_virtuel();
  //Debog::verifier("secmemP  modifier Navier_Stokes_std::corriger_derivee_impl",secmemP);
  // Correction du second membre d'apres les conditions aux limites :
  assembleur_pression_->modifier_secmem(secmemP);
  secmemP.echange_espace_virtuel();

  Debog::verifier("secmemP Navier_Stokes_std::corriger_derivee_impl",secmemP);
}

void Navier_Stokes_std_ALE::update_pressure_matrix()
{
  // BM-1Bt matrix is assembled.
  assembleur_pression_->assembler(matrice_pression_);
  solveur_pression_->reinit();

}

void Navier_Stokes_std_ALE::discretiser()
{
  Navier_Stokes_std::discretiser();
  const Discret_Thyd& dis=ref_cast(Discret_Thyd, discretisation());
  Cerr << "Mesh Velocity discretization" << finl;
  dis.discretiser_champ("vitesse", domaine_dis(), "ALEMeshVelocity","m/s", dimension,1,schema_temps().temps_courant(), ALEMeshVelocity_);
  champs_compris_.ajoute_champ(ALEMeshVelocity_);
  ALEMeshVelocity_->add_synonymous(Nom("ALEMeshVelocity"));
  Cerr << "Mesh Velocity discretization" << finl;
  dis.discretiser_champ("vitesse",domaine_dis(),"ALEMeshTotalDisplacement","m/s",dimension,1,schema_temps().temps_courant(),ALEMeshTotalDisplacement_);
  champs_compris_.ajoute_champ(ALEMeshTotalDisplacement_);
  ALEMeshTotalDisplacement_->add_synonymous(Nom("ALEMeshTotalDisplacement"));

  const Domaine_ALE& domaine_ALE=ref_cast(Domaine_ALE, probleme().domaine()) ;
  if (domaine_ALE.getMeshMotionModel() == 1)
    {
      dis.discretiser_champ("champ_elem",domaine_dis(),"ALEMeshStructuralPressure","Pa",1,1,schema_temps().temps_courant(),ALEMeshStructuralPressure_);
      champs_compris_.ajoute_champ(ALEMeshStructuralPressure_);
      ALEMeshStructuralPressure_.valeur().add_synonymous(Nom("ALEMeshStructuralPressure"));
      Cerr << "Mesh Fictitious Structural Pressure discretization" << finl;
      dis.discretiser_champ("champ_elem",domaine_dis(),"ALEMeshStructuralVonMises","Pa",1,1,schema_temps().temps_courant(),ALEMeshStructuralVonMises_);
      champs_compris_.ajoute_champ(ALEMeshStructuralVonMises_);
      ALEMeshStructuralVonMises_.valeur().add_synonymous(Nom("ALEMeshStructuralVonMises"));
      Cerr << "Mesh Fictitious Structural Von Mises discretization" << finl;
      //Rq.: "champ_sommets" would have seemed a better type for ALEMeshStructuralForces field, but it raises an error with unknown "Champ_P1_VEF" type
      //     in VEF field discretization ; to investigate ?
      dis.discretiser_champ("vitesse",domaine_dis(),"ALEMeshStructuralForces","N",dimension,1,schema_temps().temps_courant(),ALEMeshStructuralForces_);
      champs_compris_.ajoute_champ(ALEMeshStructuralForces_);
      ALEMeshStructuralForces_.valeur().add_synonymous(Nom("ALEMeshStructuralForces"));
      Cerr << "Mesh Fictitious Structural Internal Forces discretization" << finl;
    }
}

void Navier_Stokes_std_ALE::mettre_a_jour(double temps)
{
  Navier_Stokes_std::mettre_a_jour(temps);
  if(temps>0.)
    {
      const Domaine_ALE& dom_ale=ref_cast(Domaine_ALE, probleme().domaine());
      const DoubleTab& ALEMeshVelocity= dom_ale.vitesse_faces();//we access the mesh speed
      ALEMeshVelocity_->valeurs()= ALEMeshVelocity;
      double dt = schema_temps().pas_de_temps();
      DoubleTab ALEMeshVelocity_dt= ALEMeshVelocity;
      for(int dim=0; dim<dimension; dim++)
        {
          for(int i=0; i<(ALEMeshVelocity.size()/dimension); i++)
            {
              ALEMeshVelocity_dt(i, dim) *= dt;
            }
        }
      ALEMeshVelocity_dt.echange_espace_virtuel();
      ALEMeshTotalDisplacement_->valeurs() +=ALEMeshVelocity_dt;
      //ALEMeshTotalDisplacement_->valeurs().echange_espace_virtuel();
      //ALEMeshVelocity_->valeurs().echange_espace_virtuel();
      ALEMeshVelocity_->mettre_a_jour(temps);
      ALEMeshTotalDisplacement_->mettre_a_jour(temps);

      if (dom_ale.getMeshMotionModel() == 1)
        {
          const DoubleVect& meshPbPressure = dom_ale.getMeshPbPressure(); //we access the cell pressure in the fictitious structure problem for the mesh
          ALEMeshStructuralPressure_->valeurs() = meshPbPressure ;
          const DoubleVect& meshPbVonMises = dom_ale.getMeshPbVonMises(); //we access the cell von mises stress in the fictitious structure problem for the mesh
          ALEMeshStructuralVonMises_->valeurs() = meshPbVonMises ;
          const DoubleTab& meshPbForceFace = dom_ale.getMeshPbForceFace(); //we access the internal forces in the fictitious structure problem for the mesh
          ALEMeshStructuralForces_->valeurs() = meshPbForceFace ;
        }

    }

}


void Navier_Stokes_std_ALE::updateFluidForce(DoubleTab& velocity)
{
  // in case of implicit coupling with a structural code: update the fluxes (used for computing the fluid force) during implicit sub-iterations

  Domaine_ALE& dom_ale=ref_cast(Domaine_ALE, probleme().domaine());
  if(dom_ale.getCouplingMethod()) //implicit coupling case
    {
      Cout<<" Implicit coupling: Navier_Stokes_std_ALE::updateFluidForce "<<finl;
      //update diffusion operator
      DoubleTab field_value = velocity;
      field_value = 0.;
      operateur_diff().ajouter(velocity, field_value);


      //update gradient operator
      pression().mettre_a_jour(schema_temps().temps_courant());
      calculer_la_pression_en_pa();
      Op_Grad_VEF_P1B_Face& op_grad_vef = ref_cast(Op_Grad_VEF_P1B_Face, operateur_gradient().l_op_base());
      op_grad_vef.calculer_flux_bords();
    }
}
void Navier_Stokes_std_ALE::setPressureTimeN()
{

  Domaine_ALE& dom_ale=ref_cast(Domaine_ALE, probleme().domaine());
  if(dom_ale.getCouplingMethod()) //implicit coupling case

    // Equivalent of setting present = past for velocity done in for eg. Schema_Euler_Implicite::test_stationnaire.
    // Ensures sub-iterations start with the correct pressure for implicit coupling
    pression().valeurs()= getPressureTimeN();
}
bool Navier_Stokes_std_ALE::getCouplingInfoForFiltering() const
{
  const Domaine_ALE& dom_ale=ref_cast(Domaine_ALE, probleme().domaine());
  return dom_ale.getCouplingMethod();
}
