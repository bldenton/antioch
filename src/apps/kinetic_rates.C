//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014 Paul T. Bauman, Benjamin S. Kirk, Sylvain Plessis,
//                    Roy H. Stonger
// Copyright (C) 2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#include "antioch/vector_utils.h"

#include "antioch/read_reaction_set_data.h"
#include "antioch/nasa_mixture_parsing.h"
#include "antioch/chemkin_parser.h"
#include "antioch/xml_parser.h"
#include "antioch/chemical_mixture.h"
#include "antioch/reaction_set.h"
#include "antioch/nasa_evaluator.h"
#include "antioch/nasa_mixture.h"
#include "antioch/nasa_curve_fit.h"
#include "antioch/kinetics_evaluator.h"

int main( int argc, char* argv[] )
{
  std::string kinetics_filename = std::string(argv[1]);
  std::string thermo_data_filename = std::string(argv[2]);
  std::string species_data_filename = std::string(argv[3]);

  if( species_data_filename == std::string("ANTIOCH_DEFAULT") )
    species_data_filename = Antioch::DefaultFilename::chemical_mixture();

  Antioch::ChemKinParser<double> kinetics_parser( kinetics_filename, true );
  //Antioch::XMLParser<double> kinetics_parser( kinetics_filename, false );
  kinetics_parser.initialize();

  std::vector<std::string> species_names = kinetics_parser.species_list();

  Antioch::ChemicalMixture<double> chem_mixture( species_names,
                                                 false,
                                                 species_data_filename );

  Antioch::ReactionSet<double> reaction_set( chem_mixture );

  Antioch::read_reaction_set_data_chemkin<double>( kinetics_filename, true, reaction_set );
  ///Antioch::read_reaction_set_data_xml<double>( kinetics_filename, false, reaction_set );

  Antioch::NASAThermoMixture<double,Antioch::NASA9CurveFit<double> > thermo_mix( chem_mixture );

  Antioch::read_nasa_mixture_data( thermo_mix,
                                   thermo_data_filename,
                                   Antioch::CHEMKIN,
                                   //Antioch::XML,
                                   true );

  Antioch::NASAEvaluator<double,Antioch::NASA9CurveFit<double> > thermo( thermo_mix );

  const double T = 1500.0; // K
  const double P = 1.0e5; // Pa
  const Antioch::KineticsConditions<double> conditions(T);

  unsigned int n_species = chem_mixture.n_species();

  // Mass fractions
  std::vector<double> Y(n_species,0.0);
  Y[0] = 0.2;
  Y[1] = 0.8;

  const double R_mix = chem_mixture.R(Y); // get R_tot in J.kg-1.K-1

  const double rho = P/(R_mix*T); // kg.m-3

  std::vector<double> molar_densities(n_species,0.0);
  chem_mixture.molar_densities(rho,Y,molar_densities);

  std::vector<double> h_RT_minus_s_R(n_species);

  Antioch::TempCache<double> temp(T);
  thermo.h_RT_minus_s_R(temp,h_RT_minus_s_R);

  Antioch::KineticsEvaluator<double> kinetics( reaction_set, 0 );

  std::vector<double> omega_dot(n_species);

  kinetics.compute_mass_sources( T , molar_densities, h_RT_minus_s_R, omega_dot);

  return 0;
}
