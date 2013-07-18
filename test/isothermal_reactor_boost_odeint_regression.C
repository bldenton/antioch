//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
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

#include "antioch_config.h"

// C++
#include <limits>
#include <string>
#include <vector>

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/chemical_species.h"
#include "antioch/chemical_mixture.h"
#include "antioch/reaction_set.h"
#include "antioch/read_reaction_set_data_xml.h"
#include "antioch/cea_thermo.h"
#include "antioch/stirred_reactor_enum.h"
#include "antioch/isothermal_stirred_reactor.h"
#include "antioch/boost_ode_integrator.h"

#ifdef ANTIOCH_HAVE_BOOST_ODEINT
template <typename Scalar>
int tester(const std::string& input_name)
{
  using std::abs;

  std::vector<std::string> species_str_list;
  const unsigned int n_species = 5;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "O2" );
  species_str_list.push_back( "N" );
  species_str_list.push_back( "O" );
  species_str_list.push_back( "NO" );

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );
  Antioch::ReactionSet<Scalar> reaction_set( chem_mixture );
  Antioch::CEAThermoMixture<Scalar> thermo( chem_mixture );

  Antioch::read_reaction_set_data_xml<Scalar>( input_name, false, reaction_set );

  const Scalar T = 500;

  int return_flag = 0;

  Antioch::BoostODEIntegrator<Scalar,Scalar> integrator( BoostStepperType::RUNGE_KUTTA_FOURTH );

  Antioch::IsothermalStirredReactor<Scalar,Scalar> reactor( T, reaction_set, thermo, integrator, 0.0 /*example*/ );

  std::vector<Scalar> x0(n_species, 0.0);
  x0[2] = 0.5;
  x0[3] = 0.5;

  reactor.run( x0, 0.0, 1.0, 0.01 );

  return return_flag;
}

int main(int argc, char* argv[])
{
  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify reaction set XML input file." << std::endl;
      antioch_error();
    }

  int return_flag = 0;

  return_flag = tester<double>(std::string(argv[1]));

  return return_flag;
}

#else // ANTIOCH_HAVE_BOOST_ODEINT

int main()
{
  /* If Antioch was compiled without Boost ODEInt, then tell Automake
     we're skipping this test. */
  return 77;
}

#endif // ANTIOCH_HAVE_BOOST_ODEINT 
