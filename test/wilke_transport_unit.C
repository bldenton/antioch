//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014-2016 Paul T. Bauman, Benjamin S. Kirk,
//                         Sylvain Plessis, Roy H. Stonger
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// C++
#include <iostream>
#include <cmath>

// Antioch
#include "antioch_config.h"
#include "antioch/vector_utils_decl.h"

#include "antioch/default_filename.h"
#include "antioch/transport_mixture.h"
#include "antioch/stat_mech_thermo.h"
#include "antioch/cea_curve_fit.h"
#include "antioch/nasa_mixture.h"
#include "antioch/nasa_evaluator.h"
#include "antioch/nasa_mixture_parsing.h"

#include "antioch/eucken_thermal_conductivity.h"
#include "antioch/eucken_thermal_conductivity_building.h"
#include "antioch/blottner_viscosity.h"
#include "antioch/constant_lewis_diffusivity.h"
#include "antioch/kinetics_theory_thermal_conductivity.h"
#include "antioch/kinetics_theory_thermal_conductivity_building.h"

#ifdef ANTIOCH_HAVE_GSL

#include "antioch/kinetics_theory_viscosity.h"
#include "antioch/kinetics_theory_viscosity_building.h"
#include "antioch/molecular_binary_diffusion.h"
#include "antioch/gsl_spliner.h"
#include "antioch/mixture_diffusion.h"

#endif
//
#include "antioch/wilke_mixture.h"  // backward compatiblity
#include "antioch/wilke_evaluator.h"  // backward compatiblity
#include "antioch/mixture_averaged_transport_mixture.h"
#include "antioch/mixture_averaged_transport_evaluator.h"

#include "antioch/blottner_parsing.h"
#include "antioch/eucken_thermal_conductivity_building.h"
#include "antioch/constant_lewis_diffusivity_building.h"

#include "antioch/vector_utils.h"

template <typename Scalar>
int test_val( const Scalar val, const Scalar val_exact, const Scalar tol, const std::string& val_name )
{
  using std::abs;

  int return_flag = 0;

  const Scalar rel_error = abs( (val - val_exact)/val_exact);

  if( rel_error  > tol )
    {
      std::cerr << std::setprecision(20) << std::scientific
                << "Error: Mismatch in " << val_name << std::endl
		<< val_name << "    = " << val << std::endl
		<< val_name+"_exact = " << val_exact << std::endl
		<< "rel_error = " << rel_error << std::endl
		<< "tol = " << tol << std::endl;
      return_flag = 1;
    }

  return return_flag;
}

template <typename Scalar>
int tester()
{
  using std::pow;

  std::vector<std::string> species_str_list;
  const unsigned int n_species = 5;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "O2" );
  species_str_list.push_back( "N" );
  species_str_list.push_back( "O" );
  species_str_list.push_back( "NO" );

// mixture and thermo for conduction
  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );

// micro thermo
  Antioch::StatMechThermodynamics<Scalar> thermo_stat( chem_mixture );

  typedef Antioch::StatMechThermodynamics<Scalar> MicroThermo;

// macro thermo for cp (diffusion)
  Antioch::NASAThermoMixture<Scalar,Antioch::CEACurveFit<Scalar> > cea_mixture( chem_mixture );
  Antioch::read_nasa_mixture_data( cea_mixture, Antioch::DefaultFilename::thermo_data(),Antioch::ASCII, true );
  Antioch::NASAEvaluator<Scalar,Antioch::CEACurveFit<Scalar> > thermo_mix( cea_mixture );

  Antioch::TransportMixture<Scalar> tran_mixture( chem_mixture );

//
  Antioch::MixtureConductivity<Antioch::EuckenThermalConductivity<MicroThermo>,Scalar>
    k( tran_mixture );
  Antioch::build_eucken_thermal_conductivity<MicroThermo,Scalar>(k,thermo_stat);


  Antioch::MixtureViscosity<Antioch::BlottnerViscosity<Scalar>,Scalar> mu( tran_mixture );
  Antioch::read_blottner_data_ascii( mu, Antioch::DefaultFilename::blottner_data() );


// pure species set, all internally set
#ifdef ANTIOCH_HAVE_GSL
  Antioch::MixtureViscosity<Antioch::KineticsTheoryViscosity<Scalar,Antioch::GSLSpliner>,Scalar >
    ps_mu(tran_mixture);
  Antioch::build_kinetics_theory_viscosity<Scalar,Antioch::GSLSpliner>(ps_mu);

  Antioch::MixtureDiffusion<Antioch::MolecularBinaryDiffusion<Scalar,Antioch::GSLSpliner>,Scalar>
    bimol_D( tran_mixture );

#endif

  Antioch::MixtureConductivity<Antioch::KineticsTheoryThermalConductivity<MicroThermo,Scalar>,Scalar>
    ps_k( tran_mixture );

  Antioch::build_kinetics_theory_thermal_conductivity<MicroThermo,Scalar>(ps_k,thermo_stat);

//Eucken is internally set

  Antioch::MixtureDiffusion<Antioch::ConstantLewisDiffusivity<Scalar>,Scalar>
    D( tran_mixture );

  const Scalar Le = 1.4;
  Antioch::build_constant_lewis_diffusivity<Scalar>( D, Le);

// non kinetics theory
  Antioch::MixtureAveragedTransportMixture<Scalar> wilke_mixture( tran_mixture );

  Antioch::MixtureAveragedTransportEvaluator<Antioch::ConstantLewisDiffusivity<Scalar>,
                                             Antioch::BlottnerViscosity<Scalar>,
                                             Antioch::EuckenThermalConductivity<MicroThermo>,
                                             Scalar>
    wilke( wilke_mixture, D, mu, k );


// kinetics theory full
#ifdef ANTIOCH_HAVE_GSL

  Antioch::MixtureAveragedTransportEvaluator<Antioch::MolecularBinaryDiffusion<Scalar,Antioch::GSLSpliner>,
                                             Antioch::KineticsTheoryViscosity<Scalar,Antioch::GSLSpliner>,
                                             Antioch::KineticsTheoryThermalConductivity<MicroThermo,Scalar>,
                                             Scalar>
    wilke_ps_evaluator(wilke_mixture,bimol_D,ps_mu,ps_k);

#endif

  Antioch::WilkeMixture<Scalar> wilke_mix_dep( chem_mixture );

  Antioch::EuckenThermalConductivity<MicroThermo> k_eucken(thermo_stat);

  Antioch::WilkeEvaluator<Antioch::MixtureViscosity<Antioch::BlottnerViscosity<Scalar>,Scalar>,
                          Antioch::EuckenThermalConductivity<MicroThermo>,
                          Scalar> wilke_eval_dep( wilke_mix_dep, mu, k_eucken );

  int return_flag = 0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 5;

  // First check phi calculation
  // phi_s = sum_r (chi_r*(1+sqrt(mu_s/mu_r)*(Mr/Ms)^(1/4))^2)/sqrt(8*(1+Ms/Mr))
  {
    std::vector<Scalar> mu( 5 );
    mu[0] = 0.1L;
    mu[1] = 0.2L;
    mu[2] = 0.3L;
    mu[3] = 0.15L;
    mu[4] = 0.25L;

    std::vector<Scalar> chi( 5 );
    chi[0] = 0.1L;
    chi[1] = 0.2L;
    chi[2] = 0.3L;
    chi[3] = 0.15L;
    chi[4] = 0.25L;

    Scalar phi_N_exact = 0.0L;
    unsigned int N_index = 2;
    Scalar M_N = chem_mixture.M(N_index);

    for( unsigned int r = 0; r < 5; r++ )
      {
        Scalar M_r = chem_mixture.M(r);
        Scalar dummy = 1.0L + std::sqrt(mu[N_index]/mu[r])*pow( M_r/M_N, Scalar(0.25L) );
        phi_N_exact += chi[r]*dummy*dummy/std::sqrt(8.0L*( 1.0L + M_N/M_r ) );
      }

    Scalar phi_N;
    std::vector<std::vector<Scalar> > mu_mu_sqrt(mu.size());
    Antioch::init_constant(mu_mu_sqrt,mu);
    Antioch::set_zero(mu_mu_sqrt);

    wilke.compute_mu_mu_sqrt( mu, mu_mu_sqrt);
    phi_N = wilke.compute_phi( mu_mu_sqrt, chi, N_index );

    return_flag = test_val( phi_N, phi_N_exact, tol, std::string("phi") );
  }

  const Scalar P = 1e5;
  std::vector<Scalar> mass_fractions( 5, 0.2L);
  const Scalar T = 1000.0L;
  const Antioch::TempCache<Scalar> T_cache(T);

  const Scalar R_mix = chem_mixture.R(mass_fractions); // get R_tot in J.kg-1.K-1
  const Scalar rho = P/(R_mix*T); // kg.m-3
  const Scalar cp = thermo_mix.cp( T_cache, mass_fractions );

  const Scalar wilke_mu_long_double = 4.51233094078102111066e-05L;
  const Scalar wilke_k_long_double  = 8.01027375195322618301e-02L;

  Scalar wilke_mu = wilke.mu(T, mass_fractions );
  Scalar wilke_k = wilke.k(T, mass_fractions );

  return_flag = test_val( wilke_mu, wilke_mu_long_double, tol, "wilke mixture viscosity") || return_flag;
  return_flag = test_val( wilke_k, wilke_k_long_double, tol, "wilke mixture thermal conduction") || return_flag;

  std::vector<Scalar> lewis_D(5,0);

  wilke.mu_and_k(T,mass_fractions,wilke_mu,wilke_k);

  return_flag = test_val( wilke_mu, wilke_mu_long_double, tol, "wilke mixture viscosity") || return_flag;
  return_flag = test_val( wilke_k, wilke_k_long_double, tol, "wilke mixture thermal conduction") || return_flag;

  wilke.mu_and_k_and_D( T, rho, cp, mass_fractions, wilke_mu, wilke_k, lewis_D );

  Scalar D_lewis_exact = wilke_k/(rho*cp*Le);

  for(unsigned int s = 0; s < lewis_D.size(); s++)
    return_flag = test_val( lewis_D[s], D_lewis_exact, tol, "constant Lewis diffusion for species " + species_str_list[s]) || return_flag;

#if ANTIOCH_HAVE_GSL
/* \todo better the test
   Alright we need something to test, so here's the sorry version.
   It is the long double values, we're thus just verifying that
   future development won't change them.

   the smart test would be to make up molecules where the
   reduced dipole moment fall into a node (easy part, 0 works)
   and the reduced temperature falls also on a node (less easy).
   Then we will have a theoretical independant value.
*/

  const Scalar mu_kt_long_double   = 4.49877527305932602332e-05L;
  const Scalar k_kt_long_double    = 8.22050332419571328635e-02L;

  std::vector<Scalar> D_kt_long_double(5,0);
  D_kt_long_double[0] = 1.95418749838889089562e-04L;
  D_kt_long_double[1] = 1.92376915629762328034e-04L;
  D_kt_long_double[2] = 2.98006144849143296987e-04L;
  D_kt_long_double[3] = 3.08179434829991685679e-04L;
  D_kt_long_double[4] = 1.90508644119203653519e-04L;
  Scalar mu_kt, k_kt;
  std::vector<Scalar> D_kt(5,0);
  wilke_ps_evaluator.mu_and_k_and_D( T, rho, cp, mass_fractions, mu_kt, k_kt, D_kt );


  return_flag = test_val( mu_kt, mu_kt_long_double, tol, "kinetics theory viscosity") || return_flag;
  return_flag = test_val( k_kt, k_kt_long_double, tol, "kinetics theory thermal conduction") || return_flag;
  for(unsigned int s = 0; s < D_kt.size(); s++)
    return_flag = test_val( D_kt[s], D_kt_long_double[s], tol, "kinetics theory diffusion for species " + species_str_list[s]) || return_flag;

  std::vector<Scalar> D_kt_2(5,0);
  wilke_ps_evaluator.D(rho, T, mass_fractions, D_kt_2 );

  for(unsigned int s = 0; s < D_kt.size(); s++)
    return_flag = test_val( D_kt_2[s], D_kt_long_double[s], tol, "kinetics theory diffusion for species " + species_str_list[s]) || return_flag;

  // Mass flux, mass fraction test
  {
    typename Antioch::MixtureAveragedTransportEvaluator<Antioch::MolecularBinaryDiffusion<Scalar,Antioch::GSLSpliner>,
                                                        Antioch::KineticsTheoryViscosity<Scalar,Antioch::GSLSpliner>,
                                                        Antioch::KineticsTheoryThermalConductivity<MicroThermo,Scalar>,
                                                        Scalar>::DiffusivityType
      diff_type = Antioch::MixtureAveragedTransportEvaluator<Antioch::MolecularBinaryDiffusion<Scalar,Antioch::GSLSpliner>,
                                                             Antioch::KineticsTheoryViscosity<Scalar,Antioch::GSLSpliner>,
                                                             Antioch::KineticsTheoryThermalConductivity<MicroThermo,Scalar>,
                                                             Scalar>::MASS_FLUX_MASS_FRACTION;

    std::vector<Scalar> D_kt_mass_mass_long_double(5);
    D_kt_mass_mass_long_double[0] = 2.049224418207810845348e-04L;
    D_kt_mass_mass_long_double[1] = 2.065920312318615510015e-04L;
    D_kt_mass_mass_long_double[2] = 2.550809908204135170620e-04L;
    D_kt_mass_mass_long_double[3] = 2.790464132548291530901e-04L;
    D_kt_mass_mass_long_double[4] = 2.024029424511201457822e-04L;

    std::vector<Scalar> D_kt_mass_mass(5,0);
    wilke_ps_evaluator.D(rho, T, mass_fractions, D_kt_mass_mass, diff_type );

    for(unsigned int s = 0; s < D_kt.size(); s++)
      return_flag = test_val( D_kt_mass_mass[s], D_kt_mass_mass_long_double[s], tol, "kinetics theory diffusion (mass flux, mass fraction) for species " + species_str_list[s]) || return_flag;

  }

  // Mole flux, mole fraction test
  {
    typename Antioch::MixtureAveragedTransportEvaluator<Antioch::MolecularBinaryDiffusion<Scalar,Antioch::GSLSpliner>,
                                                        Antioch::KineticsTheoryViscosity<Scalar,Antioch::GSLSpliner>,
                                                        Antioch::KineticsTheoryThermalConductivity<MicroThermo,Scalar>,
                                                        Scalar>::DiffusivityType
      diff_type = Antioch::MixtureAveragedTransportEvaluator<Antioch::MolecularBinaryDiffusion<Scalar,Antioch::GSLSpliner>,
                                                             Antioch::KineticsTheoryViscosity<Scalar,Antioch::GSLSpliner>,
                                                             Antioch::KineticsTheoryThermalConductivity<MicroThermo,Scalar>,
                                                             Scalar>::MOLE_FLUX_MOLE_FRACTION;

    std::vector<Scalar> D_kt_mole_mole_long_double(5);
    D_kt_mole_mole_long_double[0] = 2.070373009368896292617e-04L;
    D_kt_mole_mole_long_double[1] = 2.083783534307470249142e-04L;
    D_kt_mole_mole_long_double[2] = 2.589403037715613871414e-04L;
    D_kt_mole_mole_long_double[3] = 2.824017881144487186270e-04L;
    D_kt_mole_mole_long_double[4] = 2.042449798230576410557e-04L;

    std::vector<Scalar> D_kt_mole_mole(5,0);
    wilke_ps_evaluator.D(rho, T, mass_fractions, D_kt_mole_mole, diff_type );

    for(unsigned int s = 0; s < D_kt.size(); s++)
      return_flag = test_val( D_kt_mole_mole[s], D_kt_mole_mole_long_double[s], tol, "kinetics theory diffusion (mass flux, mass fraction) for species " + species_str_list[s]) || return_flag;

  }



#endif

  return return_flag;
}

int main()
{
  return ( tester<double>()  ||
           tester<long double>() ||
           tester<float>()
           );
}
