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

#ifndef ANTIOCH_ISOTHERMAL_ISOBARIC_STIRRED_REACTOR_H
#define ANTIOCH_ISOTHERMAL_ISOBARIC_STIRRED_REACTOR_H

// Antioch
#include "antioch/stirred_reactor_enum.h"
#include "antioch/isothermal_stirred_reactor.h"

namespace Antioch
{
  template<typename CoeffType=double, typename StateType=CoeffType>
  class IsothermalIsobaricStirredReactor : public IsothermalStirredReactor<CoeffType,StateType>
  {

  public:

    //! Constructor.
    IsothermalIsobaricStirredReactor( const StateType& T,
                                              const ReactionSet<CoeffType>& reaction_set,
                                              const CEAThermoMixture<CoeffType>& thermo,
                                              StirredReactorTimeIntegratorBase<CoeffType,StateType>& time_integrator,
                                              const StateType example,
                                              const CoeffType pressure );
    
    virtual ~IsothermalIsobaricStirredReactor();

    template<typename VectorStateType>
    void operator()( const VectorStateType& x,
                     VectorStateType& dx_dt );

  private:

    IsothermalIsobaricStirredReactor();

  };

  /* ---------------------- Constructor/Destructor ----------------------*/
  template<typename CoeffType, typename StateType>
  inline
  IsothermalIsobaricStirredReactor<CoeffType,StateType>::IsothermalIsobaricStirredReactor
  ( const StateType& T,
    const ReactionSet<CoeffType>& reaction_set,
    const CEAThermoMixture<CoeffType>& thermo,
    StirredReactorTimeIntegratorBase<CoeffType,StateType>& time_integrator,
    const StateType example,
    const CoeffType pressure )
    : IsothermalStirredReactor<CoeffType,StateType>(T, reaction_set, thermo, time_integrator, example)
  {
    // Reset reactor_type to correct value
    this->_reactor_type = ReactorType::ISOTHERMAL_ISOBARIC;

    return;
  }

  template<typename CoeffType, typename StateType>
  inline
  IsothermalIsobaricStirredReactor<CoeffType,StateType>::~IsothermalIsobaricStirredReactor()
  {
    return;
  }

  template<typename CoeffType, typename StateType>
  template<typename VectorStateType>
  inline
  void IsothermalIsobaricStirredReactor<CoeffType,StateType>::operator()
    ( const VectorStateType& x,
      VectorStateType& dx_dt )
  {
    VectorStateType h_RT_minus_s_R = zero_clone(x);

    this->_thermo_evaluator.h_RT_minus_s_R( this->_cache, h_RT_minus_s_R );

    // Here, x is the vector of molar densities
    this->_kinetics_evaluator.compute_mole_sources( this->_cache.T, x, h_RT_minus_s_R, dx_dt );

    const unsigned int n_species = this->_kinetics_evaluator.n_species();

    dx_dt[n_species-1] = 0.0;

    for( unsigned int s = 0; s < n_species-1; s++ )
      {
        dx_dt[n_species-1] += -dx_dt[s];
      }

    return;
  }
  
} // end namespace Antioch

#endif // ANTIOCH_ISOTHERMAL_ISOBARIC_STIRRED_REACTOR_H
