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

#ifndef ANTIOCH_ISOTHERMAL_STIRRED_REACTOR_H
#define ANTIOCH_ISOTHERMAL_STIRRED_REACTOR_H

namespace Antioch
{
  template<typename CoeffType=double, typename StateType=CoeffType>
  class IsothermalStirredReactor : public StirredReactorBase
  {

  public:

    //! Constructor.
    IsothermalStirredReactor( const StateType& T,
                              const ReactionSet<CoeffType>& reaction_set,
                              const CEAThermoMixture<CoeffType>& thermo,
                              StirredReactorTimeIntegratorBase<CoeffType,StateType>& time_integrator,
                              CoeffType volume = 1.0, /* m^3 */
                              const StateType example );

    virtual ~IsothermalStirredReactor();

    template<typename VectorStateType>
    void operator()( const VectorStateType& x,
                     VectorStateType& dx_dt );

  protected:


    //! Cache temperature and derived quantities
    /*! Isothermal reactor has a fixed temperature so we optimize
        by not having to recompute T^2, ln(T), etc. */
    const TempCache<StateType> _cache;
    

  private:

    IsothermalStirredReactor();

  };

  /* ---------------------- Constructor/Destructor ----------------------*/
  template<typename CoeffType, typename StateType>
  inline
  IsothermalStirredReactor<CoeffType,StateType>::IsothermalStirredReactor
  ( const StateType& T,
    const ReactionSet<CoeffType>& reaction_set,
    const CEAThermoMixture<CoeffType>& thermo,
    StirredReactorTimeIntegratorBase<CoeffType,StateType>& time_integrator,
    CoeffType volume = 1.0, /* m^3 */
    const StateType example )
    : StirredReactorBase<CoeffType,StateType>(reaction_set, thermo, time_integrator,
                                              volume, example),
      _cache(T)
  {
    return;
  }

  template<typename CoeffType, typename StateType>
  inline
  IsothermalStirredReactor<CoeffType,StateType>::~IsothermalStirredReactor()
  {
    return;
  }

  template<typename CoeffType, typename StateType>
  template<typename VectorStateType>
  inline
  void IsothermalStirredReactor<CoeffType,StateType>::operator()
    ( const VectorStateType& x,
      VectorStateType& dx_dt )
  {
    VectorStateType h_RT_minus_s_R = zero_clone(x);

    _thermo_evaluator.h_RT_minus_s_R( _cache, h_RT_minus_s_R );

    // Here, x is the vector of molar densities
    _kinetics_evaluator.compute_mole_sources( _cache.T, x, h_RT_minus_s_R, dx_dt );

    return;
  }
  
} // end namespace Antioch

#endif // ANTIOCH_ISOTHERMAL_STIRRED_REACTOR_H
