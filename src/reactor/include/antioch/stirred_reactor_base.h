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

#ifndef ANTIOCH_STIRRED_REACTOR_BASE_H
#define ANTIOCH_STIRRED_REACTOR_BASE_H

// Antioch
#include "antioch/stirred_reactor_enum.h"
#include "antioch/chemical_mixture.h"
#include "antioch/reaction_set.h"
#include "antioch/cea_mixture.h"
#include "antioch/kinetics_evaluator.h"
#include "antioch/cea_evaluator.h"

namespace Antioch
{
  // Foward declarations
  template<typename CoeffType, typename StateType>
  class IsothermalStirredReactor;

  template<typename CoeffType, typename StateType>
  class StirredReactorTimeIntegratorBase;
  
  template<typename CoeffType=double, typename StateType=CoeffType>
  class StirredReactorBase
  {

  public:

    //! Constructor.
    StirredReactorBase( const ReactionSet<CoeffType>& reaction_set,
                        const CEAThermoMixture<CoeffType>& thermo,
                        StirredReactorTimeIntegratorBase<CoeffType,StateType>& time_integrator,
                        const StateType example,
                        CoeffType volume = 1.0 /* m^3 */ );

    virtual ~StirredReactorBase();

    template<typename VectorStateType>
    void run( const VectorStateType& x0,
              const CoeffType t0,
              const CoeffType t1,
              const CoeffType dt );
              
    void output( std::ostream& output ) const;

    //! Evaluate net production for ODE solvers
    /*! We can't make this virtual since the function is templated,
        so we manually dispatch to derived classes. */
    template<typename VectorStateType>
    void operator()( const VectorStateType& x,
                     VectorStateType& dx_dt );

  protected:

    ReactorType::ReactorType _reactor_type;

    const ReactionSet<CoeffType>& _reaction_set;

    const ChemicalMixture<CoeffType>& _chem_mixture;
    
    const CEAThermoMixture<CoeffType>& _thermo;

    StirredReactorTimeIntegratorBase<CoeffType,StateType>& _time_integrator;

    CoeffType _volume;

    KineticsEvaluator<CoeffType,StateType> _kinetics_evaluator;

    CEAEvaluator<CoeffType> _thermo_evaluator;

  private:

    StirredReactorBase();

  };

  /* ---------------------- Constructor/Destructor ----------------------*/
  template<typename CoeffType, typename StateType>
  inline
  StirredReactorBase<CoeffType,StateType>::StirredReactorBase
  ( const ReactionSet<CoeffType>& reaction_set,
    const CEAThermoMixture<CoeffType>& thermo,
    StirredReactorTimeIntegratorBase<CoeffType,StateType>& time_integrator,
    const StateType example,
    CoeffType volume )
    : _reactor_type(ReactorType::INVALID),
      _reaction_set( reaction_set ),
      _chem_mixture( reaction_set.chemical_mixture() ),
      _thermo(thermo),
      _time_integrator(time_integrator),
      _volume(volume),
      _kinetics_evaluator(_reaction_set,example),
      _thermo_evaluator(thermo)
  {
    return;
  }

  template<typename CoeffType, typename StateType>
  inline
  StirredReactorBase<CoeffType,StateType>::~StirredReactorBase()
  {
    return;
  }

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType, typename StateType>
  template<typename VectorStateType>
  inline
  void StirredReactorBase<CoeffType,StateType>::run( const VectorStateType& x0,
                                                     const CoeffType t0,
                                                     const CoeffType t1,
                                                     const CoeffType dt )
  {
    _time_integrator.integrate( x0, t0, t1, dt, (*this) );
    return;
  }

  template<typename CoeffType, typename StateType>
  template<typename VectorStateType>
  inline
  void StirredReactorBase<CoeffType,StateType>::operator()( const VectorStateType& x,
                                                            VectorStateType& dx_dt )
  {
    switch( _reactor_type )
      {
      case( ReactorType::ISOTHERMAL ):
        {
          (static_cast<IsothermalStirredReactor<CoeffType,StateType>* >(this))( x, dx_dt );
        }
        break;

      case( ReactorType::INVALID ):
        {
          std::cerr << "Error: Must use a valid, derived stirred reactor class." << std::endl;
          antioch_error();
        }
        break;

      // Wat?!
      default:
        {
          antioch_error();
        }

      } // switch( _reactor_type )
    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_STIRRED_REACTOR_BASE_H
