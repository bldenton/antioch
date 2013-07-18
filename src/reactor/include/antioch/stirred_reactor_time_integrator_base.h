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

#ifndef ANTIOCH_STIRRED_REACTOR_TIME_INTEGRATOR_BASE_H
#define ANTIOCH_STIRRED_REACTOR_TIME_INTEGRATOR_BASE_H

// Antioch
#include "antioch/stirred_reactor_enum.h"

namespace Antioch
{
  // Foward declarations
  template<typename CoeffType, typename StateType>
  class BoostODEIntergrator;

  template<typename CoeffType=double, typename StateType=CoeffType>
  class StirredReactorTimeIntegratorBase
  {
  public:

    StirredReactorTimeIntegratorBase();

    virtual ~StirredReactorTimeIntegratorBase();

    //! Principal method for integrating ODE system.
    /*! Implemented in derived classes. We can't do virtual
        templated methods, so we settle for a runtime error
        if using the base class and not a derived class. 
        Returns number of time steps taken. */
    template<typename VectorStateType>
    unsigned int integrate( const VectorStateType& x0,
                            const CoeffType t0,
                            const CoeffType t1,
                            const CoeffType dt,
                            StirredReactorBase<CoeffType,StateType>& reactor );

  protected:

    TimeIntegratorType _integrator_type;

  };

  /* ---------------------- Constructor/Destructor ----------------------*/
  template<typename CoeffType, typename StateType>
  inline
  StirredReactorTimeIntegratorBase<CoeffType,StateType>::StirredReactorTimeIntegratorBase()
    : _integrator_type(INVALID)
  {
    return;
  }

  template<typename CoeffType, typename StateType>
  inline
  StirredReactorTimeIntegratorBase<CoeffType,StateType>::~StirredReactorTimeIntegratorBase()
  {
    return;
  }

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType, typename StateType>
  template<typename VectorStateType>
  inline
  unsigned int StirredReactorTimeIntegratorBase<CoeffType,StateType>::integrate
  ( const VectorStateType& x0,
    const CoeffType t0,
    const CoeffType t1,
    const CoeffType dt,
    StirredReactorBase<CoeffType,StateType>& reactor )
  {
    unsigned int n_steps = 0;

    switch( _integrator_type )
      {
      case( TimeIntegratorType::BOOST_ODE_INTEGRATOR ):
        {
          n_steps = (static_cast<BoostODEIntergrator*>(this))->integrate(x0, t0, t1, dt, reactor);
        }
        break;

      case( TimeIntegratorType::INVALID ):
        {
          std::cerr << "Error: Must use a valid, derived time integrator class." << std::endl;
          antioch_error();
        }
        break;

      // Wat?!
      default:
        {
          antioch_error();
        }

      } // switch( _integrator_type )

    return n_steps;
  }

} // end namespace Antioch

#endif // ANTIOCH_STIRRED_REACTOR_TIME_INTEGRATOR_BASE_H
