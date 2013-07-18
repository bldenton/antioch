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

#ifndef ANTIOCH_BOOST_ODE_INTEGRATOR_H
#define ANTIOCH_BOOST_ODE_INTEGRATOR_H

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_BOOST_ODEINT

// Boost
#include "boost/numeric/odeint.hpp"

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/stirred_reactor_enum.h"
#include "antioch/stirred_reactor_base.h"
#include "antioch/stirred_reactor_time_integrator_base.h"

namespace Antioch
{

  template<typename CoeffType=double, typename StateType=CoeffType>
  class BoostODEIntegrator : public StirredReactorTimeIntegratorBase<CoeffType,StateType>
  {

  public:

    BoostODEIntegrator( BoostStepperType::BoostStepperType stepper_type );

    virtual ~BoostODEIntegrator();

    template<typename VectorStateType>
    unsigned int integrate( VectorStateType& x0,
                            CoeffType t0,
                            CoeffType t1,
                            CoeffType dt,
                            StirredReactorBase<CoeffType,StateType>& reactor );

    template<typename VectorStateType>
    void operator()( const VectorStateType& x,
                     VectorStateType& dx_dt,
                     const CoeffType t );

  protected:

    BoostStepperType::BoostStepperType _stepper_type;

    StirredReactorBase<CoeffType,StateType>* _reactor;

  private:

    BoostODEIntegrator();

  };

  /* ---------------------- Constructor/Destructor ----------------------*/
  template<typename CoeffType, typename StateType>
  inline
  BoostODEIntegrator<CoeffType,StateType>::BoostODEIntegrator
  ( BoostStepperType::BoostStepperType stepper_type )
    : StirredReactorTimeIntegratorBase<CoeffType,StateType>(),
      _stepper_type(stepper_type),
      _reactor(NULL)
  {
    this->_integrator_type = TimeIntegratorType::BOOST_ODE_INTEGRATOR;
    return;
  }
    
  template<typename CoeffType, typename StateType>
  inline
  BoostODEIntegrator<CoeffType,StateType>::~BoostODEIntegrator()
  {
    return;
  }

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType, typename StateType>
  template<typename VectorStateType>
  inline
  unsigned int BoostODEIntegrator<CoeffType,StateType>::integrate
  ( VectorStateType& x0,
    CoeffType t0,
    CoeffType t1,
    CoeffType dt,
    StirredReactorBase<CoeffType,StateType>& reactor )
  {
    _reactor = &reactor;

    unsigned int n_steps = 0;

    CoeffType abs_err = 1.0e-10;
    CoeffType rel_err = 1.0e-6;
    CoeffType a_x = 1.0;
    CoeffType a_dxdt = 1.0;

    // Shoudl error_checker be templated on StateType?
    boost::numeric::odeint::controlled_runge_kutta< boost::numeric::odeint::runge_kutta_cash_karp54<VectorStateType> > stepper( boost::numeric::odeint::default_error_checker<CoeffType>( abs_err , rel_err , a_x , a_dxdt ) );

    n_steps = boost::numeric::odeint::integrate_adaptive( stepper , (*this) , x0 , t0 , t1 , dt );

    _reactor = NULL;

    return n_steps;
  }

  template<typename CoeffType, typename StateType>
  template<typename VectorStateType>
  inline
  void BoostODEIntegrator<CoeffType,StateType>::operator()( const VectorStateType& x,
                                                            VectorStateType& dx_dt,
                                                            const CoeffType /* t */ )
  {
    antioch_assert(_reactor);

    (*_reactor)( x, dx_dt );

    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_HAVE_BOOST_ODEINT

#endif // ANTIOCH_BOOST_ODE_INTEGRATOR_H
