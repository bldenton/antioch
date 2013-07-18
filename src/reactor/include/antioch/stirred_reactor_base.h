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

namespace Antioch
{
  template<typename CoeffType=double, typename StateType=CoeffType>
  class StirredReactorBase
  {

  public:

    //! Constructor.
    StirredReactorBase( const ReactionSet<CoeffType>& reaction_set,
                        const CEAThermoMixture<CoeffType>& thermo,
                        StirredReactorTimeIntegratorBase<CoeffType,StateType>& time_integrator );

    virtual ~StirredReactorBase();

    template<typename VectorStateType>
    void run( const VectorStateType& x0,
              const CoeffType t0,
              const CoeffType t1,
              const CoeffType dt );
              
    void output( std::ostream& output ) const;

    template<typename VectorStateType>
    void operator()( const VectorStateType& x,
                     VectorStateType& dx_dt );

  protected:

    const ReactionSet<CoeffType>& _reaction_set;

    const ChemicalMixture<CoeffType>& _chem_mixture;
    
    const CEAThermoMixture<CoeffType>& _thermo;

    StirredReactorTimeIntegratorBase<CoeffType,StateType>& _time_integrator;

    KineticsEvaluator<CoeffType,StateType> _kinetics_evaluator;

    CEAEvaluator<CoeffType> _cea_evaluator;

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
    const StateType example )
    : _reaction_set( reaction_set ),
      _chem_mixture( reaction_set.chemical_mixture() ),
      _thermo(thermo),
      _time_integrator(time_integrator),
      _kinetics_evaluator(_reaction_set,example),
      _cea_evaluator(thermo)
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
    _time_integrator.integrate( x0, t0, t1, dt, (*this) )
    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_STIRRED_REACTOR_BASE_H
