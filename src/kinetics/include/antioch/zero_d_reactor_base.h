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

#ifndef ANTIOCH_ZERO_D_REACTOR_BASE_H
#define ANTIOCH_ZERO_D_REACTOR_BASE_H

namespace Antioch
{
  template<typename CoeffType=double, typename StateType=CoeffType>
  class ZeroDimensonalReactorBase
  {

  public:

    //! Constructor.
    ZeroDimensonalReactorBase( const ReactionSet<CoeffType>& reaction_set,
                               const StateType& example );

    virtual ~ZeroDimensonalReactorBase();

    template<typename VectorStateType>
    void run( const StateType& T,
              const VectorStateType& molar_densities,
              const CoeffType t0,
              const CoeffType t1,
              const CoeffType dt,
              const unsigned int estimated_n_steps = 100 );
                   
    void output( std::ostream& output ) const;

  protected:

    const ReactionSet<CoeffType>& _reaction_set;

    const ChemicalMixture<CoeffType>& _chem_mixture;

  };

  /* ---------------------- Constructor/Destructor ----------------------*/
  template<typename CoeffType, typename StateType>
  inline
  ZeroDimensonalReactorBase<CoeffType,StateType>::ZeroDimensonalReactorBase
  ( const ReactionSet<CoeffType>& reaction_set,
    const StateType& example )
    : _reaction_set( reaction_set ),
      _chem_mixture( reaction_set.chemical_mixture() )
  {
    return;
  }

  template<typename CoeffType, typename StateType>
  inline
  ZeroDimensonalReactorBase<CoeffType,StateType>::~ZeroDimensonalReactorBase()
  {
    return;
  }

  /* ------------------------- Inline Functions -------------------------*/

} // end namespace Antioch

#endif // ANTIOCH_ZERO_D_REACTOR_BASE_H
