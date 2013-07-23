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

#ifndef ANTIOCH_STIRRED_REACTOR_DATA_H
#define ANTIOCH_STIRRED_REACTOR_DATA_H

// Antioch
#include "antioch/reaction_set.h"

namespace Antioch
{

  template< typename CoeffType=double, typename VectorStateType=std::vector<CoeffType> >
  class StirredReactorData
  {
  public:

    StirredReactorData( const ReactionSet<CoeffType>& reaction_set,
                        const unsigned int n_est_steps = 100 );

    virtual ~StirredReactorData();

    void push_back_time( const CoeffType time );

    void push_back_species( const VectorStateType& species );

    void output_ascii( std::ostream& output ) const;

  protected:

    const ReactionSet<CoeffType>& _reaction_set;

    std::vector<CoeffType> _time_hist;

    std::vector<VectorStateType> _x_hist;

  };

  template<typename CoeffType, typename VectorStateType>
  inline
  StirredReactorData<CoeffType,VectorStateType>::StirredReactorData( const ReactionSet<CoeffType>& reaction_set,
                                                                     const unsigned int n_est_steps )
    : _reaction_set(reaction_set)
  {
    _time_hist.reserve( n_est_steps );
    _x_hist.reserve( n_est_steps );
    return;
  }

  template<typename CoeffType, typename VectorStateType>
  inline
  StirredReactorData<CoeffType,VectorStateType>::~StirredReactorData()
  {
    return;
  }

  template<typename CoeffType, typename VectorStateType>
  inline
  void StirredReactorData<CoeffType,VectorStateType>::push_back_time( const CoeffType time )
  {
    _time_hist.push_back(time);
    return;
  }

  template<typename CoeffType, typename VectorStateType>
  inline
  void StirredReactorData<CoeffType,VectorStateType>::push_back_species( const VectorStateType& species )
  {
    _x_hist.push_back(species);
    return;
  }

  template<typename CoeffType, typename VectorStateType>
  inline
  void StirredReactorData<CoeffType,VectorStateType>::output_ascii( std::ostream& output ) const
  {
    unsigned int n_species = _reaction_set.n_species();
    
    // Header
    output << "# species names" << std::endl;
    for( unsigned int s = 0; s < n_species; s++ )
      {
        output << _reaction_set.chemical_mixture().species_name(s) << " ";
      }
    output << std::endl;

    output << "#         t               x[s]" << std::endl;
    
    // Time history of data
    output << std::scientific << std::setprecision(10);

    for( unsigned int t = 0; t < _time_hist.size(); t++ )
      {
        output << _time_hist[t] << " ";

        for( unsigned int s = 0; s < n_species; s++ )
          {
            output << _x_hist[t][s] << " ";
          }

        output << std::endl;
      }
    
    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_STIRRED_REACTOR_DATA_H
