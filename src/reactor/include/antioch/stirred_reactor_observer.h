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

#ifndef ANTIOCH_STIRRED_REACTOR_OBSERVER_H
#define ANTIOCH_STIRRED_REACTOR_OBSERVER_H

namespace Antioch
{

  template< typename CoeffType=double, typename VectorStateType=std::vector<CoeffType> >
  class StirredReactorObserver
  {
  public:

    StirredReactorObserver( const unsigned int n_est_steps = 100 );

    virtual ~StirredReactorObserver();

    void operator()( const VectorStateType& x, CoeffType time );

    void output_ascii( std::ostream& output ) const;

  protected:

    std::vector<CoeffType> _time_hist;

    std::vector<VectorStateType> _x_hist;

  };

  /* ---------------------- Constructor/Destructor ----------------------*/
  template<typename CoeffType,  typename VectorStateType>
  inline
  StirredReactorObserver<CoeffType,VectorStateType>::StirredReactorObserver( const unsigned int n_est_steps )
  {
    _time_hist.reserve(n_est_steps);
    _x_hist.reserve(n_est_steps);

    return;
  }

  template<typename CoeffType,  typename VectorStateType>
  inline
  StirredReactorObserver<CoeffType,VectorStateType>::~StirredReactorObserver()
  {
    return;
  }

  template<typename CoeffType,  typename VectorStateType>
  inline
  void StirredReactorObserver<CoeffType,VectorStateType>::operator()( const VectorStateType& x, CoeffType time )
  {
    _time_hist.push_back(time);
    
    _x_hist.push_back( x );

    return;
  }

  template<typename CoeffType, typename VectorStateType>
  inline
  void StirredReactorObserver<CoeffType,VectorStateType>::output_ascii( std::ostream& output ) const
  {
    unsigned int n_species = _kinetics_evalutor.reaction_set().n_species();

    // Header
    output << "#         t               x[s]" << std::endl;
    
    output << std::scientific << std::setprecision(10);
    unsigned int n_steps = this->n_steps();

    for( unsigned t = 0; t < n_steps; t++ )
      {
        output << this->time(t) << " ";
        for( unsigned s = 0; s < n_species; s++ )
          {
            output << this->X(t,s) << " ";
          }
        
        output << std::endl;
      }
    
    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_STIRRED_REACTOR_OBSERVER_H
