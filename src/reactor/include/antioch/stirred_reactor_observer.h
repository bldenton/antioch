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
    std::cout << "Destructor called" << std::endl;
    return;
  }

  template<typename CoeffType,  typename VectorStateType>
  inline
  void StirredReactorObserver<CoeffType,VectorStateType>::operator()( const VectorStateType& x, CoeffType time )
  {
    
    _time_hist.push_back(time);
    
    _x_hist.push_back( x );

    std::cout << time << " " << x << std::endl;

    return;
  }

  template<typename CoeffType, typename VectorStateType>
  inline
  void StirredReactorObserver<CoeffType,VectorStateType>::output_ascii( std::ostream& output ) const
  {
    // Header
    output << "#         t               x[s]" << std::endl;
    
    output << std::scientific << std::setprecision(10);

    std::cout << "time hist size = " << _time_hist.size() << std::endl;

    for( unsigned int t = 0; t < _time_hist.size(); t++ )
      {
        output << _time_hist[t] << " " << _x_hist[t] << std::endl;
      }
    
    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_STIRRED_REACTOR_OBSERVER_H
