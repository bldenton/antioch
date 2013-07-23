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

// Antioch
#include "stirred_reactor_data.h"

namespace Antioch
{

  template< typename CoeffType=double, typename VectorStateType=std::vector<CoeffType> >
  class StirredReactorObserver
  {
  public:

    StirredReactorObserver( StirredReactorData<CoeffType,VectorStateType>& data );

    virtual ~StirredReactorObserver();

    void operator()( const VectorStateType& x, CoeffType time );

  protected:

    StirredReactorData<CoeffType,VectorStateType>& _data;

  private:

    StirredReactorObserver();

  };

  /* ---------------------- Constructor/Destructor ----------------------*/
  template<typename CoeffType,  typename VectorStateType>
  inline
  StirredReactorObserver<CoeffType,VectorStateType>::StirredReactorObserver( StirredReactorData<CoeffType,VectorStateType>& data )
    : _data(data)
  {
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
    
    _data.push_back_time(time);
    
    _data.push_back_species( x );

    return;
  }

  

} // end namespace Antioch

#endif // ANTIOCH_STIRRED_REACTOR_OBSERVER_H
