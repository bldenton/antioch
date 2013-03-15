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
//
// $Id: eigen_utils.h 37170 2013-02-19 21:40:39Z roystgnr $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_EIGEN_UTILS_H
#define ANTIOCH_EIGEN_UTILS_H

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_EIGEN

#include <Eigen/Dense>

#include "antioch/metaprogramming.h"

namespace Antioch
{

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline
_Scalar
max (const Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& in)
{
  return in.maxCoeff();
}

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
struct value_type<Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >
{
  typedef Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> 
    container_type;
  typedef _Scalar type;

  static inline
  container_type
  constant(const type& in) { return container_type::Constant(in); }
};

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline
Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
zero_clone(const Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& example)
{
  return 
  Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
  (example.rows(), example.cols()).setZero();
}

} // end namespace Antioch

#endif // ANTIOCH_HAVE_EIGEN

#endif //ANTIOCH_EIGEN_UTILS_H