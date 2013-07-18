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

#ifndef ANTIOCH_STIRRED_REACTOR_ENUM_H
#define ANTIOCH_STIRRED_REACTOR_ENUM_H

namespace Antioch
{
  namespace ReactorType
  {
    enum ReactorType{ INVALID = 0,
                      ISOTHERMAL };
  }

  namespace TimeIntegratorType
  {
    enum TimeIntegratorType{ INVALID = 0,
                             BOOST_ODE_INTEGRATOR };
  }

  namespace BoostStepperType
  {
    enum BoostStepperType{ RUNGE_KUTTA_FOURTH };
  }

} // end namespace Antioch

#endif // ANTIOCH_STIRRED_REACTOR_ENUM_H
