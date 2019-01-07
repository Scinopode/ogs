/**
 * \author Norbert Grunwald
 * \date   18.09.2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "MaterialLib/MPL/mpProperty.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;
/**
 * \class IdealGasLaw
 */
class IdealGasLaw final : public Property
{
private:
    /// A pointer to the phase object.
    Phase* _phase;
    Component* _component;

public:
    /// Constructor passing a pointer to the medium.
    IdealGasLaw(Medium* /*unused*/);
    /// Constructor passing a pointer to the phase.
    IdealGasLaw(Phase* /*p*/);
    /// Constructor passing a pointer to the component.
    IdealGasLaw(Component* /*c*/);
    /// This method overrides the base class implementation and
    /// actually computes and sets the property _value.
    PropertyDataType value(VariableArray const&) override;
    PropertyDataType dvalue(VariableArray const&, Variables const ) override;

};

}  // namespace MaterialPropertyLib
