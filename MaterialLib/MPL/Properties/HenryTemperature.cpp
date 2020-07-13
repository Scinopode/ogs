/**
 * \file
 * \author Norbert Grunwald
 * \date   Jul 08 2020
 * \brief
 *
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MaterialLib/MPL/Properties/HenryTemperature.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{
HenryTemperature::HenryTemperature(const double henryConstant,
                     const double dissolutionEnthalpy)
    : _henryConstant(henryConstant),
      _dissolutionEnthalpy(dissolutionEnthalpy){};

PropertyDataType HenryTemperature::value(VariableArray const& variable_array,
                                  ParameterLib::SpatialPosition const& pos,
                                  double const t, double const dt) const
{
    const double T = std::get<double>(
        variable_array[static_cast<int>(Variable::temperature)]);

    const double T_theta = 298.15;

    const double dH_by_R =
        _dissolutionEnthalpy / MaterialLib::PhysicalConstant::IdealGasConstant;

    // CO2-Henry coefficient depending on temperature
    const double H =
        _henryConstant * std::exp((-1.) * dH_by_R * (1. / T - 1. / T_theta));

    return H;
}

PropertyDataType HenryTemperature::dValue(VariableArray const& variable_array,
                                   Variable const primary_variable,
                                   ParameterLib::SpatialPosition const& pos,
                                   double const t, double const dt) const
{
    (void)primary_variable;
    assert((primary_variable == Variable::temperature) &&
           "HenryTemperature::dValue is implemented for "
           " derivatives with respect to temperature only.");

    if (primary_variable == Variable::temperature)
    {
        const double T = std::get<double>(
            variable_array[static_cast<int>(Variable::temperature)]);
        const double H = _component->property(PropertyType::henry_coefficient)
                         .template value<double>(variable_array, pos, t, dt);

        const double dH_by_R = _dissolutionEnthalpy /
                               MaterialLib::PhysicalConstant::IdealGasConstant;

        const double dH_dT = 1. / (T * T) * dH_by_R * H;
        return dH_dT;
    }
    else
    {
        return 0.;
    }
}

PropertyDataType HenryTemperature::d2Value(VariableArray const& variable_array,
                                    Variable const primary_variable1,
                                    Variable const primary_variable2,
                                    ParameterLib::SpatialPosition const& pos,
                                    double const t, double const dt) const
{
    OGS_FATAL("HenryTemperature::d2Value is not yet implemented.");

    return 0.;
}

}  // namespace MaterialPropertyLib
