/**
 * \file
 * \author Norbert Grunwald
 * \date   Jul 07 2020
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

#include "MaterialLib/MPL/Properties/Antoine.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{
Antoine::Antoine(const double A, const double B, const double C)
    : _A(A), _B(B), _C(C){};

PropertyDataType Antoine::value(VariableArray const& variable_array,
                                ParameterLib::SpatialPosition const& pos,
                                double const t, double const dt) const
{
    const double T = std::get<double>(
        variable_array[static_cast<int>(Variable::temperature)]);

    const double p_vap = std::pow(10., (_A - _B / (_C + T - 273.15))) * 133.322;
    return p_vap;
}

PropertyDataType Antoine::dValue(VariableArray const& variable_array,
                                 Variable const primary_variable,
                                 ParameterLib::SpatialPosition const& pos,
                                 double const t, double const dt) const
{
    (void)primary_variable;
    assert((primary_variable == Variable::temperature) &&
           "Antoine::dValue is implemented for derivatives with respect to "
           "temperature only.");

    const double T = std::get<double>(
        variable_array[static_cast<int>(Variable::temperature)]);

    auto const p_vap = std::pow(10., (_A - _B / (_C + T - 273.15))) * 133.322;
    auto const d_p_vap_dT =
        _B / (_C + T - 273.15) / (_C + T - 273.15) * std::log(10.) * p_vap;

    return d_p_vap_dT;
}

PropertyDataType Antoine::d2Value(VariableArray const& variable_array,
                                  Variable const primary_variable1,
                                  Variable const primary_variable2,
                                  ParameterLib::SpatialPosition const& pos,
                                  double const t, double const dt) const
{
    OGS_FATAL("Antoine::d2Value is not yet implemented.");

    return 0.;
}

}  // namespace MaterialPropertyLib
