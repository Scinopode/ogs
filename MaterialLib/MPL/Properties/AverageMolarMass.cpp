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

#include "MaterialLib/MPL/Properties/AverageMolarMass.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{
PropertyDataType AverageMolarMass::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    auto const numberOfComponents = _phase->numberOfComponents();
    if (numberOfComponents < 1)
    {
        return _phase->property(PropertyType::molar_mass)
            .template value<double>(variable_array, pos, t, dt);
    }

    // TODO: This should return a vector of length _phase->numberOfComponents().
    // Currently, this feature is implemented for binary phases only.
    auto const molar_fraction =
        _phase->property(PropertyType::mole_fraction)
            .template value<Eigen::Vector2d>(variable_array, pos, t, dt);

    double M = 0.;
    for (size_t c = 0; c < numberOfComponents; c++)
    {
        auto const M_zeta =
            _phase->component(c)
                .property(PropertyType::molar_mass)
                .template value<double>(variable_array, pos, t, dt);
        auto const xn_zeta = molar_fraction[c];

        M += xn_zeta * M_zeta;
    }

    // const double molar_mass = _medium->property(PropertyType::saturation)
    //                      .template value<double>(variable_array, pos, t, dt);

    return M;
}

PropertyDataType AverageMolarMass::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    (void)primary_variable;
    assert(((primary_variable == Variable::phase_pressure) ||
            (primary_variable == Variable::temperature)) &&
           "AverageMolarMass::dValue is implemented for "
           " derivatives with respect to phase_pressure or temperature only.");

    auto const numberOfComponents = _phase->numberOfComponents();
    if (numberOfComponents < 1)
    {
        return 0.;
    }
    else if (numberOfComponents > 2)
    {
        OGS_FATAL(
            "AverageMolarMass::dvalue is currently implemented two or less "
            "phase components only.");
    }
    // TODO: This should return a vector of length
    // _phase->numberOfComponents(). Currently, this feature is implemented
    // for binary phases only.
    auto const dxnC = _phase->property(PropertyType::mole_fraction)
                          .template dValue<Eigen::Vector2d>(
                              variable_array, primary_variable, pos, t, dt)[0];

    auto const M_C = _phase->component(0)
                         .property(PropertyType::molar_mass)
                         .template value<double>(variable_array, pos, t, dt);
    auto const M_W = _phase->component(1)
                         .property(PropertyType::molar_mass)
                         .template value<double>(variable_array, pos, t, dt);

    return dxnC * (M_C - M_W);

}  // namespace MaterialPropertyLib

PropertyDataType AverageMolarMass::d2Value(
    VariableArray const& variable_array, Variable const primary_variable1,
    Variable const primary_variable2, ParameterLib::SpatialPosition const& pos,
    double const t, double const dt) const
{
    OGS_FATAL("AverageMolarMass::d2Value is not yet implemented.");

    return 0.;
}

}  // namespace MaterialPropertyLib
