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

#include "MaterialLib/MPL/Properties/DaltonsLaw.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{
PropertyDataType DaltonsLaw::value(VariableArray const& variable_array,
                                   ParameterLib::SpatialPosition const& pos,
                                   double const t, double const dt) const
{
    const double pGR = std::get<double>(
        variable_array[static_cast<int>(Variable::phase_pressure)]);

    auto const numberOfComponents = _phase->numberOfComponents();
    if (numberOfComponents < 2)
    {
        return Eigen::Vector2d{1., 0.};
    }
    else if (numberOfComponents > 2)
    {
        OGS_FATAL(
            "Property \"DaltonsLaw\" for more than two phase components is not "
            "yet implemented.");
    }

    double p_vap = 0.;
    size_t evaporatingComponentIndex = -1;
    for (size_t c = 0; c < numberOfComponents; c++)
    {
        if (_phase->component(c).hasProperty(
                PropertyType::saturation_vapour_pressure))
        {
            p_vap = _phase->component(c)
                        .property(PropertyType::saturation_vapour_pressure)
                        .template value<double>(variable_array, pos, t, dt);
            evaporatingComponentIndex = c;
            break;
        }
    }

    if (evaporatingComponentIndex < 0.)
    {
        OGS_FATAL(
            "Something went wrong! Saturation vapour pressure is not given as "
            "one of the component properties or its computation didn't work.");
    }

    double xn[2];

    xn[evaporatingComponentIndex] = p_vap / pGR;
    xn[1 - evaporatingComponentIndex] = 1 - xn[evaporatingComponentIndex];

    return Eigen::Vector2d{xn[0], xn[1]};
}

PropertyDataType DaltonsLaw::dValue(VariableArray const& variable_array,
                                    Variable const primary_variable,
                                    ParameterLib::SpatialPosition const& pos,
                                    double const t, double const dt) const
{
    (void)primary_variable;
    assert(((primary_variable == Variable::phase_pressure) ||
            (primary_variable == Variable::temperature)) &&
           "DaltonsLaw::dValue is implemented for "
           " derivatives with respect to phase_pressure only.");

    const double pGR = std::get<double>(
        variable_array[static_cast<int>(Variable::phase_pressure)]);

    auto const numberOfComponents = _phase->numberOfComponents();
    if (numberOfComponents < 2)
    {
        return Eigen::Vector2d{0., 0.};
    }
    else if (numberOfComponents > 2)
    {
        OGS_FATAL(
            "Property \"DaltonsLaw\" for more than two phase components is not "
            "yet implemented.");
    }

    double p_vap = 0.;

    size_t evaporatingComponentIndex = -1;
    for (size_t c = 0; c < numberOfComponents; c++)
    {
        if (_phase->component(c).hasProperty(
                PropertyType::saturation_vapour_pressure))
        {
            p_vap = _phase->component(c)
                        .property(PropertyType::saturation_vapour_pressure)
                        .template value<double>(variable_array, pos, t, dt);
            evaporatingComponentIndex = c;
            break;
        }
    }

    if (evaporatingComponentIndex < 0.)
    {
        OGS_FATAL(
            "Something went wrong! Saturation vapour pressure is not given as "
            "one of the component properties or its computation didn't work.");
    }

    double dxn[2];

    if (primary_variable == Variable::phase_pressure)
    {
        dxn[evaporatingComponentIndex] = -p_vap / pGR / pGR;
        dxn[1 - evaporatingComponentIndex] = -dxn[evaporatingComponentIndex];
        return Eigen::Vector2d{dxn[0], dxn[1]};
    }
    else if (primary_variable == Variable::temperature)
    {
        auto const dpVap_dT =
            _phase->component(evaporatingComponentIndex)
                .property(PropertyType::saturation_vapour_pressure)
                .template dValue<double>(variable_array, Variable::temperature,
                                         pos, t, dt);

        dxn[evaporatingComponentIndex] = dpVap_dT / pGR;
        dxn[1 - evaporatingComponentIndex] = -dxn[evaporatingComponentIndex];

        return Eigen::Vector2d{dxn[0], dxn[1]};
    }
    else
    {
        OGS_FATAL(
            "Something went wrong in PropertyDataType DaltonsLaw::dValue!");
    }
}

PropertyDataType DaltonsLaw::d2Value(VariableArray const& variable_array,
                                     Variable const primary_variable1,
                                     Variable const primary_variable2,
                                     ParameterLib::SpatialPosition const& pos,
                                     double const t, double const dt) const
{
    OGS_FATAL("DaltonsLaw::d2Value is not yet implemented.");

    return 0.;
}

}  // namespace MaterialPropertyLib
