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

#include "MaterialLib/MPL/Properties/HenrysLaw.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{
PropertyDataType HenrysLaw::value(VariableArray const& variable_array,
                                  ParameterLib::SpatialPosition const& pos,
                                  double const t, double const dt) const
{
    const double pGR = std::get<double>(
        variable_array[static_cast<int>(Variable::phase_pressure)]);
    const double xnCG = std::get<double>(
        variable_array[static_cast<int>(Variable::mole_fraction)]);

    auto const numberOfComponents = _phase->numberOfComponents();
    if (numberOfComponents < 2)
    {
        return 0.;
    }
    else if (numberOfComponents > 2)
    {
        OGS_FATAL(
            "Property \"HenrysLaw\" for more than two phase components is not "
            "yet implemented.");
    }

    // find the component that is soluble and obtain its Henry-coefficient
    double H = 0.;
    for (size_t c = 0; c < numberOfComponents; c++)
    {
        if (_phase->component(c).hasProperty(PropertyType::henry_coefficient))
        {
            H = _phase->component(c)
                    .property(PropertyType::henry_coefficient)
                    .template value<double>(variable_array, pos, t, dt);
            const double cCL = H * xnCG * pGR;

            return cCL;
        }
    }

    OGS_FATAL(
        "Something went wrong! Henry-Coefficient is not given as one of "
        "the component properties or its computation didn't work.");

    return 0.;
}

PropertyDataType HenrysLaw::dValue(VariableArray const& variable_array,
                                   Variable const primary_variable,
                                   ParameterLib::SpatialPosition const& pos,
                                   double const t, double const dt) const
{
    (void)primary_variable;
    assert((primary_variable == Variable::liquid_phase_pressure) &&
           "HenrysLaw::dValue is implemented for derivatives with respect to "
           "liquid_phase_pressure only.");
    // The temperature-derivative is not implemented, since it requires a
    // property of another phase (saturation_vapuor_pressure defined in the gas
    // phase)

    auto const numberOfComponents = _phase->numberOfComponents();
    if (numberOfComponents < 2)
    {
        return 0.;
    }
    else if (numberOfComponents > 2)
    {
        OGS_FATAL(
            "Property \"HenrysLaw\" for more than two phase components is not "
            "yet implemented.");
    }

    // find the component that is soluble and obtain its Henry-coefficient
    double H = 0.;

    for (size_t c = 0; c < numberOfComponents; c++)
    {
        if (_phase->component(c).hasProperty(PropertyType::henry_coefficient))
        {
            H = _phase->component(c)
                    .property(PropertyType::henry_coefficient)
                    .template value<double>(variable_array, pos, t, dt);
            return H;
        }
    }

    OGS_FATAL(
        "Something went wrong! Henry-Coefficient is not given as one of "
        "the component properties or its computation didn't work.");
}

PropertyDataType HenrysLaw::d2Value(VariableArray const& variable_array,
                                    Variable const primary_variable1,
                                    Variable const primary_variable2,
                                    ParameterLib::SpatialPosition const& pos,
                                    double const t, double const dt) const
{
    OGS_FATAL("HenrysLaw::d2Value is not yet implemented.");

    return 0.;
}

}  // namespace MaterialPropertyLib
