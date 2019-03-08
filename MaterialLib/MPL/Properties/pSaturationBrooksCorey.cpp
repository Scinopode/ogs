/**
 * \author Norbert Grunwald
 * \date   27.06.2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MaterialLib/MPL/Properties/pSaturationBrooksCorey.h"
#include "MaterialLib/MPL/mpMedium.h"
#include "pUniversalConstants.h"

#include <algorithm>
#include <cmath>
#include <limits>


namespace MaterialPropertyLib
{
SaturationBrooksCorey::SaturationBrooksCorey(
        Medium* m,
        const double residual_liquid_saturation,
        const double residual_gas_saturation,
        const double exponent,
        const double entry_pressure)

: _medium(m),
  _residual_liquid_saturation (residual_liquid_saturation),
  _residual_gas_saturation (residual_gas_saturation),
  _exponent (exponent),
  _entry_pressure (entry_pressure) {};

/**
 */
PropertyDataType  SaturationBrooksCorey::value(VariableArray const& v)
{

    const double p_cap = std::max(1.e-14,
            getScalar(v[MaterialPropertyLib::Variables::capillary_pressure]));

    const double s_L_res =_residual_liquid_saturation;
    const double s_L_max = 1.0 - _residual_gas_saturation;
    const double lambda = _exponent;
    const double p_b = _entry_pressure;

    if (p_cap <= p_b)
        return s_L_max;

    const double s_eff = std::pow(p_b/p_cap, lambda);
    return s_eff*(s_L_max - s_L_res) + s_L_res;
}

PropertyDataType SaturationBrooksCorey::dvalue(VariableArray const& v, Variables const pv)
{

    assert((pv == Variables::capillary_pressure) && "SaturationBrooksCorey::dvalue is implemented for "
            " derivatives with respect to capillary pressure only.");
    const double p_cap = std::max(1.e-14,
            getScalar(v[MaterialPropertyLib::Variables::capillary_pressure]));

    const double s_L_res = _residual_liquid_saturation;
    const double s_L_max = 1.0 - _residual_gas_saturation;

    const double lambda = _exponent;
    const double p_b = _entry_pressure;

    return -lambda / p_cap * std::pow( p_b / p_cap, lambda) *
            (s_L_max - s_L_res);
}

PropertyDataType SaturationBrooksCorey::ddvalue(VariableArray const& v,
        Variables const pv1, Variables const pv2)
{
    assert((pv1 == Variables::capillary_pressure) &&
            (pv2 == Variables::capillary_pressure) &&
            "SaturationBrooksCorey::ddvalue is implemented for "
            " derivatives with respect to capillary pressure only.");

    const double p_cap = getScalar(
            v[MaterialPropertyLib::Variables::capillary_pressure]);

    const double s_L_res = _residual_liquid_saturation;
    const double s_L_max = 1.0 - _residual_gas_saturation;

    const double lambda = _exponent;
    const double p_b = _entry_pressure;

    return lambda * (lambda + 1) * std::pow( p_b / p_cap, lambda) /
            ( p_cap * p_cap ) * (s_L_max - s_L_res);
}



}  // MaterialPropertyLib
