/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BaseLib/ConfigTree.h"
#include "RelPermVanGenuchten.h"

namespace MaterialPropertyLib
{
std::unique_ptr<RelPermVanGenuchten> createRelPermVanGenuchten(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "RelativePermeabilityVanGenuchten");
    DBUG("Create RelativePermeabilityVanGenuchten medium property");

    auto const residual_liquid_saturation =
        //! \ogs_file_param{properties__property__RelativePermeabilityVanGenuchten__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{properties__property__RelativePermeabilityVanGenuchten__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");
    auto const min_relative_permeability_liquid =
        //! \ogs_file_param{properties__property__RelativePermeabilityVanGenuchten__minimum_relative_permeability_liquid}
        config.getConfigParameter<double>(
            "minimum_relative_permeability_liquid");
    auto const min_relative_permeability_gas =
        //! \ogs_file_param{properties__property__RelativePermeabilityVanGenuchten__minimum_relative_permeability_gas}
        config.getConfigParameter<double>("minimum_relative_permeability_gas");
    auto const exponent =
        //! \ogs_file_param{properties__property__RelativePermeabilityVanGenuchten__exponent}
        config.getConfigParameter<double>("exponent");
    if (exponent <= 0. || exponent >= 1.)
    {
        OGS_FATAL(
            "RelativePermeabilityVanGenuchten::exponent must be in the (0, 1) "
            "range.");
    }
    auto const gamma =
        //! \ogs_file_param{properties__property__RelativePermeabilityVanGenuchten__gamma}
        config.getConfigParameter<double>("gamma");
    if (gamma < 0.)
    {
        OGS_FATAL(
            "RelativePermeabilityVanGenuchten::gamma must be positive or "
            "zero.");
    }
    auto const epsilon =
        //! \ogs_file_param{properties__property__RelativePermeabilityVanGenuchten__epsilon}
        config.getConfigParameter<double>("epsilon");
    if (gamma < 0.)
    {
        OGS_FATAL(
            "RelativePermeabilityVanGenuchten::epsilon must be positive or "
            "zero.");
    }

    return std::make_unique<RelPermVanGenuchten>(
        residual_liquid_saturation,
        residual_gas_saturation,
        min_relative_permeability_liquid,
        min_relative_permeability_gas,
        exponent,
        gamma,
        epsilon);
}
}  // namespace MaterialPropertyLib
