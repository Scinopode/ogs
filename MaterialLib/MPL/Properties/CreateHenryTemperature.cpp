/**
 * \file
 * \author Norbert Grunwald
 * \date   Jul 08 2020
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "HenryTemperature.h"

namespace MaterialPropertyLib
{
std::unique_ptr<HenryTemperature> createHenryTemperature(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "HenryTemperature");
    DBUG("Create HenryTemperature medium property");

    auto const henryConstant =
        //! \ogs_file_param{properties__property__SaturationBrooksCorey__residual_liquid_saturation}
        config.getConfigParameter<double>("henry_constant");
    auto const dissolutionEnthalpy =
        //! \ogs_file_param{properties__property__SaturationBrooksCorey__residual_liquid_saturation}
        config.getConfigParameter<double>("dissolution_enthalpy");

    return std::make_unique<HenryTemperature>(henryConstant, dissolutionEnthalpy);
}
}  // namespace MaterialPropertyLib
