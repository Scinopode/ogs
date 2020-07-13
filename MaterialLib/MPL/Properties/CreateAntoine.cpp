/**
 * \file
 * \author Norbert Grunwald
 * \date   Jul 07 2020
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Antoine.h"
#include "BaseLib/ConfigTree.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Antoine> createAntoine(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "Antoine");

    DBUG("Create Antoine component property");

    auto const A =
        //! \ogs_file_param{properties__property__SaturationBrooksCorey__residual_liquid_saturation}
        config.getConfigParameter<double>("parameter_a");
    auto const B =
        //! \ogs_file_param{properties__property__SaturationBrooksCorey__residual_liquid_saturation}
        config.getConfigParameter<double>("parameter_b");
    auto const C =
        //! \ogs_file_param{properties__property__SaturationBrooksCorey__residual_liquid_saturation}
        config.getConfigParameter<double>("parameter_c");

    return std::make_unique<Antoine>(A, B, C);
}
}  // namespace MaterialPropertyLib
