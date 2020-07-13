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
#include "HenrysLaw.h"

namespace MaterialPropertyLib
{
std::unique_ptr<HenrysLaw> createHenrysLaw(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "HenrysLaw");
    DBUG("Create HenrysLaw medium property");
    return std::make_unique<HenrysLaw>();
}
}  // namespace MaterialPropertyLib
