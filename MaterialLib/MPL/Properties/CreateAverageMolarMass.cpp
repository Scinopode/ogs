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

#include "BaseLib/ConfigTree.h"
#include "AverageMolarMass.h"

namespace MaterialPropertyLib
{
std::unique_ptr<AverageMolarMass> createAverageMolarMass(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "AverageMolarMass");
    DBUG("Create AverageMolarMass medium property");
    return std::make_unique<AverageMolarMass>();
}
}  // namespace MaterialPropertyLib
