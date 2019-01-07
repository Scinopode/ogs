/**
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "MaterialLib/MPL/mpProperty.h"

namespace MaterialPropertyLib
{
/**
 * The constant property class. This property simply retrieves the stored
 * constant value. It accepts all datatypes defined in PropertyDataType
 * (currently: double, Vector, Tensor, std::string)
 */
class Constant final : public Property
{
public:
    explicit Constant(PropertyDataType const&);
};

}  // namespace MaterialPropertyLib
