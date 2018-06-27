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
#ifndef MATERIALLIB_MPL_MPCOMPONENT_H_
#define MATERIALLIB_MPL_MPCOMPONENT_H_

#include "mpProperty.h"

namespace MaterialPropertyLib
{
class Component
{
protected:
    PropertyArray _properties;
public:
    Component();
};

} //MaterialPropertyLib




#endif /* MATERIALLIB_MPL_MPCOMPONENT_H_ */
