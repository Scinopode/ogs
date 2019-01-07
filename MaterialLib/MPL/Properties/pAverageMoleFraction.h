/**
 * \author Norbert Grunwald
 * \date   12.09.2017
 * \brief  
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#ifndef MATERIALLIB_MPL_PROPERTIES_PAVERAGEMOLEFRACTION_H_
#define MATERIALLIB_MPL_PROPERTIES_PAVERAGEMOLEFRACTION_H_

#include "../mpProperty.h"

namespace MaterialPropertyLib
{

class Medium;
class Phase;
class Component;
/**
 * \class AverageMoleFraction
 * \brief A function averaging a property by mole fraction
 * \details This property is usually a phase property, it
 * computes the average of individual component properties
 * weighted by mole fraction.
 */
class AverageMoleFraction : public Property
{
private:
	/// A pointer to the phase object.
	Phase* _phase;
public:
	/// Probably useless empty constructor.
	AverageMoleFraction(){};
	/// Constructor passing a pointer to the medium.
	AverageMoleFraction(Medium*);
	/// Constructor passing a pointer to the phase.
	AverageMoleFraction(Phase*);
	/// Constructor passing a pointer to the component.
	AverageMoleFraction(Component*);
	/// This method overrides the base class implementation and
	/// actually computes and sets the property _value.
	PropertyDataType value (VariableArray const&) override;
};

}  // MaterialPropertyLib





#endif /* MATERIALLIB_MPL_PROPERTIES_PAVERAGEMOLEFRACTION_H_ */
