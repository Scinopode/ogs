/**
 * \author Norbert Grunwald
 * \date   07.09.2017
 * \brief  
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_MPL_MPMEDIUM_H_
#define MATERIALLIB_MPL_MPMEDIUM_H_

#include "BaseLib/ConfigTree.h"
#include "mpProperty.h"


#include <vector>

namespace MaterialPropertyLib
{

class Phase;
/**
 * \class Medium
 * \brief This class is for material objects on the Medium scale.
 * \details A Medium consists of an arbitrarily long vector of phases
 * and an array of properties.
 */
class Medium
{
private:
	/// The vector that holds the phases.
    std::vector<Phase*> _phases;
    /// The array that holds the maedium properties.
    PropertyArray _properties;
public:
    /// The medium is constructed based on the config tree object.
    Medium(BaseLib::ConfigTree const&);
    /// A method that parses the phase details and stores them in the
    /// private _phases member.
    void createPhases (BaseLib::ConfigTree const&);
    /// A method that parses the medium property details and stores
    /// them in the private _properties member.
    void createProperties (BaseLib::ConfigTree const&);

    /// A method that creates the default properties of the medium.
    /// Currently, these defaults is the volume fraction average.
    void createDefaultProperties(void);

    /// A get-function for a particular phase. The ul argument specifies
    /// the index within the _phases vector.
    Phase* phase (std::size_t const);
    /// A get-function for a property. The argument refers to the
    /// name of the property.
    Property* property(PropertyEnum const&);

    /// A simple get-function for retrieving the number of phases the
    /// medium consists of.
    std::size_t numberOfPhases(void);

    /// A method that prints out a summary of the constructed medium object
    /// with all its phases, components, and properties. Basically for
    /// debugging purposes only.
    void summary (void);
};

}

#endif /* MATERIALLIB_MPL_MPMEDIUM_H_ */
