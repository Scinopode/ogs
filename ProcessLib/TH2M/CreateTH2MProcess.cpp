/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateTH2MProcess.h"

#include <cassert>

#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"

#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "TH2MProcess.h"
#include "TH2MProcessData.h"

namespace ProcessLib
{
namespace TH2M
{
template <int DisplacementDim>
std::unique_ptr<Process> createTH2MProcess(
    std::string name, MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order, BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "TH2M");
    DBUG("Create TH2MProcess.");

    auto const staggered_scheme =
        //! \ogs_file_param{prj__processes__process__TH2M__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(staggered_scheme && (*staggered_scheme == "staggered"));

    // Process variable.

    //! \ogs_file_param{prj__processes__process__TH2M__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    ProcessVariable* variable_pGR;
    ProcessVariable* variable_pCap;
    ProcessVariable* variable_T;
    ProcessVariable* variable_u;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        auto per_process_variables = findProcessVariables(
            variables, pv_config,
            {//! \ogs_file_param_special{prj__processes__process__TH2M__process_variables__gas_pressure}
             "gas_pressure",
             //! \ogs_file_param_special{prj__processes__process__TH2M__process_variables__capillary_pressure}
             "capillary_pressure",
             //! \ogs_file_param_special{prj__processes__process__TH2M__process_variables__temperature}
             "temperature",
             //! \ogs_file_param_special{prj__processes__process__TH2M__process_variables__displacement}
             "displacement"});
        variable_pGR = &per_process_variables[0].get();
        variable_pCap = &per_process_variables[1].get();
        variable_T = &per_process_variables[2].get();
        variable_u = &per_process_variables[3].get();
        process_variables.push_back(std::move(per_process_variables));
    }
    else  // staggered scheme.
    {
        OGS_FATAL("A Staggered version of TH2M is not implemented yet.");
        // The following implementation is based on the
        // ThermoHydroMechanics-Process. It was expanded by a second pressure
        // (capillary pressure), but the implementation is not tested for TH2M.

        using namespace std::string_literals;
        for (auto const& variable_name :
             {"gas_pressure"s, "capillary_pressure"s, "temperature"s,
              "displacement"s})
        {
            auto per_process_variables =
                findProcessVariables(variables, pv_config, {variable_name});
            process_variables.push_back(std::move(per_process_variables));
        }
        variable_pGR = &process_variables[0][0].get();
        variable_pCap = &process_variables[1][0].get();
        variable_T = &process_variables[2][0].get();
        variable_u = &process_variables[3][0].get();
    }

    DBUG("Associate displacement with process variable '%s'.",
         variable_u->getName().c_str());

    if (variable_u->getNumberOfComponents() != DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '%s' is different "
            "from the displacement dimension: got %d, expected %d",
            variable_u->getName().c_str(),
            variable_u->getNumberOfComponents(),
            DisplacementDim);
    }

    DBUG("Associate gas (non-wetting) pressure with process variable '%s'.",
         variable_pGR->getName().c_str());
    if (variable_pGR->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Gas pressure process variable '%s' is a scalar variable but has "
            "%d components.",
            variable_pGR->getName().c_str(),
            variable_pGR->getNumberOfComponents());
    }

    DBUG("Associate capillary pressure with process variable '%s'.",
         variable_pCap->getName().c_str());
    if (variable_pCap->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Capillary pressure process variable '%s' is a scalar variable but "
            "has %d components.",
            variable_pCap->getName().c_str(),
            variable_pCap->getNumberOfComponents());
    }

    DBUG("Associate temperature with process variable '%s'.",
         variable_T->getName().c_str());
    if (variable_T->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "temperature process variable '%s' is a scalar variable but "
            "has %d components.",
            variable_T->getName().c_str(),
            variable_T->getNumberOfComponents());
    }

    auto solid_constitutive_relations =
        MaterialLib::Solids::createConstitutiveRelations<DisplacementDim>(
            parameters, local_coordinate_system, config);

    // reference temperature
    auto& reference_temperature = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__TH2M__reference_temperature}
        "reference_temperature", parameters, 1, &mesh);
    DBUG("Use '%s' as reference temperature parameter.",
         reference_temperature.name.c_str());

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__TH2M__specific_body_force}
            config.getConfigParameter<std::vector<double>>(
                "specific_body_force");
        if (b.size() != DisplacementDim)
        {
            OGS_FATAL(
                "The size of the specific body force vector does not match the "
                "displacement dimension. Vector size is %d, displacement "
                "dimension is %d",
                b.size(), DisplacementDim);
        }

        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    TH2MProcessData<DisplacementDim> process_data{
        materialIDs(mesh), std::move(media_map),
        std::move(solid_constitutive_relations), reference_temperature,
        specific_body_force};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<TH2MProcess<DisplacementDim>>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        use_monolithic_scheme);
}

template std::unique_ptr<Process> createTH2MProcess<2>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);

template std::unique_ptr<Process> createTH2MProcess<3>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);
}  // namespace TH2M
}  // namespace ProcessLib