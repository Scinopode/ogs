/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "ParameterLib/Parameter.h"

#include <iostream>

namespace ProcessLib
{
namespace TH2M
{
template <typename BMatricesType, typename ShapeMatrixTypeDisplacement,
          typename ShapeMatricesTypePressure, int DisplacementDim, int NPoints>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables())
    {
        // Initialize current time step values
        static const int kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
        sigma_eff.setZero(kelvin_vector_size);
        eps.setZero(kelvin_vector_size);
        eps_m.setZero(kelvin_vector_size);
        eps_m_prev.resize(kelvin_vector_size);

        // Previous time step values are not initialized and are set later.
        eps_prev.resize(kelvin_vector_size);
        sigma_eff_prev.resize(kelvin_vector_size);
    }

    typename ShapeMatrixTypeDisplacement::template MatrixType<
        DisplacementDim, NPoints * DisplacementDim>
        N_u_op;
    typename BMatricesType::KelvinVectorType sigma_eff, sigma_eff_prev;
    typename BMatricesType::KelvinVectorType eps, eps_prev;
    typename BMatricesType::KelvinVectorType eps_m, eps_m_prev;

    typename ShapeMatrixTypeDisplacement::NodalRowVectorType N_u;
    typename ShapeMatrixTypeDisplacement::GlobalDimNodalMatrixType dNdx_u;

    typename ShapeMatricesTypePressure::NodalRowVectorType N_p;
    typename ShapeMatricesTypePressure::GlobalDimNodalMatrixType dNdx_p;

    double saturation = std::numeric_limits<double>::quiet_NaN();
    double saturation_prev = std::numeric_limits<double>::quiet_NaN();

    // phase instrinsic densities
    double rho_GR;
    double rho_LR;
    double rho_SR;

    // solid phase pressure
    double p_SR;
    
    // real constitutent partial densities
    double rho_C_GR = std::numeric_limits<double>::quiet_NaN();
    double rho_C_GR_prev = std::numeric_limits<double>::quiet_NaN();
    double rho_W_GR = std::numeric_limits<double>::quiet_NaN();
    double rho_W_GR_prev = std::numeric_limits<double>::quiet_NaN();
    double rho_C_LR = std::numeric_limits<double>::quiet_NaN();
    double rho_C_LR_prev = std::numeric_limits<double>::quiet_NaN();
    double rho_W_LR = std::numeric_limits<double>::quiet_NaN();
    double rho_W_LR_prev = std::numeric_limits<double>::quiet_NaN();

    // phase enthalpies
    double rho_G_h_G = std::numeric_limits<double>::quiet_NaN();
    double rho_G_h_G_prev = std::numeric_limits<double>::quiet_NaN();
    double rho_L_h_L = std::numeric_limits<double>::quiet_NaN();
    double rho_L_h_L_prev = std::numeric_limits<double>::quiet_NaN();
    double rho_S_h_S = std::numeric_limits<double>::quiet_NaN();
    double rho_S_h_S_prev = std::numeric_limits<double>::quiet_NaN();

    double phi_S_p_SR = std::numeric_limits<double>::quiet_NaN();
    double phi_S_p_SR_prev = std::numeric_limits<double>::quiet_NaN();
    
    

    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;
    double integration_weight;

    void pushBackState()
    {
        eps_m_prev = eps_m;
        sigma_eff_prev = sigma_eff;
        saturation_prev = saturation;
        
        rho_G_h_G_prev = rho_G_h_G;
        rho_L_h_L_prev = rho_L_h_L;
        rho_S_h_S_prev = rho_S_h_S;

        phi_S_p_SR_prev = phi_S_p_SR;

        rho_C_GR_prev = rho_C_GR;
        rho_W_GR_prev = rho_W_GR;
        rho_C_LR_prev = rho_C_LR;
        rho_W_LR_prev = rho_W_LR;
 
        material_state_variables->pushBackState();
    }

    template <typename DisplacementVectorType>
    typename BMatricesType::KelvinMatrixType updateConstitutiveRelation(
        double const t,
        ParameterLib::SpatialPosition const& x_position,
        double const dt,
        DisplacementVectorType const& /*u*/,
        double const T)
    {
        auto&& solution = solid_material.integrateStress(
            t, x_position, dt, eps_m_prev, eps_m, sigma_eff_prev,
            *material_state_variables, T);

        if (!solution)
            OGS_FATAL("Computation of local constitutive relation failed.");

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma_eff, material_state_variables, C) = std::move(*solution);

        return C;
    }

    template <typename DisplacementVectorType>
    typename BMatricesType::KelvinMatrixType updateConstitutiveRelationThermal(
        double const t,
        ParameterLib::SpatialPosition const& x_position,
        double const dt,
        DisplacementVectorType const& /*u*/,
        double const T,
        double const thermal_strain)
    {
        auto const& identity2 = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>::identity2;

        // assume isotropic thermal expansion
        eps_m.noalias() = eps - thermal_strain * identity2;
        auto&& solution = solid_material.integrateStress(
            t, x_position, dt, eps_m_prev, eps_m, sigma_eff_prev,
            *material_state_variables, T);

        if (!solution)
            OGS_FATAL("Computation of local constitutive relation failed.");

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma_eff, material_state_variables, C) = std::move(*solution);

        return C;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace TH2M
}  // namespace ProcessLib
