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

#include <iostream>
#include "TH2MFEM.h"

#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"

#include "MaterialLib/MPL/Components/GetThermalExpansivity.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/Utils/FormEffectiveThermalConductivity.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"

#include "MathLib/KelvinVector.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
namespace TH2M
{
namespace MPL = MaterialPropertyLib;

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                   IntegrationMethod, DisplacementDim>::
    TH2MLocalAssembler(MeshLib::Element const& e,
                       std::size_t const /*local_matrix_size*/,
                       bool const is_axially_symmetric,
                       unsigned const integration_order,
                       TH2MProcessData<DisplacementDim>& process_data)
    : _process_data(process_data),
      _integration_method(integration_order),
      _element(e),
      _is_axially_symmetric(is_axially_symmetric)
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _secondary_data.N_u.resize(n_integration_points);

    auto const shape_matrices_u =
        initShapeMatrices<ShapeFunctionDisplacement,
                          ShapeMatricesTypeDisplacement, IntegrationMethod,
                          DisplacementDim>(e, is_axially_symmetric,
                                           _integration_method);

    auto const shape_matrices_p =
        initShapeMatrices<ShapeFunctionPressure, ShapeMatricesTypePressure,
                          IntegrationMethod, DisplacementDim>(
            e, is_axially_symmetric, _integration_method);

    auto const& solid_material =
        MaterialLib::Solids::selectSolidConstitutiveRelation(
            _process_data.solid_materials, _process_data.material_ids,
            e.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        _ip_data.emplace_back(solid_material);
        auto& ip_data = _ip_data[ip];
        auto const& sm_u = shape_matrices_u[ip];
        ip_data.integration_weight =
            _integration_method.getWeightedPoint(ip).getWeight() *
            sm_u.integralMeasure * sm_u.detJ;

        ip_data.N_u_op = ShapeMatricesTypeDisplacement::template MatrixType<
            DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                      displacement_size);
        for (int i = 0; i < DisplacementDim; ++i)
            ip_data.N_u_op
                .template block<1, displacement_size / DisplacementDim>(
                    i, i * displacement_size / DisplacementDim)
                .noalias() = sm_u.N;

        ip_data.N_u = sm_u.N;
        ip_data.dNdx_u = sm_u.dNdx;

        ip_data.N_p = shape_matrices_p[ip].N;
        ip_data.dNdx_p = shape_matrices_p[ip].dNdx;

        _secondary_data.N_u[ip] = shape_matrices_u[ip].N;
    }
}

// Assembles the local Jacobian matrix. So far, the linearisation of HT part is
// not considered as that in HT process.
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                        IntegrationMethod, DisplacementDim>::
    assembleWithJacobian(double const t, double const dt,
                         std::vector<double> const& local_x,
                         std::vector<double> const& local_xdot,
                         const double /*dxdot_dx*/, const double /*dx_dx*/,
                         std::vector<double>& /*local_M_data*/,
                         std::vector<double>& /*local_K_data*/,
                         std::vector<double>& local_rhs_data,
                         std::vector<double>& local_Jac_data)
{
    auto const matrix_size = gas_pressure_size + capillary_pressure_size +
                             temperature_size + displacement_size;

    assert(local_x.size() == matrix_size);

    auto pGR =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            gas_pressure_size> const>(local_x.data() + gas_pressure_index,
                                      gas_pressure_size);

    auto pCap =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            capillary_pressure_size> const>(
            local_x.data() + capillary_pressure_index, capillary_pressure_size);

    auto T = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
        temperature_size> const>(local_x.data() + temperature_index,
                                 temperature_size);

    auto u =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + displacement_index,
                                      displacement_size);

    auto pGR_dot =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            gas_pressure_size> const>(local_xdot.data() + gas_pressure_index,
                                      gas_pressure_size);

    auto pCap_dot =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            capillary_pressure_size> const>(
            local_xdot.data() + capillary_pressure_index,
            capillary_pressure_size);

    auto T_dot =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            temperature_size> const>(local_xdot.data() + temperature_index,
                                     temperature_size);

    auto u_dot =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_xdot.data() + displacement_index,
                                      displacement_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            matrix_size, matrix_size>>(local_Jac_data, matrix_size,
                                       matrix_size);

    auto local_rhs =
        MathLib::createZeroedVector<typename ShapeMatricesTypeDisplacement::
                                        template VectorType<matrix_size>>(
            local_rhs_data, matrix_size);

    // gas pressure equation
    //   - mass submatrices
//    typename ShapeMatricesTypePressure::NodalMatrixType MGpG;
//    MGpG.setZero(gas_pressure_size, gas_pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType MGpG =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(gas_pressure_size,
                                                         gas_pressure_size);

//    typename ShapeMatricesTypePressure::NodalMatrixType MGpC;
//    MGpC.setZero(gas_pressure_size, capillary_pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType MGpC =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(gas_pressure_size,
                                                         capillary_pressure_size);

//    typename ShapeMatricesTypePressure::NodalMatrixType MGT;
//    MGT.setZero(gas_pressure_size, temperature_size);

    typename ShapeMatricesTypePressure::NodalMatrixType MGT =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(gas_pressure_size,
                                                         temperature_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        gas_pressure_size, displacement_size>
        MGu = ShapeMatricesTypeDisplacement::template MatrixType<
            gas_pressure_size, displacement_size>::Zero(gas_pressure_size,
                                                    displacement_size);


    //  - laplace matrix
//    typename ShapeMatricesTypePressure::NodalMatrixType LGpG;
//    LGpG.setZero(gas_pressure_size, gas_pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType LGpG =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(gas_pressure_size,
                                                         gas_pressure_size);


    //  - rhs vector
//    typename ShapeMatricesTypePressure::NodalVectorType fG;
//    fG.setZero(gas_pressure_size);

    typename ShapeMatricesTypePressure::template MatrixType<
        gas_pressure_size, 1>
        fG = ShapeMatricesTypePressure::template MatrixType<
            gas_pressure_size, 1>::Zero(gas_pressure_size,
                                                    1);

    // capillary pressure equation
    //  - mass submatrices
//    typename ShapeMatricesTypePressure::NodalMatrixType MLpG;
//    MLpG.setZero(capillary_pressure_size, gas_pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType MLpG =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(capillary_pressure_size,
                                                         gas_pressure_size);


//    typename ShapeMatricesTypePressure::NodalMatrixType MLpC;
//    MLpC.setZero(capillary_pressure_size, capillary_pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType MLpC =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(capillary_pressure_size,
                                                         capillary_pressure_size);


//    typename ShapeMatricesTypePressure::NodalMatrixType MLT;
//    MLT.setZero(capillary_pressure_size, temperature_size);

    typename ShapeMatricesTypePressure::NodalMatrixType MLT =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(capillary_pressure_size,
                                                         temperature_size);
 
//    typename ShapeMatricesTypePressure::NodalMatrixType MLu;
//    MLu.setZero(capillary_pressure_size, displacement_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        capillary_pressure_size, displacement_size>
        MLu = ShapeMatricesTypeDisplacement::template MatrixType<
            capillary_pressure_size, displacement_size>::Zero(capillary_pressure_size,
                                                    displacement_size);


    //  - laplace submatrices
//    typename ShapeMatricesTypePressure::NodalMatrixType LLpG;
//    LLpG.setZero(capillary_pressure_size, gas_pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType LLpG =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(capillary_pressure_size,
                                                         gas_pressure_size);


//    typename ShapeMatricesTypePressure::NodalMatrixType LLpC;
//    LLpC.setZero(capillary_pressure_size, capillary_pressure_size);

typename ShapeMatricesTypePressure::NodalMatrixType LLpC =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(capillary_pressure_size,
                                                         capillary_pressure_size);


    //  - rhs vector
//    typename ShapeMatricesTypePressure::NodalVectorType fL;
//    fL.setZero(capillary_pressure_size);


    typename ShapeMatricesTypePressure::template MatrixType<
        capillary_pressure_size, 1>
        fL = ShapeMatricesTypePressure::template MatrixType<
            capillary_pressure_size, 1>::Zero(capillary_pressure_size,
                                                    1);


    // temperature equation
    //  - mass submatrices
//    typename ShapeMatricesTypePressure::NodalMatrixType MTpG;
//    MTpG.setZero(temperature_size, gas_pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType MTpG =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(temperature_size,
                                                         gas_pressure_size);


//    typename ShapeMatricesTypePressure::NodalMatrixType MTpC;
//    MTpC.setZero(temperature_size, capillary_pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType MTpC =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(temperature_size,
                                                         capillary_pressure_size);

    //typename ShapeMatricesTypePressure::NodalMatrixType MTT;
    //MTT.setZero(temperature_size, temperature_size);

    typename ShapeMatricesTypePressure::NodalMatrixType MTT =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(temperature_size,
                                                         temperature_size);

    //  - advection submatrices
//    typename ShapeMatricesTypePressure::NodalMatrixType ATpG;
//    ATpG.setZero(temperature_size, gas_pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType ATpG =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(temperature_size,
                                                         gas_pressure_size);

//    typename ShapeMatricesTypePressure::NodalMatrixType ATpC;
//    ATpC.setZero(temperature_size, capillary_pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType ATpC =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(temperature_size,
                                                         capillary_pressure_size);

//    typename ShapeMatricesTypePressure::NodalMatrixType ATT;
//    ATT.setZero(temperature_size, temperature_size);

    typename ShapeMatricesTypePressure::NodalMatrixType ATT =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(temperature_size,
                                                         temperature_size);

    //  - laplace submatrix
//    typename ShapeMatricesTypePressure::NodalMatrixType LTT;
//    LTT.setZero(temperature_size, temperature_size);

    typename ShapeMatricesTypePressure::NodalMatrixType LTT =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(temperature_size,
                                                         temperature_size);

    //  - rhs vector
//    typename ShapeMatricesTypePressure::NodalVectorType fT;
//    fT.setZero(temperature_size);

    typename ShapeMatricesTypePressure::template MatrixType<
        temperature_size, 1>
        fT = ShapeMatricesTypePressure::template MatrixType<
            temperature_size, 1>::Zero(temperature_size,
                                                    1);



    // displacement equation
    //  - stiffness submatrices
//    typename ShapeMatricesTypePressure::NodalMatrixType KUpG;
//    KUpG.setZero(displacement_size, gas_pressure_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        displacement_size, gas_pressure_size>
        KUpG = ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size, gas_pressure_size>::Zero(displacement_size,
                                                    gas_pressure_size);


//    typename ShapeMatricesTypePressure::NodalMatrixType KUpC;
//    KUpC.setZero(displacement_size, capillary_pressure_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        displacement_size, capillary_pressure_size>
        KUpC = ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size, capillary_pressure_size>::Zero(displacement_size,
                                                    capillary_pressure_size);

    //  - rhs vector
//    typename ShapeMatricesTypePressure::NodalVectorType fU;
//    fU.setZero(displacement_size);

//typename ShapeMatricesTypeDisplacement::NodalVectorType fU =
//        ShapeMatricesTypeDisplacement::NodalVectorType::Zero(displacement_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        displacement_size, 1>
        fU = ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size, 1>::Zero(displacement_size,
                                                    1);


    // pointer-matrices to the jacobian matrix
    auto JGpG = local_Jac.template block<gas_pressure_size, gas_pressure_size>(
        gas_pressure_index, gas_pressure_index);
    auto JGpC =
        local_Jac.template block<gas_pressure_size, capillary_pressure_size>(
            gas_pressure_index, capillary_pressure_index);
    auto JGT = local_Jac.template block<gas_pressure_size, temperature_size>(
        gas_pressure_index, temperature_index);
    auto JGu = local_Jac.template block<gas_pressure_size, displacement_size>(
        gas_pressure_index, displacement_index);

    auto JLpG =
        local_Jac.template block<capillary_pressure_size, gas_pressure_size>(
            capillary_pressure_index, gas_pressure_index);
    auto JLpC =
        local_Jac
            .template block<capillary_pressure_size, capillary_pressure_size>(
                capillary_pressure_index, capillary_pressure_index);
    auto JLT =
        local_Jac.template block<capillary_pressure_size, temperature_size>(
            capillary_pressure_index, temperature_index);
    auto JLu =
        local_Jac.template block<capillary_pressure_size, displacement_size>(
            capillary_pressure_index, displacement_index);

    auto JTpG = local_Jac.template block<temperature_size, gas_pressure_size>(
        temperature_index, gas_pressure_index);
    auto JTpC =
        local_Jac.template block<temperature_size, capillary_pressure_size>(
            temperature_index, capillary_pressure_index);
    auto JTT = local_Jac.template block<temperature_size, temperature_size>(
        temperature_index, temperature_index);

    auto JUpG = local_Jac.template block<displacement_size, gas_pressure_size>(
        displacement_index, gas_pressure_index);
    auto JUpC =
        local_Jac.template block<displacement_size, capillary_pressure_size>(
            displacement_index, capillary_pressure_index);
    auto JUu = local_Jac.template block<displacement_size, displacement_size>(
        displacement_index, displacement_index);

    // pointer-matrices to the residuum column-matrix
    auto rG = local_rhs.template segment<gas_pressure_size>(gas_pressure_index);
    auto rL = local_rhs.template segment<capillary_pressure_size>(
        capillary_pressure_index);
    auto rT = local_rhs.template segment<temperature_size>(temperature_index);
    auto rU = local_rhs.template segment<displacement_size>(displacement_index);

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    auto const& medium = *_process_data.media_map->getMedium(_element.getID());
    auto const& liquid_phase = medium.phase("AqueousLiquid");
    auto const& gas_phase = medium.phase("Gas");
    auto const& solid_phase = medium.phase("Solid");

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);

        auto const& Np = _ip_data[ip].N_p;
        auto const& NT = Np;
        auto const& Nu = _ip_data[ip].N_u;

        auto const& NpT = Np.transpose();
        auto const& NTT = NT.transpose();

        auto const& gradNp = _ip_data[ip].dNdx_p;
        auto const& gradNT = gradNp;
        auto const& gradNu = _ip_data[ip].dNdx_u;

        auto const& gradNpT = gradNp.transpose();
        auto const& gradNTT = gradNT.transpose();

        auto const& Nu_op = _ip_data[ip].N_u_op;
        auto const& w = _ip_data[ip].integration_weight;

        auto const& m = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>::identity2;

        auto const mT = m.transpose();

        auto const x_coord =
            interpolateXCoordinate<ShapeFunctionDisplacement,
                                   ShapeMatricesTypeDisplacement>(_element, Nu);

        auto const Bu =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                gradNu, Nu, x_coord, _is_axially_symmetric);

        auto const BuT = Bu.transpose();

        auto& eps = _ip_data[ip].eps;
        auto const& sigma_eff = _ip_data[ip].sigma_eff;
        double const T0 = _process_data.reference_temperature(t, pos)[0];
#define nDEBUG_TH2M

        auto const T_int_pt = NT.dot(T);
        auto const pGR_int_pt = Np.dot(pGR);
        auto const pCap_int_pt = Np.dot(pCap);
        auto const pLR_int_pt = pGR_int_pt - pCap_int_pt;

#ifdef DEBUG_TH2M
        std::cout << "-----------------\n";
        std::cout << "--- unknowns: ---\n";
        std::cout << "pGR: " << pGR << "\n";
        std::cout << "pCap: " << pCap << "\n";
        std::cout << "T: " << T << "\n";
        std::cout << "--------------------\n";

        std::cout << "---------------------\n";
        std::cout << "--- unknowns(IP): ---\n";
        std::cout << "pGR_int_pt: " << pGR_int_pt << "\n";
        std::cout << "pCap_int_pt: " << pCap_int_pt << "\n";
        std::cout << "T_int_pt: " << T_int_pt << "\n";
        std::cout << "--------------------\n";
#endif

#ifdef DEBUG_TH2M
        std::cout << "*************************************\n";
        std::cout << " Shape matrices: \n";
        std::cout << " --------------- \n";
        std::cout << " Np:\n" << Np << "\n";
        std::cout << " --------------- \n";
        std::cout << " Nu:\n" << Nu << "\n";
        std::cout << " --------------- \n";
        std::cout << " Nu_op:\n" << Nu_op << "\n";
        std::cout << " --------------- \n";
        std::cout << " gradNp:\n" << gradNp << "\n";
        std::cout << " --------------- \n";
        std::cout << " Bu:\n" << Bu << "\n";
        std::cout << " --------------- \n";
        std::cout << "*************************************\n";
        std::cout << " Process variables: \n";
        std::cout << " --------------- \n";
        std::cout << " Bu:\n" << Bu << "\n";
        std::cout << " --------------- \n";

        std::cout << " Calculate constitutive parameters. \n";

#endif


        MPL::VariableArray vars;
        vars[static_cast<int>(MPL::Variable::temperature)] = T_int_pt;
        vars[static_cast<int>(MPL::Variable::gas_phase_pressure)] = pGR_int_pt;
        vars[static_cast<int>(MPL::Variable::capillary_pressure)] = pCap_int_pt;
        vars[static_cast<int>(MPL::Variable::liquid_phase_pressure)] =
            pLR_int_pt;

        // Material properties
        //  - solid phase properties
        auto const beta_p_SR =
            solid_phase.property(MPL::PropertyType::compressibility)
                .template value<double>(vars, pos, t);

        auto const beta_T_SR =
            solid_phase.property(MPL::PropertyType::thermal_expansivity)
                .template value<double>(vars, pos, t);

        auto const rho_SR = solid_phase.property(MPL::PropertyType::density)
                                .template value<double>(vars, pos, t);

        auto const c_p_S =
            solid_phase.property(MPL::PropertyType::specific_heat_capacity)
                .template value<double>(vars, pos, t);

        //  - gas phase properties
        auto const beta_p_GR =
            gas_phase.property(MPL::PropertyType::compressibility)
                .template value<double>(vars, pos, t);

        auto const beta_T_GR =
            gas_phase.property(MPL::PropertyType::thermal_expansivity)
                .template value<double>(vars, pos, t);

        auto const mu_GR = gas_phase.property(MPL::PropertyType::viscosity)
                               .template value<double>(vars, pos, t);

        auto const rho_GR = gas_phase.property(MPL::PropertyType::density)
                                .template value<double>(vars, pos, t);

        auto const c_p_G =
            gas_phase.property(MPL::PropertyType::specific_heat_capacity)
                .template value<double>(vars, pos, t);

        //  - liquid phase properties
        auto const beta_p_LR =
            liquid_phase.property(MPL::PropertyType::compressibility)
                .template value<double>(vars, pos, t);

        auto const beta_T_LR =
            liquid_phase.property(MPL::PropertyType::thermal_expansivity)
                .template value<double>(vars, pos, t);

        auto const mu_LR = liquid_phase.property(MPL::PropertyType::viscosity)
                               .template value<double>(vars, pos, t);

        auto const rho_LR = liquid_phase.property(MPL::PropertyType::density)
                                .template value<double>(vars, pos, t);

        auto const c_p_L =
            liquid_phase.property(MPL::PropertyType::specific_heat_capacity)
                .template value<double>(vars, pos, t);

        //  - medium properties
        auto const k_S = MPL::formEigenTensor<DisplacementDim>(
            medium.property(MPL::PropertyType::permeability)
                .value(vars, pos, t));

        auto const s_L = medium.property(MPL::PropertyType::saturation)
                             .template value<double>(vars, pos, t);

        auto const s_G = 1. - s_L;

        auto const dsLdPc =
            medium.property(MPL::PropertyType::saturation)
                .template dValue<double>(
                    vars, MPL::Variable::capillary_pressure, pos, t);

        auto const alpha_B =
            medium.property(MPL::PropertyType::biot_coefficient)
                .template value<double>(vars, pos, t);

        auto const k_rel =
            medium.property(MPL::PropertyType::relative_permeability)
                .template value<MPL::Pair>(vars, pos, t);

        auto const k_rel_L = k_rel[0];
        auto const k_rel_G = k_rel[1];

        auto const& b = _process_data.specific_body_force;

        auto const rho = rho_GR + rho_LR + rho_SR;
        auto const rho_c_p = rho_GR * c_p_G + rho_LR * c_p_L + rho_SR * c_p_S;

        auto const lambda = MPL::formEigenTensor<DisplacementDim>(
            medium.property(MPL::PropertyType::specific_heat_capacity)
                .template value<double>(vars, pos, t));

        auto const k_over_mu_G = k_S * k_rel_G / mu_GR;
        auto const k_over_mu_L = k_S * k_rel_L / mu_LR;

        auto const phi = medium.property(MPL::PropertyType::porosity)
                             .template value<double>(vars, pos, t);

        auto const phi_G = s_G * phi;
        auto const phi_L = s_L * phi;
        auto const phi_S = 1. - phi;

#ifdef DEBUG_TH2M
        std::cout << " assemble sub-matrices: \n";
#endif

        // coefficient matrices
        //  - gas pressure equation
        MGpG.noalias() +=
            (NpT * s_G * (phi * beta_p_GR + (alpha_B - phi) * beta_p_SR) * Np) *
            w;
        MGpC.noalias() += 
            (NpT *
             (s_G * (alpha_B - phi) * beta_p_SR * (s_L + pCap_int_pt * dsLdPc) +
              phi * dsLdPc) *
             Np) *
            w;
        MGT.noalias() +=
            (NpT * s_G * (phi * beta_T_GR + (alpha_B - phi) * beta_T_SR) * Np) *
            w;

       MGu.noalias() += (NpT * mT * Bu).eval() * s_G * alpha_B * w;


        LGpG.noalias() += (gradNpT * k_over_mu_G * gradNp) * w;
        fG.noalias() += (gradNpT * rho_GR * k_over_mu_G * b) * w;

        //  - liquid pressure equation
        MLpG.noalias() +=
            (NpT * s_L * (phi * beta_p_LR + (alpha_B - phi) * beta_p_SR) * Np) *
            w;
        MLpC.noalias() +=
            (NpT *
             (s_L * (alpha_B - phi) * beta_p_SR * (s_L + pCap_int_pt * dsLdPc) +
              phi_L * beta_p_LR - phi * dsLdPc) *
             Np) *
            w;
        MLT.noalias() +=
            (NpT * s_L * (phi * beta_T_LR + (alpha_B - phi) * beta_T_SR) * Np) *
            w;
        MLu.noalias() += (NpT * s_L * alpha_B * mT * Bu) * w;
        LLpG.noalias() += (gradNpT * k_over_mu_L * gradNp) * w;
        LLpC.noalias() += (gradNpT * k_over_mu_L * gradNp) * w;
        fL.noalias() += (gradNpT * rho_LR * k_over_mu_L * b) * w;

        // darcy-velocities
        auto const w_GS = -k_over_mu_G * (gradNp * pGR - rho_GR * b);

        auto const w_LS =
            -k_over_mu_L * (gradNp * pGR - gradNp * pCap - rho_GR * b);

#ifdef DEBUG_TH2M
        std::cout << "--------------------\n";
        std::cout << "--- velocities:  ---\n";
        std::cout << "w_GS: " << w_GS << "\n";
        std::cout << "w_LS: " << w_LS << "\n";
        std::cout << "--------------------\n";
#endif

        //  - temperature equation
        MTpG.noalias() +=
            (NTT * (phi_G * beta_T_GR + phi_L * beta_T_LR + phi_S * beta_T_SR) *
             T_int_pt * NT) *
            w;
        MTpC.noalias() +=
            (NTT *
             ((phi_L * beta_T_LR +
               phi_S * beta_T_SR * (s_L + pCap_int_pt * dsLdPc) * T_int_pt) +
              phi * pCap_int_pt * dsLdPc) *
             NT) *
            w;
        MTT.noalias() += (NTT * rho_c_p * NT) * w;
        ATpG.noalias() +=
            (NTT *
             (beta_T_GR * w_GS.transpose() + beta_T_LR * w_LS.transpose()) *
             gradNT) *
            w;
        ATpC.noalias() += (NTT * beta_T_LR * w_LS.transpose() * gradNT) * w;
        ATT.noalias() += (NTT *
                          (rho_LR * c_p_L * w_LS.transpose() +
                           rho_GR * c_p_G * w_GS.transpose()) *
                          gradNT) *
                         w;
        LTT.noalias() += (gradNTT * lambda * gradNT) * w;

        //  - displacement equation
        KUpG.noalias() += (BuT * alpha_B * m * Np) * w;
        KUpC.noalias() += (BuT * alpha_B * s_L * m * Np) * w;

#ifdef DEBUG_TH2M
        std::cout << " Bu " << Bu.size() << " " << Bu.rows() << " " << Bu.cols() << "\n";
        std::cout << " -------------------------------------------------------\n";
        std::cout << " BuT " << BuT.size() << " " << BuT.rows() << " " << BuT.cols() << "\n";
        std::cout << " Nu_op " << Nu_op << "\n";
        std::cout << " Np " << Np.size() << " " << Np.rows() << " " << Np.cols() << "\n";
        std::cout << " NpT " << NpT.size() << " " << NpT.rows() << " " << NpT.cols() << "\n";
        std::cout << " NpT " << NpT << "\n";
        std::cout << " m " << m << "\n";
        std::cout << " m " << m.size() << " " << m.rows() << " " << m.cols() << "\n";
        std::cout << " mT " << mT.size() << " " << mT.rows() << " " << mT.cols() << "\n";


        std::cout << " MGu " << MGu << "\n" << MGu.rows() << " " << MGu.cols() << "\n";
        std::cout << " b " << b << "\n";
		OGS_FATAL("_______________________");

#endif

   fU.noalias() += (BuT * sigma_eff - Nu_op.transpose() * rho * b).eval() *w;
     
        // TODO (Wenqing) : Change dT to time step wise increment
        double const delta_T(T_int_pt - T0);
        double const thermal_strain = beta_T_SR * delta_T;
        //
        // displacement equation, displacement part
        //
        eps.noalias() = Bu * u;
        auto C = _ip_data[ip].updateConstitutiveRelationThermal(
            t, pos, dt, u, _process_data.reference_temperature(t, pos)[0],
            thermal_strain);

        JUu.noalias() += BuT * C * Bu * w;

    }

    JGpG.noalias() = MGpG / dt + LGpG;

#ifdef DEBUG_TH2M
        std::cout << " JGpG:\n" << "\n";
        std::cout << JGpG << "\n";
#endif

    JGpC.noalias() = -MGpC / dt;
    JGT.noalias() = -MGT / dt;
    JGu.noalias() = MGu / dt;

    JLpG.noalias() = MLpG / dt + LLpG;
    JLpC.noalias() = -MLpC / dt - LLpC;
    JLT.noalias() = -MLT / dt;
    JLu.noalias() = MLu / dt;

    JTpG.noalias() = -MTpG / dt - ATpG;
    JTpC.noalias() = MTpC / dt + ATpC;
    JTT.noalias() = MTT / dt + ATT + LTT;

    JUpG.noalias() = -KUpG;
    JUpC.noalias() = -KUpC;

    rG.noalias() = -1*(MGpG * pGR_dot - MGpC * pCap_dot - MGT * T_dot +
                   MGu * u_dot + LGpG * pGR - fG);
    rL.noalias() = -1*(MLpG * pGR_dot - MLpC * pCap_dot - MLT * T_dot +
                   MLu * u_dot + LLpG * pGR - LLpC * pCap - fL);
    rT.noalias() = -1*(-MTpG * pGR_dot + MTpC * pCap_dot + MTT * T_dot +
                   (ATT + LTT) * T - ATpG * pGR + ATpC * pCap + fT);
    rU.noalias() = -1*(fU - KUpG * pGR + KUpC * pCap);

#ifdef DEBUG_TH2M
    std::cout << "--------------------------\n";
    std::cout << local_Jac << "\n";
    std::cout << "--------------------------\n";
    std::cout << local_rhs << "\n";
#endif

}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const&
TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                   IntegrationMethod, DisplacementDim>::
    getIntPtDarcyVelocityGas(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const
{
    auto const num_intpts = _ip_data.size();

    constexpr int process_id = 0;  // monolithic scheme;
    auto const indices =
        NumLib::getIndices(_element.getID(), *dof_table[process_id]);
    assert(!indices.empty());
    auto const local_x = x[process_id]->get(indices);

    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, num_intpts);

    auto pGR =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            gas_pressure_size> const>(local_x.data() + gas_pressure_index,
                                      gas_pressure_size);
    auto T = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
        temperature_size> const>(local_x.data() + temperature_index,
                                 temperature_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    auto const& medium = _process_data.media_map->getMedium(_element.getID());
    auto const& gas_phase = medium->phase("NonAqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");
    MPL::VariableArray vars;

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);

        auto const& N_p = _ip_data[ip].N_p;

        vars[static_cast<int>(MPL::Variable::temperature)] =
            N_p * T;  // N_p = N_T
        vars[static_cast<int>(MPL::Variable::phase_pressure)] = N_p * pGR;

        auto const viscosity = gas_phase.property(MPL::PropertyType::viscosity)
                                   .template value<double>(vars, pos, t);
        GlobalDimMatrixType K_over_mu =
            MPL::formEigenTensor<DisplacementDim>(
                solid_phase.property(MPL::PropertyType::permeability)
                    .value(vars, pos, t)) /
            viscosity;

        auto const fluid_density =
            gas_phase.property(MPL::PropertyType::density)
                .template value<double>(vars, pos, t);
        auto const& b = _process_data.specific_body_force;

        // Compute the velocity
        auto const& dNdx_p = _ip_data[ip].dNdx_p;
        cache_matrix.col(ip).noalias() =
            -K_over_mu * dNdx_p * pGR + K_over_mu * fluid_density * b;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const&
TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                   IntegrationMethod, DisplacementDim>::
    getIntPtDarcyVelocityLiquid(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const
{
    auto const num_intpts = _ip_data.size();

    constexpr int process_id = 0;  // monolithic scheme;
    auto const indices =
        NumLib::getIndices(_element.getID(), *dof_table[process_id]);
    assert(!indices.empty());
    auto const local_x = x[process_id]->get(indices);

    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, num_intpts);

    auto pGR =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            gas_pressure_size> const>(local_x.data() + gas_pressure_index,
                                      gas_pressure_size);
    auto pCap =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            capillary_pressure_size> const>(
            local_x.data() + capillary_pressure_index, capillary_pressure_size);
    auto pLR = pGR - pCap;
    auto T = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
        temperature_size> const>(local_x.data() + temperature_index,
                                 temperature_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    auto const& medium = _process_data.media_map->getMedium(_element.getID());
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");
    MPL::VariableArray vars;

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);

        auto const& N_p = _ip_data[ip].N_p;

        vars[static_cast<int>(MPL::Variable::temperature)] =
            N_p * T;  // N_p = N_T
        vars[static_cast<int>(MPL::Variable::phase_pressure)] = N_p * pLR;

        auto const viscosity =
            liquid_phase.property(MPL::PropertyType::viscosity)
                .template value<double>(vars, pos, t);
        GlobalDimMatrixType K_over_mu =
            MPL::formEigenTensor<DisplacementDim>(
                solid_phase.property(MPL::PropertyType::permeability)
                    .value(vars, pos, t)) /
            viscosity;

        auto const fluid_density =
            liquid_phase.property(MPL::PropertyType::density)
                .template value<double>(vars, pos, t);
        auto const& b = _process_data.specific_body_force;

        // Compute the velocity
        auto const& dNdx_p = _ip_data[ip].dNdx_p;
        cache_matrix.col(ip).noalias() =
            -K_over_mu * dNdx_p * pLR + K_over_mu * fluid_density * b;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                        IntegrationMethod, DisplacementDim>::
    postNonLinearSolverConcrete(std::vector<double> const& local_x,
                                double const t, double const dt,
                                bool const use_monolithic_scheme)
{
    DBUG(
        "Warning(TODO): postNonLinearSolverConcrete is not fully "
        "configurated for TH2M!");

    const int displacement_offset =
        use_monolithic_scheme ? displacement_index : 0;

    auto u =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + displacement_offset,
                                      displacement_size);

    auto T = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
        temperature_size> const>(local_x.data() + temperature_index,
                                 temperature_size);
    auto pGR =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            gas_pressure_size> const>(local_x.data() + gas_pressure_index,
                                      gas_pressure_size);

    auto pCap =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            capillary_pressure_size> const>(
            local_x.data() + capillary_pressure_index, capillary_pressure_size);

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());
    auto const& medium = _process_data.media_map->getMedium(_element.getID());
    auto const& solid_phase = medium->phase("Solid");
    MPL::VariableArray vars;

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);
        auto const& N_u = _ip_data[ip].N_u;
        auto const& N_p = _ip_data[ip].N_p;
        auto const& N_T = N_p;
        auto const& dNdx_u = _ip_data[ip].dNdx_u;

        auto const x_coord =
            interpolateXCoordinate<ShapeFunctionDisplacement,
                                   ShapeMatricesTypeDisplacement>(_element,
                                                                  N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        double const T0 = _process_data.reference_temperature(t, pos)[0];

        double const T_int_pt = N_T * T;
        vars[static_cast<int>(MPL::Variable::temperature)] = T_int_pt;
        vars[static_cast<int>(MPL::Variable::gas_phase_pressure)] = N_p * pGR;
        vars[static_cast<int>(MPL::Variable::capillary_pressure)] = N_p * pCap;

        auto const solid_linear_thermal_expansion_coefficient =
            solid_phase.property(MPL::PropertyType::thermal_expansivity)
                .template value<double>(vars, pos, t);

        double const delta_T(T_int_pt - T0);
        double const thermal_strain =
            solid_linear_thermal_expansion_coefficient * delta_T;

        auto& eps = _ip_data[ip].eps;
        eps.noalias() = B * u;

        _ip_data[ip].updateConstitutiveRelationThermal(
            t, pos, dt, u, _process_data.reference_temperature(t, pos)[0],
            thermal_strain);
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                        IntegrationMethod, DisplacementDim>::
    computeSecondaryVariableConcrete(double const /*t*/,
                                     std::vector<double> const& local_x)
{
    auto pGR =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            gas_pressure_size> const>(local_x.data() + gas_pressure_index,
                                      gas_pressure_size);
    auto pCap =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            capillary_pressure_size> const>(
            local_x.data() + capillary_pressure_index, capillary_pressure_size);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(_element, _is_axially_symmetric, pGR,
                         *_process_data.gas_pressure_interpolated);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(_element, _is_axially_symmetric, pCap,
                         *_process_data.capillary_pressure_interpolated);

    auto T = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
        temperature_size> const>(local_x.data() + temperature_index,
                                 temperature_size);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(_element, _is_axially_symmetric, T,
                         *_process_data.temperature_interpolated);
}

}  // namespace TH2M
}  // namespace ProcessLib
