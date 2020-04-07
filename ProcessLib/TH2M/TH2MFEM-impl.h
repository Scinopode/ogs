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

#include <iomanip>
#include <iostream>

#include "MaterialLib/MPL/Components/GetThermalExpansivity.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/Utils/FormEffectiveThermalConductivity.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"
#include "TH2MFEM.h"

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
      _is_axially_symmetric(is_axially_symmetric),
      _liquid_pressure(
          std::vector<double>(_integration_method.getNumberOfPoints())),
      _liquid_density(
          std::vector<double>(_integration_method.getNumberOfPoints())),
      _gas_density(
          std::vector<double>(_integration_method.getNumberOfPoints())),
      _porosity(std::vector<double>(_integration_method.getNumberOfPoints())),
      _saturation(std::vector<double>(_integration_method.getNumberOfPoints()))
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
void TH2MLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::assemble(double const t, double const dt,
                               std::vector<double> const& local_x,
                               std::vector<double> const& /*local_xdot*/,
                               std::vector<double>& local_M_data,
                               std::vector<double>& local_K_data,
                               std::vector<double>& local_rhs_data)
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

    // pointer to local_M_data vector
    auto M = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            matrix_size, matrix_size>>(local_M_data, matrix_size, matrix_size);

    // pointer to local_K_data vector
    auto K = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            matrix_size, matrix_size>>(local_K_data, matrix_size, matrix_size);

    // pointer to local_K_data vector
    auto f =
        MathLib::createZeroedVector<typename ShapeMatricesTypeDisplacement::
                                        template VectorType<matrix_size>>(
            local_rhs_data, matrix_size);

    // pointer-matrices to the mass matrix - gas phase equation
    auto MGpG = M.template block<gas_pressure_size, gas_pressure_size>(
        gas_pressure_index, gas_pressure_index);
    auto MGpC = M.template block<gas_pressure_size, capillary_pressure_size>(
        gas_pressure_index, capillary_pressure_index);
    auto MGT = M.template block<gas_pressure_size, temperature_size>(
        gas_pressure_index, temperature_index);
    auto MGu = M.template block<gas_pressure_size, displacement_size>(
        gas_pressure_index, displacement_index);

    // pointer-matrices to the mass matrix - liquid phase equation
    auto MLpG = M.template block<capillary_pressure_size, gas_pressure_size>(
        capillary_pressure_index, gas_pressure_index);
    auto MLpC =
        M.template block<capillary_pressure_size, capillary_pressure_size>(
            capillary_pressure_index, capillary_pressure_index);
    auto MLT = M.template block<capillary_pressure_size, temperature_size>(
        capillary_pressure_index, temperature_index);
    auto MLu = M.template block<capillary_pressure_size, displacement_size>(
        capillary_pressure_index, displacement_index);

    // pointer-matrices to the mass matrix - temperature equation
    auto MTpG = M.template block<temperature_size, gas_pressure_size>(
        temperature_index, gas_pressure_index);
    auto MTpC = M.template block<temperature_size, capillary_pressure_size>(
        temperature_index, capillary_pressure_index);
    auto MTT = M.template block<temperature_size, temperature_size>(
        temperature_index, temperature_index);

    // pointer-matrices to the stiffness matrix - gas phase equation
    auto KGpG = K.template block<gas_pressure_size, gas_pressure_size>(
        gas_pressure_index, gas_pressure_index);

    // pointer-matrices to the stiffness matrix - liquid phase equation
    auto KLpG = K.template block<capillary_pressure_size, gas_pressure_size>(
        capillary_pressure_index, gas_pressure_index);
    auto KLpC =
        K.template block<capillary_pressure_size, capillary_pressure_size>(
            capillary_pressure_index, capillary_pressure_index);

    // pointer-matrices to the stiffness matrix - temperature equation
    auto KTpG = K.template block<temperature_size, gas_pressure_size>(
        temperature_index, gas_pressure_index);
    auto KTpC = K.template block<temperature_size, capillary_pressure_size>(
        temperature_index, capillary_pressure_index);
    auto KTT = K.template block<temperature_size, temperature_size>(
        temperature_index, temperature_index);

    // pointer-matrices to the stiffness matrix - displacement equation
    auto KUpG = K.template block<displacement_size, gas_pressure_size>(
        displacement_index, gas_pressure_index);
    auto KUpC = K.template block<displacement_size, capillary_pressure_size>(
        displacement_index, capillary_pressure_index);
    auto KUU = K.template block<displacement_size, displacement_size>(
        displacement_index, displacement_index);

    // pointer-vectors to the right hand side terms - gas phase equation
    auto fG = f.template segment<gas_pressure_size>(gas_pressure_index);
    // pointer-vectors to the right hand side terms - liquid phase equation
    auto fL =
        f.template segment<capillary_pressure_size>(capillary_pressure_index);
    // pointer-vectors to the right hand side terms - displacement equation
    auto fU = f.template segment<displacement_size>(displacement_index);

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
        auto& ip_data = _ip_data[ip];

        auto const& Np = ip_data.N_p;
        auto const& NT = Np;
        auto const& Nu = ip_data.N_u;

        auto const& NpT = Np.transpose();
        auto const& NTT = NT.transpose();

        auto const& gradNp = ip_data.dNdx_p;
        auto const& gradNT = gradNp;
        auto const& gradNu = ip_data.dNdx_u;

        auto const& gradNpT = gradNp.transpose();
        auto const& gradNTT = gradNT.transpose();

        auto const& Nu_op = ip_data.N_u_op;
        auto const& w = ip_data.integration_weight;

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

        auto& eps = ip_data.eps;
        auto const& sigma_eff = _ip_data[ip].sigma_eff;
        double const T0 = _process_data.reference_temperature(t, pos)[0];

        auto const T_int_pt = NT.dot(T);
        auto const pGR_int_pt = Np.dot(pGR);
        auto const pCap_int_pt = Np.dot(pCap);
        auto const pLR_int_pt = pGR_int_pt - pCap_int_pt;

#define nDEBUG_OUTPUT
#ifdef DEBUG_OUTPUT
#define UNKNOWNS
#define SHAPE_MATRICES
#define MATERIAL_PROPERTIES
#endif

#ifdef UNKNOWNS
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

#ifdef SHAPE_MATRICES
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
#endif

        MPL::VariableArray vars;
        vars[static_cast<int>(MPL::Variable::temperature)] = T_int_pt;
        vars[static_cast<int>(MPL::Variable::phase_pressure)] = pGR_int_pt;
        vars[static_cast<int>(MPL::Variable::capillary_pressure)] = pCap_int_pt;
        vars[static_cast<int>(MPL::Variable::liquid_phase_pressure)] =
            pLR_int_pt;

        // Henry constants for carbon dioxide and hydrogen
        const double H_theta_CO2 = 3.3e-4;              // mol/m^3/Pa
        const double dH_solve_CO2 = -19954.7102835678;  // J/mol

        const double H_theta_H2 = 7.8e-6;             // mol/m^3/Pa
        const double dH_solve_H2 = -4406.6651876212;  // J/mol
        const double R = 8.31446261815324;            // J/mol/K
        const double T_theta = 298.15;                // K

        const double dH_by_R = dH_solve_CO2 / R;

        // CO2-Henry coefficient depending on T
        const double H_C =
            H_theta_CO2 *
            std::exp((-1.) * dH_by_R * (1. / T_int_pt - 1. / T_theta));
        const double dHc_dT = 1. / ((T_int_pt) * (T_int_pt)) * dH_by_R * H_C;

        // coefficients for Antoine equation
        const double A = 8.07131;
        const double B = 1730.63;
        const double C = 233.426;

        // saturation vapour pressure (Antoine equation)
        const double p_vap = std::pow(10., A - B / (C + T_int_pt));
        const double ln10 = 2.30258509299;
        const double dp_vap_dT =
            B * ln10 / ((C + T_int_pt) * (C + T_int_pt)) * p_vap;

        // partial pressure of constituents in gas phase
        const double p_W_GR = p_vap;
        const double p_C_GR = pGR_int_pt - p_W_GR;  // Daltons law

        // molar fraction of gas phase constituents
        const double xn_W_G = p_W_GR / pGR_int_pt;
        const double xn_CO2_G = 1. - xn_W_G;
        const double xn_H2_G = 1. - xn_W_G;
        const double xn_C_G = xn_CO2_G;

        // molar masses and average molar mass
        const double M_W = .01801528;    // kg/mol
        const double M_CO2 = .0440095;   // kg/mol
        const double M_H2 = 0.00201948;  // kg/mol

        const double M_G = xn_C_G * M_CO2 + xn_W_G * M_W;

        // density of the gas phase (ideal gas law)
        const double rho_GR = pGR_int_pt * M_G / R / T_int_pt;

        // mass fractions of gas phase constituents
        const double xm_W_G = xn_W_G * M_W / M_G;
        const double xm_C_G = 1. - xm_W_G;

        // partial densities of gas phase constituents
        const double rho_W_GR = xm_W_G * pGR_int_pt;
        const double rho_C_GR = xm_C_G * pGR_int_pt;

        // concentration of dissolved gas in water
        const double c_C_L = rho_C_GR * H_C;

        // Vapour-Liquid-Equilibrium:
        const double beta_pLR = 1.0e-8;
        const double beta_TLR = -1.0e-4;
        const double beta_CLR = 1.0e-5;

        const double T_ref = 290.;
        const double p_ref = 1.0e5;
        const double rho_LR_ref = 999.0;

        // water component partial density liquid phase
        const double rho_W_LR =
            rho_LR_ref * (1. + beta_pLR * (pLR_int_pt - p_ref) +
                          beta_TLR * (T_int_pt - T_ref));

        // liquid phase density
        const double rho_LR =
            rho_LR_ref * (1. + beta_pLR * (pLR_int_pt - p_ref) +
                          beta_TLR * (T_int_pt - T_ref) + beta_CLR * c_C_L);

        // partial density of dissolved gas component in water
        const double rho_C_LR = rho_LR - rho_W_LR;

        // mass fractions of liquid phase constituents
        const double xm_W_L = rho_W_LR / rho_W_LR;
        const double xm_C_L = 1. - xm_W_L;

        // constituent compressibility of the gas phase
        const double beta_C_pGR = 1. / p_C_GR;
        const double beta_W_pGR = 1. / p_W_GR;

        // constituent thermal expansivity of the gas phase
        const double beta_C_TGR = -beta_C_pGR * dp_vap_dT - 1 / T_int_pt;
        const double beta_W_TGR = beta_W_pGR * dp_vap_dT - 1 / T_int_pt;

        // constituent compressibility of the liquid phase
        const double beta_C_pLR = beta_CLR * H_C;
        const double beta_C_TLR =
            beta_CLR * (dHc_dT * p_C_GR - dp_vap_dT * H_C);

        // constituent thermal expansivity of the liquid phase
        const double beta_W_pLR = beta_pLR;
        const double beta_W_TLR = beta_TLR;

        const double beta_pGR = 1. / pGR_int_pt;
        const double beta_TGR = -1. / T_int_pt;

        const double dxm_C_G_dpGR = xm_C_G * (beta_C_pGR - beta_pGR);
        const double dxm_W_G_dpGR = xm_W_G * (beta_W_pGR - beta_pGR);
        const double dxm_C_L_dpLR = xm_C_L * (beta_C_pLR - beta_pLR);
        const double dxm_W_L_dpLR = xm_W_L * (beta_W_pLR - beta_pLR);

        const double dxm_C_G_dT = xm_C_G * (beta_C_TGR - beta_TGR);
        const double dxm_C_L_dT = xm_C_L * (beta_C_TLR - beta_TLR);
        const double dxm_W_G_dT = xm_W_G * (beta_W_TGR - beta_TGR);
        const double dxm_W_L_dT = xm_W_L * (beta_W_TLR - beta_TLR);

        //  - solid phase properties
        const double beta_pS = 1e-10;
        const double beta_TSR = -1e-7;
        const double rho_ref_SR = 2500.;

        const double c_pS = 1000.;
        const double lambda_SR = 0.5;

        //  - gas phase properties
        const double mu_GR = 1.e-5;
        const double c_pG = 1000.;

        const double lambda_GR = 0.5;

        //  - liquid phase properties
        const double mu_LR = 1.e-3;
        const double c_pL = 4000.;

        const double lambda_LR = 0.5;

        //  - medium properties
        auto const k_S = MPL::formEigenTensor<DisplacementDim>(
            medium.property(MPL::PropertyType::permeability)
                .value(vars, pos, t, dt));

        auto const s_L = medium.property(MPL::PropertyType::saturation)
                             .template value<double>(vars, pos, t, dt);

        vars[static_cast<int>(MPL::Variable::liquid_saturation)] = s_L;

        auto const s_G = 1. - s_L;

        auto const dsLdPc =
            medium.property(MPL::PropertyType::saturation)
                .template dValue<double>(
                    vars, MPL::Variable::capillary_pressure, pos, t, dt);

        auto const alpha_B =
            medium.property(MPL::PropertyType::biot_coefficient)
                .template value<double>(vars, pos, t, dt);

        auto const beta_pSR = (1. - alpha_B) * beta_pS;

        auto const k_rel =
            medium.property(MPL::PropertyType::relative_permeability)
                .template value<Eigen::Vector2d>(vars, pos, t, dt);

        auto const k_rel_L = k_rel[0];
        auto const k_rel_G = k_rel[1];

        auto const& b = _process_data.specific_body_force;

        auto const k_over_mu_G = k_S * k_rel_G / mu_GR;
        auto const k_over_mu_L = k_S * k_rel_L / mu_LR;

        auto const phi = medium.property(MPL::PropertyType::porosity)
                             .template value<double>(vars, pos, t, dt);

        auto const phi_G = s_G * phi;
        auto const phi_L = s_L * phi;
        auto const phi_S = 1. - phi;

        // TODO (Grunwald) replace effective thermal conductivity by a more
        // sophisticated law, maybe even allow the law to be chosen in the
        // project file as medium property
        auto const lambda = MPL::formEigenTensor<DisplacementDim>(
            phi_S * lambda_SR + phi_L * lambda_LR + phi_G * lambda_GR);

        // TODO (Wenqing) : Change dT to time step wise increment
        double const delta_T(T_int_pt - T0);
        double const thermal_strain = beta_TSR * delta_T;
        auto const rho_SR = rho_ref_SR * (1 - 3 * thermal_strain);

        auto const rho = rho_GR + rho_LR + rho_SR;
        // TODO: change back to rho_SR
        auto const rho_c_p = phi_G * rho_GR * c_pG + phi_L * rho_LR * c_pL +
                             phi_S * rho_ref_SR * c_pS;

        // update secondary variables. TODO: Refactoring potential!!

        _liquid_pressure[ip] = pLR_int_pt;
        _liquid_density[ip] = rho_LR;
        _gas_density[ip] = rho_GR;
        _porosity[ip] = phi;
        _saturation[ip] = s_L;

        GlobalDimVectorType const w_GS =
            k_over_mu_G * rho_GR * b - k_over_mu_G * gradNp * pGR;

        GlobalDimVectorType const w_LS = k_over_mu_L * gradNp * pCap +
                                         k_over_mu_L * rho_GR * b -
                                         k_over_mu_L * gradNp * pGR;

#ifdef MATERIAL_PROPERTIES
        std::cout << "######################################################\n";
        std::cout << "#    Material properties:\n";
        std::cout << "#----------------------------------------------------#\n";
        std::cout << "#         rho_GR:  " << rho_GR << "\n";
        std::cout << "#         rho_LR:  " << rho_LR << "\n";
        std::cout << "#         rho_SR:  " << rho_SR << "\n";
        std::cout << "#            rho:  " << rho << "\n";
        std::cout << "#----------------------------------------------------#\n";
        std::cout << "#          c_p_G:  " << c_pG << "\n";
        std::cout << "#          c_p_L:  " << c_pL << "\n";
        std::cout << "#          c_p_S:  " << c_pS << "\n";
        std::cout << "#        rho_c_p:  " << rho_c_p << "\n";
        std::cout << "#----------------------------------------------------#\n";
        std::cout << "#      beta_p_GR:  " << beta_p_GR << "\n";
        std::cout << "#      beta_p_LR:  " << beta_p_LR << "\n";
        std::cout << "#      beta_p_SR:  " << beta_p_SR << "\n";
        std::cout << "#      beta_T_GR:  " << beta_T_GR << "\n";
        std::cout << "#      beta_T_LR:  " << beta_T_LR << "\n";
        std::cout << "#      beta_T_SR:  " << beta_T_SR << "\n";
        std::cout << "#----------------------------------------------------#\n";
        std::cout << "#            k_S:\n" << k_S << "\n";
        std::cout << "#----------------------------------------------------#\n";
        std::cout << "#        k_rel_L:  " << k_rel_L << "\n";
        std::cout << "#        k_rel_G:  " << k_rel_G << "\n";
        std::cout << "#          mu_GR:  " << mu_GR << "\n";
        std::cout << "#          mu_LR:  " << mu_LR << "\n";
        std::cout << "#    k_over_mu_G:  " << k_over_mu_G << "\n";
        std::cout << "#    k_over_mu_L:  " << k_over_mu_L << "\n";
        std::cout << "#----------------------------------------------------#\n";
        std::cout << "#            s_L:  " << s_L << "\n";
        std::cout << "#            s_G:  " << s_G << "\n";
        std::cout << "#         dsLdPc:  " << dsLdPc << "\n";
        std::cout << "#----------------------------------------------------#\n";
        std::cout << "#        alpha_B:  " << alpha_B << "\n";
        std::cout << "#----------------------------------------------------#\n";
        std::cout << "#         lambda:\n" << lambda << "\n";
        std::cout << "#----------------------------------------------------#\n";
        std::cout << "#            phi:  " << phi << "\n";
        std::cout << "#          phi_G:  " << phi_G << "\n";
        std::cout << "#          phi_L:  " << phi_L << "\n";
        std::cout << "#          phi_S:  " << phi_S << "\n";
        std::cout << "######################################################\n";
        std::cout << "#  Darcy-Velocities: "
                  << "\n";
        std::cout << "#----------------------------------------------------#\n";
        std::cout << "#           w_LS:\n" << w_LS << "\n";
        std::cout << "#           w_GS:\n" << w_GS << "\n";
        std::cout << "######################################################\n";
#endif

        // coefficient matrices
        //  - gas pressure equation
        MGpG.noalias() +=
            (NpT * s_G * (phi * beta_pGR + (alpha_B - phi) * beta_pSR) * Np) *
            w;
        MGpC.noalias() -=
            (NpT *
             (s_G * (alpha_B - phi) * beta_pSR * (s_L + pCap_int_pt * dsLdPc) +
              phi * dsLdPc) *
             Np) *
            w;
        MGT.noalias() -=
            (NpT * s_G * (phi * beta_TGR + (alpha_B - phi) * beta_TSR) * Np) *
            w;
        MGu.noalias() += NpT * mT * Bu * s_G * alpha_B * w;
        KGpG.noalias() += (gradNpT * k_over_mu_G * gradNp) * w;
        fG.noalias() += (gradNpT * rho_GR * k_over_mu_G * b) * w;
        //  - liquid pressure equation
        MLpG.noalias() +=
            (NpT * s_L * (phi * beta_pLR + (alpha_B - phi) * beta_pSR) * Np) *
            w;
        MLpC.noalias() -=
            (NpT *
             (s_L * (alpha_B - phi) * beta_pSR * (s_L + pCap_int_pt * dsLdPc) +
              phi_L * beta_pLR - phi * dsLdPc) *
             Np) *
            w;
        MLT.noalias() -=
            (NpT * s_L * (phi * beta_TLR + (alpha_B - phi) * beta_TSR) * Np) *
            w;
        MLu.noalias() += (NpT * s_L * alpha_B * mT * Bu) * w;
        KLpG.noalias() += (gradNpT * k_over_mu_L * gradNp) * w;
        KLpC.noalias() -= (gradNpT * k_over_mu_L * gradNp) * w;
        fL.noalias() += (gradNpT * rho_LR * k_over_mu_L * b) * w;

        //  - temperature equation
        MTpG.noalias() -=
            (NTT * (phi_G * beta_TGR + phi_L * beta_TLR + phi_S * beta_TSR) *
             T_int_pt * NT) *
            w;

        MTpC.noalias() +=
            (NTT *
             ((phi_L * beta_TLR +
               phi_S * beta_TSR * (s_L + pCap_int_pt * dsLdPc) * T_int_pt) +
              phi * pCap_int_pt * dsLdPc) *
             NT) *
            w;

        MTT.noalias() += (NTT * rho_c_p * NT) * w;
        KTpG.noalias() -= (NTT * beta_TGR * w_GS.transpose() * gradNT +
                           NTT * beta_TLR * w_LS.transpose() * gradNT) *
                          w;
        KTpC.noalias() += (NTT * beta_TLR * w_LS.transpose() * gradNT) * w;

        // ATT
        KTT.noalias() += (NTT * w_LS.transpose() * gradNT * rho_LR * c_pL +
                          NTT * w_GS.transpose() * gradNT * rho_GR * c_pG) *
                         w;
        // LTT
        KTT.noalias() += (gradNTT * lambda * gradNT) * w;

        //  - displacement equation
        KUpG.noalias() -= (BuT * alpha_B * m * Np) * w;
        KUpC.noalias() += (BuT * alpha_B * s_L * m * Np) * w;

        //
        // displacement equation, displacement part
        //
        eps.noalias() = Bu * u;
        ip_data.updateConstitutiveRelationThermal(
            t, pos, dt, u, _process_data.reference_temperature(t, pos)[0],
            thermal_strain);
        fU.noalias() -= (BuT * sigma_eff - Nu_op.transpose() * rho * b) * w;

#ifdef DEBUG_OUTPUT
        std::cout << "------------------------------------------------------\n";
        std::cout << " MGpG:\n" << MGpG << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " MGpC:\n" << MGpC << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " MGT:\n" << MGT << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " KGpG (LGpG):\n" << KGpG << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " fG:\n" << fG << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " MLpG:\n" << MLpG << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " MLpC:\n" << MLpC << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " MLT:\n" << MLT << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " MLu:\n" << MLu << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " KLpG (LLpG):\n" << KLpG << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " KLpC (LLpC):\n" << KLpC << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " fL:\n" << fL << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " MTpG:\n" << MTpG << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " MTpC:\n" << MTpC << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " MTT:\n" << MTT << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " KTpG (ATpG):\n" << KTpG << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " KTpC (ATpC):\n" << KTpC << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " KTT (ATT):\n" << KTT << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " KTT (ATT+LTT):\n" << KTT << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " KUpG:\n" << KUpG << "\n";
        std::cout << " KUpC:\n" << KUpC << "\n";
        std::cout << " KUU:\n" << KUU << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " fU:\n" << fU << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " eps:\n" << eps << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " C:\n" << C << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " thermal_strain:\n" << thermal_strain << "\n";
        std::cout << "------------------------------------------------------\n";
        std::cout << " sigma_eff:\n" << sigma_eff << "\n";
        std::cout << "------------------------------------------------------\n";
        OGS_FATAL("Stop intended.");
#endif
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
    OGS_FATAL("TH2M-assembleWithJaCobian is not implemented!");
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

    auto const& medium = *_process_data.media_map->getMedium(_element.getID());
    auto const& gas_phase = medium.phase("Gas");

    MPL::VariableArray vars;

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);

        auto const& N_p = _ip_data[ip].N_p;

        vars[static_cast<int>(MPL::Variable::temperature)] =
            N_p.dot(T);  // N_p = N_T
        vars[static_cast<int>(MPL::Variable::phase_pressure)] = N_p.dot(pGR);

        // TODO (naumov) Temporary value not used by current material models.
        // Need extension of secondary variables interface.
        double const dt = std::numeric_limits<double>::quiet_NaN();

        auto const mu_GR = gas_phase.property(MPL::PropertyType::viscosity)
                               .template value<double>(vars, pos, t, dt);

        GlobalDimMatrixType k_S = MPL::formEigenTensor<DisplacementDim>(
            medium.property(MPL::PropertyType::permeability)
                .value(vars, pos, t, dt));

        auto const k_rel =
            medium.property(MPL::PropertyType::relative_permeability)
                .template value<Eigen::Vector2d>(vars, pos, t, dt)[1];

        auto const k_over_mu = k_S * k_rel / mu_GR;

        auto const rho_GR = gas_phase.property(MPL::PropertyType::density)
                                .template value<double>(vars, pos, t, dt);
        auto const& b = _process_data.specific_body_force;

        // Compute the velocity
        auto const& dNdx_p = _ip_data[ip].dNdx_p;
        cache_matrix.col(ip).noalias() =
            -k_over_mu * dNdx_p * pGR + k_over_mu * rho_GR * b;
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

    auto const& medium = *_process_data.media_map->getMedium(_element.getID());
    auto const& liquid_phase = medium.phase("AqueousLiquid");

    MPL::VariableArray vars;

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);

        auto const& N_p = _ip_data[ip].N_p;

        vars[static_cast<int>(MPL::Variable::temperature)] = N_p.dot(T);
        vars[static_cast<int>(MPL::Variable::phase_pressure)] = N_p.dot(pLR);

        // TODO (naumov) Temporary value not used by current material models.
        // Need extension of secondary variables interface.
        double const dt = std::numeric_limits<double>::quiet_NaN();

        auto const mu_LR = liquid_phase.property(MPL::PropertyType::viscosity)
                               .template value<double>(vars, pos, t, dt);
        GlobalDimMatrixType k_S = MPL::formEigenTensor<DisplacementDim>(
            medium.property(MPL::PropertyType::permeability)
                .value(vars, pos, t, dt));

        auto const k_rel_L =
            medium.property(MPL::PropertyType::relative_permeability)
                .template value<Eigen::Vector2d>(vars, pos, t, dt)[0];

        auto const k_over_mu = k_S * k_rel_L / mu_LR;

        auto const rho_LR = liquid_phase.property(MPL::PropertyType::density)
                                .template value<double>(vars, pos, t, dt);
        auto const& b = _process_data.specific_body_force;

        // Compute the velocity
        auto const& dNdx_p = _ip_data[ip].dNdx_p;
        cache_matrix.col(ip).noalias() =
            -k_over_mu * dNdx_p * pLR + k_over_mu * rho_LR * b;
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

    auto const& medium = *_process_data.media_map->getMedium(_element.getID());
    auto const& liquid_phase = medium.phase("AqueousLiquid");
    auto const& gas_phase = medium.phase("Gas");
    auto const& solid_phase = medium.phase("Solid");

    MPL::VariableArray vars;

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);
        auto const& Nu = _ip_data[ip].N_u;
        auto const& Np = _ip_data[ip].N_p;
        auto const& NT = Np;
        auto const& dNdx_u = _ip_data[ip].dNdx_u;

        auto const x_coord =
            interpolateXCoordinate<ShapeFunctionDisplacement,
                                   ShapeMatricesTypeDisplacement>(_element, Nu);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, Nu, x_coord, _is_axially_symmetric);

        double const T0 = _process_data.reference_temperature(t, pos)[0];

        auto const T_int_pt = NT.dot(T);
        auto const pGR_int_pt = Np.dot(pGR);
        auto const pCap_int_pt = Np.dot(pCap);
        auto const pLR_int_pt = pGR_int_pt - pCap_int_pt;

        vars[static_cast<int>(MPL::Variable::temperature)] = T_int_pt;
        vars[static_cast<int>(MPL::Variable::phase_pressure)] = pGR_int_pt;
        vars[static_cast<int>(MPL::Variable::capillary_pressure)] = pCap_int_pt;

        auto const solid_linear_thermal_expansion_coefficient =
            solid_phase.property(MPL::PropertyType::thermal_expansivity)
                .template value<double>(vars, pos, t, dt);

        double const delta_T(T_int_pt - T0);
        double const thermal_strain =
            solid_linear_thermal_expansion_coefficient * delta_T;

        auto& eps = _ip_data[ip].eps;
        eps.noalias() = B * u;

        _ip_data[ip].updateConstitutiveRelationThermal(
            t, pos, dt, u, _process_data.reference_temperature(t, pos)[0],
            thermal_strain);

        auto const rho_LR = liquid_phase.property(MPL::PropertyType::density)
                                .template value<double>(vars, pos, t, dt);

        auto const phi = medium.property(MPL::PropertyType::porosity)
                             .template value<double>(vars, pos, t, dt);

        auto const rho_GR = gas_phase.property(MPL::PropertyType::density)
                                .template value<double>(vars, pos, t, dt);

        auto const s_L = medium.property(MPL::PropertyType::saturation)
                             .template value<double>(vars, pos, t, dt);

        _liquid_pressure[ip] = pLR_int_pt;
        _liquid_density[ip] = rho_LR;
        _gas_density[ip] = rho_GR;
        _porosity[ip] = phi;
        _saturation[ip] = s_L;
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