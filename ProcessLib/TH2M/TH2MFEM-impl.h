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

void VLE(
    const bool phase_transition, const double T, const double pGR,
    const double pCap, double& rho_GR /*gas phase density*/,
    double& rho_C_GR /*C constituent partial density in gas phase*/,
    double& rho_W_GR /*W constituent partial density in gas phase*/,
    double& rho_LR /*liquid phase density*/,
    double& rho_C_LR /*C constituent partial density in liquid phase*/,
    double& rho_W_LR /*W constituent partial density in liquid phase*/,
    double& xm_C_G /*C constituent gas phase mass fraction*/,
    double& xm_W_G /*W constituent gas phase mass fraction*/,
    double& xm_C_L /*C constituent liquid phase mass fraction*/,
    double& xm_W_L /*W constituent liquid phase mass fraction*/,
    double& d_xm_C_G_d_pGR /*C mass fraction derivative w.r.t. pGR in G phase*/,
    double& d_xm_W_G_d_pGR /*W mass fraction derivative w.r.t. pGR in G phase*/,
    double& d_xm_C_L_d_pLR /*C mass fraction derivative w.r.t. pLR in L phase*/,
    double& d_xm_W_L_d_pLR /*W mass fraction derivative w.r.t. pLR in L phase*/,
    double& d_xm_C_G_d_T /*C mass fraction derivative w.r.t. T in G phase*/,
    double& d_xm_W_G_d_T /*W mass fraction derivative w.r.t. T in G phase*/,
    double& d_xm_C_L_d_T /*C mass fraction derivative w.r.t. T in L phase*/,
    double& d_xm_W_L_d_T /*W mass fraction derivative w.r.t. T in L phase*/
)
{
    // Hard-coded VLE-properties for the moment. Will be changed later and
    // incorporated into MPL-structure. Henry constants for carbon dioxide
    // and hydrogen
    const double H_theta_CO2 = 3.3e-4;              // mol/m^3/Pa
    const double dH_solve_CO2 = -19954.7102835678;  // J/mol

    const double H_theta_H2 = 7.8e-6;             // mol/m^3/Pa
    const double dH_solve_H2 = -4406.6651876212;  // J/mol
    const double R = 8.31446261815324;            // J/mol/K
    const double T_theta = 298.15;                // K

    const double dH_by_R = dH_solve_CO2 / R;

    // CO2-Henry coefficient depending on T
    const double H_C =
        phase_transition
            ? H_theta_CO2 * std::exp((-1.) * dH_by_R * (1. / T - 1. / T_theta))
            : 0.;
    const double dHc_dT = 1. / (T * T) * dH_by_R * H_C;

    // coefficients for Antoine equation
    const double A = 8.07131;
    const double B = 1730.63;
    const double C = 233.426;

    const auto pLR = pGR - pCap;

    // saturation vapour pressure (Antoine equation)
    const double theta = T - 273.15;
    const double p_vap =
        (phase_transition ? std::pow(10., A - B / (C + theta)) : 0.) *
        133.322;  // converted from Torr to Pascal
    const double ln10 = 2.30258509299;
    const double dp_vap_dT = B * ln10 / ((C + theta) * (C + theta)) * p_vap;

    // partial pressure of constituents in gas phase
    const double p_W_GR = p_vap;
    const double p_C_GR = pGR - p_W_GR;  // Daltons law

    // molar fraction of gas phase constituents
    const double xn_W_G = p_W_GR / pGR;
    const double xn_CO2_G = 1. - xn_W_G;
    const double xn_H2_G = 1. - xn_W_G;
    const double xn_C_G = xn_CO2_G;

    const double d_xn_C_G_d_pGR = xn_W_G / pGR;
    const double d_xn_C_G_d_T = -dp_vap_dT / pGR;

    // molar masses and average molar mass
    const double M_W = .01801528;    // kg/mol
    const double M_CO2 = .0440095;   // kg/mol
    const double M_H2 = 0.00201948;  // kg/mol

    const double M_C = M_CO2;
    const double M_G = xn_C_G * M_C + xn_W_G * M_W;

    // density of the gas phase (ideal gas law)
    rho_GR = pGR * M_G / R / T;
    const double d_rho_GR_d_pGR = M_C / R / T;
    const double d_rho_GR_d_T = dp_vap_dT * (M_W - M_C) / R / T - rho_GR / T;

    // mass fractions of gas phase constituents
    xm_W_G = xn_W_G * M_W / M_G;
    xm_C_G = 1. - xm_W_G;

    // partial densities of gas phase constituents
    rho_W_GR = xm_W_G * rho_GR;
    const double d_rho_W_GR_d_pGR = 0.;
    const double d_rho_W_GR_d_T = rho_W_GR * (1. / p_W_GR * dp_vap_dT - 1 / T);
    rho_C_GR = xm_C_G * rho_GR;
    const double d_rho_C_GR_d_pGR = d_rho_GR_d_pGR;
    const double d_rho_C_GR_d_T = d_rho_GR_d_T - d_rho_W_GR_d_T;

    // concentration of dissolved gas in water
    const double c_C_L = p_C_GR * H_C;
    const double d_c_C_L_d_T = dHc_dT * xn_C_G * pGR - H_C * dp_vap_dT;

    // Liquid-EOS, slopes of hyperplane
    const double beta_pLR = 1.0e-8;
    const double beta_TLR = -1.0e-4;
    const double beta_CLR = 1.0e-5;

    const double T_ref = 290.;
    const double p_ref = 1.0e5;
    const double rho_LR_ref = 999.0;

    // water component partial density liquid phase
    rho_W_LR =
        rho_LR_ref * (1. + beta_pLR * (pLR - p_ref) + beta_TLR * (T - T_ref));
    const double d_rho_W_LR_d_pLR = rho_LR_ref * beta_pLR;
    const double d_rho_W_LR_d_T = rho_LR_ref * beta_TLR;

    // liquid phase density
    rho_LR = rho_LR_ref * (1. + beta_pLR * (pLR - p_ref) +
                           beta_TLR * (T - T_ref) + beta_CLR * c_C_L);

    const double d_rho_LR_d_pLR = rho_LR_ref * (beta_pLR + beta_CLR * H_C);
    const double d_rho_LR_d_T =
        rho_LR_ref * (beta_TLR + beta_CLR * d_c_C_L_d_T);

    // partial density of dissolved gas component in water
    rho_C_LR = rho_LR - rho_W_LR;
    const double d_rho_C_LR_d_pLR =
        rho_LR_ref * beta_CLR * H_C * (d_xn_C_G_d_pGR * pGR + xn_C_G);
    const double d_rho_C_LR_d_T =
        rho_LR_ref * beta_CLR *
        (dHc_dT * xn_C_G * pGR + H_C * d_xn_C_G_d_T * pGR);

    // mass fractions of liquid phase constituents
    xm_W_L = rho_W_LR / rho_LR;
    xm_C_L = 1. - xm_W_L;

    // constituent compressibility of the gas phase
    const double beta_C_pGR = (p_C_GR != 0) ? 1. / p_C_GR : 0.;
    const double beta_W_pGR = (p_W_GR != 0) ? 1. / p_W_GR : 0.;

    // constituent thermal expansivity of the gas phase
    const double beta_C_TGR = -beta_C_pGR * dp_vap_dT - 1 / T;
    const double beta_W_TGR = beta_W_pGR * dp_vap_dT - 1 / T;

    // constituent compressibility of the liquid phase
    const double beta_C_pLR = beta_CLR * H_C;
    const double beta_C_TLR = beta_CLR * (dHc_dT * p_C_GR - dp_vap_dT * H_C);

    // constituent thermal expansivity of the liquid phase
    const double beta_W_pLR = beta_pLR;
    const double beta_W_TLR = beta_TLR;

    const double beta_pGR = 1. / pGR;
    const double beta_TGR = -1. / T;

    // dxm_C_G_dpGR = xm_C_G * (beta_C_pGR - beta_pGR);
    // dxm_W_G_dpGR = xm_W_G * (beta_W_pGR - beta_pGR);
    // dxm_C_L_dpLR = xm_C_L * (beta_C_pLR - beta_pLR);
    // dxm_W_L_dpLR = xm_W_L * (beta_W_pLR - beta_pLR);

    // dxm_C_G_dT = xm_C_G * (beta_C_TGR - beta_TGR);
    // dxm_C_L_dT = xm_C_L * (beta_C_TLR - beta_TLR);
    // dxm_W_G_dT = xm_W_G * (beta_W_TGR - beta_TGR);
    // dxm_W_L_dT = xm_W_L * (beta_W_TLR - beta_TLR);

    if (phase_transition)
    {
    d_xm_C_G_d_pGR = 1. / rho_GR * (d_rho_C_GR_d_pGR - xm_C_G * d_rho_GR_d_pGR);
    d_xm_W_G_d_pGR = 1. / rho_GR * (d_rho_W_GR_d_pGR - xm_W_G * d_rho_GR_d_pGR);
    d_xm_C_L_d_pLR = 1. / rho_LR * (d_rho_C_LR_d_pLR - xm_C_L * d_rho_LR_d_pLR);
    d_xm_W_L_d_pLR = 1. / rho_LR * (d_rho_W_LR_d_pLR - xm_W_L * d_rho_LR_d_pLR);

    d_xm_C_G_d_T = 1. / rho_GR * (d_rho_C_GR_d_T - xm_C_G * d_rho_GR_d_T);
    d_xm_W_G_d_T = 1. / rho_GR * (d_rho_W_GR_d_T - xm_W_G * d_rho_GR_d_T);
    d_xm_C_L_d_T = 1. / rho_LR * (d_rho_C_LR_d_T - xm_C_L * d_rho_LR_d_T);
    d_xm_W_L_d_T = 1. / rho_LR * (d_rho_W_LR_d_T - xm_W_L * d_rho_LR_d_T);
    }
    else
    {
    d_xm_C_G_d_pGR = 0.;
    d_xm_W_G_d_pGR = 0.;
    d_xm_C_L_d_pLR = 0.;
    d_xm_W_L_d_pLR = 0.;

    d_xm_C_G_d_T = 0.;
    d_xm_W_G_d_T = 0.;
    d_xm_C_L_d_T = 0.;
    d_xm_W_L_d_T = 0.;
    }
    
    // std::cout << std::setprecision(16);

    // std::cout << T << " " << pGR << " " << pCap << " " << pLR << " " <<
    // rho_GR
    //           << " " << rho_C_GR << " " << rho_W_GR << " " << xm_C_G << " "
    //           << d_xm_C_G_d_pGR << " " << d_xm_C_G_d_T << " " << xm_W_G << "
    //           "
    //           << d_xm_W_G_d_pGR << " " << d_xm_W_G_d_T << " " << rho_LR << "
    //           "
    //           << rho_C_LR << " " << rho_W_LR << " " << xm_C_L << " "
    //           << d_xm_C_L_d_pLR << " " << d_xm_C_L_d_T << " " << xm_W_L << "
    //           "
    //           << d_xm_W_L_d_pLR << " " << d_xm_W_L_d_T << "\n";

#define nVLE_DEBUG_OUTPUT
#ifdef VLE_DEBUG_OUTPUT
    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "######################################################\n";
    std::cout << "#    VLE-Model:                                      #\n";
    std::cout << "#----------------------------------------------------#\n";
    std::cout << "# Phase transition is " << (phase_transition ? "ON" : "OFF")
              << "\n";
    std::cout << "#----------------------------------------------------#\n";
    std::cout << "# gas phase pressure: " << pGR << "\n";
    std::cout << "# liquid phase pressure: " << pLR << "\n";
    std::cout << "# capillary pressure: " << pCap << "\n";
    std::cout << "# temperature: " << T << "\n";
    std::cout << "#----------------------------------------------------#\n";
    std::cout << "# dH_by_R: " << dH_by_R << "\n";
    std::cout << "# Henry-coefficient: " << H_C << "\n";
    std::cout << "# derivative dHc_dT: " << dHc_dT << "\n";
    std::cout << "# water vapour pressure: " << p_vap << "\n";
    std::cout << "# derivative dp_vap_dT: " << dp_vap_dT << "\n";
    std::cout << "#----------------------------------------------------#\n";
    std::cout << "# partial pressures:\n";
    std::cout << "# W: " << p_W_GR << "\n";
    std::cout << "# C: " << p_C_GR << "\n";
    std::cout << "#----------------------------------------------------#\n";
    std::cout << "gas phase density: " << rho_GR << "\n ";
    std::cout << "partial densities gas phase:\n";
    std::cout << "W: " << rho_W_GR << "\n";
    std::cout << "C: " << rho_C_GR << "\n";
    std::cout << "molar fractions of gas phase:\n";
    std::cout << "W: " << xn_W_G << "\n";
    std::cout << "C: " << xn_C_G << "\n";
    std::cout << "mass fractions gas phase:\n";
    std::cout << "W: " << xm_W_G << "\n";
    std::cout << "C: " << xm_C_G << "\n";
    std::cout << "average molar mass: " << M_G << "\n ";
    std::cout << "#----------------------------------------------------#\n";
    std::cout << "liquid phase density: " << rho_LR << "\n";
    std::cout << "partial densities liquid phase:\n";
    std::cout << "W: " << rho_W_LR << "\n";
    std::cout << "C: " << rho_C_LR << "\n";
    std::cout << "dissolved gas concentration: " << c_C_L << "\n";
    std::cout << "mass fractions liquid phase:\n";
    std::cout << "W: " << xm_W_L << "\n";
    std::cout << "C: " << xm_C_L << "\n";
    std::cout << "#----------------------------------------------------#\n";
    std::cout << "Gas phase compressibilities:\n";
    std::cout << "beta_pGR: " << beta_pGR << "\n";
    std::cout << "beta_TGR: " << beta_TGR << "\n";
    std::cout << "Gas phase component compressibilities:\n";
    std::cout << "W: beta_W_pGR: " << beta_W_pGR << "\n";
    std::cout << "W: beta_W_TGR: " << beta_W_TGR << "\n";
    std::cout << "C: beta_C_pGR: " << beta_C_pGR << "\n";
    std::cout << "C: beta_C_TGR: " << beta_C_TGR << "\n";
    std::cout << "#----------------------------------------------------#\n";
    std::cout << "Liquid phase compressibilities:\n";
    std::cout << "beta_pLR: " << beta_pLR << "\n";
    std::cout << "beta_TLR: " << beta_TLR << "\n";
    std::cout << "beta_CLR: " << beta_CLR << "\n";
    std::cout << "Liquid phase component compressibilities:\n";
    std::cout << "W: beta_W_pLR: " << beta_W_pLR << "\n";
    std::cout << "W: beta_W_TLR: " << beta_W_TLR << "\n";
    std::cout << "C: beta_C_pLR: " << beta_C_pLR << "\n";
    std::cout << "C: beta_C_TLR: " << beta_C_TLR << "\n";
    std::cout << "#----------------------------------------------------#\n";
    std::cout << "derivatives of mass fraction:\n";
    std::cout << "#----------------------------------------------------#\n";
    std::cout << "gas phase:\n";
    std::cout << "W: dxm_W_G/dpGR: " << dxm_W_G_dpGR << "\n";
    std::cout << "W: dxm_W_G/dT: " << dxm_W_G_dT << "\n";
    std::cout << "C: dxm_C_G/dpGR: " << dxm_C_G_dpGR << "\n";
    std::cout << "C: dxm_C_G/dT: " << dxm_C_G_dT << "\n";
    std::cout << "#----------------------------------------------------#\n";
    std::cout << "liquid phase:\n";
    std::cout << "W: dxm_W_L/dpLR: " << dxm_W_L_dpLR << "\n";
    std::cout << "W: dxm_W_L/dT: " << dxm_W_L_dT << "\n";
    std::cout << "C: dxm_C_L/dpLR: " << dxm_C_L_dpLR << "\n";
    std::cout << "C: dxm_C_L/dT: " << dxm_C_L_dT << "\n";
    std::cout << "#----------------------------------------------------#\n";
    std::cout << "######################################################\n";
    OGS_FATAL(".");
#endif
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                        IntegrationMethod, DisplacementDim>::
    getConstitutiveVariables(std::vector<double> const& local_x, double const t,
                             double const dt)
{
    auto const matrix_size = gas_pressure_size + capillary_pressure_size +
                             temperature_size + displacement_size;

    assert(local_x.size() == matrix_size);

    auto gas_pressure =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            gas_pressure_size> const>(local_x.data() + gas_pressure_index,
                                      gas_pressure_size);

    auto capillary_pressure =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            capillary_pressure_size> const>(
            local_x.data() + capillary_pressure_index, capillary_pressure_size);

    auto temperature =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            temperature_size> const>(local_x.data() + temperature_index,
                                     temperature_size);

    auto displacement =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + displacement_index,
                                      displacement_size);

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    auto const& medium = *_process_data.media_map->getMedium(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);
        auto& ip_data = _ip_data[ip];

        auto const& Np = ip_data.N_p;
        auto const& NT = Np;
        auto const& Nu = ip_data.N_u;
        auto const& gradNu = ip_data.dNdx_u;

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

        auto const T = NT.dot(temperature);
        auto const pGR = Np.dot(gas_pressure);
        auto const pCap = Np.dot(capillary_pressure);
        auto const pLR = pGR - pCap;

        double div_u = mT * Bu * displacement;

        MPL::VariableArray vars;
        vars[static_cast<int>(MPL::Variable::temperature)] = T;
        vars[static_cast<int>(MPL::Variable::phase_pressure)] = pGR;
        vars[static_cast<int>(MPL::Variable::capillary_pressure)] = pCap;
        vars[static_cast<int>(MPL::Variable::liquid_phase_pressure)] = pLR;

        ip_data.saturation =
            medium.property(MPL::PropertyType::saturation)
                .template value<double>(
                    vars, pos, t, std::numeric_limits<double>::quiet_NaN());

        const double beta_pS = 1e-10;

        auto const phi = medium.property(MPL::PropertyType::porosity)
                             .template value<double>(vars, pos, t, dt);

        auto const s_L = ip_data.saturation;
        auto const phi_G = (1. - s_L) * phi;
        auto const phi_L = s_L * phi;
        auto const phi_S = 1. - phi;

        const double p_FR =
            (1. - ip_data.saturation) * pGR + ip_data.saturation * pLR;

        double const T0 = _process_data.reference_temperature(t, pos)[0];
        double const delta_T(T - T0);
        const double beta_T_SR = -1e-7;
        const double rho_ref_SR = 2500.;

        double const thermal_strain = beta_T_SR * delta_T;

        const double p_SR =
            p_FR - 1. / (beta_pS * phi_S) * (div_u - thermal_strain);
        auto const rho_SR = rho_ref_SR * (1. - 3. * thermal_strain);
        ip_data.phi_S_p_SR = phi_S * p_SR;

        double rho_GR;
        double rho_LR;
        double rho_C_GR;
        double rho_W_GR;
        double rho_C_LR;
        double rho_W_LR;

        double xm_C_G;
        double xm_W_G;
        double xm_C_L;
        double xm_W_L;

        double dummy;
        bool phase_transition = _process_data.phase_transition;

        VLE(phase_transition, T, pGR, pCap, rho_GR, rho_C_GR, rho_W_GR, rho_LR,
            rho_C_LR, rho_W_LR, xm_C_G, xm_W_G, xm_C_L, xm_W_L, dummy, dummy,
            dummy, dummy, dummy, dummy, dummy, dummy);

        ip_data.rho_C_GR = rho_C_GR;
        ip_data.rho_W_GR = rho_W_GR;
        ip_data.rho_C_LR = rho_C_LR;
        ip_data.rho_W_LR = rho_W_LR;
        
        const double c_p_S = 1000.;
        const auto c_p_C_G = 1000.;  // J/kg/K
        const auto c_p_W_G = 1914.;  // J/kg/K
        const auto h_C_G = c_p_C_G * T;
        const auto h_W_G = c_p_W_G * T;
        const auto h_W_L = 112654.;  // J/kg
        const auto h_C_L = 648584.;  // J/kg

        // phase specific enthalpies
        const auto h_G = xm_C_G * h_C_G + xm_W_G * h_W_G;
        const auto h_L = xm_W_L * h_W_L + xm_C_L * h_C_L;
        const auto h_S = c_p_S * T;

        // phase enthalpies
        ip_data.rho_G_h_G = phi_G * rho_GR * h_G;
        ip_data.rho_L_h_L = phi_L * rho_LR * h_L;
        ip_data.rho_S_h_S = phi_S * rho_SR * h_S;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                        IntegrationMethod, DisplacementDim>::
    setInitialConditionsConcrete(std::vector<double> const& local_x,
                                 double const t)
{
    auto const matrix_size = gas_pressure_size + capillary_pressure_size +
                             temperature_size + displacement_size;

    assert(local_x.size() == matrix_size);

    auto gas_pressure =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            gas_pressure_size> const>(local_x.data() + gas_pressure_index,
                                      gas_pressure_size);

    auto capillary_pressure =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            capillary_pressure_size> const>(
            local_x.data() + capillary_pressure_index, capillary_pressure_size);

    auto temperature =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            temperature_size> const>(local_x.data() + temperature_index,
                                     temperature_size);

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    getConstitutiveVariables(local_x, t,
                             std::numeric_limits<double>::quiet_NaN());

    auto const& medium = *_process_data.media_map->getMedium(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);
        auto& ip_data = _ip_data[ip];

        auto const& Np = ip_data.N_p;
        auto const& NT = Np;

        auto const T = NT.dot(temperature);
        auto const pGR = Np.dot(gas_pressure);
        auto const pCap = Np.dot(capillary_pressure);
        auto const pLR = pGR - pCap;

        MPL::VariableArray vars;
        vars[static_cast<int>(MPL::Variable::temperature)] = T;
        vars[static_cast<int>(MPL::Variable::phase_pressure)] = pGR;
        vars[static_cast<int>(MPL::Variable::capillary_pressure)] = pCap;
        vars[static_cast<int>(MPL::Variable::liquid_phase_pressure)] = pLR;

        ip_data.pushBackState();
        // ip_data.saturation_prev =
        //     medium.property(MPL::PropertyType::saturation)
        //         .template value<double>(
        //             vars, pos, t, std::numeric_limits<double>::quiet_NaN());

        // auto const phi =
        //     medium.property(MPL::PropertyType::porosity)
        //         .template value<double>(
        //             vars, pos, t, std::numeric_limits<double>::quiet_NaN());

        // double const T0 = _process_data.reference_temperature(t, pos)[0];
        // double const delta_T(T - T0);
        // const double beta_T_SR = -1e-7;
        // const double rho_ref_SR = 2500.;

        // double const thermal_strain = beta_T_SR * delta_T;
        // auto const rho_SR = rho_ref_SR * (1. - 3. * thermal_strain);

        // const double s_L = ip_data.saturation_prev;
        // const double s_G = (1. - s_L);
        // const double phi_L = s_L * phi;
        // const double phi_G = s_G * phi;
        // const double phi_S = 1. - phi;

        // double rho_GR;
        // double rho_LR;
        // double rho_C_GR;
        // double rho_W_GR;
        // double rho_C_LR;
        // double rho_W_LR;

        // double xm_C_G;
        // double xm_W_G;
        // double xm_C_L;
        // double xm_W_L;

        // double dummy;
        // bool phase_transition = _process_data.phase_transition;

        // VLE(phase_transition, T, pGR, pCap, rho_GR, rho_C_GR, rho_W_GR,
        // rho_LR,
        //     rho_C_LR, rho_W_LR, xm_C_G, xm_W_G, xm_C_L, xm_W_L, dummy, dummy,
        //     dummy, dummy, dummy, dummy, dummy, dummy);

        // const double c_p_S = 1000.;
        // const auto c_p_C_G = 1000.;  // J/kg/K
        // const auto c_p_W_G = 1914.;  // J/kg/K
        // const auto h_C_G = c_p_C_G * T;
        // const auto h_W_G = c_p_W_G * T;
        // const auto h_W_L = 112654.;  // J/kg
        // const auto h_C_L = 648584.;  // J/kg

        // const auto h_G = xm_C_G * h_C_G + xm_W_G * h_W_G;
        // const auto h_L = xm_W_L * h_W_L + xm_C_L * h_C_L;
        // const auto h_S = c_p_S * T;

        // ip_data.rho_C_GR_prev =
        //     rho_C_GR;  // C constituent partial density in gas phase
        // ip_data.rho_W_GR_prev =
        //     rho_W_GR;  // W constituent partial density in gas phase
        // ip_data.rho_C_LR_prev =
        //     rho_C_LR;  // C constituent partial density in liquid phase
        // ip_data.rho_W_LR_prev =
        //     rho_W_LR;  // W constituent partial density in liquid phase

        // ip_data.rho_G_h_G_prev = phi_G * rho_GR * h_G;
        // ip_data.rho_L_h_L_prev = phi_L * rho_LR * h_L;
        // ip_data.rho_S_h_S_prev = phi_S * rho_SR * h_S;
    }
}

// Assembles the local Jacobian matrix.
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void TH2MLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::assemble(double const t, double const dt,
                               std::vector<double> const& local_x,
                               std::vector<double> const& local_xdot,
                               std::vector<double>& local_M_data,
                               std::vector<double>& local_K_data,
                               std::vector<double>& local_rhs_data)
{
    auto const matrix_size = gas_pressure_size + capillary_pressure_size +
                             temperature_size + displacement_size;

    assert(local_x.size() == matrix_size);

    auto gas_pressure =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            gas_pressure_size> const>(local_x.data() + gas_pressure_index,
                                      gas_pressure_size);

    auto capillary_pressure =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            capillary_pressure_size> const>(
            local_x.data() + capillary_pressure_index, capillary_pressure_size);

    auto temperature =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            temperature_size> const>(local_x.data() + temperature_index,
                                     temperature_size);

    auto displacement =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + displacement_index,
                                      displacement_size);

    auto gas_pressure_dot =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            gas_pressure_size> const>(local_xdot.data() + gas_pressure_index,
                                      gas_pressure_size);
    auto capillary_pressure_dot =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            capillary_pressure_size> const>(
            local_xdot.data() + capillary_pressure_index,
            capillary_pressure_size);

    auto temperature_dot =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            temperature_size> const>(local_xdot.data() + temperature_index,
                                     temperature_size);
    auto displacement_dot =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_xdot.data() + displacement_index,
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

    // component-formulation
    // W - liquid phase main component
    // C - gas phase main component
    // pointer-matrices to the mass matrix - C component equation
    auto MCpG = M.template block<C_size, gas_pressure_size>(C_index,
                                                            gas_pressure_index);
    auto MCpC = M.template block<C_size, capillary_pressure_size>(
        C_index, capillary_pressure_index);
    auto MCT =
        M.template block<C_size, temperature_size>(C_index, temperature_index);
    auto MCu = M.template block<C_size, displacement_size>(C_index,
                                                           displacement_index);

    // pointer-matrices to the stiffness matrix - C component equation
    auto LCpG = K.template block<C_size, gas_pressure_size>(C_index,
                                                            gas_pressure_index);
    auto LCpC = K.template block<C_size, capillary_pressure_size>(
        C_index, capillary_pressure_index);
    auto LCT =
        K.template block<C_size, temperature_size>(C_index, temperature_index);

    // pointer-matrices to the mass matrix - W component equation
    auto MWpG = M.template block<W_size, gas_pressure_size>(W_index,
                                                            gas_pressure_index);
    auto MWpC = M.template block<W_size, capillary_pressure_size>(
        W_index, capillary_pressure_index);
    auto MWT =
        M.template block<W_size, temperature_size>(W_index, temperature_index);
    auto MWu = M.template block<W_size, displacement_size>(W_index,
                                                           displacement_index);

    // pointer-matrices to the stiffness matrix - W component equation
    auto LWpG = K.template block<W_size, gas_pressure_size>(W_index,
                                                            gas_pressure_index);
    auto LWpC = K.template block<W_size, capillary_pressure_size>(
        W_index, capillary_pressure_index);
    auto LWT =
        K.template block<W_size, temperature_size>(C_index, temperature_index);

    // pointer-matrices to the mass matrix - temperature equation
    auto MTpG = M.template block<temperature_size, gas_pressure_size>(
        temperature_index, gas_pressure_index);
    auto MTpC = M.template block<temperature_size, capillary_pressure_size>(
        temperature_index, capillary_pressure_index);
    auto MTT = M.template block<temperature_size, temperature_size>(
        temperature_index, temperature_index);
    auto MTu = M.template block<temperature_size, displacement_size>(
        temperature_index, displacement_index);

    // pointer-matrices to the stiffness matrix - liquid phase equation
    auto KLpG = K.template block<capillary_pressure_size, gas_pressure_size>(
        capillary_pressure_index, gas_pressure_index);
    auto KLpC =
        K.template block<capillary_pressure_size, capillary_pressure_size>(
            capillary_pressure_index, capillary_pressure_index);

    // pointer-matrices to the stiffness matrix - temperature equation
    auto ATpG = K.template block<temperature_size, gas_pressure_size>(
        temperature_index, gas_pressure_index);
    auto ATpC = K.template block<temperature_size, capillary_pressure_size>(
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

    // pointer-vectors to the right hand side terms - C-component equation
    auto fC = f.template segment<C_size>(C_index);
    // pointer-vectors to the right hand side terms - W-component equation
    auto fW = f.template segment<W_size>(W_index);
    // pointer-vectors to the right hand side terms - temperature equation
    auto fT = f.template segment<temperature_size>(temperature_index);
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

    getConstitutiveVariables(local_x, t, dt);

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
        auto const& sigma_eff = ip_data.sigma_eff;

        auto& s_L = ip_data.saturation;

        double const T0 = _process_data.reference_temperature(t, pos)[0];

        auto const T = NT.dot(temperature);
        auto const pGR = Np.dot(gas_pressure);
        auto const pCap = Np.dot(capillary_pressure);
        auto const pLR = pGR - pCap;

        auto const gradpGR = gradNp * gas_pressure;
        auto const gradpCap = gradNp * capillary_pressure;

        auto const T_dot = NT.dot(temperature_dot);
        auto const pGR_dot = Np.dot(gas_pressure_dot);
        auto const pCap_dot = Np.dot(capillary_pressure_dot);

        double div_u = mT * Bu * displacement;
        double div_u_dot = mT * Bu * displacement_dot;

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
        std::cout << "pGR: " << pGR << "\n";
        std::cout << "pCap: " << pCap << "\n";
        std::cout << "T: " << T << "\n";
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

        //   Constitutive properties

        MPL::VariableArray vars;
        vars[static_cast<int>(MPL::Variable::temperature)] = T;
        vars[static_cast<int>(MPL::Variable::phase_pressure)] = pGR;
        vars[static_cast<int>(MPL::Variable::capillary_pressure)] = pCap;
        vars[static_cast<int>(MPL::Variable::liquid_phase_pressure)] = pLR;

        bool phase_transition = _process_data.phase_transition;

        auto& rho_GR = ip_data.rho_GR;  // gas phase density
        auto& rho_C_GR =
            ip_data.rho_C_GR;  // C constituent partial density in gas phase
        auto& rho_W_GR =
            ip_data.rho_W_GR;  // W constituent partial density in gas phase
        auto& rho_LR = ip_data.rho_LR;  // liquid phase density
        auto& rho_C_LR =
            ip_data.rho_C_LR;  // C constituent partial density in liquid phase
        auto& rho_W_LR =
            ip_data.rho_W_LR;  // W constituent partial density in liquid phase

        double xm_C_G;        // C constituent gas phase mass fraction
        double xm_W_G;        // W constituent gas phase mass fraction
        double xm_C_L;        // C constituent liquid phase mass fraction
        double xm_W_L;        // W constituent liquid phase mass fraction
        double dxm_C_G_dpGR;  // C mass fraction derivative w.r.t. pGR in G
                              // phase
        double dxm_W_G_dpGR;  // W mass fraction derivative w.r.t. pGR in G
                              // phase
        double dxm_C_L_dpLR;  // C mass fraction derivative w.r.t. pLR in L
                              // phase
        double dxm_W_L_dpLR;  // W mass fraction derivative w.r.t. pLR in L
                              // phase
        double dxm_C_G_dT;    // C mass fraction derivative w.r.t. T in G phase
        double dxm_W_G_dT;    // W mass fraction derivative w.r.t. T in G phase
        double dxm_C_L_dT;    // C mass fraction derivative w.r.t. T in L phase
        double dxm_W_L_dT;    // W mass fraction derivative w.r.t. T in L phase

        // std::cout
        //     << "T pGR pCap pLR "
        //     << "rho_GR rho_C_GR rho_W_GR xm_C_G d_xm_C_G_d_pGR d_xm_C_G_d_T "
        //     << "xm_W_G d_xm_W_G_d_pGR d_xm_W_G_d_T rho_LR rho_C_LR rho_W_LR "
        //     << "xm_C_L d_xm_C_L_d_pLR d_xm_C_L_d_T xm_W_L d_xm_W_L_d_pLR "
        //     << "d_xm_W_L_d_T\n";

        // double _temperature;
        // double _gas_pressure = 100000.;
        // double _capillary_pressure = 1000.;

        // for (_temperature = 290.; _temperature < 370; _temperature += 5.)
        //     VLE(phase_transition, _temperature, _gas_pressure,
        //         _capillary_pressure, rho_GR, rho_C_GR, rho_W_GR, rho_LR,
        //         rho_C_LR, rho_W_LR, xm_C_G, xm_W_G, xm_C_L, xm_W_L,
        //         dxm_C_G_dpGR, dxm_W_G_dpGR, dxm_C_L_dpLR, dxm_W_L_dpLR,
        //         dxm_C_G_dT, dxm_W_G_dT, dxm_C_L_dT, dxm_W_L_dT);

        // _temperature = 300.;
        // for (_gas_pressure = 100000.; _gas_pressure < 200000;
        //      _gas_pressure += 5000.)
        //     VLE(phase_transition, _temperature, _gas_pressure,
        //         _capillary_pressure, rho_GR, rho_C_GR, rho_W_GR, rho_LR,
        //         rho_C_LR, rho_W_LR, xm_C_G, xm_W_G, xm_C_L, xm_W_L,
        //         dxm_C_G_dpGR, dxm_W_G_dpGR, dxm_C_L_dpLR, dxm_W_L_dpLR,
        //         dxm_C_G_dT, dxm_W_G_dT, dxm_C_L_dT, dxm_W_L_dT);

        // _gas_pressure = 100000;
        // for (_capillary_pressure = 1000.; _capillary_pressure < 2000;
        //      _capillary_pressure += 50.)
        //     VLE(phase_transition, _temperature, _gas_pressure,
        //         _capillary_pressure, rho_GR, rho_C_GR, rho_W_GR, rho_LR,
        //         rho_C_LR, rho_W_LR, xm_C_G, xm_W_G, xm_C_L, xm_W_L,
        //         dxm_C_G_dpGR, dxm_W_G_dpGR, dxm_C_L_dpLR, dxm_W_L_dpLR,
        //         dxm_C_G_dT, dxm_W_G_dT, dxm_C_L_dT, dxm_W_L_dT);

        // OGS_FATAL("Stop");

        VLE(phase_transition, T, pGR, pCap, rho_GR, rho_C_GR, rho_W_GR, rho_LR,
            rho_C_LR, rho_W_LR, xm_C_G, xm_W_G, xm_C_L, xm_W_L, dxm_C_G_dpGR,
            dxm_W_G_dpGR, dxm_C_L_dpLR, dxm_W_L_dpLR, dxm_C_G_dT, dxm_W_G_dT,
            dxm_C_L_dT, dxm_W_L_dT);

        //  - solid phase properties
        const double beta_pS = 1e-10;
        const double beta_T_SR = -1e-7;
        const double rho_ref_SR = 2500.;

        const double c_p_S = 1000.;
        const double lambda_SR = 0.5;

        //  - gas phase properties
        const double mu_GR = 1.e-5;
        const double c_p_G = 1000.;

        const double lambda_GR = 0.5;

        const auto c_p_C_G = 1000.;  // J/kg/K
        const auto c_p_W_G = 1914.;  // J/kg/K
        const auto h_C_G = c_p_C_G * T;
        const auto h_W_G = c_p_W_G * T;

        //  - liquid phase properties
        const double mu_LR = 1.e-3;
        const double c_p_L = 4000.;

        const double lambda_LR = 0.5;
        // constant liquid constituent enthalpies ftw!
        const auto h_W_L = 112654.;  // J/kg
        const auto h_C_L = 648584.;  // J/kg

        const auto I =
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Identity();
        const double sD_G = 0.1;
        const double sD_L = 0.1;

        const auto D_C_G = sD_G * I;
        const auto D_W_G = sD_G * I;
        const auto D_C_L = sD_L * I;
        const auto D_W_L = sD_L * I;

        //  - medium properties
        auto const k_S = MPL::formEigenTensor<DisplacementDim>(
            medium.property(MPL::PropertyType::permeability)
                .value(vars, pos, t, dt));

        const auto s_L_dot = (s_L - ip_data.saturation_prev) / dt;

        vars[static_cast<int>(MPL::Variable::liquid_saturation)] = s_L;

        auto const s_G = 1. - s_L;

        auto const alpha_B =
            medium.property(MPL::PropertyType::biot_coefficient)
                .template value<double>(vars, pos, t, dt);

        auto const beta_p_SR = (1. - alpha_B) * beta_pS;

        auto const k_rel =
            medium.property(MPL::PropertyType::relative_permeability)
                .template value<Eigen::Vector2d>(vars, pos, t, dt);

        auto const k_rel_L = k_rel[0];
        auto const k_rel_G = k_rel[1];

        auto const& b = _process_data.specific_body_force;

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
        double const delta_T(T - T0);
        double const thermal_strain = beta_T_SR * delta_T;
        auto const rho_SR = rho_ref_SR * (1. - 3. * thermal_strain);

        auto const rho = rho_GR + rho_LR + rho_SR;
        // TODO: change back to rho_SR
        // auto const rho_c_p = phi_G * rho_GR * c_p_G + phi_L * rho_LR * c_p_L
        // +
        //                      phi_S * rho_ref_SR * c_p_S;

        // update secondary variables. TODO: Refactoring potential!!

        _liquid_pressure[ip] = pLR;
        _liquid_density[ip] = rho_LR;
        _gas_density[ip] = rho_GR;
        _porosity[ip] = phi;
        _saturation[ip] = s_L;

        // abbreviations
        const double rho_C_FR = s_G * rho_C_GR + s_L * rho_C_LR;
        const double rho_W_FR = s_G * rho_W_GR + s_L * rho_W_LR;

        // pore fluid pressure
        const double p_FR = s_G * pGR + s_L * pLR;
        const double p_FR_dot = pGR_dot - (s_L_dot * pCap + pCap_dot * s_L);

        // porosity derivative
        auto const phi_dot = (alpha_B - phi) * (div_u_dot - beta_T_SR * T_dot +
                                                beta_p_SR * p_FR_dot);

        // solid phase pressure
        const double p_SR =
            p_FR - 1. / (beta_pS * phi_S) * (div_u - thermal_strain);

        // const double p_SR_dot =
        //     p_FR_dot -
        //     1. / (beta_pS * phi_S) * (div_u_dot - beta_T_SR * T_dot) -
        //     phi_dot / phi_S * (p_FR - p_SR);

        auto& phi_S_p_SR = ip_data.phi_S_p_SR;
        phi_S_p_SR = phi_S * p_SR;

        const double phi_S_p_SR_dot =
            (phi_S_p_SR - ip_data.phi_S_p_SR_prev) / dt;

        // phase specific enthalpies
        const auto h_G = xm_C_G * h_C_G + xm_W_G * h_W_G;
        const auto h_L = xm_W_L * h_W_L + xm_C_L * h_C_L;
        const auto h_S = c_p_S * T;

        // phase enthalpies
        const double rho_G_h_G_dot =
            (ip_data.rho_G_h_G - ip_data.rho_G_h_G_prev) / dt;
        const double rho_L_h_L_dot =
            (ip_data.rho_L_h_L - ip_data.rho_L_h_L_prev) / dt;
        const double rho_S_h_S_dot =
            (ip_data.rho_S_h_S - ip_data.rho_S_h_S_prev) / dt;

        const double rho_C_GR_dot = (rho_C_GR - ip_data.rho_C_GR_prev) / dt;
        const double rho_C_LR_dot = (rho_C_LR - ip_data.rho_C_LR_prev) / dt;
        const double rho_W_GR_dot = (rho_W_GR - ip_data.rho_W_GR_prev) / dt;
        const double rho_W_LR_dot = (rho_W_LR - ip_data.rho_W_LR_prev) / dt;

        const auto rho_h_eff =
            ip_data.rho_G_h_G + ip_data.rho_L_h_L + ip_data.rho_S_h_S;

        const auto k_over_mu_G = k_S * k_rel_G / mu_GR;
        const auto k_over_mu_L = k_S * k_rel_L / mu_LR;

        GlobalDimVectorType const w_GS =
            k_over_mu_G * rho_GR * b - k_over_mu_G * gradpGR;

        GlobalDimVectorType const w_LS = k_over_mu_L * gradpCap +
                                         k_over_mu_L * rho_GR * b -
                                         k_over_mu_L * gradpGR;


#ifndef MATERIAL_PROPERTIES
        std::cout << "######################################################\n";
        std::cout << "#    Material properties:\n";
        std::cout << "#----------------------------------------------------#\n";
        std::cout << "#         rho_GR:  " << rho_GR << "\n";
        std::cout << "#       rho_C_GR:  " << rho_C_GR << "\n";
        std::cout << "#       rho_W_GR:  " << rho_W_GR << "\n";
        std::cout << "#    .    .    .    .    .    .    .    .    .    .  #\n";
        std::cout << "#   rho_C_GR_dot:  " << rho_C_GR_dot << "\n";
        std::cout << "#   rho_W_GR_dot:  " << rho_W_GR_dot << "\n";
        std::cout << "#    .    .    .    .    .    .    .    .    .    .  #\n";
        std::cout << "#         xm_C_G:  " << xm_C_G << "\n";
        std::cout << "#         xm_W_G:  " << xm_W_G << "\n";
        std::cout << "#    .    .    .    .    .    .    .    .    .    .  #\n";
        std::cout << "#  d_xm_C_G_dpGR:  " << dxm_C_G_dpGR << "\n";
        std::cout << "#  d_xm_W_G_dpGR:  " << dxm_W_G_dpGR << "\n";
        std::cout << "#    d_xm_C_G_dT:  " << dxm_C_G_dT << "\n";
        std::cout << "#    d_xm_W_G_dT:  " << dxm_W_G_dT << "\n";
        std::cout << "#----------------------------------------------------#\n";
        std::cout << "#         rho_LR:  " << rho_LR << "\n";
        std::cout << "#       rho_C_LR:  " << rho_C_GR << "\n";
        std::cout << "#       rho_W_LR:  " << rho_W_GR << "\n";
        std::cout << "#    .    .    .    .    .    .    .    .    .    .  #\n";
        std::cout << "#   rho_C_LR_dot:  " << rho_C_LR_dot << "\n";
        std::cout << "#   rho_W_LR_dot:  " << rho_W_LR_dot << "\n";
        std::cout << "#    .    .    .    .    .    .    .    .    .    .  #\n";
        std::cout << "#         xm_C_L:  " << xm_C_G << "\n";
        std::cout << "#         xm_W_L:  " << xm_W_G << "\n";
        std::cout << "#    .    .    .    .    .    .    .    .    .    .  #\n";
        std::cout << "#  d_xm_C_L_dpLR:  " << dxm_C_L_dpLR << "\n";
        std::cout << "#  d_xm_W_L_dpLR:  " << dxm_W_L_dpLR << "\n";
        std::cout << "#    d_xm_C_L_dT:  " << dxm_C_L_dT << "\n";
        std::cout << "#    d_xm_W_L_dT:  " << dxm_W_L_dT << "\n";
        std::cout << "#----------------------------------------------------#\n";
        std::cout << "#         rho_SR:  " << rho_SR << "\n";
        std::cout << "#            rho:  " << rho << "\n";
        std::cout << "#----------------------------------------------------#\n";
        std::cout << "#          c_p_G:  " << c_p_G << "\n";
        std::cout << "#          c_p_L:  " << c_p_L << "\n";
        std::cout << "#          c_p_S:  " << c_p_S << "\n";
        std::cout << "#----------------------------------------------------#\n";
        std::cout << "#      beta_p_SR:  " << beta_p_SR << "\n";
        std::cout << "#      beta_T_SR:  " << beta_T_SR << "\n";
        std::cout << "#----------------------------------------------------#\n";
        std::cout << "#            h_G:  " << h_G << "\n";
        std::cout << "#            h_L:  " << h_L << "\n";
        std::cout << "#            h_S:  " << h_S << "\n";
        std::cout << "#    .    .    .    .    .    .    .    .    .    .  #\n";
        std::cout << "#      rho_h_eff:  " << rho_h_eff << "\n";
        std::cout << "#  rho_G_h_G:  " << ip_data.rho_G_h_G << "\n";
        std::cout << "#  rho_L_h_L:  " << ip_data.rho_L_h_L << "\n";
        std::cout << "#  rho_S_h_S:  " << ip_data.rho_S_h_S << "\n";
        std::cout << "#    .    .    .    .    .    .    .    .    .    .  #\n";
        std::cout << "#  rho_G_h_G_dot:  " << rho_G_h_G_dot << "\n";
        std::cout << "#  rho_L_h_L_dot:  " << rho_L_h_L_dot << "\n";
        std::cout << "#  rho_S_h_S_dot:  " << rho_S_h_S_dot << "\n";
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
        OGS_FATAL("It has to stop!");
#endif

        // coefficient matrices
        // C-component equation
        MCpG.noalias() += NpT * rho_C_FR * (alpha_B - phi) * beta_p_SR * Np * w;

        MCpC.noalias() -=
            NpT * rho_C_FR * (alpha_B - phi) * beta_p_SR * s_L * Np * w;

        MCT.noalias() -= NpT * rho_C_FR * (alpha_B - phi) * beta_T_SR * Np * w;

        MCu.noalias() += NpT * rho_C_FR * alpha_B * mT * Bu * w;

        const auto advection_C_G = rho_C_GR * k_over_mu_G;
        const auto advection_C_L = rho_C_LR * k_over_mu_L;
        const auto diffusion_C_G_p = phi_G * rho_GR * D_C_G * dxm_C_L_dpLR;
        const auto diffusion_C_L_p = phi_L * rho_LR * D_C_L * dxm_C_G_dpGR;
        const auto diffusion_C_G_T = phi_G * rho_GR * D_C_G * dxm_C_G_dT;
        const auto diffusion_C_L_T = phi_L * rho_LR * D_C_L * dxm_C_L_dT;

        const auto advection_C = advection_C_G + advection_C_L;
        const auto diffusion_C_p = diffusion_C_G_p + diffusion_C_L_p;
        const auto diffusion_C_T = diffusion_C_G_T + diffusion_C_L_T;

        LCpG.noalias() += gradNpT * (advection_C + diffusion_C_p) * gradNp * w;

        LCpC.noalias() -=
            gradNpT * (advection_C_L + diffusion_C_L_p) * gradNp * w;

        LCT.noalias() += gradNpT * (diffusion_C_T)*gradNp * w;

        fC.noalias() +=
            gradNpT * (advection_C_G * rho_GR + advection_C_L * rho_LR) * b * w;

        fC.noalias() -= NpT *
                        (phi * (rho_C_LR - rho_C_GR) -
                         rho_C_FR * pCap * (alpha_B - phi) * beta_p_SR) *
                        s_L_dot * w;

        fC.noalias() -=
            NpT * phi * (s_G * rho_C_GR_dot + s_L * rho_C_LR_dot) * w;

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        // W-component equation
        MWpG.noalias() += NpT * rho_W_FR * (alpha_B - phi) * beta_p_SR * Np * w;

        MWpC.noalias() -=
            NpT * rho_W_FR * (alpha_B - phi) * beta_p_SR * s_L * Np * w;

        MWT.noalias() -= NpT * rho_W_FR * (alpha_B - phi) * beta_T_SR * Np * w;

        MWu.noalias() += NpT * rho_W_FR * alpha_B * mT * Bu * w;

        const auto advection_W_G = rho_W_GR * k_over_mu_G;
        const auto advection_W_L = rho_W_LR * k_over_mu_L;
        const auto diffusion_W_G_p = phi_G * rho_GR * D_W_G * dxm_W_L_dpLR;
        const auto diffusion_W_L_p = phi_L * rho_LR * D_W_L * dxm_W_G_dpGR;
        const auto diffusion_W_G_T = phi_G * rho_GR * D_W_G * dxm_W_G_dT;
        const auto diffusion_W_L_T = phi_L * rho_LR * D_W_L * dxm_W_L_dT;

        const auto advection_W = advection_W_G + advection_W_L;
        const auto diffusion_W_p = diffusion_W_G_p + diffusion_W_L_p;
        const auto diffusion_W_T = diffusion_W_G_T + diffusion_W_L_T;

        LWpG.noalias() += gradNpT * (advection_W + diffusion_W_p) * gradNp * w;

        LWpC.noalias() -=
            gradNpT * (advection_W_L + diffusion_W_L_p) * gradNp * w;

        LWT.noalias() += gradNpT * (diffusion_W_T)*gradNp * w;

        fW.noalias() +=
            gradNpT * (advection_W_G * rho_GR + advection_W_L * rho_LR) * b * w;

        fW.noalias() -= NpT *
                        (phi * (rho_W_LR - rho_W_GR) -
                         rho_W_FR * pCap * (alpha_B - phi) * beta_p_SR) *
                        s_L_dot * w;

        fW.noalias() -=
            NpT * phi * (s_G * rho_W_GR_dot + s_L * rho_W_LR_dot) * w;

        //  - temperature equation
        MTpG.noalias() -= NTT * phi * NT * w;
        MTpC.noalias() += NTT * phi_L * NT * w;

        MTu.noalias() += NTT * rho_h_eff * mT * Bu * w;

        // MTT.noalias() += NTT * 0. * NT * w;

        // ATpG.noalias() -= NTT *
        //                   ((1. - beta_TGR * T) * w_GS.transpose() +
        //                    (1. - beta_TLR * T) * w_LS.transpose()) *
        //                   gradNT * w;
        // ATpC.noalias() +=
        //     NTT * (1. - beta_TLR * T) * w_LS.transpose() * gradNT * w;

        // // ATT
        // KTT.noalias() += NTT *
        //                  (rho_GR * c_p_G * w_GS.transpose() +
        //                   rho_LR * c_p_L * w_LS.transpose()) *
        //                  gradNT * w;
        // LTT

        KTT.noalias() += gradNTT * lambda * gradNT * w;

        // fT
        fT.noalias() -=
            NTT *
            (rho_G_h_G_dot + rho_L_h_L_dot + rho_S_h_S_dot +
             phi * pCap * s_L_dot - p_FR * phi_dot - phi_S_p_SR_dot) *
            w;

        fT.noalias() +=
            NTT * (rho_GR * w_GS.transpose() + rho_LR * w_LS.transpose()) * b *
            w;

        fT.noalias() +=
            gradNTT * (rho_GR * h_G * w_GS + rho_LR * h_L * w_LS) * w;

        //  - displacement equation
        KUpG.noalias() -= (BuT * alpha_B * m * Np) * w;
        KUpC.noalias() += (BuT * alpha_B * s_L * m * Np) * w;

        //        eps.noalias() = Bu * displacement;
        ip_data.updateConstitutiveRelationThermal(
            t, pos, dt, displacement,
            _process_data.reference_temperature(t, pos)[0], thermal_strain);
        fU.noalias() += (BuT * sigma_eff - Nu_op.transpose() * rho * b) * w;

#ifdef DEBUG_OUTPUT
        std::cout << "------------------------------------------------------\n";
        std::cout << " MCpG:\n" << MCpG << "\n";

        std::cout << " MCpG \n" << MCpG  << "\n";
        std::cout << " MCpC\n" << MCpC << "\n";
        std::cout << " MCT\n" << MCT << "\n";
        std::cout << " MCu\n" << MCu << "\n";
        std::cout << " advection_C_G\n" << advection_C_G << "\n";
        std::cout << " advection_C_L\n" << advection_C_L << "\n";
        std::cout << " diffusion_C_G_p\n" << diffusion_C_G_p << "\n";
        std::cout << " diffusion_C_L_p\n" << diffusion_C_L_p << "\n";
        std::cout << " diffusion_C_G_T\n" << diffusion_C_G_T << "\n";
        std::cout << " diffusion_C_L_T\n" << diffusion_C_L_T << "\n";
        std::cout << " advection_C\n" << advection_C << "\n";
        std::cout << " diffusion_C_p\n" << diffusion_C_p << "\n";
        std::cout << " diffusion_C_T\n" << diffusion_C_T << "\n";
        std::cout << " LCpG\n" << LCpG << "\n";
        std::cout << " LCpC\n" << LCpC << "\n";
        std::cout << " LCT\n" << LCT << "\n";
        std::cout << " fC\n" << fC << "\n";
        std::cout << " MWpG\n" << MWpG << "\n";
        std::cout << " MWpC\n" << MWpC << "\n";
        std::cout << " MWT\n" << MWT << "\n";
        std::cout << " MWu\n" << MWu << "\n";
        std::cout << " advection_W_G\n" << advection_W_G << "\n";
        std::cout << " advection_W_L\n" << advection_W_L << "\n";
        std::cout << " diffusion_W_G_p\n" << diffusion_W_G_p << "\n";
        std::cout << " diffusion_W_L_p\n" << diffusion_W_L_p << "\n";
        std::cout << " diffusion_W_G_T\n" << diffusion_W_G_T << "\n";
        std::cout << " diffusion_W_L_T\n" << diffusion_W_L_T << "\n";
        std::cout << " advection_W\n" << advection_W << "\n";
        std::cout << " diffusion_W_p\n" << diffusion_W_p << "\n";
        std::cout << " diffusion_W_T\n" << diffusion_W_T << "\n";
        std::cout << " LWpG\n" << LWpG << "\n";
        std::cout << " LWpC\n" << LWpC << "\n";
        std::cout << " LWT\n" << LWT << "\n";
        std::cout << " fW\n" << fW << "\n";
        std::cout << " MTpG\n" << MTpG << "\n";
        std::cout << " MTpC\n" << MTpC << "\n";
        std::cout << " MTu\n" << MTu << "\n";
        std::cout << " KTT\n" << KTT << "\n";
        std::cout << " fT\n" << fT << "\n";
        std::cout << " KUpG\n" << KUpG << "\n";
        std::cout << " KUpC\n" << KUpC << "\n";
        std::cout << " fU\n" << fU << "\n";

        OGS_FATAL("Stop intended.");
#endif
    }
}

// Assembles the local Jacobian matrix.
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

    auto displacement =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + displacement_offset,
                                      displacement_size);

    auto temperature =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            temperature_size> const>(local_x.data() + temperature_index,
                                     temperature_size);
    auto gas_pressure =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            gas_pressure_size> const>(local_x.data() + gas_pressure_index,
                                      gas_pressure_size);

    auto capillary_pressure =
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

        auto const T = NT.dot(temperature);
        auto const pGR = Np.dot(gas_pressure);
        auto const pCap = Np.dot(capillary_pressure);
        auto const pLR = pGR - pCap;

        vars[static_cast<int>(MPL::Variable::temperature)] = T;
        vars[static_cast<int>(MPL::Variable::phase_pressure)] = pGR;
        vars[static_cast<int>(MPL::Variable::capillary_pressure)] = pCap;

        auto const solid_linear_thermal_expansion_coefficient =
            solid_phase.property(MPL::PropertyType::thermal_expansivity)
                .template value<double>(vars, pos, t, dt);

        double const delta_T(T - T0);
        double const thermal_strain =
            solid_linear_thermal_expansion_coefficient * delta_T;

        auto& eps = _ip_data[ip].eps;
        eps.noalias() = B * displacement;

        _ip_data[ip].updateConstitutiveRelationThermal(
            t, pos, dt, displacement,
            _process_data.reference_temperature(t, pos)[0], thermal_strain);

        auto const rho_LR = liquid_phase.property(MPL::PropertyType::density)
                                .template value<double>(vars, pos, t, dt);

        auto const phi = medium.property(MPL::PropertyType::porosity)
                             .template value<double>(vars, pos, t, dt);

        auto const rho_GR = gas_phase.property(MPL::PropertyType::density)
                                .template value<double>(vars, pos, t, dt);

        auto const s_L = medium.property(MPL::PropertyType::saturation)
                             .template value<double>(vars, pos, t, dt);

        _liquid_pressure[ip] = pLR;
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
    auto gas_pressure =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            gas_pressure_size> const>(local_x.data() + gas_pressure_index,
                                      gas_pressure_size);
    auto capillary_pressure =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            capillary_pressure_size> const>(
            local_x.data() + capillary_pressure_index, capillary_pressure_size);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(_element, _is_axially_symmetric, gas_pressure,
                         *_process_data.gas_pressure_interpolated);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(_element, _is_axially_symmetric, capillary_pressure,
                         *_process_data.capillary_pressure_interpolated);

    auto temperature =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            temperature_size> const>(local_x.data() + temperature_index,
                                     temperature_size);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(_element, _is_axially_symmetric, temperature,
                         *_process_data.temperature_interpolated);
}

}  // namespace TH2M
}  // namespace ProcessLib