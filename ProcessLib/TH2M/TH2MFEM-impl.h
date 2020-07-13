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

#define PRINT(var) \
    std::cout << std::setprecision(16) << #var << ": " << var << "\n";

#define PRINT2(bool, text, var) \
    if (bool)                   \
        std::cout << std::setprecision(16) << #text << ": " << var << "\n";

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
      _saturation(std::vector<double>(_integration_method.getNumberOfPoints())),
      _mole_fraction_gas(
          std::vector<double>(_integration_method.getNumberOfPoints())),
      _mole_fraction_liquid(
          std::vector<double>(_integration_method.getNumberOfPoints())),
      _mass_fraction_gas(
          std::vector<double>(_integration_method.getNumberOfPoints())),
      _mass_fraction_liquid(
          std::vector<double>(_integration_method.getNumberOfPoints()))
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
        d_xm_C_G_d_pGR =
            1. / rho_GR * (d_rho_C_GR_d_pGR - xm_C_G * d_rho_GR_d_pGR);
        d_xm_W_G_d_pGR =
            1. / rho_GR * (d_rho_W_GR_d_pGR - xm_W_G * d_rho_GR_d_pGR);
        d_xm_C_L_d_pLR =
            1. / rho_LR * (d_rho_C_LR_d_pLR - xm_C_L * d_rho_LR_d_pLR);
        d_xm_W_L_d_pLR =
            1. / rho_LR * (d_rho_W_LR_d_pLR - xm_W_L * d_rho_LR_d_pLR);

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

        double const T = NT.dot(temperature);
        double const pGR = Np.dot(gas_pressure);
        double const pCap = Np.dot(capillary_pressure);
        double const pLR = pGR - pCap;

        double div_u = mT * Bu * displacement;

        MPL::VariableArray vars;
        vars[static_cast<int>(MPL::Variable::temperature)] = T;
        vars[static_cast<int>(MPL::Variable::phase_pressure)] = pGR;
        vars[static_cast<int>(MPL::Variable::capillary_pressure)] = pCap;
        vars[static_cast<int>(MPL::Variable::liquid_phase_pressure)] = pLR;

        auto const s_L =
            medium.property(MPL::PropertyType::saturation)
                .template value<double>(
                    vars, pos, t, std::numeric_limits<double>::quiet_NaN());

        vars[static_cast<int>(MPL::Variable::liquid_saturation)] = s_L;

        auto const phi = medium.property(MPL::PropertyType::porosity)
                             .template value<double>(vars, pos, t, dt);

        ip_data.phi = phi;
        ip_data.saturation = s_L;

        auto const phi_G = (1. - s_L) * phi;
        auto const phi_L = s_L * phi;
        auto const phi_S = 1. - phi;

        const double p_FR =
            (1. - ip_data.saturation) * pGR + ip_data.saturation * pLR;

        double const T0 = _process_data.reference_temperature(t, pos)[0];
        double const delta_T(T - T0);

        ip_data.beta_pS =
            solid_phase.property(MPL::PropertyType::compressibility)
                .template value<double>(
                    vars, pos, t, std::numeric_limits<double>::quiet_NaN());
        ip_data.beta_T_SR =
            solid_phase.property(MPL::PropertyType::thermal_expansivity)
                .template value<double>(
                    vars, pos, t, std::numeric_limits<double>::quiet_NaN());
        const auto rho_ref_SR =
            solid_phase.property(MPL::PropertyType::density)
                .template value<double>(
                    vars, pos, t, std::numeric_limits<double>::quiet_NaN());

        ip_data.c_p_S =
            solid_phase.property(MPL::PropertyType::specific_heat_capacity)
                .template value<double>(vars, pos, t, dt);

        ip_data.lambdaSR =
            solid_phase.property(MPL::PropertyType::thermal_conductivity)
                .template value<double>(vars, pos, t, dt);

        ip_data.alpha_B = medium.property(MPL::PropertyType::biot_coefficient)
                              .template value<double>(vars, pos, t, dt);

        ip_data.beta_p_SR = (1. - ip_data.alpha_B) * ip_data.beta_pS;

        ip_data.thermal_strain = ip_data.beta_T_SR * delta_T;

        // cf. Eq. 79
        ip_data.p_SR = p_FR - 1. / (ip_data.beta_pS * phi_S) *
                                  (div_u - ip_data.thermal_strain);

        // const double p_SR_dot =
        //     p_FR_dot -
        //     1. / (beta_pS * phi_S) * (div_u_dot - beta_T_SR * T_dot) -
        //     phi_dot / phi_S * (p_FR - p_SR);

        // cf. WW: THM
        const auto rhoSR = rho_ref_SR * (1. - 3. * ip_data.thermal_strain);

        // TODO: check if needed!
        // ip_data.phi_S_p_SR = phi_S * p_SR;
        ip_data.phi_S_p_SR = 0.;

        ip_data.k_S = MPL::formEigenTensor<DisplacementDim>(
            medium.property(MPL::PropertyType::permeability)
                .value(vars, pos, t, dt));

        auto const k_rel =
            medium.property(MPL::PropertyType::relative_permeability)
                .template value<Eigen::Vector2d>(vars, pos, t, dt);

        ip_data.kRel = k_rel;

        // #####################################################

        auto const nComponentsGas = gas_phase.numberOfComponents();
        auto const nComponentsLiquid = liquid_phase.numberOfComponents();
        // PRINT(nComponentsGas);

        // PRINT(pGR);
        // PRINT(pCap);
        // PRINT(pLR);
        // PRINT(T);

        Eigen::Vector2d xnG = {1., 0.};

        if (nComponentsGas == 2)
        {
            xnG = gas_phase.property(MPL::PropertyType::mole_fraction)
                      .template value<Eigen::Vector2d>(vars, pos, t, dt);
        }

        // PRINT(xnG[0]);
        // PRINT(xnG[1]);

        ip_data.xnCG = xnG[0];

        vars[static_cast<int>(MPL::Variable::mole_fraction)] = xnG[0];

        auto const pCGR = xnG[0] * pGR;
        auto const pWGR = xnG[1] * pGR;

        // PRINT(pCGR);
        // PRINT(pWGR);

        auto rhoGR = gas_phase.property(MPL::PropertyType::density)
                         .template value<double>(vars, pos, t, dt);
        // PRINT(rhoGR);

        auto drhoGR_dpGR =
            gas_phase.property(MPL::PropertyType::density)
                .template dValue<double>(vars, MPL::Variable::phase_pressure,
                                         pos, t, dt);
        // PRINT(drhoGR_dpGR);

        auto beta_pGR = 1. / rhoGR * drhoGR_dpGR;
        // PRINT(beta_pGR);

        auto drhoGR_dT = gas_phase.property(MPL::PropertyType::density)
                             .template dValue<double>(
                                 vars, MPL::Variable::temperature, pos, t, dt);
        // PRINT(drhoGR_dT);

        auto beta_TGR = 1. / rhoGR * drhoGR_dT;
        // PRINT(beta_TGR);

        auto pWVap = 0.;
        auto dpWVap_dT = 0.;

        Eigen::Vector2d dxm_dpGR = {0., 0.};
        Eigen::Vector2d xmG = {1., 0.};

        if (nComponentsGas == 2)
        {
            pWVap = gas_phase.component(1)
                        .property(MPL::PropertyType::saturation_vapour_pressure)
                        .template value<double>(vars, pos, t, dt);
            dpWVap_dT =
                gas_phase.component(1)
                    .property(MPL::PropertyType::saturation_vapour_pressure)
                    .template dValue<double>(vars, MPL::Variable::temperature,
                                             pos, t, dt);

            auto const MG = gas_phase.property(MPL::PropertyType::molar_mass)
                                .template value<double>(vars, pos, t, dt);
            // PRINT(MG);
            xmG[0] = xnG[0] *
                     gas_phase.component(0)
                         .property(MPL::PropertyType::molar_mass)
                         .template value<double>(vars, pos, t, dt) /
                     MG;
            xmG[1] = 1. - xmG[0];
        }

        auto xmCG = xmG[0];
        auto xmWG = xmG[1];

        // PRINT(xmG[0]);
        // PRINT(xmG[1]);

        auto const dxmWG_dpGR = xmG[1] * beta_pGR;
        auto const dxmCG_dpGR = -dxmWG_dpGR;
        // PRINT(dxmWG_dpGR);
        // PRINT(dxmCG_dpGR);

        // PRINT(pWVap);
        // PRINT(dpWVap_dT);

        auto rhoCGR = xmG[0] * rhoGR;
        auto rhoWGR = xmG[1] * rhoGR;
        // PRINT(rhoCGR);
        // PRINT(rhoWGR);

        auto drhoWGR_dT = 0.;
        if (pWGR != 0)
        {
            drhoWGR_dT = rhoWGR * (1. / pWGR * dpWVap_dT - 1. / T);
        }
        // PRINT(drhoWGR_dT);

        auto const dxmWG_dT = 1. / rhoGR * drhoWGR_dT - xmG[1] * beta_TGR;
        auto const dxmCG_dT = -dxmWG_dT;
        // PRINT(dxmCG_dT);
        // PRINT(dxmWG_dT);

        auto H = 0.;
        auto dH_dT = 0.;
        if (nComponentsLiquid == 2)
        {
            H = liquid_phase.component(0)
                    .property(MPL::PropertyType::henry_coefficient)
                    .template value<double>(vars, pos, t, dt);
            dH_dT = liquid_phase.component(0)
                        .property(MPL::PropertyType::henry_coefficient)
                        .template dValue<double>(
                            vars, MPL::Variable::temperature, pos, t, dt);
        }
        // PRINT(H);
        // PRINT(dH_dT);

        auto cCL = 0.;
        auto dcCL_dpLR = 0.;
        if (liquid_phase.hasProperty(MPL::PropertyType::concentration))
        {
            cCL = liquid_phase.property(MPL::PropertyType::concentration)
                      .template value<double>(vars, pos, t,
                                              dt);  // in mol*m^(-3)
            dcCL_dpLR =
                liquid_phase.property(MPL::PropertyType::concentration)
                    .template dValue<double>(
                        vars, MPL::Variable::liquid_phase_pressure, pos, t,
                        dt);  // in mol*m^(-3)
        }
        // PRINT(cCL);
        // PRINT(dcCL_dpLR);

        auto dcCL_dT = dH_dT * xnG[0] * pGR - H * dpWVap_dT;
        // PRINT(dcCL_dT);

        vars[static_cast<int>(MPL::Variable::concentration)] = cCL;

        auto rhoLR =
            liquid_phase.property(MPL::PropertyType::density)
                .template value<double>(vars, pos, t, dt);  // in mol*m^(-3)

        vars[static_cast<int>(MPL::Variable::concentration)] = 0.;
        auto rhoWLR = liquid_phase.property(MPL::PropertyType::density)
                          .template value<double>(vars, pos, t, dt);

        auto rhoCLR = rhoLR - rhoWLR;

        // PRINT(rhoLR);
        // PRINT(rhoWLR);
        // PRINT(rhoCLR);

        const auto xmWL = rhoWLR / rhoLR;
        const auto xmCL = 1. - xmWL;

        // PRINT(xmWL);
        // PRINT(xmCL);

        // HACK!! int the linear EOS of the liquid, the third 'independent'
        // variable is actually pressure-dependent, thus not independent.
        // The exact derivative would reed
        // drho_dpLR=rho_ref*(betaP+betaC*dC_dpLR). In order to get the
        // parameters of the eos (rho_ref, betaP, and betaC), I firstly
        // calculate the pressure derivative and then the concenctration
        // derivative only.
        const auto rho_ref_betaP =
            liquid_phase.property(MPL::PropertyType::density)
                .template dValue<double>(
                    vars, MPL::Variable::liquid_phase_pressure, pos, t, dt);
        // PRINT(rho_ref_betaP);

        const auto rho_ref_betaC =
            liquid_phase.property(MPL::PropertyType::density)
                .template dValue<double>(vars, MPL::Variable::concentration,
                                         pos, t, dt);
        // PRINT(rho_ref_betaC);
        const auto drhoLR_dpLR = rho_ref_betaP + rho_ref_betaC * H;
        // PRINT(drhoLR_dpLR);

        const auto rho_ref_betaT =
            liquid_phase.property(MPL::PropertyType::density)
                .template dValue<double>(vars, MPL::Variable::temperature, pos,
                                         t, dt);
        // PRINT(rho_ref_betaT);
        const auto drhoLR_dT = rho_ref_betaT + rho_ref_betaC * dcCL_dT;
        // PRINT(drhoLR_dT);

        auto drhoWLR_dpLR = rho_ref_betaP;
        // PRINT(drhoWLR_dpLR);

        const auto drhoWLR_dT = rho_ref_betaT;
        // PRINT(drhoWLR_dT);

        const auto dxmWL_dpLR =
            1. / rhoLR * (drhoWLR_dpLR - xmWL * drhoLR_dpLR);
        const auto dxmCL_dpLR = -dxmWL_dpLR;
        // PRINT(dxmWL_dpLR);
        // PRINT(dxmCL_dpLR);
        auto const dxmWL_dT = 1. / rhoLR * (drhoWLR_dT - xmWL * drhoLR_dT);
        auto const dxmCL_dT = -dxmWL_dT;
        // PRINT(dxmWL_dT);
        // PRINT(dxmCL_dT);

        const double c_p_S =
            solid_phase.property(MPL::PropertyType::specific_heat_capacity)
                .template value<double>(vars, pos, t, dt);
        auto h_G = 0.;
        auto h_L = 0.;
        const auto h_S = c_p_S * T;

        if (nComponentsGas == 2)
        {
            auto c_p_C_G =
                gas_phase.component(0)
                    .property(MPL::PropertyType::specific_heat_capacity)
                    .template value<double>(vars, pos, t, dt);
            auto c_p_W_G =
                gas_phase.component(1)
                    .property(MPL::PropertyType::specific_heat_capacity)
                    .template value<double>(vars, pos, t, dt);

            auto h_C_G = c_p_C_G * T;
            auto h_W_G = c_p_W_G * T;
            h_G = xmCG * h_C_G + xmWG * h_W_G;
        }
        else
        {
            h_G = gas_phase.property(MPL::PropertyType::specific_heat_capacity)
                      .template value<double>(vars, pos, t, dt) *
                  T;
        }
        if (nComponentsLiquid == 2)
        {
            auto c_p_C_L =
                liquid_phase.component(0)
                    .property(MPL::PropertyType::specific_heat_capacity)
                    .template value<double>(vars, pos, t, dt);
            auto c_p_W_L =
                liquid_phase.component(1)
                    .property(MPL::PropertyType::specific_heat_capacity)
                    .template value<double>(vars, pos, t, dt);
            auto h_C_L = c_p_C_L * T;
            auto h_W_L = c_p_W_L * T;
            h_L = xmCL * h_C_L + xmWL * h_W_L;
        }
        else
        {
            h_L =
                liquid_phase.property(MPL::PropertyType::specific_heat_capacity)
                    .template value<double>(vars, pos, t, dt) *
                T;
        }

        ip_data.rho_GR = rhoGR;
        ip_data.rho_LR = rhoLR;

        ip_data.rho_C_GR = rhoCGR;
        ip_data.rho_W_GR = rhoWGR;
        ip_data.rho_C_LR = rhoCLR;
        ip_data.rho_W_LR = rhoWLR;

        ip_data.xmCG = xmCG;
        ip_data.xmWG = xmWG;
        ip_data.xmCL = xmCL;
        ip_data.xmWL = xmWL;

        ip_data.dxmCG_dpGR = dxmCG_dpGR;
        ip_data.dxmWG_dpGR = dxmWG_dpGR;
        ip_data.dxmCL_dpLR = dxmCL_dpLR;
        ip_data.dxmWL_dpLR = dxmWL_dpLR;
        ip_data.dxmCG_dT = dxmCG_dT;
        ip_data.dxmWG_dT = dxmWG_dT;
        ip_data.dxmCL_dT = dxmCL_dT;
        ip_data.dxmWL_dT = dxmWL_dT;

        ip_data.h_G = h_G;
        ip_data.h_L = h_L;
        ip_data.h_S = h_S;

        ip_data.rho_G_h_G = phi_G * rhoGR * h_G;
        ip_data.rho_L_h_L = phi_L * rhoLR * h_L;
        ip_data.rho_S_h_S = phi_S * rhoSR * h_S;
        ip_data.muGR = gas_phase.property(MPL::PropertyType::viscosity)
                           .template value<double>(vars, pos, t, dt);

        ip_data.lambdaGR =
            gas_phase.property(MPL::PropertyType::thermal_conductivity)
                .template value<double>(vars, pos, t, dt);

        ip_data.muLR = liquid_phase.property(MPL::PropertyType::viscosity)
                           .template value<double>(vars, pos, t, dt);

        ip_data.lambdaLR =
            liquid_phase.property(MPL::PropertyType::thermal_conductivity)
                .template value<double>(vars, pos, t, dt);
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
        // pos.setIntegrationPoint(ip);
        auto& ip_data = _ip_data[ip];

        // auto const& Np = ip_data.N_p;
        // auto const& NT = Np;

        // auto const T = NT.dot(temperature);
        // auto const pGR = Np.dot(gas_pressure);
        // auto const pCap = Np.dot(capillary_pressure);
        // auto const pLR = pGR - pCap;

        // MPL::VariableArray vars;
        // vars[static_cast<int>(MPL::Variable::temperature)] = T;
        // vars[static_cast<int>(MPL::Variable::phase_pressure)] = pGR;
        // vars[static_cast<int>(MPL::Variable::capillary_pressure)] = pCap;
        // vars[static_cast<int>(MPL::Variable::liquid_phase_pressure)] =
        // pLR;

        ip_data.pushBackState();
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
        K.template block<W_size, temperature_size>(W_index, temperature_index);

    // pointer-matrices to the mass matrix - temperature equation
    auto MTpG = M.template block<temperature_size, gas_pressure_size>(
        temperature_index, gas_pressure_index);
    auto MTpC = M.template block<temperature_size, capillary_pressure_size>(
        temperature_index, capillary_pressure_index);
    auto MTT = M.template block<temperature_size, temperature_size>(
        temperature_index, temperature_index);
    auto MTu = M.template block<temperature_size, displacement_size>(
        temperature_index, displacement_index);

    // // pointer-matrices to the stiffness matrix - liquid phase equation
    // auto KLpG = K.template block<capillary_pressure_size,
    // gas_pressure_size>(
    //     capillary_pressure_index, gas_pressure_index);
    // auto KLpC =
    //     K.template block<capillary_pressure_size,
    //     capillary_pressure_size>(
    //         capillary_pressure_index, capillary_pressure_index);

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

        auto const& NpT = Np.transpose().eval();
        auto const& NTT = NT.transpose().eval();

        auto const& gradNp = ip_data.dNdx_p;
        auto const& gradNT = gradNp;
        auto const& gradNu = ip_data.dNdx_u;

        auto const& gradNpT = gradNp.transpose().eval();
        auto const& gradNTT = gradNT.transpose().eval();

        auto const& Nu_op = ip_data.N_u_op;
        auto const& w = ip_data.integration_weight;

        auto const& m = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>::identity2;

        auto const mT = m.transpose().eval();

        auto const x_coord =
            interpolateXCoordinate<ShapeFunctionDisplacement,
                                   ShapeMatricesTypeDisplacement>(_element, Nu);

        auto const Bu =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                gradNu, Nu, x_coord, _is_axially_symmetric);

        auto const BuT = Bu.transpose().eval();

        auto& eps = ip_data.eps;
        auto const& sigma_eff = ip_data.sigma_eff;

        double const T0 = _process_data.reference_temperature(t, pos)[0];

        double const T = NT.dot(temperature);
        double const pGR = Np.dot(gas_pressure);
        double const pCap = Np.dot(capillary_pressure);
        double const pLR = pGR - pCap;

        GlobalDimVectorType const gradpGR = gradNp * gas_pressure;
        GlobalDimVectorType const gradpCap = gradNp * capillary_pressure;

        double const T_dot = NT.dot(temperature_dot);
        double const pGR_dot = Np.dot(gas_pressure_dot);
        double const pCap_dot = Np.dot(capillary_pressure_dot);

        double const div_u = mT * Bu * displacement;
        double const div_u_dot = mT * Bu * displacement_dot;

#define nDEBUG_OUTPUT
#ifdef DEBUG_OUTPUT
#define UNKNOWNS
#define SHAPE_MATRICES
#define MATERIAL_PROPERTIES
#endif

        bool output = ((ip == 0) && (_element.getID() == 0)) ? 0 : 0;

#ifndef UNKNOWNS
        if (output)
        {
            std::cout << "#################################################"
                         "#####\n";
            std::cout << "--- unknowns: ---\n";
            std::cout << "pGR: " << gas_pressure << "\n";
            std::cout << "pCap: " << capillary_pressure << "\n";
            std::cout << "T: " << temperature << "\n";
            std::cout << "u: " << displacement << "\n";
            std::cout << "--------------------\n";
            std::cout << "pGR_dot: " << gas_pressure_dot << "\n";
            std::cout << "pCap_dot: " << capillary_pressure_dot << "\n";
            std::cout << "T_dot: " << temperature_dot << "\n";
            std::cout << "u_dot: " << displacement_dot << "\n";
            std::cout << "--- unknowns (IP): ---\n";
            std::cout << "pGR: " << pGR << "\n";
            std::cout << "pCap: " << pCap << "\n";
            std::cout << "T: " << T << "\n";
            std::cout << "--------------------\n";
            std::cout << "pGR_dot: " << pGR_dot << "\n";
            std::cout << "pCap_dot: " << pCap_dot << "\n";
            std::cout << "T_dot: " << T_dot << "\n";
        }
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
        // MPL::VariableArray vars;
        // vars[static_cast<int>(MPL::Variable::temperature)] = T;
        // vars[static_cast<int>(MPL::Variable::phase_pressure)] = pGR;
        // vars[static_cast<int>(MPL::Variable::capillary_pressure)] = pCap;
        // vars[static_cast<int>(MPL::Variable::liquid_phase_pressure)] = pLR;

        // std::cout <<
        // "#####################################################\n";
        // {
        //     auto const nComponents = gas_phase.numberOfComponents();
        //     PRINT(nComponents);

        //     PRINT(pGR);
        //     PRINT(pCap);
        //     PRINT(pLR);
        //     PRINT(T);

        //     auto const xnG =
        //         gas_phase.property(MPL::PropertyType::mole_fraction)
        //             .template value<Eigen::Vector2d>(vars, pos, t, dt);
        //     PRINT(xnG[0]);
        //     PRINT(xnG[1]);

        //     vars[static_cast<int>(MPL::Variable::mole_fraction)] = xnG[0];

        //     auto const MG = gas_phase.property(MPL::PropertyType::molar_mass)
        //                         .template value<double>(vars, pos, t, dt);
        //     PRINT(MG);

        //     auto const pCGR = xnG[0] * pGR;
        //     auto const pWGR = xnG[1] * pGR;

        //     PRINT(pCGR);
        //     PRINT(pWGR);

        //     auto rhoGR = gas_phase.property(MPL::PropertyType::density)
        //                      .template value<double>(vars, pos, t, dt);
        //     PRINT(rhoGR);

        //     auto drhoGR_dpGR =
        //         gas_phase.property(MPL::PropertyType::density)
        //             .template dValue<double>(
        //                 vars, MPL::Variable::phase_pressure, pos, t, dt);
        //     PRINT(drhoGR_dpGR);

        //     auto beta_pGR = 1. / rhoGR * drhoGR_dpGR;
        //     PRINT(beta_pGR);

        //     auto drhoGR_dT =
        //         gas_phase.property(MPL::PropertyType::density)
        //             .template dValue<double>(vars,
        //             MPL::Variable::temperature,
        //                                      pos, t, dt);
        //     PRINT(drhoGR_dT);

        //     auto beta_TGR = 1. / rhoGR * drhoGR_dT;
        //     PRINT(beta_TGR);

        //     auto pWVap = 0.;
        //     auto dpWVap_dT = 0.;

        //     Eigen::Vector2d dxm_dpGR = {0., 0.};
        //     Eigen::Vector2d xmG = {1., 0.};

        //     if (nComponents == 2)
        //     {
        //         pWVap =
        //             gas_phase.component(1)
        //                 .property(MPL::PropertyType::saturation_vapour_pressure)
        //                 .template value<double>(vars, pos, t, dt);
        //         dpWVap_dT =
        //             gas_phase.component(1)
        //                 .property(MPL::PropertyType::saturation_vapour_pressure)
        //                 .template dValue<double>(
        //                     vars, MPL::Variable::temperature, pos, t, dt);
        //     }

        //     if (nComponents == 2)
        //     {
        //         xmG[0] = xnG[0] *
        //                  gas_phase.component(0)
        //                      .property(MPL::PropertyType::molar_mass)
        //                      .template value<double>(vars, pos, t, dt) /
        //                  MG;
        //         xmG[1] = 1. - xmG[0];
        //     }
        //     else
        //     {
        //         xmG[0] = 1.;
        //         xmG[1] = 0.;
        //     }

        //     PRINT(xmG[0]);
        //     PRINT(xmG[1]);

        //     auto const dxmWG_dpGR = xmG[1] * beta_pGR;
        //     auto const dxmCG_dpGR = -dxmWG_dpGR;
        //     PRINT(dxmWG_dpGR);
        //     PRINT(dxmCG_dpGR);

        //     PRINT(pWVap);
        //     PRINT(dpWVap_dT);

        //     auto rhoCGR = xmG[0] * rhoGR;
        //     auto rhoWGR = xmG[1] * rhoGR;
        //     PRINT(rhoCGR);
        //     PRINT(rhoWGR);

        //     auto drhoWGR_dT = 0.;
        //     if (pWGR != 0)
        //     {
        //         drhoWGR_dT = rhoWGR * (1. / pWGR * dpWVap_dT - 1. / T);
        //     }
        //     PRINT(drhoWGR_dT);

        //     auto const dxmWG_dT = 1. / rhoGR * drhoWGR_dT - xmG[1] *
        //     beta_TGR; auto const dxmCG_dT = -dxmWG_dT; PRINT(dxmCG_dT);
        //     PRINT(dxmWG_dT);

        //     auto H = 0.;
        //     auto dH_dT = 0.;
        //     if (liquid_phase.numberOfComponents() == 2)
        //     {
        //         H = liquid_phase.component(0)
        //                 .property(MPL::PropertyType::henry_coefficient)
        //                 .template value<double>(vars, pos, t, dt);
        //         dH_dT = liquid_phase.component(0)
        //                     .property(MPL::PropertyType::henry_coefficient)
        //                     .template dValue<double>(
        //                         vars, MPL::Variable::temperature, pos, t,
        //                         dt);
        //     }
        //     PRINT(H);
        //     PRINT(dH_dT);

        //     auto cCL = 0.;
        //     auto dcCL_dpLR = 0.;
        //     if (liquid_phase.hasProperty(MPL::PropertyType::concentration))
        //     {
        //         cCL = liquid_phase.property(MPL::PropertyType::concentration)
        //                   .template value<double>(vars, pos, t,
        //                                           dt);  // in mol*m^(-3)
        //         dcCL_dpLR =
        //             liquid_phase.property(MPL::PropertyType::concentration)
        //                 .template dValue<double>(
        //                     vars, MPL::Variable::liquid_phase_pressure, pos,
        //                     t, dt);  // in mol*m^(-3)
        //     }
        //     PRINT(cCL);
        //     PRINT(dcCL_dpLR);

        //     auto dcCL_dT = dH_dT * xnG[0] * pGR - H * dpWVap_dT;
        //     PRINT(dcCL_dT);

        //     vars[static_cast<int>(MPL::Variable::concentration)] = cCL;

        //     auto rhoLR = liquid_phase.property(MPL::PropertyType::density)
        //                      .template value<double>(vars, pos, t,
        //                                              dt);  // in mol*m^(-3)

        //     vars[static_cast<int>(MPL::Variable::concentration)] = 0.;
        //     auto rhoWLR = liquid_phase.property(MPL::PropertyType::density)
        //                       .template value<double>(vars, pos, t, dt);

        //     auto rhoCLR = rhoLR - rhoWLR;

        //     PRINT(rhoLR);
        //     PRINT(rhoWLR);
        //     PRINT(rhoCLR);

        //     const auto xmWL = rhoWLR / rhoLR;
        //     const auto xmCL = 1. - xmWL;

        //     PRINT(xmWL);
        //     PRINT(xmCL);

        //     // HACK!! int the linear EOS of the liquid, the third
        //     // 'independent' variable is actually pressure-dependent, thus
        //     // not independent. The exact derivative would reed
        //     // drho_dpLR=rho_ref*(betaP+betaC*dC_dpLR). In order to get the
        //     // parameters of the eos (rho_ref, betaP, and betaC), I firstly
        //     // calculate the pressure derivative and then the concenctration
        //     // derivative only.
        //     const auto rho_ref_betaP =
        //         liquid_phase.property(MPL::PropertyType::density)
        //             .template dValue<double>(
        //                 vars, MPL::Variable::liquid_phase_pressure, pos, t,
        //                 dt);
        //     PRINT(rho_ref_betaP);
        //     const auto rho_ref_betaC =
        //         liquid_phase.property(MPL::PropertyType::density)
        //             .template dValue<double>(vars,
        //             MPL::Variable::concentration,
        //                                      pos, t, dt);
        //     PRINT(rho_ref_betaC);

        //     const auto drhoLR_dpLR = rho_ref_betaP + rho_ref_betaC * H;
        //     PRINT(drhoLR_dpLR);

        //     const auto rho_ref_betaT =
        //         liquid_phase.property(MPL::PropertyType::density)
        //             .template dValue<double>(vars,
        //             MPL::Variable::temperature,
        //                                      pos, t, dt);
        //     PRINT(rho_ref_betaT);

        //     const auto drhoLR_dT = rho_ref_betaT + rho_ref_betaC * dcCL_dT;
        //     PRINT(drhoLR_dT);

        //     auto drhoWLR_dpLR = rho_ref_betaP;
        //     PRINT(drhoWLR_dpLR);

        //     const auto drhoWLR_dT = rho_ref_betaT;
        //     PRINT(drhoWLR_dT);

        //     const auto dxmWL_dpLR =
        //         1. / rhoLR * (drhoWLR_dpLR - xmWL * drhoLR_dpLR);
        //     const auto dxmCL_dpLR = -dxmWL_dpLR;
        //     PRINT(dxmWL_dpLR);
        //     PRINT(dxmCL_dpLR);

        //     auto const dxmWL_dT = 1. / rhoLR * (drhoWLR_dT - xmWL *
        //     drhoLR_dT); auto const dxmCL_dT = -dxmWL_dT; PRINT(dxmWL_dT);
        //     PRINT(dxmCL_dT);

        //     OGS_FATAL("stop.");
        // }

        bool phase_transition = _process_data.phase_transition;

        auto& rho_GR = ip_data.rho_GR;  // gas phase density
        auto& rho_C_GR = ip_data.rho_C_GR;
        auto& rho_W_GR = ip_data.rho_W_GR;

        auto& rho_LR = ip_data.rho_LR;  // liquid phase density
        auto& rho_W_LR = ip_data.rho_W_LR;
        auto& rho_C_LR = ip_data.rho_C_LR;

        auto& xm_C_G = ip_data.xmCG;  // C constituent gas phase mass fraction
        auto& xm_W_G = ip_data.xmWG;  // W constituent gas phase mass fraction

        auto& xm_C_L =
            ip_data.xmCL;  // C constituent liquid phase mass fraction
        auto& xm_W_L =
            ip_data.xmWL;  // W constituent liquid phase mass fraction

        auto& dxm_C_G_dpGR = ip_data.dxmCG_dpGR;
        auto& dxm_W_G_dpGR = ip_data.dxmWG_dpGR;
        auto& dxm_C_L_dpLR = ip_data.dxmCL_dpLR;
        auto& dxm_W_L_dpLR = ip_data.dxmWL_dpLR;
        auto& dxm_C_G_dT = ip_data.dxmCG_dT;
        auto& dxm_W_G_dT = ip_data.dxmWG_dT;
        auto& dxm_C_L_dT = ip_data.dxmCL_dT;
        auto& dxm_W_L_dT = ip_data.dxmWL_dT;

        //  - solid phase properties
        auto& beta_pS = ip_data.beta_pS;
        auto& beta_T_SR = ip_data.beta_T_SR;
        auto& c_p_S = ip_data.c_p_S;
        auto& lambda_SR = ip_data.lambdaSR;

        //  - gas phase properties
        auto& mu_GR = ip_data.muGR;
        auto& lambda_GR = ip_data.lambdaGR;

        //  - liquid phase properties
        auto& mu_LR = ip_data.muLR;
        auto& lambda_LR = ip_data.lambdaLR;

        // constant liquid constituent enthalpies ftw!
        // const auto h_W_L = 112654.;  // J/kg
        // const auto h_C_L = 648584.;  // J/kg

        const auto I =
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Identity();
        const double sD_G = 0.1;  // or whatever
        const double sD_L = 0.1;

        const auto D_C_G = (sD_G * I).eval();  // I know, I know...
        const auto D_W_G = (sD_G * I).eval();
        const auto D_C_L = (sD_L * I).eval();
        const auto D_W_L = (sD_L * I).eval();

        //  - medium properties
        MPL::VariableArray vars;
        // auto const k_S = MPL::formEigenTensor<DisplacementDim>(
        //     medium.property(MPL::PropertyType::permeability)
        //         .value(vars, pos, t, dt));
        auto& k_S = ip_data.k_S;

        auto& s_L = ip_data.saturation;
        auto const s_G = 1. - s_L;
        const auto s_L_dot = (s_L - ip_data.saturation_prev) / dt;

        // vars[static_cast<int>(MPL::Variable::liquid_saturation)] = s_L;

        auto& alpha_B = ip_data.alpha_B;
        auto& beta_p_SR = ip_data.beta_p_SR;

        auto const k_rel = ip_data.kRel;
        auto const k_rel_L = k_rel[0];
        auto const k_rel_G = k_rel[1];

        auto const& b = _process_data.specific_body_force;

        // pore fluid pressure
        const double p_FR = s_G * pGR + s_L * pLR;
        const double p_FR_dot = pGR_dot - (s_L_dot * pCap + pCap_dot * s_L);

        auto& phi = ip_data.phi;
        // porosity derivative
        auto const phi_dot = (alpha_B - phi) * (div_u_dot - beta_T_SR * T_dot +
                                                beta_p_SR * p_FR_dot);
        // cf. Eq. 86
        // porosity update
        //      phi = phi_dot * dt + ip_data.phi_prev;

        auto const phi_G = s_G * phi;
        auto const phi_L = s_L * phi;
        auto const phi_S = 1. - phi;

        // TODO (Grunwald) replace effective thermal conductivity by a more
        // sophisticated law, maybe even allow the law to be chosen in the
        // project file as medium property
        auto const lambda = MPL::formEigenTensor<DisplacementDim>(
            phi_S * lambda_SR + phi_L * lambda_LR + phi_G * lambda_GR);

        // TODO (Wenqing) : Change dT to time step wise increment
        // double const delta_T(T - T0);
        // double const thermal_strain = beta_T_SR * delta_T;
        // prevent rho_SR and p_SR being killed by crazy T

        auto& rho_SR = ip_data.rho_SR;
        auto const rho = phi_G * rho_GR + phi_L * rho_LR + phi_S * rho_SR;

        // TODO: change back to rho_SR
        // auto const rho_c_p = phi_G * rho_GR * c_p_G + phi_L * rho_LR *
        // c_p_L
        // +
        //                      phi_S * rho_ref_SR * c_p_S;

        // update secondary variables. TODO: Refactoring potential!!
        _liquid_pressure[ip] = pLR;
        _liquid_density[ip] = rho_LR;
        _gas_density[ip] = rho_GR;
        _porosity[ip] = phi;
        _saturation[ip] = s_L;
        _mole_fraction_gas[ip] = ip_data.xnCG;
        _mass_fraction_gas[ip] = ip_data.xmCG;
        _mass_fraction_liquid[ip] = ip_data.xmWL;

        // abbreviations
        const double rho_C_FR = s_G * rho_C_GR + s_L * rho_C_LR;
        const double rho_W_FR = s_G * rho_W_GR + s_L * rho_W_LR;

        // solid phase pressure
        auto& p_SR = ip_data.p_SR;
        // const double p_SR_dot =
        //     p_FR_dot -
        //     1. / (beta_pS * phi_S) * (div_u_dot - beta_T_SR * T_dot) -
        //     phi_dot / phi_S * (p_FR - p_SR);

        auto& phi_S_p_SR = ip_data.phi_S_p_SR;
        //        phi_S_p_SR = phi_S * p_SR;

        const double phi_S_p_SR_dot =
            (phi_S_p_SR - ip_data.phi_S_p_SR_prev) / dt;

        // phase specific enthalpies
        auto& h_G = ip_data.h_G;
        auto& h_L = ip_data.h_L;
        auto& h_S = ip_data.h_S;

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
        // cf. Eq 101

        const auto k_over_mu_G = k_S * k_rel_G / mu_GR;
        const auto k_over_mu_L = k_S * k_rel_L / mu_LR;

        GlobalDimVectorType const w_GS =
            k_over_mu_G * rho_GR * b - k_over_mu_G * gradpGR;

        GlobalDimVectorType const w_LS = k_over_mu_L * gradpCap +
                                         k_over_mu_L * rho_GR * b -
                                         k_over_mu_L * gradpGR;

#ifndef MATERIAL_PROPERTIES
        if (output)
        {
            std::cout << "#################################################"
                         "#####\n";
            std::cout << "#    Material properties:\n";
            std::cout << "#------------------------------------------------"
                         "----#\n";
            std::cout << "#         rho_GR:  " << rho_GR << "\n";
            std::cout << "#       rho_C_GR:  " << rho_C_GR << "\n";
            std::cout << "#       rho_W_GR:  " << rho_W_GR << "\n";
            std::cout << "#    .    .    .    .    .    .    .    .    .   "
                         " .  #\n";
            std::cout << "#   rho_C_GR_dot:  " << rho_C_GR_dot << "\n";
            std::cout << "#   rho_W_GR_dot:  " << rho_W_GR_dot << "\n";
            std::cout << "#    .    .    .    .    .    .    .    .    .   "
                         " .  #\n";
            std::cout << "       rho_C_FR : " << rho_C_FR << "\n";
            std::cout << "#    .    .    .    .    .    .    .    .    .   "
                         " .  #\n";
            std::cout << "#         xm_C_G:  " << xm_C_G << "\n";
            std::cout << "#         xm_W_G:  " << xm_W_G << "\n";
            std::cout << "#    .    .    .    .    .    .    .    .    .   "
                         " .  #\n";
            std::cout << "#  d_xm_C_G_dpGR:  " << dxm_C_G_dpGR << "\n";
            std::cout << "#  d_xm_W_G_dpGR:  " << dxm_W_G_dpGR << "\n";
            std::cout << "#    d_xm_C_G_dT:  " << dxm_C_G_dT << "\n";
            std::cout << "#    d_xm_W_G_dT:  " << dxm_W_G_dT << "\n";
            std::cout << "#------------------------------------------------"
                         "----#\n";
            std::cout << "#         rho_LR:  " << rho_LR << "\n";
            std::cout << "#       rho_C_LR:  " << rho_C_LR << "\n";
            std::cout << "#       rho_W_LR:  " << rho_W_LR << "\n";
            std::cout << "#    .    .    .    .    .    .    .    .    .   "
                         " .  #\n";
            std::cout << "#   rho_C_LR_dot:  " << rho_C_LR_dot << "\n";
            std::cout << "#   rho_W_LR_dot:  " << rho_W_LR_dot << "\n";
            std::cout << "#    .    .    .    .    .    .    .    .    .   "
                         " .  #\n";
            std::cout << "       rho_W_FR : " << rho_W_FR << "\n";
            std::cout << "#    .    .    .    .    .    .    .    .    .   "
                         " .  #\n";
            std::cout << "#         xm_C_L:  " << xm_C_L << "\n";
            std::cout << "#         xm_W_L:  " << xm_W_L << "\n";
            std::cout << "#    .    .    .    .    .    .    .    .    .   "
                         " .  #\n";
            std::cout << "#  d_xm_C_L_dpLR:  " << dxm_C_L_dpLR << "\n";
            std::cout << "#  d_xm_W_L_dpLR:  " << dxm_W_L_dpLR << "\n";
            std::cout << "#    d_xm_C_L_dT:  " << dxm_C_L_dT << "\n";
            std::cout << "#    d_xm_W_L_dT:  " << dxm_W_L_dT << "\n";
            std::cout << "#------------------------------------------------"
                         "----#\n";
            std::cout << "#         rho_SR:  " << rho_SR << "\n";
            std::cout << "#            rho:  " << rho << "\n";
            std::cout << "#------------------------------------------------"
                         "----#\n";
            // std::cout << "#          c_p_G:  " << c_p_G << "\n";
            // std::cout << "#          c_p_L:  " << c_p_L << "\n";
            // std::cout << "#          c_p_S:  " << c_p_S << "\n";
            std::cout << "#------------------------------------------------"
                         "----#\n";
            std::cout << "#      beta_p_SR:  " << beta_p_SR << "\n";
            std::cout << "#      beta_T_SR:  " << beta_T_SR << "\n";
            std::cout << "#------------------------------------------------"
                         "----#\n";
            std::cout << "#            h_G:  " << h_G << "\n";
            std::cout << "#            h_L:  " << h_L << "\n";
            std::cout << "#            h_S:  " << h_S << "\n";
            std::cout << "#    .    .    .    .    .    .    .    .    .   "
                         " .  #\n";
            // std::cout << "#          h_C_G:  " << h_C_G << "\n";
            // std::cout << "#          h_W_G:  " << h_W_G << "\n";
            // std::cout << "#          h_C_L:  " << h_C_L << "\n";
            // std::cout << "#          h_W_L:  " << h_W_L << "\n";
            std::cout << "#    .    .    .    .    .    .    .    .    .   "
                         " .  #\n";
            std::cout << "#      rho_h_eff:  " << rho_h_eff << "\n";
            std::cout << "#      rho_G_h_G:  " << ip_data.rho_G_h_G << "\n";
            std::cout << "#      rho_L_h_L:  " << ip_data.rho_L_h_L << "\n";
            std::cout << "#      rho_S_h_S:  " << ip_data.rho_S_h_S << "\n";
            std::cout << "#    .    .    .    .    .    .    .    .    .   "
                         " .  #\n";
            std::cout << "#  rho_G_h_G_dot:  " << rho_G_h_G_dot << "\n";
            std::cout << "#  rho_L_h_L_dot:  " << rho_L_h_L_dot << "\n";
            std::cout << "#  rho_S_h_S_dot:  " << rho_S_h_S_dot << "\n";
            std::cout << "#------------------------------------------------"
                         "----#\n";
            std::cout << "#            k_S:\n" << k_S << "\n";
            std::cout << "#------------------------------------------------"
                         "----#\n";
            std::cout << "#        k_rel_L:  " << k_rel_L << "\n";
            std::cout << "#        k_rel_G:  " << k_rel_G << "\n";
            std::cout << "#          mu_GR:  " << mu_GR << "\n";
            std::cout << "#          mu_LR:  " << mu_LR << "\n";
            std::cout << "#    k_over_mu_G:\n" << k_over_mu_G << "\n";
            std::cout << "#    k_over_mu_L:\n" << k_over_mu_L << "\n";
            std::cout << "#------------------------------------------------"
                         "----#\n";
            std::cout << "#            s_L:  " << s_L << "\n";
            std::cout << "#            s_G:  " << s_G << "\n";
            std::cout << "#        s_L_dot:  " << s_L_dot << "\n";
            std::cout << "#------------------------------------------------"
                         "----#\n";
            std::cout << "#        alpha_B:  " << alpha_B << "\n";
            std::cout << "#        beta_p_SR:  " << beta_p_SR << "\n";
            std::cout << "#        beta_T_SR:  " << beta_T_SR << "\n";
            std::cout << "#------------------------------------------------"
                         "----#\n";
            std::cout << "#         lambda:\n" << lambda << "\n";
            std::cout << "#------------------------------------------------"
                         "----#\n";
            std::cout << "#            phi:  " << phi << "\n";
            std::cout << "#          phi_G:  " << phi_G << "\n";
            std::cout << "#          phi_L:  " << phi_L << "\n";
            std::cout << "#          phi_S:  " << phi_S << "\n";
            std::cout << "#################################################"
                         "#####\n";
            std::cout << "#  Darcy-Velocities: "
                      << "\n";
            std::cout << "#------------------------------------------------"
                         "----#\n";
            std::cout << "#           w_LS:\n" << w_LS << "\n";
            std::cout << "#           w_GS:\n" << w_GS << "\n";
            std::cout << "#################################################"
                         "#####\n";
        }
        // OGS_FATAL("It has to stop!");
#endif
        PRINT2(output, --------------------------------------, 0);
        PRINT2(output, Component C equation, 0);

        // coefficient matrices
        // C-component equation
        MCpG.noalias() += NpT * rho_C_FR * (alpha_B - phi) * beta_p_SR * Np * w;
        PRINT2(output, MCpG, rho_C_FR * (alpha_B - phi) * beta_p_SR);
        // std::cout << MCpG << "\n";
        MCpC.noalias() -=
            NpT * rho_C_FR * (alpha_B - phi) * beta_p_SR * s_L * Np * w;

        if (pCap_dot != 0.)  // avoid division by Zero
        {
            MCpC.noalias() += NpT *
                              (phi * (rho_C_LR - rho_C_GR) -
                               rho_C_FR * pCap * (alpha_B - phi) * beta_p_SR) *
                              s_L_dot / pCap_dot * Np * w;
        }

        PRINT2(output, MCpC, rho_C_FR * (alpha_B - phi) * beta_p_SR * s_L);
        // std::cout << MCpC << "\n";

        MCT.noalias() -= NpT * rho_C_FR * (alpha_B - phi) * beta_T_SR * Np * w;
        PRINT2(output, MCT, rho_C_FR * (alpha_B - phi) * beta_T_SR);

        MCu.noalias() += NpT * rho_C_FR * alpha_B * mT * Bu * w;
        PRINT2(output, MCu, rho_C_FR * alpha_B);

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
        PRINT2(output, LCpG, (advection_C + diffusion_C_p));

        LCpC.noalias() -=
            gradNpT * (advection_C_L + diffusion_C_L_p) * gradNp * w;
        PRINT2(output, LCpC, (advection_C_L + diffusion_C_L_p));

        LCT.noalias() += gradNpT * (diffusion_C_T)*gradNp * w;
        PRINT2(output, LCT, (diffusion_C_T));

        fC.noalias() +=
            gradNpT * (advection_C_G * rho_GR + advection_C_L * rho_LR) * b * w;
        PRINT2(output, fC_I,
               (advection_C_G * rho_GR + advection_C_L * rho_LR) * b);

        // fC.noalias() -= NpT *
        //                 (phi * (rho_C_LR - rho_C_GR) -
        //                  rho_C_FR * pCap * (alpha_B - phi) * beta_p_SR) *
        //                 s_L_dot * w;

        PRINT2(output, fC_II,
               (phi * (rho_C_LR - rho_C_GR) -
                rho_C_FR * pCap * (alpha_B - phi) * beta_p_SR) *
                   s_L_dot);
        fC.noalias() -=
            NpT * phi * (s_G * rho_C_GR_dot + s_L * rho_C_LR_dot) * w;
        PRINT2(output, fC_III, phi * (s_G * rho_C_GR_dot + s_L * rho_C_LR_dot));

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        // W-component equation
        PRINT2(output, --------------------------------------, 0);
        PRINT2(output, Component W equation, 0);

        MWpG.noalias() += NpT * rho_W_FR * (alpha_B - phi) * beta_p_SR * Np * w;
        PRINT2(output, MWpG, rho_W_FR * (alpha_B - phi) * beta_p_SR * Np);
        // std::cout << MWpG << "\n";
        MWpC.noalias() -=
            NpT * rho_W_FR * (alpha_B - phi) * beta_p_SR * s_L * Np * w;
        PRINT2(output, MWpC, rho_W_FR * (alpha_B - phi) * beta_p_SR * s_L);
        // std::cout << MWpC << "\n";
        // OGS_FATAL("sdfdasf");

        MWT.noalias() -= NpT * rho_W_FR * (alpha_B - phi) * beta_T_SR * Np * w;
        PRINT2(output, MWT, rho_W_FR * (alpha_B - phi) * beta_T_SR);

        MWu.noalias() += NpT * rho_W_FR * alpha_B * mT * Bu * w;
        PRINT2(output, MWu, rho_W_FR * alpha_B);

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
        PRINT2(output, LWpG, (advection_W + diffusion_W_p));

        LWpC.noalias() -=
            gradNpT * (advection_W_L + diffusion_W_L_p) * gradNp * w;
        PRINT2(output, LWpC, (advection_W_L + diffusion_W_L_p));

        LWT.noalias() += gradNpT * (diffusion_W_T)*gradNp * w;
        PRINT2(output, LWT, (diffusion_W_T));

        fW.noalias() +=
            gradNpT * (advection_W_G * rho_GR + advection_W_L * rho_LR) * b * w;
        PRINT2(output, fW_III,
               (advection_W_G * rho_GR + advection_W_L * rho_LR) * b);

        fW.noalias() -= NpT *
                        (phi * (rho_W_LR - rho_W_GR) -
                         rho_W_FR * pCap * (alpha_B - phi) * beta_p_SR) *
                        s_L_dot * w;
        PRINT2(output, fW_III,
               (phi * (rho_W_LR - rho_W_GR) -
                rho_W_FR * pCap * (alpha_B - phi) * beta_p_SR) *
                   s_L_dot);

        fW.noalias() -=
            NpT * phi * (s_G * rho_W_GR_dot + s_L * rho_W_LR_dot) * w;
        PRINT2(output, fW_III, phi * (s_G * rho_W_GR_dot + s_L * rho_W_LR_dot));

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        //  - temperature equation
        PRINT2(output, --------------------------------------, 0);
        PRINT2(output, Temperature equation, 0);

        MTpG.noalias() -= NTT * phi * NT * w;
        PRINT2(output, MTpG, phi);
        MTpC.noalias() += NTT * phi_L * NT * w;
        PRINT2(output, MTpC, phi_L);

        MTu.noalias() += NTT * rho_h_eff * mT * Bu * w;
        PRINT2(output, MTu, rho_h_eff);

        KTT.noalias() += gradNTT * lambda * gradNT * w;
        PRINT2(output, KTT, gradNTT * lambda * gradNT);

        // fT
        fT.noalias() -=
            NTT * (rho_G_h_G_dot + rho_L_h_L_dot + rho_S_h_S_dot) * w;
        PRINT2(output, fT_I, (rho_G_h_G_dot + rho_L_h_L_dot + rho_S_h_S_dot));
        fT.noalias() -=
            NTT * (phi * pCap * s_L_dot - p_FR * phi_dot - phi_S_p_SR_dot) * w;
        PRINT2(output, fT_II,
               (phi * pCap * s_L_dot - p_FR * phi_dot - phi_S_p_SR_dot));
        fT.noalias() +=
            NTT * (rho_GR * w_GS.transpose() + rho_LR * w_LS.transpose()) * b *
            w;
        PRINT2(output, fT_III,
               (rho_GR * w_GS.transpose() + rho_LR * w_LS.transpose()) * b);
        fT.noalias() +=
            gradNTT * (rho_GR * h_G * w_GS + rho_LR * h_L * w_LS) * w;
        PRINT2(output, fT_IV, (rho_GR * h_G * w_GS + rho_LR * h_L * w_LS));

        //  - displacement equation
        PRINT2(output, --------------------------------------, 0);
        PRINT2(output, Displacement equation, 0);

        KUpG.noalias() -= (BuT * alpha_B * m * Np) * w;
        PRINT2(output, KUpG, (BuT * alpha_B * m * Np));
        KUpC.noalias() += (BuT * alpha_B * s_L * m * Np) * w;
        PRINT2(output, KUpC, (BuT * alpha_B * s_L * m * Np));

        eps.noalias() = Bu * displacement;
        ip_data.updateConstitutiveRelationThermal(
            t, pos, dt, displacement,
            _process_data.reference_temperature(t, pos)[0]);

        fU.noalias() -= (BuT * sigma_eff - Nu_op.transpose() * rho * b) * w;
        PRINT2(output, fU, (BuT * sigma_eff - Nu_op.transpose() * rho * b));

        if (output)
        {
            PRINT(sigma_eff);
            PRINT(eps);
        }

        if (0 && (_process_data.apply_mass_lumping))
        {
            for (unsigned row = 0; row < MCpG.cols(); row++)
            {
                for (unsigned column = 0; column < MCpG.cols(); column++)
                {
                    if (row != column)
                    {
                        // std::cout << "row: " << row << " column: " <<
                        // column
                        //           << "\n";

                        // std::cout << MCpG(row, row) << " " << MCpG(row,
                        // column)
                        //           << "\n";

                        MCpG(row, row) += MCpG(row, column);
                        MCpG(row, column) = 0.0;

                        // std::cout << "MCpCG:\n";
                        // std::cout << MCpC(row, row) << " " << MCpC(row,
                        // column)
                        //           << "\n";

                        MCpC(row, row) += MCpC(row, column);
                        MCpC(row, column) = 0.0;

                        // std::cout << "MWpG:\n";
                        // std::cout << MWpG(row, row) << " " << MWpG(row,
                        // column)
                        //           << "\n";

                        MWpG(row, row) += MWpG(row, column);
                        MWpG(row, column) = 0.0;

                        // std::cout << "MWpC:\n";
                        // std::cout << MWpC(row, row) << " " << MWpC(row,
                        // column)
                        //           << "\n";

                        MWpC(row, row) += MWpC(row, column);
                        MWpC(row, column) = 0.0;
                    }
                }
            }
            // auto MpG = M.template block<C_size, gas_pressure_size>(
            //     C_index, gas_pressure_index);
            // MpG = MpG.colwise().sum().eval().asDiagonal();
            // auto MpC = M.template block<W_size, capillary_pressure_size>(
            //     W_index, capillary_pressure_index);
            // MpC = MpC.colwise().sum().eval().asDiagonal();
        }

        // std::cout << "M (post_lumping):\n" << M << "\n";

        // OGS_FATAL("So.");

#ifdef DEBUG_OUTPUT
        std::cout << "------------------------------------------------------\n";
        std::cout << " MCpG:\n" << MCpG << "\n";

        std::cout << " MCpG \n" << MCpG << "\n";
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

        OGS_FATAL("Intentional stop.");
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

        // TODO (naumov) Temporary value not used by current material
        // models. Need extension of secondary variables interface.
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

        // TODO (naumov) Temporary value not used by current material
        // models. Need extension of secondary variables interface.
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
        auto& ip_data = _ip_data[ip];

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
            _process_data.reference_temperature(t, pos)[0]);

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
        _mole_fraction_gas[ip] = ip_data.xnCG;
        _mass_fraction_gas[ip] = ip_data.xmCG;
        _mass_fraction_liquid[ip] = ip_data.xmWL;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void TH2MLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                        IntegrationMethod, DisplacementDim>::
    computeSecondaryVariableConcrete(double const t,
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

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    // TODO: get dt!
    getConstitutiveVariables(local_x, t, 0);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = _ip_data[ip];

        auto const& Np = ip_data.N_p;

        double const pGR = Np.dot(gas_pressure);
        double const pCap = Np.dot(capillary_pressure);
        double const pLR = pGR - pCap;

        _liquid_pressure[ip] = pLR;
        _liquid_density[ip] = ip_data.rho_LR;
        _gas_density[ip] = ip_data.rho_GR;
        _saturation[ip] = ip_data.saturation;
        _mole_fraction_gas[ip] = ip_data.xnCG;
        _mass_fraction_gas[ip] = ip_data.xmCG;
        _mass_fraction_liquid[ip] = ip_data.xmWL;
    }
}

}  // namespace TH2M
}  // namespace ProcessLib