AddTest(
    NAME TH2M_T
    PATH TH2M/numerical_jacobian/T
    EXECUTABLE ogs
    EXECUTABLE_ARGS T.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB th2m_t_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB th2m_t_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 0
    GLOB th2m_t_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 0
    GLOB th2m_t_pcs_0_ts_*.vtu saturation saturation 1e-15 0
    GLOB th2m_t_pcs_0_ts_*.vtu temperature temperature 1e-15 0
    GLOB th2m_t_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB th2m_t_pcs_0_ts_*.vtu porosity porosity 1e-15 0
    GLOB th2m_t_pcs_0_ts_*.vtu sigma sigma 5e-10 0
)

AddTest(
    NAME TH2M_H_Gas_compressible
    PATH TH2M/numerical_jacobian/H_Gas_compressible
    EXECUTABLE ogs
    EXECUTABLE_ARGS H_Gas_compressible.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB th2m_h_g_comp_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB th2m_h_g_comp_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 0
    GLOB th2m_h_g_comp_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 0
    GLOB th2m_h_g_comp_pcs_0_ts_*.vtu saturation saturation 1e-15 0
    GLOB th2m_h_g_comp_pcs_0_ts_*.vtu temperature temperature 1e-15 0
    GLOB th2m_h_g_comp_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB th2m_h_g_comp_pcs_0_ts_*.vtu porosity porosity 1e-15 0
    GLOB th2m_h_g_comp_pcs_0_ts_*.vtu sigma sigma 5e-10 0
)

AddTest(
    NAME TH2M_H_Gas_incompressible
    PATH TH2M/numerical_jacobian/H_Gas_incompressible
    EXECUTABLE ogs
    EXECUTABLE_ARGS H_Gas_incompressible.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB th2m_h_g_incomp_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB th2m_h_g_incomp_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 0
    GLOB th2m_h_g_incomp_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 0
    GLOB th2m_h_g_incomp_pcs_0_ts_*.vtu saturation saturation 1e-15 0
    GLOB th2m_h_g_incomp_pcs_0_ts_*.vtu temperature temperature 1e-15 0
    GLOB th2m_h_g_incomp_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB th2m_h_g_incomp_pcs_0_ts_*.vtu porosity porosity 1e-15 0
    GLOB th2m_h_g_incomp_pcs_0_ts_*.vtu sigma sigma 5e-10 0
)

AddTest(
    NAME TH2M_HT_dirichlet_left
    PATH TH2M/numerical_jacobian/HT
    EXECUTABLE ogs
    EXECUTABLE_ARGS HT.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB th2m_ht_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB th2m_ht_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 0
    GLOB th2m_ht_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 0
    GLOB th2m_ht_pcs_0_ts_*.vtu saturation saturation 1e-15 0
    GLOB th2m_ht_pcs_0_ts_*.vtu temperature temperature 1e-15 0
    GLOB th2m_ht_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB th2m_ht_pcs_0_ts_*.vtu porosity porosity 1e-15 0
    GLOB th2m_ht_pcs_0_ts_*.vtu sigma sigma 5e-10 0
)

AddTest(
    NAME TH2M_H2T
    PATH TH2M/numerical_jacobian/H2T
    EXECUTABLE ogs
    EXECUTABLE_ARGS H2T.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB th2m_h2t_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB th2m_h2t_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 0
    GLOB th2m_h2t_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 0
    GLOB th2m_h2t_pcs_0_ts_*.vtu saturation saturation 1e-15 0
    GLOB th2m_h2t_pcs_0_ts_*.vtu temperature temperature 1e-15 0
    GLOB th2m_h2t_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB th2m_h2t_pcs_0_ts_*.vtu porosity porosity 1e-15 0
    GLOB th2m_h2t_pcs_0_ts_*.vtu sigma sigma 5e-10 0
)

AddTest(
    NAME TH2M_M_simple
    PATH TH2M/numerical_jacobian/M
    EXECUTABLE ogs
    EXECUTABLE_ARGS M.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB th2m_m_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB th2m_m_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 0
    GLOB th2m_m_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 0
    GLOB th2m_m_pcs_0_ts_*.vtu saturation saturation 1e-15 0
    GLOB th2m_m_pcs_0_ts_*.vtu temperature temperature 1e-15 0
    GLOB th2m_m_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB th2m_m_pcs_0_ts_*.vtu porosity porosity 1e-15 0
    GLOB th2m_m_pcs_0_ts_*.vtu sigma sigma 5e-10 0
)

AddTest(
    NAME TH2M_HM_linear
    PATH TH2M/numerical_jacobian/HM
    EXECUTABLE ogs
    EXECUTABLE_ARGS HM.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB th2m_hm_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB th2m_hm_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 0
    GLOB th2m_hm_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 0
    GLOB th2m_hm_pcs_0_ts_*.vtu saturation saturation 1e-15 0
    GLOB th2m_hm_pcs_0_ts_*.vtu temperature temperature 1e-15 0
    GLOB th2m_hm_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB th2m_hm_pcs_0_ts_*.vtu porosity porosity 1e-15 0
    GLOB th2m_hm_pcs_0_ts_*.vtu sigma sigma 5e-10 0
)