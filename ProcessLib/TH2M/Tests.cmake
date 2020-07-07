AddTest(
    NAME TH2M_T
    PATH TH2M/T
    EXECUTABLE ogs
    EXECUTABLE_ARGS T_unitSquare.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB th2m_T_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB th2m_T_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 0
    GLOB th2m_T_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 0
    GLOB th2m_T_pcs_0_ts_*.vtu saturation saturation 1e-15 0
    GLOB th2m_T_pcs_0_ts_*.vtu temperature temperature 1e-15 0
    GLOB th2m_T_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB th2m_T_pcs_0_ts_*.vtu porosity porosity 1e-15 0
    GLOB th2m_T_pcs_0_ts_*.vtu sigma sigma 5e-10 0
)

AddTest(
    NAME TH2M_TH
    PATH TH2M/TH
    EXECUTABLE ogs
    EXECUTABLE_ARGS TH_unitSquare.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB th2m_TH_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB th2m_TH_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 0
    GLOB th2m_TH_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 0
    GLOB th2m_TH_pcs_0_ts_*.vtu saturation saturation 1e-15 0
    GLOB th2m_TH_pcs_0_ts_*.vtu temperature temperature 1e-15 0
    GLOB th2m_TH_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB th2m_TH_pcs_0_ts_*.vtu porosity porosity 1e-15 0
    GLOB th2m_TH_pcs_0_ts_*.vtu sigma sigma 5e-10 0
)


AddTest(
    NAME TH2M_M
    PATH TH2M/M
    EXECUTABLE ogs
    EXECUTABLE_ARGS M_unitSquare.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB th2m_M_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB th2m_M_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 0
    GLOB th2m_M_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 0
    GLOB th2m_M_pcs_0_ts_*.vtu saturation saturation 1e-15 0
    GLOB th2m_M_pcs_0_ts_*.vtu temperature temperature 1e-15 0
    GLOB th2m_M_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB th2m_M_pcs_0_ts_*.vtu porosity porosity 1e-15 0
    GLOB th2m_M_pcs_0_ts_*.vtu sigma sigma 5e-10 0
)

AddTest(
    NAME TH2M_H2_McWhorter
    PATH TH2M/H2_McWhorter
    EXECUTABLE ogs
    EXECUTABLE_ARGS H2_unitSquare.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB th2m_H_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu saturation saturation 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu temperature temperature 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu porosity porosity 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu sigma sigma 5e-10 0
)

AddTest(
    NAME TH2M_TH2M_Liakopoulos
    PATH TH2M/TH2M_Liakopoulos
    EXECUTABLE ogs
    EXECUTABLE_ARGS TH2M_Liakopoulos.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB th2m_liakopoulos_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB th2m_liakopoulos_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 0
    GLOB th2m_liakopoulos_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 0
    GLOB th2m_liakopoulos_pcs_0_ts_*.vtu saturation saturation 1e-15 0
    GLOB th2m_liakopoulos_pcs_0_ts_*.vtu temperature temperature 1e-15 0
    GLOB th2m_liakopoulos_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB th2m_liakopoulos_pcs_0_ts_*.vtu porosity porosity 1e-15 0
    GLOB th2m_liakopoulos_pcs_0_ts_*.vtu sigma sigma 5e-10 0
)

AddTest(
    NAME TH2M_TM
    PATH TH2M/TM
    EXECUTABLE ogs
    EXECUTABLE_ARGS TM_unitSquare.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB th2m_H_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu saturation saturation 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu temperature temperature 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu porosity porosity 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu sigma sigma 5e-10 0
)


AddTest(
    NAME TH2M_H2_desaturation
    PATH TH2M/H2_desaturation
    EXECUTABLE ogs
    EXECUTABLE_ARGS H2_unitSquare_1e2.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB th2m_H2_desaturation_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB th2m_H2_desaturation_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 0
    GLOB th2m_H2_desaturation_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 0
    GLOB th2m_H2_desaturation_pcs_0_ts_*.vtu saturation saturation 1e-15 0
    GLOB th2m_H2_desaturation_pcs_0_ts_*.vtu temperature temperature 1e-15 0
    GLOB th2m_H2_desaturation_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB th2m_H2_desaturation_pcs_0_ts_*.vtu porosity porosity 1e-15 0
    GLOB th2m_H2_desaturation_pcs_0_ts_*.vtu sigma sigma 5e-10 0
)

AddTest(
    NAME TH2M_H2_desaturation_mono
    PATH TH2M/H2_desaturation_mono
    EXECUTABLE ogs
    EXECUTABLE_ARGS H2_unitSquare.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB th2m_H2_desaturation_mono_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB th2m_H2_desaturation_mono_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 0
    GLOB th2m_H2_desaturation_mono_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 0
    GLOB th2m_H2_desaturation_mono_pcs_0_ts_*.vtu saturation saturation 1e-15 0
    GLOB th2m_H2_desaturation_mono_pcs_0_ts_*.vtu temperature temperature 1e-15 0
    GLOB th2m_H2_desaturation_mono_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB th2m_H2_desaturation_mono_pcs_0_ts_*.vtu porosity porosity 1e-15 0
    GLOB th2m_H2_desaturation_mono_pcs_0_ts_*.vtu sigma sigma 5e-10 0
)

AddTest(
    NAME TH2M_TH_TM
    PATH TH2M/TH_TM
    EXECUTABLE ogs
    EXECUTABLE_ARGS THM_unitSquare.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB th2m_H_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu saturation saturation 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu temperature temperature 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu porosity porosity 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu sigma sigma 5e-10 0
)

AddTest(
    NAME TH2M_H_L
    PATH TH2M/H_L
    EXECUTABLE ogs
    EXECUTABLE_ARGS H_unitSquare.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB th2m_H_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu saturation saturation 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu temperature temperature 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu porosity porosity 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu sigma sigma 5e-10 0
)

AddTest(
    NAME TH2M_TH2M_Liakopoulos_Richards
    PATH TH2M/TH2M_Liakopoulos_Richards
    EXECUTABLE ogs
    EXECUTABLE_ARGS TH2M_Liakopoulos_Richards.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB th2m_liakopoulos_richards_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB th2m_liakopoulos_richards_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 0
    GLOB th2m_liakopoulos_richards_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 0
    GLOB th2m_liakopoulos_richards_pcs_0_ts_*.vtu saturation saturation 1e-15 0
    GLOB th2m_liakopoulos_richards_pcs_0_ts_*.vtu temperature temperature 1e-15 0
    GLOB th2m_liakopoulos_richards_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB th2m_liakopoulos_richards_pcs_0_ts_*.vtu porosity porosity 1e-15 0
    GLOB th2m_liakopoulos_richards_pcs_0_ts_*.vtu sigma sigma 5e-10 0
)

AddTest(
    NAME TH2M_H2_monoelemental
    PATH TH2M/H2_monoelemental
    EXECUTABLE ogs
    EXECUTABLE_ARGS H2_unitSquare.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB th2m_H_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu saturation saturation 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu temperature temperature 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu porosity porosity 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu sigma sigma 5e-10 0
)

AddTest(
    NAME TH2M_H_G
    PATH TH2M/H_G
    EXECUTABLE ogs
    EXECUTABLE_ARGS H_unitSquare.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB th2m_H_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu gas_pressure gas_pressure 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu capillary_pressure capillary_pressure 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu saturation saturation 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu temperature temperature 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu porosity porosity 1e-15 0
    GLOB th2m_H_pcs_0_ts_*.vtu sigma sigma 5e-10 0
)

