#ifndef INITIALIZE_FUNCTIONS_SOLVERLOOP_H_
#define INITIALIZE_FUNCTIONS_SOLVERLOOP_H_

void initialize_functions_solverloop(){
  if (FUNCTION_F == 1) {    
       free_energy               = function_F_01_free_energy;
       dc_dmu                    = function_F_01_dc_dmu;
       c_mu                      = function_F_01_c_mu;
       Mu                        = function_F_01_Mu;
       dpsi                      = function_F_01_dpsi;
//        function_A                = function_F_01_function_A;
       function_B                = function_F_01_function_B;
       function_C                = function_F_01_function_C;
       init_propertymatrices     = function_F_01_init_propertymatrices;
  }
  if ((FUNCTION_F == 2)) {
    thermo_phase = (long*)malloc(NUMPHASES*sizeof(long));
    long a, b;
    for (a=0; a<NUMPHASES; a++) {
      for (b=0; b<NUM_THERMO_PHASES; b++) {
        if (strcmp(phase_map[a],Phases_tdb[b])==0) {
          thermo_phase[a] = b;
          break;
        }
      }
    }
    free_energy               = function_F_02_free_energy;
    dc_dmu                    = function_F_02_dc_dmu;
    c_mu                      = function_F_02_c_mu;
    Mu                        = function_F_02_Mu;
    dpsi                      = function_F_02_dpsi;
  }
  if ((FUNCTION_F == 3)) {
    thermo_phase = (long*)malloc(NUMPHASES*sizeof(long));
    long a, b;
    for (a=0; a<NUMPHASES; a++) {
      for (b=0; b<NUM_THERMO_PHASES; b++) {
        if (strcmp(phase_map[a],Phases_tdb[b])==0) {
          thermo_phase[a] = b;
          break;
        }
      }
    }
    free_energy               = function_F_03_free_energy;
    dc_dmu                    = function_F_03_dc_dmu;
    c_mu                      = function_F_03_c_mu;
    Mu                        = function_F_03_Mu;
    dpsi                      = function_F_03_dpsi;
    function_A                = function_F_03_function_A;
    function_B                = function_F_03_function_B;
    function_C                = function_F_03_function_C;
    init_propertymatrices     = function_F_03_init_propertymatrices;
  }
  if ((FUNCTION_F == 4)) {
    thermo_phase = (long*)malloc(NUMPHASES*sizeof(long));
    long a, b;
    for (a=0; a<NUMPHASES; a++) {
      for (b=0; b<NUM_THERMO_PHASES; b++) {
        if (strcmp(phase_map[a],Phases_tdb[b])==0) {
          thermo_phase[a] = b;
          break;
        }
      }
    }
    free_energy               = function_F_04_free_energy;
    dc_dmu                    = function_F_04_dc_dmu;
    c_mu                      = function_F_04_c_mu;
    Mu                        = function_F_04_Mu;
    dpsi                      = function_F_04_dpsi;
    function_A                = function_F_04_function_A;
    function_B                = function_F_04_function_B;
    function_C                = function_F_04_function_C;
    init_propertymatrices     = function_F_04_init_propertymatrices;
  }
   if ((FUNCTION_F == 5)) {
//     thermo_phase = (long*)malloc(NUMPHASES*sizeof(long));
//     long a, b;
//     for (a=0; a<NUMPHASES; a++) {
//       for (b=0; b<NUM_THERMO_PHASES; b++) {
//         if (strcmp(phase_map[a],Phases_tdb[b])==0) {
//           thermo_phase[a] = b;
//           break;
//         }
//       }
//     }
//     free_energy               = function_F_04_free_energy;
//     dc_dmu                    = function_F_04_dc_dmu;
//     c_mu                      = function_F_04_c_mu;
//     Mu                        = function_F_04_Mu;
//     dpsi                      = function_F_04_dpsi;
//     function_A                = function_F_04_function_A;
//     function_B                = function_F_04_function_B;
//     function_C                = function_F_04_function_C;
//     init_propertymatrices     = function_F_04_init_propertymatrices;
     dpsi                      = function_F_05_dpsi;
  }
  if ((FUNCTION_F == 6)) {
    thermo_phase = (long*)malloc(NUMPHASES*sizeof(long));
    long a, b;
    for (a=0; a < NUMPHASES; a++) {
      for (b=0; b < NUM_THERMO_PHASES; b++) {
        if (strcmp(phase_map[a],Phases_tdb[b])==0) {
          thermo_phase[a] = b;
          break;
        }
      }
    }
    free_energy               = function_F_06_free_energy;
    dc_dmu                    = function_F_06_dc_dmu;
    c_mu                      = function_F_06_c_mu;
    Mu                        = function_F_06_Mu;
    dpsi                      = function_F_06_dpsi;
    function_A                = function_F_06_function_A;
    function_B                = function_F_06_function_B;
    function_C                = function_F_06_function_C;
    init_propertymatrices     = function_F_06_init_propertymatrices;
  }
  
  if (FUNCTION_W == 1) {
    dwdphi        = function_W_01_dwdphi;
    dwdphi_smooth = function_W_01_dwdphi_smooth;
    OBSTACLE = 1;
    WELL = 0;
  }
  if (FUNCTION_W == 2) {
    dwdphi        = function_W_02_dwdphi;
    dwdphi_smooth = function_W_02_dwdphi_smooth;
    OBSTACLE = 0;
    WELL = 1;
  }
  if(FUNCTION_ANISOTROPY == 0) {
    dAdphi               = function_A_00_dAdphi;
    divdAdgradphi        = function_A_00_divdAdgradphi;
    dAdphi_smooth        = function_A_00_dAdphi;
    divdAdgradphi_smooth = function_A_00_divdAdgradphi;
    ANISOTROPY = 0;
  }
  if(FUNCTION_ANISOTROPY == 1) {
    dAdphi               = function_A_01_dAdphi;
    divdAdgradphi        = function_A_01_divdAdgradphi;
    dAdphi_smooth        = function_A_01_dAdphi_smooth;
    divdAdgradphi_smooth = function_A_01_divdAdgradphi_smooth;
    ANISOTROPY = 1;
  }

  if (fab != NULL) {
    dAdq        = anisotropy_gen_dAdq;
    function_ac = anisotropy_gen_function_ac;
  } else {
    if (FOLD == 4) {
      dAdq        = anisotropy_01_dAdq;
      function_ac = anisotropy_01_function_ac;
    }
    if (FOLD == 2) {
      dAdq        = anisotropy_02_dAdq;
      function_ac = anisotropy_02_function_ac;
    }
  }
  
  if (DIMENSION==2) {
    calculate_gradients_phasefield                          = calculate_gradients_phasefield_2D;
    calculate_gradients_concentration                       = calculate_gradients_concentration_2D;
    calculate_fluxes_concentration                          = calculate_fluxes_concentration_2D;
    calculate_divergence_concentration                      = calculate_divergence_concentration_2D;
    calculate_divergence_concentration_smooth               = calculate_divergence_concentration_smooth_2D;
    calculate_divergence_concentration_smooth_concentration = calculate_divergence_concentration_smooth_concentration_2D;
    calculate_divergence_phasefield                         = calculate_divergence_phasefield_2D;
    calculate_divergence_phasefield_smooth                  = calculate_divergence_phasefield_smooth_2D;
    calculate_divergence_stress                             = calculate_divergence_stress_2D;
    calculate_gradients_stress                              = calculate_gradients_stress_2D;
    calculate_streaming_row_forward                         = calculate_streaming_row_forward_2D;
    calculate_streaming_row_back                            = calculate_streaming_row_back_2D;
    calculate_streaming_col_forward                         = calculate_streaming_col_forward_2D;
    calculate_streaming_col_back                            = calculate_streaming_col_back_2D;
    compute_field_variables_LBM                             = compute_field_variables_LBM_2D;
    calculate_distribution_function                         = calculate_distribution_function_2D;
    solverloop_concentration                                = solverloop_concentration_tdb;
    solverloop_phasefield                                   = solverloop_phasefield_tdb;
    writetofile_mpi                                         = writetofile_mpi2D;
    writetofile_mpi_binary                                  = writetofile_mpi2D_binary;
    writetofile_mpi_hdf5                                    = writetofile_mpi2D_hdf5;
    readfromfile_mpi_hdf5                                   = readfromfile_mpi2D_hdf5;
    readfromfile_mpi                                        = readfromfile_mpi2D;
    readfromfile_mpi_binary                                 = readfromfile_mpi2D_binary;
    compute_error                                           = compute_error_2D;
    df_elast                                                = df_elast_2D;
  } else {
    calculate_gradients_phasefield                          = calculate_gradients_phasefield_3D;
    calculate_gradients_concentration                       = calculate_gradients_concentration_3D;
    calculate_fluxes_concentration                          = calculate_fluxes_concentration_3D;
    calculate_divergence_concentration                      = calculate_divergence_concentration_3D;
    calculate_divergence_concentration_smooth               = calculate_divergence_concentration_smooth_3D;
    calculate_divergence_concentration_smooth_concentration = calculate_divergence_concentration_smooth_concentration_3D;
    calculate_divergence_phasefield                         = calculate_divergence_phasefield_3D;
    calculate_divergence_phasefield_smooth                  = calculate_divergence_phasefield_smooth_3D;
    calculate_divergence_stress                             = calculate_divergence_stress_3D;
    calculate_gradients_stress                              = calculate_gradients_stress_3D;
    calculate_streaming_row_forward                         = calculate_streaming_row_forward_3D;
    calculate_streaming_row_back                            = calculate_streaming_row_back_3D;
    calculate_streaming_col_forward                         = calculate_streaming_col_forward_3D;
    calculate_streaming_col_back                            = calculate_streaming_col_back_3D;
    compute_field_variables_LBM                             = compute_field_variables_LBM_3D;
    calculate_distribution_function                         = calculate_distribution_function_3D;
    solverloop_concentration                                = solverloop_concentration_tdb;
    solverloop_phasefield                                   = solverloop_phasefield_tdb;
    writetofile_mpi                                         = writetofile_mpi3D;
    writetofile_mpi_binary                                  = writetofile_mpi3D_binary;
    writetofile_mpi_hdf5                                    = writetofile_mpi3D_hdf5;
    readfromfile_mpi_hdf5                                   = readfromfile_mpi3D_hdf5;
    readfromfile_mpi_binary                                 = readfromfile_mpi3D_binary; 
    compute_error                                           = compute_error_3D;
    df_elast                                                = df_elast_3D;
  }
}
#endif
