      LOGICAL L_FCORR_WITHZ, L_FCORR, L_UPD_FCORR, L_PERIODIC_FCORR      
      INTEGER ndtupdfcorr, fcorr_period
      REAL fcorr_withz(NZP1)
      REAL fcorr_saltz(NZP1)   ! DIDI: add for qflux
      REAL fcorr_twod
      REAL tanual(NZP1),sanual(NZP1)
      REAL tclm(NZP1,12),sclm(NZP1,12)
      CHARACTER*40 fcorrin_file
      common /fcorr_in/ L_FCORR_WITHZ,L_FCORR,
     &     ndtupdfcorr, fcorrin_file,L_UPD_FCORR, L_PERIODIC_FCORR,
     &     fcorr_period
