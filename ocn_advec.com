      INTEGER nmodeadv(2),modeadv(maxmodeadv,2)
      REAL advection(maxmodeadv,2)
      LOGICAL L_ADVECT
      CHARACTER*40 advect_file
      INTEGER ncid_advec
      common/ ocn advec / nmodeadv,modeadv,
     &     advection,L_ADVECT,advect_file,ncid_advec

      LOGICAL L_RELAX_SST,L_RELAX_CALCONLY
      REAL relax_sst,SST0,fcorr
      REAL relax_sst_in(ny)

      common/ ocn_relax / relax_sst,relax_sst_in,SST0,fcorr,
     &   L_RELAX_SST,L_RELAX_CALCONLY
