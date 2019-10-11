      INTEGER jerlov
      common/ proc_pars / jerlov
      
      INTEGER ncid_paras
      LOGICAL L_JERLOV
      common/ l_proc_pars / ncid_paras,
     $     l_jerlov
      
