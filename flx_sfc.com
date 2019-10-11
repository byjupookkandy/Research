      REAL sflux(NSFLXS,5,0:NJDT),VAF(NSFLXSP2)
      common/ flx_sfc   / sflux,VAF

