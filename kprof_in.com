      REAL ustar,B0,buoy(NZP1tmax),B0sol
      common/ kprof_in  / ustar,B0,buoy,B0sol

