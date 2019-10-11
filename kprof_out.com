      INTEGER kmix
      REAL hmix,difm(0:NZtmax),
     $          difs(0:NZtmax),ghat(NZtmax)
      common/ kprof out / hmix,difm,difs,ghat,kmix

