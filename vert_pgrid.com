      REAL DMAX,zm(NZP1),hm(NZP1),dm(0:NZ)
      common/ vert pgrid/ DMAX,zm,hm,dm

