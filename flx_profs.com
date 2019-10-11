      REAL wU(0:NZtmax,NVP1),wX(0:NZtmax,NSP1),
     +     wXNT(0:NZtmax,NSCLR),
     +     wtI(0:NZ,0:NJDT),wsI(0:NZ,0:NJDT) 
      common/ flx profs / wU,wX,
     +     wXNT,
     +     wtI,wsI

