      REAL CP(0:NZP1tmax),rho(0:NZP1tmax),rhoh2o,
     &     rhob,talpha(0:NZP1tmax),sbeta(0:NZP1tmax),
     &     Sref,SSref,epsw  
      LOGICAL L_SSRef
      common/ ocn_paras / CP,rho,rhoh2o,rhob,   
     &     talpha,sbeta,Sref,SSref,epsw,L_SSref

