      REAL dlat,dlon,rlat,rlon,f
      LOGICAL L_REGGRID
      common/ location  / dlat,dlon,rlat,rlon,f,L_REGGRID

