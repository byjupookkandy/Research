      REAL dtsec,time,startt
      INTEGER ntime,nstart,nend,NRECS,finalt
      common/ times / dtsec,time,startt,finalt,ntime,nstart,nend,NRECS

