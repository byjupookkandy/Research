#!/bin/bash
cdo ymonsub forcing.mon.nc -ymonmean forcing.mon.nc ncep.force.mon.Anom.nc

/opt/local/bin/gfortran-mp-4.8 -I/opt/local/include -L/opt/local/lib -lnetcdf -lnetcdff -fdefault-real-8 convert_6hr2mon_atmforc.f
./a.out

cdo -r settaxis,1901-01-15,00:00,1mon -setcalendar,360_day noise.force.mon.nc test.nc
mv test.nc noise.force.mon.nc
cdo ymonsub noise.force.mon.nc -ymonmean noise.force.mon.nc noise.force.mon.Anom.nc
