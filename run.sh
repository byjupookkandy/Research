#!/bin/sh
make clean
make forced
./KPP_ocean
rm -f *.o

#convert 6hrly 500 year kpp output to monthly, and then anomalies
#/opt/local/bin/gfortran-mp-4.8 -I/opt/local/include -L/opt/local/lib -lnetcdf -lnetcdff -fdefault-real-8 convert_6hr2mon.f
#./a.out
cdo -r settaxis,1901-01-15,00:00,1mon -setcalendar,360_day KPPocean_0009.nc kppdummy.nc
cdo selyear,2001/2500 kppdummy.nc kpp.mon.500yr.nc
cdo ymonsub kpp.mon.500yr.nc -ymonmean kpp.mon.500yr.nc kpp.mon.Anom.500yr.nc
rm -f KPPocean_0009.nc kppdummy.nc

cdo ymonmean kpp.mon.500yr.nc kpp.mon.clim.500yr.nc
cdo ymonmean noise.force.mon.nc noise.force.mon.clim.nc
cdo ymonmean forcing.mon.nc ncep.forcing.mon.clim.nc
