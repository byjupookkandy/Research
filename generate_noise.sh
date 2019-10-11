#!/bin/bash
lon=100
lat=-40

data=HadISST_sst_1870-2011_Anom.nc
ferret << hadsst
set mem/size=250
use "/Users/byjup/Data/sst/${data}"
save/keep_axisnames/file="hadlSSTa.nc"/x=${lon}/y=${lat}/clob hadsst
exit
hadsst

data=TS_clim_at_modelevls_12mnths.nc
ferret << initcond
set mem/size=250
use "../extract_data/${data}"
save/keep_axisnames/file="TS.Clim.mon.WOA09.nc"/x=${lon}/y=${lat}/clob temp,sal
exit
initcond

data=ncep_monthly_forcing.nc
ferret << monthly
set mem/size=250
use "../extract_data/${data}"
save/keep_axisnames/file="forcing.mon.nc"/x=${lon}/y=${lat}/l=1:804/clob solar_in,nsolar_in,taux_in,tauy_in,pminuse_in
exit
monthly

data=ncep-ncar-4xdaily.Clim.nc
ferret << sixhrly
set mem/size=250
use "../extract_data/${data}"
save/keep_axisnames/file="forcing.6hrly.clim.nc"/x=${lon}/y=${lat}/l=10:1449/clob solar_in,nsolar_in,taux_in,tauy_in,pminuse_in
exit
sixhrly

/opt/local/bin/gfortran-mp-4.8 -I/opt/local/include -L/opt/local/lib -lnetcdf -lnetcdff -fdefault-real-8 generate_6hrly_forcing.f 
./a.out

cdo ymonsub forcing.mon.nc -ymonmean forcing.mon.nc ncep.force.mon.Anom.nc

/opt/local/bin/gfortran-mp-4.8 -I/opt/local/include -L/opt/local/lib -lnetcdf -lnetcdff -fdefault-real-8 convert_6hr2mon_atmforc.f
./a.out

cdo -r settaxis,1901-01-15,00:00,1mon -setcalendar,360_day noise.force.mon.nc test.nc
mv test.nc noise.force.mon.nc
cdo ymonsub noise.force.mon.nc -ymonmean noise.force.mon.nc noise.force.mon.Anom.nc
