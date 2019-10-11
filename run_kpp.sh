#!/bin/bash
lon=100
lat=-40

#data=HadISST_sst_clim.nc
#ferret << hadsstclm
#set mem/size=200
#set data/order=xyt "/Users/byjup/Data/sst/${data}"
#save/file="hadlSST_clim.nc"/x=${lon}/y=${lat}/clob hadsst
#exit
#hadsstclm

#data=HadISST_sst_1870-2011_Anom.nc
#ferret << hadsst
#set mem/size=200
#use "/Users/byjup/Data/sst/${data}"
#save/keep_axisnames/file="hadlSSTa.nc"/x=${lon}/y=${lat}/clob hadsst
#exit
#hadsst

#data=TS_clim_at_modelevls_12mnths.nc
#ferret << initcond
#set mem/size=200
#use "/Users/byjup/kpp_model/extract_data/${data}"
#save/keep_axisnames/file="TS.Clim.mon.WOA09.nc"/x=${lon}/y=${lat}/clob temp,sal
#exit
#initcond

#data=ncep_monthly_forcing.nc
#ferret << monthly
#set mem/size=200
#use "/Users/byjup/kpp_model/extract_data/${data}"
#save/keep_axisnames/file="forcing.mon.nc"/x=${lon}/y=${lat}/l=1:804/clob solar_in,nsolar_in,taux_in,tauy_in,pminuse_in
#exit
#monthly

#data=ncep-ncar-4xdaily.Clim.nc
#ferret << sixhrly
#set mem/size=200
#use "/Users/byjup/kpp_model/extract_data/${data}"
#save/keep_axisnames/file="forcing.6hrly.clim.nc"/x=${lon}/y=${lat}/l=10:1449/clob solar_in,nsolar_in,taux_in,tauy_in,pminuse_in
#exit
#sixhrly

#generate noise forces
/opt/local/bin/gfortran-mp-4.8 -I/opt/local/include -L/opt/local/lib -lnetcdf -lnetcdff -fdefault-real-8 generate_6hrly_forcing.f 
./a.out
# compute monthly anomalies
cdo ymonsub forcing.mon.nc -ymonmean forcing.mon.nc ncep.force.mon.Anom.nc
#convert 6hrly noise forces to monthly
/opt/local/bin/gfortran-mp-4.8 -I/opt/local/include -L/opt/local/lib -lnetcdf -lnetcdff -fdefault-real-8 convert_6hr2mon_atmforc.f
./a.out
#compute the anomalies
cdo -r settaxis,1901-01-15,00:00,1mon -setcalendar,360_day noise.force.mon.nc test.nc
mv test.nc noise.force.mon.nc
cdo ymonsub noise.force.mon.nc -ymonmean noise.force.mon.nc noise.force.mon.Anom.nc
# KPP single point model simulation
make clean
make forced
./KPP_ocean
rm -f *.o ferret.jnl*

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
~     
