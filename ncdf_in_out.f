	SUBROUTINE initialise(temp,sal,lati,longi)
!********************************************************
!*** This program reads the initial condition temp and salt profile
!***** edited by Byju Pookkandy on 20-aug-2015 
!********************************************************
!******* MAIN PROGRAM ************************
C     This is part of the netCDF package.
        implicit none
        include 'netcdf.inc'
	include 'parameter.inc'

c	real :: startcp, finish !just to know the time taken for the run
C     This is the name of the data file we will read.
        character*(*) FILE_INIT
	parameter (FILE_INIT='TS.Clim.mon.WOA09.nc')
        integer ncid
     
C     Local variable
        integer NDIMS,NTIM
        parameter (NDIMS = 4,NTIM = 12)
        character*(*) LVL_NAME, LAT_NAME, LON_NAME, REC_NAME
        parameter (LVL_NAME = 'Z')
        parameter (LAT_NAME = 'LATITUDE', LON_NAME = 'LONGITUDE')
        parameter (REC_NAME = 'TIME')
     
C     The start and count arrays will tell the netCDF library where to
C     read our data.
        integer start(NDIMS), count(NDIMS)  ! for temperature and salinity
 
C     In addition to the latitude and longitude dimensions, we will also
C     create latitude and longitude variables which will hold the actual
C     latitudes and longitudes. Since they hold data about the
C     coordinate system, the netCDF term for these is: "coordinate
C     variables."
        real lati,longi,deps(NZP1),tims(1)
        integer lon_varid, lat_varid, lvl_varid, rec_varid
     
C     We will read surface temperature and pressure fields. In netCDF
C     terminology these are called "variables."
        character*(*)  TEMP_NAME,SAL_NAME
        parameter (TEMP_NAME='TEMP',SAL_NAME='SAL')
	integer temp_varid,sal_varid
	     
C     Error handling.
        integer retval

       real temp(NX,NY,NZP1,NTIM)
       real sal(NX,NY,NZP1,NTIM)

!**************************
C     Open the file to read temperature data
!**************************
        retval = nf_open(FILE_INIT, nf_nowrite, ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)
	write(*,*)'file opened = ',FILE_INIT 
C     Get the varids of the latitude and longitude coordinate variables.
        retval = nf_inq_varid(ncid, LON_NAME, lon_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_inq_varid(ncid, LAT_NAME, lat_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
C     Read the latitude and longitude data.
        retval = nf_get_var_double(ncid, lat_varid, lati)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_get_var_double(ncid, lon_varid, longi)
        if (retval .ne. nf_noerr) call handle_err(retval)
	
! get the varid and read depth
        retval = nf_inq_varid(ncid,LVL_NAME,lvl_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_get_var_double(ncid, lvl_varid, deps)
        if (retval .ne. nf_noerr) call handle_err(retval)

! get the varid and read time
	retval = nf_inq_varid(ncid,REC_NAME,rec_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_get_var_double(ncid, rec_varid, tims)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     Get the varids of temperature netCDF variables.
        retval = nf_inq_varid(ncid,TEMP_NAME, temp_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
C     Read the temperature data from the file, one
        count = (/NX, NY, NZP1,NTIM/)
	start = (/1,1,1,1/)
	retval = nf_get_vara_double(ncid,temp_varid,start,count,temp)
        if (retval .ne. nf_noerr) call handle_err(retval)

C	open and get the varid of salinity data
	retval = nf_inq_varid(ncid,SAL_NAME, sal_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
C     Read the salinity data from the file, one
C     record at a time.
        count = (/NX, NY, NZP1,NTIM/)
        start = (/1,1,1,1/)
        retval = nf_get_vara_double(ncid,sal_varid,start,count,sal)
        if (retval .ne. nf_noerr) call handle_err(retval)

	retval = nf_close(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

	write(*,1001)longi,lati
1001	format('1D model is initialised with Temperature,Salinity',/,
     .     'profile at longitude=',f6.2,' and latitide=',f6.2)        
!****************************     
C     If we got this far, everything worked as expected. Yipee!
!********************************************************************
!c	call cpu_time(finish)
!c        print '("Time = ",f10.3," minutss.")',(finish-startcp)/60
	return 
        end subroutine

!********************************************************************
!********************************************************************
!********************************************************************

	SUBROUTINE forcing(taux,tauy,solar,nsolar,pme)
!********************************************************
!***** this program reads the input forcing flux 
!*****  Byju Pookkandy on 20-Aug-2015 
!********************************************************
!******* MAIN PROGRAM ************************
C     This is part of the netCDF package.
        implicit none
        include 'netcdf.inc'
	include 'parameter.inc'

c	real :: startcp, finish !just to know the time taken for the run
C     This is the name of the data file we will read.
        character*(*) FILE_FORC
        parameter (FILE_FORC='noise.force.6hrly.nc')
        integer ncid
     
C     Local variable
        integer NDIMS,nt
	parameter (nt=totim*12)
        parameter (NDIMS = 3)
        character*(*) LAT_NAME, LON_NAME, REC_NAME
        parameter (LAT_NAME = 'latitude', LON_NAME = 'longitude')
        parameter (REC_NAME = 'time')
     
C     The start and count arrays will tell the netCDF library where to
C     read our data.
       integer start(NDIMS), count(NDIMS)
 
C     In addition to the latitude and longitude dimensions, we will also
C     create latitude and longitude variables which will hold the actual
C     latitudes and longitudes. Since they hold data about the
C     coordinate system, the netCDF term for these is: "coordinate
C     variables."
        real lats,lons,tims(nt)
        integer lon_varid, lat_varid, rec_varid
     
C     We will read surface temperature and pressure fields. In netCDF
C     terminology these are called "variables."
	character*(*) TAUX_NAME,TAUY_NAME
	parameter (TAUX_NAME='taux_in',TAUY_NAME='tauy_in')
	character*(*) SOLAR_NAME,NSOLAR_NAME,PME_NAME
	parameter (SOLAR_NAME='solar_in',NSOLAR_NAME='nsolar_in')
	parameter (PME_NAME='pme_in')
        integer taux_varid,tauy_varid,solar_varid,nsolar_varid
	integer pme_varid
	     
C     Error handling.
        integer retval
       real taux(NX,NY,nt)
       real tauy(NX,NY,nt)
       real solar(NX,NY,nt)
       real nsolar(NX,NY,nt)
       real pme(NX,NY,nt)

!	write(*,*)'started reading forces'

!**************************
C     Open the file to read temperature data
!**************************
	retval = nf_open(FILE_FORC, nf_nowrite, ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)
       write(*,*)'file opened = ',FILE_FORC
C     Get the varids of the latitude and longitude coordinate variables.
        retval = nf_inq_varid(ncid, LON_NAME, lon_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_inq_varid(ncid, LAT_NAME, lat_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
C     Read the latitude and longitude data.
        retval = nf_get_var_double(ncid, lat_varid, lats)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_get_var_double(ncid, lon_varid, lons)
        if (retval .ne. nf_noerr) call handle_err(retval)
	
! get the varid and read time
	retval = nf_inq_varid(ncid,REC_NAME,rec_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_get_var_double(ncid, rec_varid, tims)
        if (retval .ne. nf_noerr) call handle_err(retval)

	write(*,*)'Starts reading input forcing flux data'

!Get the var ids of forcing flux TAUX
        retval = nf_inq_varid(ncid,TAUX_NAME, taux_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        count = (/NX, NY, nt/)
        start = (/1,1,1/)
        retval = nf_get_vara_double(ncid,taux_varid,start,count,taux)
        if (retval .ne. nf_noerr) call handle_err(retval)

!Get the var ids of forcing flux TAUY
	retval = nf_inq_varid(ncid,TAUY_NAME, tauy_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_get_vara_double(ncid,tauy_varid,start,count,tauy)
        if (retval .ne. nf_noerr) call handle_err(retval)

!Get the var ids of forcing flux SOLAR
	retval = nf_inq_varid(ncid,SOLAR_NAME, solar_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_get_vara_double(ncid,solar_varid,start,count,solar)
        if (retval .ne. nf_noerr) call handle_err(retval)

!Get the var ids of forcing flux NSOLAR
        retval = nf_inq_varid(ncid,NSOLAR_NAME, nsolar_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_get_vara_double(ncid,nsolar_varid,start,
     .              count,nsolar)
        if (retval .ne. nf_noerr) call handle_err(retval)

!Get the var ids of forcing flux P-E
        retval = nf_inq_varid(ncid,PME_NAME, pme_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_get_vara_double(ncid,pme_varid,start,count,pme)
        if (retval .ne. nf_noerr) call handle_err(retval)

        retval = nf_close(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)
!****************************     
C     If we got this far, everything worked as expected. Yipee!
!********************************************************************
c	call cpu_time(finish)
c        print '("Time = ",f10.3," minutss.")',(finish-startcp)/60
	return 
        end
!********************************************************************
!********************************************************************

        subroutine handle_err(errcode)
        implicit none
        include 'netcdf.inc'
        integer errcode
     
        print *, 'Error: ', nf_strerror(errcode)
        stop 2
        end

!********************************************************************
!********************************************************************
!Now routine for writing the output
	subroutine write_output(uvel,vvel,tout,sout,mixd,
     +                                    nmm,lats,lons,FILE_OUT)

!********************************************************
!***** this program reads the input forcing flux 
!*****  Byju Pookkandy on 20-Aug-2015 
!********************************************************
!******* MAIN PROGRAM ************************
C     This is part of the netCDF package.
        implicit none
        include 'netcdf.inc'
        include 'parameter.inc'
        include 'vert_pgrid.com'
        include 'times.com'

        real :: startcp, finish !just to know the time taken for the run
C     This is the name of the data file we will read.
        character*(*) FILE_OUT
!        parameter (FILE_OUT='kpp_output.nc')
        integer ncid,i

C     Local variable
        integer NDIMS,nmm
        parameter (NDIMS = 4)
        character*(*) LAT_NAME, LON_NAME, REC_NAME,LVL_NAME
        parameter (LAT_NAME = 'latitude', LON_NAME = 'longitude')
        parameter (REC_NAME = 'time',LVL_NAME = 'z')

        character*(*) LAT_UNITS, LON_UNITS
        character*(*) REC_UNITS,LVL_UNITS
        character*(*) UNITS
        parameter (UNITS = 'units')
        parameter (LAT_UNITS = 'degrees_north')
        parameter (LON_UNITS = 'degrees_east')
        parameter (REC_UNITS = 'days')
        parameter (LVL_UNITS = 'm')
        character*(*) VEL_UNITS, TEMP_UNITS
        character*(*) SALT_UNITS,MLD_UNITS
        parameter (VEL_UNITS = 'm/s',MLD_UNITS='m')
        parameter (TEMP_UNITS='oC',SALT_UNITS='PSU')
        character*(*) RHO_UNITS
        parameter ( RHO_UNITS='kg/m3')

C     The start and count arrays will tell the netCDF library where to
C     read our data.
        integer start(NDIMS), count(NDIMS)
        integer starth(3), counth(3)   ! this is for hmix

C     In addition to the latitude and longitude dimensions, we will also
C     create latitude and longitude variables which will hold the actual
C     latitudes and longitudes. Since they hold data about the
C     coordinate system, the netCDF term for these is: "coordinate
C     variables."
       real lats,lons,tims(nmm),lvl(NZP1)
       integer lon_varid, lat_varid, rec_varid, lvl_varid,z_varid
       integer lat_dimid,lon_dimid,rec_dimid,lvl_dimid,z_dimid

C     We will read surface temperature and pressure fields. In netCDF
C     terminology these are called "variables."
        character*(*) UVEL_NAME,VVEL_NAME
        parameter (UVEL_NAME='uvel',VVEL_NAME='vvel')
        character*(*) TEMP_NAME,SALT_NAME
        parameter (TEMP_NAME='temp',SALT_NAME='sal')
        character*(*) MLD_NAME,RHO_NAME
	parameter (MLD_NAME='hmix',RHO_NAME='rho')
        integer uvel_varid,vvel_varid,temp_varid,salt_varid
        integer mld_varid,rho_varid
        
C     Error handling.
       integer retval
       integer dimids(NDIMS),dimidh(3)

! index variable
       integer k,l

       real uvel(NX,NY,NZP1,nmm)
       real vvel(NX,NY,NZP1,nmm)
       real tout(NX,NY,NZP1,nmm)
       real sout(NX,NY,NZP1,nmm)
       real mixd(NX,NY,nmm)

! local data rho=f(sal,temp)
       real,allocatable, dimension(:,:,:,:) ::  den
       real Sig0,Sig
       logical KapFlg
       allocate(den(NX,NY,NZP1,nmm))

       KapFlg=.False.

       lvl=zm
       do i=1,nmm
        tims(i)=i
       end do
       write(*,*)'nmm=',nmm
! compute density 
       do l=1,nmm
        do k=1,NZP1
          call Sig80(sout(1,1,k,l),tout(1,1,k,l),0.,
     +                                       KapFlg,Sig0,Sig)
          den(1,1,k,l)=1000.+Sig0
        end do
       end do 
! end computing density

        write(*,*)'*** Now starts writing output netcdf file'
!        call cpu_time(startcp)
!        write(*,*)'called cputime=',startcp
!       Create the file
        retval=nf_create(FILE_OUT,ior(nf_clobber,nf_64bit_offset),ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     Define the dimensions.
        retval = nf_def_dim(ncid, LAT_NAME, NY, lat_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, LON_NAME, NX, lon_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, LVL_NAME, NZP1, lvl_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, REC_NAME, NF_UNLIMITED, rec_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     Define the coordinate variables.
        retval = nf_def_var(ncid, LAT_NAME, NF_DOUBLE, 1, lat_dimid,
     +     lat_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, LON_NAME, NF_DOUBLE, 1, lon_dimid,
     +     lon_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, LVL_NAME, NF_DOUBLE, 1, lvl_dimid,
     +     lvl_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, REC_NAME, NF_DOUBLE, 1, rec_dimid,
     +     rec_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     Assign units attributes to coordinate variables.
        retval = nf_put_att_text(ncid, lat_varid, UNITS, len(LAT_UNITS),
     +     LAT_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, lon_varid, UNITS, len(LON_UNITS),
     +     LON_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, lvl_varid, UNITS, len(LVL_UNITS),
     +     LVL_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, rec_varid, UNITS, len(REC_UNITS),
     +     REC_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     The dimids array is used to pass the dimids of the dimensions of
C     the netCDF variables.
        dimids = (/lon_dimid,lat_dimid,lvl_dimid,rec_dimid/)
!       Define netcdf variable for Mixed layer depth (MLD)
        retval = nf_def_var(ncid, UVEL_NAME, NF_DOUBLE, 4, dimids,
     +     uvel_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, VVEL_NAME, NF_DOUBLE, 4, dimids,
     +     vvel_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, TEMP_NAME, NF_DOUBLE, 4, dimids,
     +     temp_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, SALT_NAME, NF_DOUBLE, 4, dimids,
     +     salt_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, RHO_NAME, NF_DOUBLE, 4, dimids,
     +     rho_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        dimidh = (/lon_dimid,lat_dimid,rec_dimid/)
        retval = nf_def_var(ncid, MLD_NAME, NF_DOUBLE, 3, 
     +          dimidh, mld_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

!       Assign units attributes to the netCDF variables.
        retval = nf_put_att_text(ncid, uvel_varid, UNITS,len(VEL_UNITS),
     +     VEL_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, vvel_varid, UNITS,len(VEL_UNITS),
     +     VEL_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, temp_varid,UNITS,len(temp_UNITS),
     +     TEMP_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, salt_varid,UNITS,len(SALT_UNITS),
     +     SALT_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, rho_varid,UNITS,len(RHO_UNITS),
     +     RHO_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, mld_varid, UNITS, len(MLD_UNITS),
     +     MLD_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
C     End define mode.
        retval = nf_enddef(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)
C     Write the coordinate variable data. This will put the latitudes
C     ,longitudes,depth and time of our data grid into the netCDF file.
        retval = nf_put_var_double(ncid, lat_varid, lats)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_double(ncid, lon_varid, lons)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_double(ncid, lvl_varid, lvl)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_double(ncid, rec_varid, tims)
        if (retval .ne. nf_noerr) call handle_err(retval)

        count = (/NX,NY,NZP1,nmm/)
        start = (/1,1,1,1/)
        retval = nf_put_vara_double(ncid, uvel_varid, start, count,uvel)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_vara_double(ncid, vvel_varid, start, count,vvel)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_vara_double(ncid, temp_varid, start, count,tout)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_vara_double(ncid, salt_varid, start, count,sout)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_vara_double(ncid, rho_varid, start, count,den)
        if (retval .ne. nf_noerr) call handle_err(retval)

	counth = (/NX,NY,nmm/)
        starth = (/1,1,1/)
        retval = nf_put_vara_double(ncid, mld_varid, starth, 
     +                                           counth,mixd)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     Close the file. This causes netCDF to flush all buffers and make
C     sure your data are really written to disk.
        retval = nf_close(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        write(*,*)'*** SUCCESS writing output file  ', FILE_OUT
!        call cpu_time(finish)
!        print '("Time = ",f10.3," minutss.")',(finish-startcp)/60
        deallocate(den)
        return
	end
!********** End of creating netcdf output ***************************

!********************************************************************
!********************************************************************
!Now routine for writing the output forces
	subroutine write_forces(taux_in,tauy_in,solar_in,nsolar_in,
     +                   pme_in,noise,lats,lons,FILE_OUT)

!********************************************************
!***** this program writes the noise generated forcing 
!*****  Byju Pookkandy on 09-Nov-2015 
!********************************************************
!******* MAIN PROGRAM ************************
C     This is part of the netCDF package.
        implicit none
        include 'netcdf.inc'
        include 'parameter.inc'
        include 'vert_pgrid.com'
        include 'times.com'

        real :: startcp, finish !just to know the time taken for the run
C     This is the name of the data file we will read.
        character*(*) FILE_OUT
        integer ncid,i

C     Local variable
        integer NDIMS
        parameter (NDIMS = 3)
        character*(*) LAT_NAME, LON_NAME, REC_NAME,LVL_NAME
        parameter (LAT_NAME = 'latitude', LON_NAME = 'longitude')
        parameter (REC_NAME = 'time',LVL_NAME = 'z')

        character*(*) LAT_UNITS, LON_UNITS
        character*(*) REC_UNITS
        character*(*) UNITS
        parameter (UNITS = 'units')
        parameter (LAT_UNITS = 'degrees_north')
        parameter (LON_UNITS = 'degrees_east')
        parameter (REC_UNITS = 'days')
        character*(*) TAU_UNITS, NHF_UNITS
        character*(*) PME_UNITS,NOS_UNITS
        parameter (TAU_UNITS = 'N/m2',NOS_UNITS='random')
        parameter (NHF_UNITS='W/m2',PME_UNITS='mm/s')

C     The start and count arrays will tell the netCDF library where to
C     read our data.
        integer start(NDIMS), count(NDIMS)

C     In addition to the latitude and longitude dimensions, we will also
C     create latitude and longitude variables which will hold the actual
C     latitudes and longitudes. Since they hold data about the
C     coordinate system, the netCDF term for these is: "coordinate
C     variables."
       real lats,lons,tims(NRECS)
       integer lon_varid, lat_varid, rec_varid
       integer lat_dimid,lon_dimid,rec_dimid

C     We will read surface temperature and pressure fields. In netCDF
C     terminology these are called "variables."
        character*(*) TAUX_NAME,TAUY_NAME
        parameter (TAUX_NAME='taux_in',TAUY_NAME='tauy_in')
        character*(*) SWF_NAME,LWF_NAME
        parameter (SWF_NAME='solar_in',LWF_NAME='nsolar_in')
        character*(*) PME_NAME,NOS_NAME
	parameter (PME_NAME='pme_in',NOS_NAME='noise')
        integer taux_varid,tauy_varid,solar_varid,nsolar_varid
        integer pme_varid,nos_varid
        
C     Error handling.
       integer retval
       integer dimids(NDIMS)

! index variable
       integer k,l

       real taux_in(NX,NY,NRECS)
       real tauy_in(NX,NY,NRECS)
       real solar_in(NX,NY,NRECS)
       real nsolar_in(NX,NY,NRECS)
       real pme_in(NX,NY,NRECS)
       real noise(NX,NY,NRECS)

       do i=1,NRECS
        tims(i)=i
       end do

        write(*,*)'*** Now starts writing output netcdf file'
!        call cpu_time(startcp)
!        write(*,*)'called cputime=',startcp
!       Create the file
        retval=nf_create(FILE_OUT,ior(nf_clobber,nf_64bit_offset),ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     Define the dimensions.
        retval = nf_def_dim(ncid, LAT_NAME, NY, lat_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, LON_NAME, NX, lon_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, REC_NAME, NRECS, rec_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     Define the coordinate variables.
        retval = nf_def_var(ncid, LAT_NAME, NF_DOUBLE, 1, lat_dimid,
     +     lat_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, LON_NAME, NF_DOUBLE, 1, lon_dimid,
     +     lon_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, REC_NAME, NF_DOUBLE, 1, rec_dimid,
     +     rec_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     Assign units attributes to coordinate variables.
        retval = nf_put_att_text(ncid, lat_varid, UNITS, len(LAT_UNITS),
     +     LAT_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, lon_varid, UNITS, len(LON_UNITS),
     +     LON_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, rec_varid, UNITS, len(REC_UNITS),
     +     REC_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     The dimids array is used to pass the dimids of the dimensions of
C     the netCDF variables.
        dimids = (/lon_dimid,lat_dimid,rec_dimid/)
!       Define netcdf variable for Mixed layer depth (MLD)
        retval = nf_def_var(ncid, TAUX_NAME, NF_DOUBLE, 3, dimids,
     +     taux_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, TAUY_NAME, NF_DOUBLE, 3, dimids,
     +     tauy_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, SWF_NAME, NF_DOUBLE, 3, dimids,
     +     solar_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, LWF_NAME, NF_DOUBLE, 3, dimids,
     +     nsolar_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, PME_NAME, NF_DOUBLE, 3, dimids,
     +     pme_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, NOS_NAME, NF_DOUBLE, 3, 
     +          dimids, nos_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

!       Assign units attributes to the netCDF variables.
        retval = nf_put_att_text(ncid, taux_varid, UNITS,len(TAU_UNITS),
     +     TAU_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, tauy_varid, UNITS,len(TAU_UNITS),
     +     TAU_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, solar_varid,UNITS,len(NHF_UNITS),
     +     NHF_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid,nsolar_varid,UNITS,len(NHF_UNITS),
     +     NHF_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, pme_varid,UNITS,len(PME_UNITS),
     +     PME_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, nos_varid, UNITS, len(NOS_UNITS),
     +     NOS_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
C     End define mode.
        retval = nf_enddef(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)
C     Write the coordinate variable data. This will put the latitudes
C     ,longitudes,depth and time of our data grid into the netCDF file.
        retval = nf_put_var_double(ncid, lat_varid, lats)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_double(ncid, lon_varid, lons)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_double(ncid, rec_varid, tims)
        if (retval .ne. nf_noerr) call handle_err(retval)

        count = (/NX,NY,NRECS/)
        start = (/1,1,1/)
        retval = nf_put_vara_double(ncid, taux_varid, start, count,
     +                                                     taux_in)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_vara_double(ncid, tauy_varid, start, count,
     +                                                     tauy_in)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_vara_double(ncid, solar_varid, start, count,
     +                                                     solar_in)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_vara_double(ncid, nsolar_varid, start, count,
     +                                                     nsolar_in)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_vara_double(ncid, pme_varid, start, count,
     +                                                     pme_in)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_vara_double(ncid, nos_varid, start, count,
     +                                                     noise)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     Close the file. This causes netCDF to flush all buffers and make
C     sure your data are really written to disk.
        retval = nf_close(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        write(*,*)'*** SUCCESS writing output file  ', FILE_OUT
!        call cpu_time(finish)
!        print '("Time = ",f10.3," minutss.")',(finish-startcp)/60
        return
	end
!********** End of writing netcdf output forces *********************

!********** creating flux correction netcdf output ******************

        subroutine write_fcor(qcorr,scorr,lats,lons,FILE_OUT)

!********************************************************
!***** this program reads the input forcing flux 
!*****  Byju Pookkandy on 20-Aug-2015 
!********************************************************
!******* MAIN PROGRAM ************************
C     This is part of the netCDF package.
        implicit none
        include 'netcdf.inc'
        include 'parameter.inc'
        include 'vert_pgrid.com'

        real :: startcp, finish !just to know the time taken for the run
C     This is the name of the data file we will read.
        character(len=40) FILE_OUT
!        parameter (FILE_OUT='fcorr_1.nc')
        integer ncid,i

C     Local variable
        integer NDIMS
        parameter (NDIMS = 4)
        character*(*) LAT_NAME, LON_NAME, REC_NAME,LVL_NAME
        parameter (LAT_NAME = 'latitude', LON_NAME = 'longitude')
        parameter (REC_NAME = 'time',LVL_NAME = 'z')

        character*(*) LAT_UNITS, LON_UNITS
        character*(*) REC_UNITS,LVL_UNITS
        character*(*) UNITS
        parameter (UNITS = 'units')
        parameter (LAT_UNITS = 'degrees_north')
        parameter (LON_UNITS = 'degrees_east')
        parameter (REC_UNITS = 'days')
        parameter (LVL_UNITS = 'm')
        character*(*) TEMP_UNITS,SALT_UNITS
        parameter (TEMP_UNITS='W/m3',SALT_UNITS='mm/day.m3')

C     The start and count arrays will tell the netCDF library where to
C     read our data.
        integer start(NDIMS), count(NDIMS)

C     In addition to the latitude and longitude dimensions, we will also
C     create latitude and longitude variables which will hold the actual
C     latitudes and longitudes. Since they hold data about the
C     coordinate system, the netCDF term for these is: "coordinate
C     variables."
       real lats,lons,tims(12),lvl(NZP1)
       integer lon_varid, lat_varid, rec_varid, lvl_varid,z_varid
       integer lat_dimid,lon_dimid,rec_dimid,lvl_dimid,z_dimid

C     We will read surface temperature and pressure fields. In netCDF
C     terminology these are called "variables."
        character*(*) TEMP_NAME,SALT_NAME
        parameter (TEMP_NAME='qcorr',SALT_NAME='scorr')
        integer temp_varid,salt_varid

C     Error handling.
       integer retval
       integer dimids(NDIMS),l,nm
       parameter(nm=12)   ! length of time in flux correction
       real qcorr(NX,NY,NZP1,nm)
       real scorr(NX,NY,NZP1,nm)

       lvl=zm
       do i=1,nm
        tims(i)=i
       end do
        write(*,*)'*** Now starts writing FCORR file'
!        call cpu_time(startcp)
!        write(*,*)'called cputime=',startcp
!       Create the file
        retval=nf_create(FILE_OUT,ior(nf_clobber,nf_64bit_offset),ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     Define the dimensions.
        retval = nf_def_dim(ncid, LAT_NAME, NY, lat_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, LON_NAME, NX, lon_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, LVL_NAME, NZP1, lvl_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, REC_NAME, nm, rec_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     Define the coordinate variables.
        retval = nf_def_var(ncid, LAT_NAME, NF_DOUBLE, 1, lat_dimid,
     +     lat_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, LON_NAME, NF_DOUBLE, 1, lon_dimid,
     +     lon_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, LVL_NAME, NF_DOUBLE, 1, lvl_dimid,
     +     lvl_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, REC_NAME, NF_DOUBLE, 1, rec_dimid,
     +     rec_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     Assign units attributes to coordinate variables.
        retval = nf_put_att_text(ncid, lat_varid, UNITS, len(LAT_UNITS),
     +     LAT_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, lon_varid, UNITS, len(LON_UNITS),
     +     LON_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, lvl_varid, UNITS, len(LVL_UNITS),
     +     LVL_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, rec_varid, UNITS, len(REC_UNITS),
     +     REC_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     The dimids array is used to pass the dimids of the dimensions of
C     the netCDF variables.
        dimids = (/lon_dimid,lat_dimid,lvl_dimid,rec_dimid/)
!       Define netcdf variable for Mixed layer depth (MLD)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, TEMP_NAME, NF_DOUBLE, 4, dimids,
     +     temp_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, SALT_NAME, NF_DOUBLE, 4, dimids,
     +     salt_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

!       Assign units attributes to the netCDF variables.
        retval = nf_put_att_text(ncid, temp_varid,UNITS,len(temp_UNITS),
     +     TEMP_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, salt_varid,UNITS,len(SALT_UNITS),
     +     SALT_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     End define mode.
        retval = nf_enddef(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     Write the coordinate variable data. This will put the latitudes
C     ,longitudes,depth and time of our data grid into the netCDF file.
        retval = nf_put_var_double(ncid, lat_varid, lats)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_double(ncid, lon_varid, lons)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_double(ncid, lvl_varid, lvl)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_double(ncid, rec_varid, tims)
        if (retval .ne. nf_noerr) call handle_err(retval)

        count = (/NX,NY,NZP1,nm/)
        start = (/1,1,1,1/)
        retval = nf_put_vara_double(ncid, temp_varid, start, 
     +                                           count,qcorr)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_vara_double(ncid, salt_varid, start,
     +                                           count,scorr)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     Close the file. This causes netCDF to flush all buffers and make
C     sure your data are really written to disk.
        retval = nf_close(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        write(*,*)'*** SUCCESS writing fcorr file  ', FILE_OUT
!        call cpu_time(finish)
!        print '("Time = ",f10.3," minutss.")',(finish-startcp)/60
        return
        end
!********** End of creating fcorrection output ***************************

!************************************************************
!***** reading input flux correction with Z for 12 months ***
!************************************************************
        SUBROUTINE read_fcorrwithz(FILE_INIT,qcorr_withz,scorr_withz)
!********************************************************
!*** This program reads the flux correction data from the file
!***** edited by Byju Pookkandy on 27-sep-2015 
!********************************************************
!******* MAIN PROGRAM ************************
C     This is part of the netCDF package.
        implicit none
        include 'netcdf.inc'
        include 'parameter.inc'

c       real :: startcp, finish !just to know the time taken for the run
C     This is the name of the data file we will read.
        character(len=40) FILE_INIT
!        parameter (FILE_INIT='fcorr_1.nc')
        integer ncid

C     Local variable
        integer NDIMS,NTIM
        parameter (NDIMS = 4,NTIM = 12)
        character*(*) LVL_NAME, LAT_NAME, LON_NAME, REC_NAME
        parameter (LVL_NAME = 'z')
        parameter (LAT_NAME = 'latitude', LON_NAME = 'longitude')
        parameter (REC_NAME = 'time')

C     The start and count arrays will tell the netCDF library where to
C     read our data.
        integer start(NDIMS), count(NDIMS)  ! for temperature and salinity

C     In addition to the latitude and longitude dimensions, we will also
C     create latitude and longitude variables which will hold the actual
C     latitudes and longitudes. Since they hold data about the
C     coordinate system, the netCDF term for these is: "coordinate
C     variables."
        real lati,longi,deps(NZP1),tims(12)
        integer lon_varid, lat_varid, lvl_varid, rec_varid

C     We will read surface temperature and pressure fields. In netCDF
C     terminology these are called "variables."
        character*(*)  TEMP_NAME,SAL_NAME
        parameter (TEMP_NAME='qcorr',SAL_NAME='scorr')
        integer temp_varid,sal_varid
        
C     Error handling.
        integer retval

       real qcorr_withz(NX,NY,NZP1,NTIM)
       real scorr_withz(NX,NY,NZP1,NTIM)

!**************************
C     Open the file to read temperature data
!**************************
        retval = nf_open(FILE_INIT, nf_nowrite, ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        write(*,*)'file opened = ',FILE_INIT
C     Get the varids of the latitude and longitude coordinate variables.
        retval = nf_inq_varid(ncid, LON_NAME, lon_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_inq_varid(ncid, LAT_NAME, lat_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
C     Read the latitude and longitude data.
        retval = nf_get_var_double(ncid, lat_varid, lati)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_get_var_double(ncid, lon_varid, longi)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
! get the varid and read depth
        retval = nf_inq_varid(ncid,LVL_NAME,lvl_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_get_var_double(ncid, lvl_varid, deps)
        if (retval .ne. nf_noerr) call handle_err(retval)

! get the varid and read time
        retval = nf_inq_varid(ncid,REC_NAME,rec_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_get_var_double(ncid, rec_varid, tims)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     Get the varids of temperature netCDF variables.
        retval = nf_inq_varid(ncid,TEMP_NAME, temp_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
C     Read the temperature data from the file, one
        count = (/NX, NY, NZP1,NTIM/)
        start = (/1,1,1,1/)
        retval = nf_get_vara_double(ncid,temp_varid,start,count,
     +                                               qcorr_withz)
        if (retval .ne. nf_noerr) call handle_err(retval)

C       open and get the varid of salinity data
        retval = nf_inq_varid(ncid,SAL_NAME, sal_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
C     Read the salinity data from the file, one
C     record at a time.
        count = (/NX, NY, NZP1,NTIM/)
        start = (/1,1,1,1/)
        retval = nf_get_vara_double(ncid,sal_varid,start,count,
     +                                               scorr_withz)
        if (retval .ne. nf_noerr) call handle_err(retval)

        retval = nf_close(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

!****************************     
C     If we got this far, everything worked as expected. Yipee!
!********************************************************************
!c      call cpu_time(finish)
!c        print '("Time = ",f10.3," minutss.")',(finish-startcp)/60
        return
        end subroutine
