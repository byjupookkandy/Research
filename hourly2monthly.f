        subroutine monthlymean(uvel,vvel,tout,sout,mixd
!this program convert the 6 hourly KPP model output to monethly means
!byju pookkandy 
        implicit none

        integer, parameter :: NX=1,NY=1,NRECS=720000, NZP1=40
        integer, parameter :: mons=6000
        integer i,j,nmm,cas
        real lats,lons,lvl(NZP1)

        real,allocatable,dimension(:,:,:,:):: uvel
        real,allocatable,dimension(:,:,:,:):: vvel
        real,allocatable,dimension(:,:,:,:):: tout
        real,allocatable,dimension(:,:,:,:):: sout
        real,allocatable,dimension(:,:,:,:):: den
        real,allocatable,dimension(:,:,:):: mixd
! dummy variable
        real,allocatable,dimension(:,:,:,:):: var
        real,allocatable,dimension(:,:):: dat

        real,allocatable,dimension(:,:):: uvel_mon
        real,allocatable,dimension(:,:):: vvel_mon
        real,allocatable,dimension(:,:):: tout_mon
        real,allocatable,dimension(:,:):: sout_mon
        real,allocatable,dimension(:,:):: den_mon
        real,allocatable,dimension(:):: mixd_mon

        allocate(uvel(NX,NY,NZP1,NRECS))
        allocate(vvel(NX,NY,NZP1,NRECS))
        allocate(tout(NX,NY,NZP1,NRECS))
        allocate(sout(NX,NY,NZP1,NRECS))
        allocate(den(NX,NY,NZP1,NRECS))
        allocate(mixd(NX,NY,NRECS))

        allocate(uvel_mon(NZP1,mons))
        allocate(vvel_mon(NZP1,mons))
        allocate(tout_mon(NZP1,mons))
        allocate(sout_mon(NZP1,mons))
        allocate(den_mon(NZP1,mons))
        allocate(mixd_mon(mons))

!dummy variable allocation
        allocate(var(NX,NY,NZP1,NRECS))
        allocate(dat(NZP1,mons))

        call read_input(uvel,vvel,tout,sout,mixd,den,lats,lons,lvl)

       ! write(*,'(10 f10.2)')lvl
        do cas=1,5
         if (cas.eq.1) then
           do i=1,NZP1
                var(1,1,i,:) = uvel(1,1,i,:)
           end do
         else if (cas.eq.2) then
           do i=1,NZP1       
                var(1,1,i,:) = vvel(1,1,i,:)
           end do
         else  if (cas.eq.3) then
           do i=1,NZP1
                var(1,1,i,:) = tout(1,1,i,:)
           end do
         else if (cas.eq.4) then
           do i=1,NZP1
                var(1,1,i,:) = sout(1,1,i,:)
           end do
         else 
           do i=1,NZP1 
                var(1,1,i,:) = den(1,1,i,:)
           end do
        end if

         do i=1,NZP1
           nmm=0
           do j=1,NRECS,120
              nmm=nmm+1
              dat(i,nmm)=sum(var(1,1,i,j:j+119))/120.
           end do
         end do
         if (cas.eq.1) uvel_mon = dat
         if (cas.eq.2) vvel_mon = dat
         if (cas.eq.3) tout_mon = dat
         if (cas.eq.4) sout_mon = dat
         if (cas.eq.5) den_mon = dat
        end do

        nmm=0
        do j=1,NRECS,120
           nmm=nmm+1
           mixd_mon(nmm)=sum(mixd(1,1,j:j+119))/120.
        end do
        call write_output(uvel_mon,vvel_mon,tout_mon,sout_mon,
     +                           mixd_mon,den_mon,lats,lons,lvl)

       
        deallocate(uvel)
        deallocate(vvel)
        deallocate(tout)
        deallocate(sout)
        deallocate(den)
        deallocate(mixd)

        deallocate(uvel_mon)
        deallocate(vvel_mon)
        deallocate(tout_mon)
        deallocate(sout_mon)
        deallocate(den_mon)
        deallocate(mixd_mon)

!dummy variable allocation
        deallocate(var)
        deallocate(dat)

        stop
        end

!***************************************************************


        SUBROUTINE read_input(uvel,vvel,tout,sout,mixd,den,
     +                                        lats,lons,lvl)

C     This is part of the netCDF package.
        implicit none
        include 'netcdf.inc'

c	real :: startcp, finish !just to know the time taken for the run
C     This is the name of the data file we will read.
C     This is the name of the data file we will read.
        character*(*) FILE_IN
        parameter (FILE_IN='KPPocean_0009.nc')
        integer ncid,i

C     Local variable
        integer, parameter :: NX=1, NY=1, NZP1=40
        integer, parameter ::  NDIMS=4, NRECS=720000
        character*(*) LAT_NAME, LON_NAME, REC_NAME,LVL_NAME
        parameter (LAT_NAME = 'latitude', LON_NAME = 'longitude')
        parameter (REC_NAME = 'time',LVL_NAME = 'z')

C     The start and count arrays will tell the netCDF library where to
C     read our data.
        integer start(NDIMS), count(NDIMS)
        integer starth(3), counth(3)   ! this is for hmix

        real lats,lons,tims(NRECS),lvl(NZP1)
        integer lon_varid, lat_varid, rec_varid, lvl_varid

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

        real uvel(NX,NY,NZP1,NRECS)
        real vvel(NX,NY,NZP1,NRECS)
        real tout(NX,NY,NZP1,NRECS)
        real sout(NX,NY,NZP1,NRECS)
        real mixd(NX,NY,NRECS)
        real den(NX,NY,NZP1,NRECS)

!**************************
C     Open the file to read temperature data
!**************************
        retval = nf_open(FILE_IN, nf_nowrite, ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)
	    write(*,*)'file opened = ',FILE_IN
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
	
! get the varid and read depth
        retval = nf_inq_varid(ncid,LVL_NAME,lvl_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_get_var_double(ncid, lvl_varid, lvl)
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
        count = (/NX, NY, NZP1,NRECS/)
	    start = (/1,1,1,1/)
	    retval = nf_get_vara_double(ncid,temp_varid,start,count,tout)
        if (retval .ne. nf_noerr) call handle_err(retval)

C	open and get the varid of salinity data
	retval = nf_inq_varid(ncid,SALT_NAME, salt_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_get_vara_double(ncid,salt_varid,start,count,sout)
        if (retval .ne. nf_noerr) call handle_err(retval)

C	open and get the varid of density data
        retval = nf_inq_varid(ncid,RHO_NAME, rho_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_get_vara_double(ncid,rho_varid,start,count,den)
        if (retval .ne. nf_noerr) call handle_err(retval)

C	open and get the varid of uvel data
        retval = nf_inq_varid(ncid,UVEL_NAME, uvel_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_get_vara_double(ncid,uvel_varid,start,count,uvel)
        if (retval .ne. nf_noerr) call handle_err(retval)

C	open and get the varid of vvel data
        retval = nf_inq_varid(ncid,VVEL_NAME, vvel_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_get_vara_double(ncid,vvel_varid,start,count,vvel)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     Get the varids of temperature netCDF variables.
        retval = nf_inq_varid(ncid,MLD_NAME, mld_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
C     Read the temperature data from the file, one
        counth = (/NX, NY, NRECS/)
        starth = (/1,1,1/)
        retval = nf_get_vara_double(ncid,mld_varid,starth,counth,mixd)
        if (retval .ne. nf_noerr) call handle_err(retval)

        retval = nf_close(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

	return 
        end subroutine

!********************************************************************
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
!********************************************************************
!Now routine for writing the output
        subroutine write_output(uvel,vvel,tout,sout,mixd,den,
     +                                        lats,lons,lvl)

!********************************************************
!***** this program reads the input forcing flux 
!*****  Byju Pookkandy on 20-Aug-2015 
!********************************************************
!******* MAIN PROGRAM ************************
C     This is part of the netCDF package.
        implicit none
        include 'netcdf.inc'

        integer, parameter :: NX=1, NY=1, NZP1=40
        integer, parameter :: NDIMS=4, NRECS=6000

        real :: startcp, finish !just to know the time taken for the run
C     This is the name of the data file we will read.
        character*(*) FILE_OUT
        parameter (FILE_OUT='kpp_monthly_500yr.nc')
        integer ncid,i

C     Local variable
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
       real lats,lons,tims(NRECS),lvl(NZP1)
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

       real uvel(NX,NY,NZP1,NRECS)
       real vvel(NX,NY,NZP1,NRECS)
       real tout(NX,NY,NZP1,NRECS)
       real sout(NX,NY,NZP1,NRECS)
       real mixd(NX,NY,NRECS)
       real den(NX,NY,NZP1,NRECS)


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
        retval = nf_def_dim(ncid, REC_NAME, NRECS, rec_dimid)
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

        count = (/NX,NY,NZP1,NRECS/)
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

	    counth = (/NX,NY,NRECS/)
        starth = (/1,1,1/)
        retval = nf_put_vara_double(ncid, mld_varid, starth, 
     +                                           counth,mixd)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     Close the file. This causes netCDF to flush all buffers and make
C     sure your data are really written to disk.
        retval = nf_close(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        write(*,*)'*** SUCCESS writing output file  ', FILE_OUT
        return
	end
!********** End of creating netcdf output ***************************
