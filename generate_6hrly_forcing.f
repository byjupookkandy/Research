      PROGRAM random_forcing
*******************************************************************
!tis porgram generate noise forcing AR(1) out of monthly variance (inter annual)
!gf -I/opt/local/include -L/opt/local/lib -lnetcdf -lnetcdff -fdefault-real-8 generate_forcing.f
! input netcdf files : NCEP.forcing.x160y30.mon.nc (1948-2014 monthly)
!                    : ncep-4xday.Clim_x161y31.nc (6 hourly climatology)
!                    : Need to provide lag-1 correlation values for each variable.
! Output data        : noise_forcing.x160y30.6hrly.nc
******************************************************************** 
      implicit none

      include 'parameter.inc'

      integer i,j,k,ijk,nyy,cas,l
      real mon(12,67),monmean(67),variance(12)
      real lag1c, lat,lon
      real, allocatable :: var(:)
      integer it1, it2
! Local variables
      real,allocatable,dimension(:,:,:):: taux
      real,allocatable,dimension(:,:,:):: tauy
      real,allocatable,dimension(:,:,:):: solar
      real,allocatable,dimension(:,:,:):: nsolar
      real,allocatable,dimension(:,:,:):: pme

!output variables
      real,allocatable,dimension(:):: taux_in
      real,allocatable,dimension(:):: taux_noise
      real,allocatable,dimension(:):: tauy_in
      real,allocatable,dimension(:):: tauy_noise
      real,allocatable,dimension(:):: solar_in
      real,allocatable,dimension(:):: solar_noise
      real,allocatable,dimension(:):: nsolar_in
      real,allocatable,dimension(:):: nsolar_noise
      real,allocatable,dimension(:):: pme_in
      real,allocatable,dimension(:):: pme_noise
      real,allocatable,dimension(:):: dat     ! dummy variable
      real,allocatable,dimension(:):: noise
      real,allocatable,dimension(:):: dat_noise   

      integer, parameter :: nclim=1440
      integer, parameter :: nrecs=804
      real :: txm(nclim),tym(nclim),swm(nclim),nswm(nclim),pem(nclim)
      real :: clm(nclim)  ! dummy variable for the above

      allocate(taux_in(totim*12))
      allocate(tauy_in(totim*12))
      allocate(solar_in(totim*12))
      allocate(nsolar_in(totim*12))
      allocate(pme_in(totim*12))
      allocate(dat(totim*12))
      allocate(noise(120))
      allocate(taux_noise(totim*12))
      allocate(tauy_noise(totim*12))
      allocate(solar_noise(totim*12))
      allocate(nsolar_noise(totim*12))
      allocate(pme_noise(totim*12))
      allocate(dat_noise(totim*12))

      allocate(taux(NX,NY,nrecs))
      allocate(tauy(NX,NY,nrecs))
      allocate(solar(NX,NY,nrecs))
      allocate(nsolar(NX,NY,nrecs))
      allocate(pme(NX,NY,nrecs))

! read the 6hourly climatology
      call read_climseries(txm,tym,swm,nswm,pem)
! read monthly time series data
      call read_timeseries(taux,tauy,solar,nsolar,pme)

      allocate(var(size(solar(1,1,:)))) !allocate memory for the dummy variable: monthly

      do cas=1,5  ! 5 cases for 5 variable, assigning the data to the dummy variable
      if (cas.eq.1) then
        var(:)=solar(1,1,:)
        call lag1correlation(var,nrecs,lag1c)
        clm = swm
        write(*,*)'solar lag1c=',lag1c
      end if
      if (cas.eq.2) then
        var(:)=nsolar(1,1,:)
        call lag1correlation(var,nrecs,lag1c)
        clm = nswm
        write(*,*)'nsolar lag1c=',lag1c
      end if
      if (cas.eq.3) then
        var(:)=taux(1,1,:)
        clm = txm
        call lag1correlation(var,nrecs,lag1c)
        write(*,*)'taux lag1c=',lag1c
      end if
      if (cas.eq.4) then
        var(:)=tauy(1,1,:)
        call lag1correlation(var,nrecs,lag1c)
        clm = tym
        write(*,*)'tauy lag1c=',lag1c
      end if
      if (cas.eq.5) then
        var(:)=pme(1,1,:)
        clm = pem
        call lag1correlation(var,nrecs,lag1c)
        write(*,*)'pme lag1c=',lag1c
      end if

! now segregate months
      do i=1,12
         k=0
         do j=i,nrecs,12
            k=k+1
            mon(i,k)=var(j)
         end do
      end do

      nyy=size(mon(1,:))

! computing the monthly (interannual) variance
      do i=1,12
        monmean(i)=sum(mon(i,:))/nyy
! compute the monthly variance of the noise, monthly variance of the
! real data scaled by the lag-1 auto-correlation
        variance(i)=(sum((mon(i,:)-monmean(i))**2)/(nyy-1)) 
     +                                       *(1.-lag1c**2)
      end do

! call the subroutine that generate 500yrsX6hrly random noise (0 mean, 1 stdv)
      k=1
      do i=1,totim/120
                ijk = i*120-120
        nyy=0    ! index for 1 yr 6hrly time step (1440)
        do l=1,12
        it2=0
        call randomnormal(noise)
          do j=ijk+1,ijk+120
          nyy=nyy+1
          it2=it2+1
          it1=nyy-1; if (it1 .eq. 0) it1 = 1440
          if (k.eq.1) then  ! only for first time step
            dat_noise(k)=sqrt(variance(l))*noise(it2)
            dat(k)      =dat_noise(k)+clm(nyy)
            k=k+1
          else
            if (cas.eq.3.or.cas.eq.5) then  !taux and emp
               dat_noise(k)=lag1c*dat_noise(k-1)+
     +                   sqrt(5*variance(l))*noise(it2)
               dat(k) = dat_noise(k)+clm(nyy)
               k=k+1
            else if (cas.eq.4) then  ! tauy
                dat_noise(k)=lag1c*dat_noise(k-1)+
     +                   sqrt(25*variance(l))*noise(it2)
                dat(k) = dat_noise(k)+clm(nyy)
                k=k+1
            else if (cas.eq.1) then         ! solar
                dat_noise(k)=lag1c*dat_noise(k-1)+
     +               sqrt(variance(l))*noise(it2)
                dat(k) = dat_noise(k)+clm(nyy)
                k=k+1
            else
                dat_noise(k)=lag1c*dat_noise(k-1)+
     +               sqrt(variance(l))*noise(it2)
                dat(k) = dat_noise(k)+clm(nyy)
                k=k+1
            end if
          end if
          end do
        end do
      end do
      if (cas.eq.1) then
         solar_in=dat
         solar_noise=dat_noise
      end if
      if (cas.eq.2) then
         nsolar_in=dat
         nsolar_noise=dat_noise
      end if
      if (cas.eq.3) then
         taux_in=dat
         taux_noise=dat_noise
      end if
      if (cas.eq.4) then
         tauy_in=dat
         tauy_noise=dat_noise
      end if
      if (cas.eq.5) then
         pme_in=dat
         pme_noise=dat_noise
      end if
      end do
! now write the output in to netcdf format
      CALL write_forces(taux_in,taux_noise,tauy_in,tauy_noise,
     +                    solar_in,solar_noise,nsolar_in,nsolar_noise,
     +                   pme_in,pme_noise)
      deallocate(taux_in)
      deallocate(tauy_in)
      deallocate(solar_in)
      deallocate(nsolar_in)
      deallocate(pme_in)
      deallocate(dat)
      deallocate(var)
      deallocate(noise)
      deallocate(taux_noise)
      deallocate(tauy_noise)
      deallocate(solar_noise)
      deallocate(nsolar_noise)
      deallocate(pme_noise)
      deallocate(dat_noise)

      deallocate(taux)
      deallocate(tauy)
      deallocate(solar)
      deallocate(nsolar)
      deallocate(pme)

      stop
      end        


      SUBROUTINE read_timeseries(taux,tauy,solar,nsolar,pme)
!********************************************************
!***** this program reads the input forcing flux 
!*****  Byju Pookkandy on 20-Aug-2015 
!********************************************************
!******* MAIN PROGRAM ************************
C     This is part of the netCDF package.
      implicit none
      include 'netcdf.inc'
      include 'parameter.inc'

c       real :: startcp, finish !just to know the time taken for the run
C     This is the name of the data file we will read.
      character*(*) FILE_FORCING
      parameter (FILE_FORCING='forcing.mon.nc')
      integer ncid

C     Local variable
      integer :: NDIMS,nt
      parameter(NDIMS=3,nt=804)

C     The start and count arrays will tell the netCDF library where to
C     read our data.
      integer start(NDIMS), count(NDIMS)

C     We will read surface temperature and pressure fields. In netCDF
C     terminology these are called "variables."
      character*(*) TAUX_NAME,TAUY_NAME
      parameter (TAUX_NAME='TAUX_IN',TAUY_NAME='TAUY_IN')
      character*(*) SOLAR_NAME,NSOLAR_NAME,PME_NAME
      parameter (SOLAR_NAME='SOLAR_IN',NSOLAR_NAME='NSOLAR_IN')
      parameter (PME_NAME='PMINUSE_IN')
      integer taux_varid,tauy_varid,solar_varid,nsolar_varid
      integer pme_varid
        
C     Error handling.
      integer retval
      real taux(NX,NY,nt)
      real tauy(NX,NY,nt)
      real solar(NX,NY,nt)
      real nsolar(NX,NY,nt)
      real pme(NX,NY,nt)

!       write(*,*)'started reading forces'
!**************************
C     Open the file to read temperature data
!**************************
      retval = nf_open(FILE_FORCING, nf_nowrite, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(*,*)'file opened = ',FILE_FORCING

      write(*,*)'Starts reading input forcing flux data'

      count = (/NX,NY,nt/)
      start = (/1,1,1/)

!Get the var ids of forcing flux TAUX
      retval = nf_inq_varid(ncid,TAUX_NAME, taux_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
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
!**************************     
C     If we got this far, everything worked as expected. Yipee!
!********************************************************************
c       call cpu_time(finish)
c        print '("Time = ",f10.3," minutss.")',(finish-startcp)/60
      return
      end
!********************************************************************
!********************************************************************
      SUBROUTINE read_climseries(taux,tauy,solar,nsolar,pme)
!********************************************************
!***** this program reads the input forcing flux 
!*****  Byju Pookkandy on 20-Aug-2015 
!********************************************************
!******* MAIN PROGRAM ************************
C     This is part of the netCDF package.
        implicit none
        include 'netcdf.inc'
        include 'parameter.inc'

c       real :: startcp, finish !just to know the time taken for the run
C     This is the name of the data file we will read.
        character*(*) FILE_FORC
        parameter(FILE_FORC='forcing.6hrly.clim.nc')
        integer ncid

C     Local variable
        integer NDIMS,nt
        parameter (nt=1440, NDIMS = 3)

C     The start and count arrays will tell the netCDF library where to
C     read our data.
        integer start(NDIMS), count(NDIMS)

C     We will read surface temperature and pressure fields. In netCDF
C     terminology these are called "variables."
        character*(*) TAUX_NAME,TAUY_NAME
        parameter (TAUX_NAME='TAUX_IN',TAUY_NAME='TAUY_IN')
        character*(*) SOLAR_NAME,NSOLAR_NAME,PME_NAME
        parameter (SOLAR_NAME='SOLAR_IN',NSOLAR_NAME='NSOLAR_IN')
        parameter (PME_NAME='PMINUSE_IN')
        integer taux_varid,tauy_varid,solar_varid,nsolar_varid
        integer pme_varid
        
C     Error handling.
        integer retval
       real taux(NX,NY,nt)
       real tauy(NX,NY,nt)
       real solar(NX,NY,nt)
       real nsolar(NX,NY,nt)
       real pme(NX,NY,nt)

!       write(*,*)'started reading forces'
!**************************
C     Open the file to read temperature data
!**************************
        retval = nf_open(FILE_FORC, nf_nowrite, ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        write(*,*)'file opened = ',FILE_FORC

        write(*,*)'Starts reading input forcing flux data'

        count = (/NX, NY, nt/)
        start = (/1,1,1/)

!Get the var ids of forcing flux TAUX
        retval = nf_inq_varid(ncid,TAUX_NAME, taux_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
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
c       call cpu_time(finish)
c        print '("Time = ",f10.3," minutss.")',(finish-startcp)/60
        return
        end
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
!Now routine for writing the output forces
        subroutine write_forces(taux_in,taux_noise,tauy_in,tauy_noise,
     +                    solar_in,solar_noise,nsolar_in,nsolar_noise,
     +                   pme_in,pme_noise)

!********************************************************
!***** this program writes the noise generated forcing 
!*****  Byju Pookkandy on 09-Nov-2015 
!********************************************************
!******* MAIN PROGRAM ************************
C     This is part of the netCDF package.
        implicit none
        include 'netcdf.inc'
        include 'parameter.inc'

        real :: startcp, finish !just to know the time taken for the run
C     This is the name of the data file we will read.
        character*(*) FILE_OUT
        parameter(FILE_OUT='noise.force.6hrly.nc')
        integer ncid,i

C     Local variable
        integer NDIMS , ntim
        parameter (NDIMS = 3, ntim=totim*12)
        character*(*) LAT_NAME, LON_NAME, REC_NAME,LVL_NAME
        parameter (LAT_NAME = 'latitude', LON_NAME = 'longitude')
        parameter (REC_NAME = 'time',LVL_NAME = 'z')

        character*(*) LAT_UNITS, LON_UNITS
        character*(*) REC_UNITS
        character*(*) UNITS
        parameter (UNITS = 'units')
        parameter (LAT_UNITS = 'degrees_north')
        parameter (LON_UNITS = 'degrees_east')
        parameter (REC_UNITS = 'hours since 1901-01-01 00:00 UTC')
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
        real lats,lons,tims(ntim)
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
        integer txn_varid,tyn_varid,son_varid,nsn_varid
        integer pmn_varid

C     Error handling.
        integer retval
        integer dimids(NDIMS)

! index variable
        integer k,l

        real taux_in(NX,NY,ntim)
        real tauy_in(NX,NY,ntim)
        real solar_in(NX,NY,ntim)
        real nsolar_in(NX,NY,ntim)
        real pme_in(NX,NY,ntim)

        real taux_noise(NX,NY,ntim)
        real tauy_noise(NX,NY,ntim)
        real solar_noise(NX,NY,ntim)
        real nsolar_noise(NX,NY,ntim)
        real pme_noise(NX,NY,ntim)

        do i=1,ntim
        tims(i)=i
        end do
        lats=1 ; lons=1
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
        retval = nf_def_dim(ncid, REC_NAME, NF_UNLIMITED, rec_dimid)
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
        retval = nf_def_var(ncid, TAUX_NAME, NF_DOUBLE, 3, dimids,
     +     taux_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, 'taux_noise', NF_DOUBLE, 3, dimids,
     +     txn_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, TAUY_NAME, NF_DOUBLE, 3, dimids,
     +     tauy_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, 'tauy_noise', NF_DOUBLE, 3, dimids,
     +     tyn_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, SWF_NAME, NF_DOUBLE, 3, dimids,
     +     solar_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, 'solar_noise', NF_DOUBLE, 3, dimids,
     +     son_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, LWF_NAME, NF_DOUBLE, 3, dimids,
     +     nsolar_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, 'nsolar_noise', NF_DOUBLE, 3, dimids,
     +     nsn_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, PME_NAME, NF_DOUBLE, 3, dimids,
     +     pme_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, 'pme_noise', NF_DOUBLE, 3, dimids,
     +     pmn_varid)
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

        count = (/NX,NY,ntim/)
        start = (/1,1,1/)
        retval = nf_put_vara_double(ncid, taux_varid, start, count,
     +                                                     taux_in)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_vara_double(ncid, txn_varid, start, count,
     +                                                   taux_noise)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_vara_double(ncid, tauy_varid, start, count,
     +                                                     tauy_in)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_vara_double(ncid, tyn_varid, start, count,
     +                                                  tauy_noise)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_vara_double(ncid, solar_varid, start, count,
     +                                                     solar_in)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_vara_double(ncid, son_varid, start, count,
     +                                                  solar_noise)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_vara_double(ncid, nsolar_varid, start, count,
     +                                                     nsolar_in)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_vara_double(ncid, nsn_varid, start, count,
     +                                                nsolar_noise)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_vara_double(ncid, pme_varid, start, count,
     +                                                     pme_in)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_vara_double(ncid, pmn_varid, start, count,
     +                                                 pme_noise)
        if (retval .ne. nf_noerr) call handle_err(retval)

C     Close the file. This causes netCDF to flush all buffers and make
C     sure your data are really written to disk.
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      write(*,*)'*** SUCCESS writing output file  ', FILE_OUT
!       call cpu_time(finish)
!       print '("Time = ",f10.3," minutss.")',(finish-startcp)/60
       return
      end
!********** End of writing netcdf output forces *********************
!---------------------------------------------------------
!generate random noise (normal distribution with mean=0
!and standard deviation = 1
!---------------------------------------------------------
       subroutine randomnormal(noise)

      implicit none
      include 'parameter.inc'

! this code is actaully from the following web site
! "http://rosettacode.org/wiki/Random_numbers"

       INTEGER :: i
       integer, parameter :: sz=120
       REAL    :: noise(sz), pi, temp_noise, mean = 0.0, sd = 1.

       pi = 4.0*ATAN(1.0)
!       CALL RANDOM_SEED()
       CALL RANDOM_NUMBER(noise) ! Uniform distribution

! Now convert to normal distribution
      DO i = 1, sz-1, 2
         temp_noise = sd * SQRT(-2.0*LOG(noise(i))) *
     +                COS(2*pi*noise(i+1)) + mean
         noise(i+1) = sd * SQRT(-2.0*LOG(noise(i))) *
     +          SIN(2*pi*noise(i+1)) + mean
         noise(i) = temp_noise
      END DO
      return
      end

!---------------------------------------------------------
!compute lag 1 autocorrelation
!---------------------------------------------------------

      subroutine lag1correlation(varg,nmm,lag1c)
!this code actually from the following link
!https://vibrationdata.wordpress.com/2012/11/17/autocorrelation-cross-correlation/
C********************************************************

      implicit none
!      include 'times.com'
      include 'parameter.inc'
      integer :: nmm
      real :: varg(nmm)
      real, ALLOCATABLE :: XC(:,:)
      real DT,DF,DUR, SR
      real MAXA,MAXT, lag1c
      real lag0var
C
      INTEGER I,IJK,J,K
      INTEGER N2,NH,NUM
      INTEGER stat_alloc

      NUM=size(varg)
      DUR=NUM-1
      DT=DUR/FLOAT(NUM-1)
      SR=1/DT

      NH=NINT(NUM/2.)
      DF=1/(NUM*DT)
      N2=2*NUM

      ALLOCATE(XC(N2-1,2))
!      allocate(varg(NUM))

! initialise the lead/lag time and correlation
      DO I=1,N2-1
         XC(I,1)=-DUR+(I-1)*DT    ! lead/lag time
         XC(I,2)=0.               ! correlation
      END DO

! lead time variance     
      DO I=1,NUM
        IJK=I+NUM-1
        DO J=1,NUM
            IF((I+J-1).GT.NUM)THEN
                EXIT
            ELSE
                XC(IJK,2)=XC(IJK,2)+varg(J)*varg(I+J-1)
            ENDIF
        ENDDO
      ENDDO

!lag time variance
      DO I=1,NUM-1
        IJK=NUM-I
        DO J=1,NUM
            IF((I+J)>NUM)THEN
                EXIT
            ELSE
               XC(IJK,2)=XC(IJK,2)+varg(J)*varg(I+J)
            ENDIF
        ENDDO
      ENDDO

      lag0var = XC(NUM,2)  ! variance? at lag 0
      XC(:,2)=XC(:,2)/lag0var ! Normalise with lag0 variance
      lag1c = XC(NUM+1,2)     ! lag 1 correlation

      DEALLOCATE(XC,stat=stat_alloc)
!      DEALLOCATE(varg,stat=stat_alloc)

      return
      end
