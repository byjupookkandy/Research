      PROGRAM ocn_model_1D
**************************************************************************
* 1D ocean model using the kpp mixing scheme of 
* Large et al, with his interface between the ocean and kpp scheme
* calls his subroutines
* init_ocn : to initialize the ocean model
* ocn_step : to update the model
*
* Also uses his parameter namelists and common blocks
* There maybe a lot of unnecessary variables and parameters here, but until
* we get something working it seems foolish to try and strip too much out.
* Written 19 August 2015
*
*  BYJU POOKKANDY
*************************************************************************

      IMPLICIT NONE

	include 'parameter.inc'
	include 'timocn.com'
	include 'times.com'
	include 'kprof_out.com'
	include 'restore.com'
	include 'fcorr_in.com'
	include 'ocn_paras.com'
	include 'ocn_state.com'

	real    ::  U(NZP1,NVEL), X(NZP1,2),lats,lons
	integer :: i,k,rep,mon,l,iterf,startf   ! loop variables
	real    :: dscale,diff(12),std,stdold
        integer :: er,mz,mk
        real    :: startcp, finish !just to know the time taken for the run
! Local variables
        real,allocatable,dimension(:,:,:,:):: temp
        real,allocatable,dimension(:,:,:,:):: sal
        real,allocatable,dimension(:,:,:):: taux
        real,allocatable,dimension(:,:,:):: tauy
        real,allocatable,dimension(:,:,:):: solar
        real,allocatable,dimension(:,:,:):: nsolar
        real,allocatable,dimension(:,:,:):: pme

        real,allocatable,dimension(:):: noise,nsolar_in,pme_in
        real,allocatable,dimension(:):: taux_in,tauy_in,solar_in
        real :: taux_var,tauy_var,solar_var,nsolar_var,pme_var
        real :: mean

! model output
        real,dimension(NX,NY,NZP1,120):: uvel
        real,dimension(NX,NY,NZP1,120):: vvel
        real,dimension(NX,NY,NZP1,120):: tout
        real,dimension(NX,NY,NZP1,120):: sout
        real,dimension(NX,NY,120):: mixd

        real,allocatable,dimension(:,:,:,:):: uvel_mon
        real,allocatable,dimension(:,:,:,:):: vvel_mon
        real,allocatable,dimension(:,:,:,:):: tout_mon
        real,allocatable,dimension(:,:,:,:):: sout_mon
        real,allocatable,dimension(:,:,:):: mixd_mon

! heat and salt flux correction storage
        real qcorr_withz(NX,NY,NZP1,12),scorr_withz(NX,NY,NZP1,12)
        real qcorr(NZP1,12),scorr(NZP1,12)
        real Qold(NZP1,12),Sold(NZP1,12) ! temporary storage for flux correction

        character(len=40) :: output_file, fcorr_file
        integer flen

        integer :: srno,runlength(50) ! objects in timstep.nml

        integer nmm
        real nyy
        real, allocatable :: var(:)
        real :: lag1c
        
        COMMON /save_fcorr_withz/ fcorr_withz , fcorr_saltz
        COMMON /clim_in/ sal_clim, temp_clim

! 1440 is 360 days x 4 times daily (6hrly)
        allocate(temp(NX,NY,NZP1,12))
        allocate(sal(NX,NY,NZP1,12))
	allocate(taux(NX,NY,totim*12))
        allocate(tauy(NX,NY,totim*12))
        allocate(solar(NX,NY,totim*12))
        allocate(nsolar(NX,NY,totim*12))
        allocate(pme(NX,NY,totim*12))


      CALL initialise(temp,sal,lats,lons) 
      write(*,*) '1D KPP model initialized with, ',NZP1,' levels'
!set observed annual mean temp/sal for restoring the fileds
      do i=1,NZP1
         temp_clim(i)=sum(temp(1,1,i,:))/12  ! obs annual mean temp
         sal_clim(i) =sum(sal(1,1,i,:))/12   ! obs annal mean sal 
      end do

      call init_env(dscale,lats)
      write(*,*)'model levels are set'
      CALL forcing(taux,tauy,solar,nsolar,pme)
      where(solar(1,1,:) .lt. 0.0 ) solar(1,1,:) = 0.0 !BP-added to avoide -ve noise
      where(nsolar(1,1,:) .gt. 0.0 ) nsolar(1,1,:) = 0.0
      WRITE(*,*) ' Called input forcing flux data'
! finish computing noise variance

      dtsec = 21600.  ! time period of input forcing
      ndtocn= 6.    ! number of ocean time step per update to the forcing
      dto   = dtsec/float(ndtocn)  ! ocean model time period
      ndtupdfcorr = INT(30*ndtocn*86400/dtsec)      ! time step to update flux correction
                               ! 30days x ndtocn x 4daily
      stdold=1 ! to check model drift  
      er = 0   ! error count for model drift
      
      open(10,file='timestep.nml')
      read(10,*)
      read(10,*)
        call cpu_time(startcp)
      do iterf = 1,50  ! iteration loop for FCORR
        !write(*,*)'called cputime=',startcp
       read(10,*,end=1001) srno,runlength(iterf)
      end do
1001  er=iterf-1
      do iterf=1,er
       write(*,*)
       write(*,*)'fcor loop=', iterf
       write(*,*)'run length =',runlength(iterf),'years'
       NRECS = 1440*runlength(iterf)

! initialise the ocean fields
       U(:,1)=0
       U(:,2)=0
       X(:,1)=temp(1,1,:,1)
       X(:,2)=sal(1,1,:,1)

       Sref  = (X(1,2)+X(NZP1,2))/2.
       SSref = Sref
       Ssurf = X(1,2)+Sref
       uref = U(1,1)
       vref = U(1,2)
       Tref = X(1,1)

       nyy = NRECS/1440
       nmm = nyy*12
!output allocate; NRECS=model output times steps
       allocate(uvel_mon(NX,NY,NZP1,nmm))
       allocate(vvel_mon(NX,NY,NZP1,nmm))
       allocate(tout_mon(NX,NY,NZP1,nmm))
       allocate(sout_mon(NX,NY,NZP1,nmm))
       allocate(mixd_mon(NX,NY,nmm))

       call init_ocn(U,X)

       if (iterf.eq.1) then   
          L_FCORR_WITHZ=.False. ! free run
       else 
          L_FCORR_WITHZ=.True.  ! add flux corr in model levels
       end if

       if (L_FCORR_WITHZ) then
          fcorr_file='Fcorr'
          flen=INDEX(fcorr_file,' ')-1
          write(fcorr_file(flen+1:flen+1),'(a)') '_'
          write(fcorr_file(flen+2:flen+5),'(i4.4)') iterf-1
          write(fcorr_file(flen+6:flen+8),'(3A)') '.nc'
! previous run fcorr .nc file
        call read_fcorrwithz(fcorr_file,qcorr_withz,scorr_withz)
          do l=1,12
            Qold(:,l) = qcorr_withz(1,1,:,l)
            Sold(:,l) = scorr_withz(1,1,:,l)
          end do
       else
          Qold(:,:) = 0.
          Sold(:,:) = 0.
       end if
       
! call normal random distribution -> noise forcing
       rep=1
       k=0; mk=0;  ! k is the index for model output time step
       write(*,*)'***entering into ocstep loop******'
!      write(*,*)'NRECS=',X(:,1)
       DO ntime=1,NRECS*ndtocn
!        write(*,*)'time step=',ntime
        if (MOD(ntime-1,12*ndtupdfcorr).eq.0) mon=0 ! repeat fcorr after 12
                                                                     ! months
        if (MOD(ntime-1,ndtupdfcorr).eq.0) then     ! update fcorr 
           mon=mon+1
           fcorr_withz(:) = qcorr_withz(1,1,:,mon)
           fcorr_saltz(:) = scorr_withz(1,1,:,mon)
        end if
!        if (MOD(ntime-1,1440*ndtocn).eq.0) rep=0 
        if (MOD(ntime-1,ndtocn).eq.0) then
           CALL fluxes(taux(NX,NY,rep),tauy(NX,NY,rep), ! testing with 4xdaily forcing
     +              solar(NX,NY,rep),nsolar(NX,NY,rep),
     +              pme(NX,NY,rep))
           rep=rep+1
        end if

! main ocean integration subroutine 
        CALL ocnstep(U,X)
! store the data for each forcing time step
        if (MOD(ntime,ndtocn).eq.0) then
           k=k+1
           uvel(NX,NY,:,k)=U(:,1)
           vvel(NX,NY,:,k)=U(:,2)
           tout(NX,NY,:,k)=X(:,1)
           sout(NX,NY,:,k)=X(:,2)
           mixd(NX,NY,k)=hmix
        end if
        if (MOD(ntime,ndtupdfcorr).eq.0) then
           mk=mk+1
           do mz=1,NZP1
              uvel_mon(NX,NY,mz,mk)=sum(uvel(NX,NY,mz,:))/120
              vvel_mon(NX,NY,mz,mk)=sum(vvel(NX,NY,mz,:))/120
              tout_mon(NX,NY,mz,mk)=sum(tout(NX,NY,mz,:))/120
              sout_mon(NX,NY,mz,mk)=sum(sout(NX,NY,mz,:))/120
              uvel(NX,NY,mz,:)=0.; vvel(NX,NY,mz,:)=0.;
              tout(NX,NY,mz,:)=0.; sout(NX,NY,mz,:)=0.;
           end do
              mixd_mon(NX,NY,mk)=sum(mixd(NX,NY,:))/120
              mixd(NX,NY,:)=0.; k=0;
!           uvel=0; vvel=0; tout=0; sout=0; mixd=0; k=0;
        end if
       ENDDO   ! end of ocean step loop
       write(*,*)'***ocstep loop completed******'
!      WRITE(*,*) 'KPP: Successful termination of the model ',
!     +     'integration'
!**********************************************************
! compute flux correction
       call day2monclm(nyy,nmm,tout_mon,sout_mon,
     +                  tclm,sclm,tanual,sanual) ! get the simulated output into
                                        ! long term monthly average

! compute the flux corr based on observed and model climatology
       call fcorrection(temp,sal,tclm,sclm,tanual,sanual,Qold,Sold,
     +                                       iterf,qcorr,scorr)
      
! create the flux correction output at each iteration
       fcorr_file='Fcorr'
       flen=INDEX(fcorr_file,' ')-1
       write(fcorr_file(flen+1:flen+1),'(a)') '_'
       write(fcorr_file(flen+2:flen+5),'(i4.4)') iterf
       write(fcorr_file(flen+6:flen+8),'(3A)') '.nc'
       call write_fcor(qcorr,scorr,lats,lons,fcorr_file)

! create the kpp output file name and write instantaneous output
!       output_file='KPPocean'
!       flen=INDEX(output_file,' ')-1
!       write(output_file(flen+1:flen+1),'(a)') '_'
!       write(output_file(flen+2:flen+5),'(i4.4)') iterf
!       write(output_file(flen+6:flen+8),'(3A)') '.nc'
!       CALL write_output(uvel,vvel,tout,sout,mixd,lats,lons,
!    +                                                 output_file) 

! create the kpp output file name and write monthly mean output
       output_file='KPPocean'
       flen=INDEX(output_file,' ')-1
       write(output_file(flen+1:flen+1),'(a)') '_'
       write(output_file(flen+2:flen+5),'(i4.4)') iterf
       write(output_file(flen+6:flen+8),'(3A)') '.nc'
       CALL write_output(uvel_mon,vvel_mon,tout_mon,sout_mon,mixd_mon,
     +                                    nmm,lats,lons,output_file)
! Testing for  sst convergence
       diff(:) = (tclm(1,:)-temp(1,1,1,:))**2
       std = sqrt(sum(diff)/12)
       write(*,*)'stdv SST=',std
!       if (std.le.0.2) go to 1001

!       if (std .gt. stdold) er=er+1
       stdold=std
!       if (er.ge.3) then
!         write(*,*)'*************************'
!         write(*,*)'error in model simulation'
!         write(*,*)'it has a tendency to drift more'
!         write(*,*)'*************************'
!         go to 1001
!       end if
!test for fcorr convergence
!       diff(:) = (Qold(1,:)-qcorr(1,:))**2
!       std = sqrt(sum(diff)/12)
!       write(*,*)'diff max=',std
!       if (std.le.1.5E-2) go to 1001
      call cpu_time(finish)
      print '("Time = ",f15.3," minuts.")',(finish-startcp)/60

!de-allocate output fields
       deallocate(uvel_mon)
       deallocate(vvel_mon)
       deallocate(tout_mon)
       deallocate(sout_mon)
       deallocate(mixd_mon) 

      end do   ! end of fcorr iteration loop

      WRITE(*,*) 'KPP: Successful termination of the model ',
     +     'integration'

      deallocate(temp)
      deallocate(sal)
      deallocate(taux)
      deallocate(tauy)
      deallocate(solar)
      deallocate(nsolar)
      deallocate(pme)

      END      !PROGRAM ocn_model_1D 
!************************************************************************
!************************************************************************
!****************routines to compute flux correction*********************

        subroutine day2monclm(nyy,nmm,tout_mon,sout_mon,
     +                                       tclm,sclm,tanual,sanual)
!here we convert the daily output to monthly climatology
!Also computes the annual mean at the end of the simulation
        implicit none

        include 'parameter.inc'
        include 'fcorr_in.com'
        include 'times.com'

        integer nmm
        real nyy ! totla number of years and months 

! input variables
        real tout_mon(NX,NY,NZP1,nmm),sout_mon(NX,NY,NZP1,nmm)

! local and out variables
        real tmon(NZP1,nmm),smon(NZP1,nmm)
        real stemp(NZP1),ssalt(NZP1)

        integer i,j,k
        integer dt

! monthly to monthly climatology
      do j=1,12
         stemp(:)=0.
         ssalt(:)=0.
         do k=j,nmm,12
            stemp(:)=stemp(:)+tout_mon(1,1,:,k)
            ssalt(:)=ssalt(:)+sout_mon(1,1,:,k)
         end do
         tclm(:,j)=stemp(:)/nyy
         sclm(:,j)=ssalt(:)/nyy
      end do
! compute annual mean for the last simulations
      finalt = INT(nyy-nyy/10)*12+1
      dt     = (nmm-finalt)+1

      do k=1,NZP1
         tanual(k)=sum(tout_mon(1,1,k,finalt:nmm))/dt
         sanual(k)=sum(sout_mon(1,1,k,finalt:nmm))/dt
      end do
      write(*,*)'tanual=',tanual(1),tanual(40)
     
      return
      end

!***********************************************************
! the following routines compute the flux correction for the model
! output relative to the departure from an observed profile

        subroutine fcorrection(temp,sal,tclm,sclm,tanual,sanual,
     +                              Qold,Sold,n,qcorr,scorr)

        include 'parameter.inc'
        include 'vert_pgrid.com'
        include 'timocn.com'
        include 'times.com'
        include 'fcorr_in.com'

! local variable
        real qcorr(NZP1,12),scorr(NZP1,12)
        real Qqs(NZP1,12),Qss(NZP1,12)
        real qanual(NZP1),qsanual(NZP1)
        real cp,rho,tmean,smean
        real Sig0,Sig
! DIDI 
        real xdt, xdtc
        integer it1,it2

        logical KapFlg

! input variable
        real sal(NX,NY,NZP1,12),temp(NX,NY,NZP1,12) 
        integer dt,n
! previous flux correction 
        real Qold(NZP1,12),Sold(NZP1,12)

        KapFlg=.False.
       
! compute the seasonal mean flux correction
!=========================================
! compute the seasonal mean flux correction
!=========================================
        dt = 30*86400
        do l=1,12
            it2 = l+1;    if (it2 .gt. 12) it2 = 1
            it1 = l-1;    if (it1 .eq. 0) it1 = 12
         do i=1,NZP1
! specific heat capacity - observed
           cp=CPSW(sal(NX,NY,i,l),temp(NX,NY,i,l),-zm(i))
! compute the density -observed
           call Sig80(sal(NX,NY,i,l),temp(NX,NY,i,l),0.,
     +                                       KapFlg,Sig0,Sig)
           rho=1000.+Sig0
! heat flux correction in W/m3 (seasonal)
!       if ((NRECS/1440).le.5.) then    ! no seasonal correction beyond...
!           if (n.le.2) then                ! XX years
! next step xdtc/xdt computes the tendency of the variable for the corresponding
! month compared to the month before and after
           xdtc = (temp(NX,NY,i,it2)-temp(NX,NY,i,it1)) / 2
           xdt = (tclm(i,it2)-tclm(i,it1)) / 2
           Qqs(i,l)  =rho*cp*(xdtc-xdt)/dt

           xdtc = (sal(NX,NY,i,it2)-sal(NX,NY,i,it1)) / 2
           xdt = (sclm(i,it2)-sclm(i,it1)) / 2
           Qss(i,l)  =(xdtc-xdt)/dt

!          else
!          Qqs(i,l)  =0.
!           Qss(i,l)  =0.
!          end if

         end do
        end do

!no seasonal corrections at some bottom levels
       Qqs(10:NZP1,:) = 0.0
!       Qss = 0.0
       Qss(2:NZP1,:) = 0.0 !BP added salinity seasonal correction at the surface only

         qcorr = Qold + Qqs
         scorr = Sold + Qss
! eliminate the seasonal fcorr residual forcing by making annual mean 
! seasonal fcorr = 0
        do i=1,NZP1
           qcorr(i,:)=qcorr(i,:)-sum(Qqs(i,:))/12
           scorr(i,:)=scorr(i,:)-sum(Qss(i,:))/12
        end do
!        write(*,*)'surface Qanual=',sum(qcorr(1,:))/12
         
! compute the annual mean flux correction
        dt = 86400*NRECS*360/1440
        do i=1,NZP1
           tmean=sum(temp(1,1,i,:))/12  ! obs annual mean temp
           smean=sum(sal(1,1,i,:))/12   ! obs annal mean sal
           cp=CPSW(smean,tmean,-zm(i))
           call Sig80(smean,tmean,0.,KapFlg,Sig0,Sig)
           rho=1000.+Sig0
           qanual(i)=rho*cp*(tmean-tanual(i))/dt
           qcorr(i,:)= qcorr(i,:)+qanual(i)
           qsanual(i)=(smean-sanual(i))/dt
           scorr(i,:)= scorr(i,:) + qsanual(i)
        end do
!  DIDI:no flux correction in last level
        qcorr(NZP1,:) = 0.0
        scorr(NZP1,:) = 0.0

        return
        end
