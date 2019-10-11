*********************************************************************
      subroutine vmix(U,X,hmixn,kmixn)
c  Interface between 1-d model and vertical mixing

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'parameter.inc'
      include 'constants.com'
      include 'times.com'
      include 'vert_pgrid.com'
      include 'location.com'
      include 'proc_swit.com'
      include 'proc_pars.com'
      include 'ocn_state.com'
      include 'ocn_paras.com'
      include 'flx_sfc.com'
      include 'flx_profs.com'
      include 'kprof_in.com'
      include 'kprof_out.com'
      include 'dble_diff.com'

c inputs including those from parameter.inc
      real U(nzp1,nvel),X(nzp1,nsclr)

      real hmixn          ! boundary layer depth (m)
      integer kmixn       ! hmixn index 

c local
      real Shsq(nzp1)    ! (local velocity shear)^2       (m/s)^2
      real dVsq(nzp1)    ! (velocity shear re sfc)^2      (m/s)^2
      real dbloc(nz)     ! local delta buoyancy            (m/s^2)
      real Ritop(nz)     ! numerator of bulk Richardson Number (m/s)^2
      real alphaDT(nz)   ! alpha * DT across interfaces
      real betaDS(nz)    ! beta  * DS across interfaces
      real epsilon,epsln
      real alpha,beta,exppr
      real sigma,sigma0
      real cpsw
      real tau
      real zref,wz,bref
      integer k,n,kl
      real del
      real dlimit,vlimit
c     logical LKPP       ! kpp boundary layer mixing switch
      integer jwtype
      integer jerl(12)
c          month  1   2   3   4   5   6   7   8   9   10  11  12
      data jerl / 2 , 2 , 2 , 3 , 3 , 3 , 4 , 4 , 4 , 4 , 3 , 2 /


      epsilon = 0.1
      epsln   = 1.e-20
      jerlov=1
      Sice=4.
c                        ensure that bottom layer isn't zero thickness
      hm(nzp1) = AMAX1(hm(nzp1),epsln)
c                        find the jerlov water type for this month 

      if(jerlov.lt.1) then
        write(nuerr,*) 'Make sure you have got time in days'
      else
         jwtype = jerlov
      endif

c calculate density of fresh water and brine in surface layer
      alpha = 1.
      beta  = 1.
      exppr = 0.0
      call ABK80(0.0,X(1,1),-zm(1),alpha,beta,exppr,sigma0,sigma)
      rhoh2o    = 1000. + sigma0
      call ABK80(Sice,X(1,1),-zm(1),alpha,beta,exppr,sigma0,sigma)
      rhob      = 1000. + sigma0
 
c     calculate temperature and salt contributions of buoyancy gradients  
c      calculate buoyancy profile (m/s**2) on gridlevels
      do k=1,nzp1 
         call ABK80(X(k,2)+Sref,X(k,1),-zm(k),
     +        alpha,beta,exppr,sigma0,sigma)
         rho(k)    = 1000. + sigma0
         CP(k)     = CPSW(X(k,2)+Sref,X(k,1),-zm(k))
         talpha(k) = alpha
         sbeta(k)  = beta
         buoy(k)   = - grav * sigma0 / 1000.
      end do

         rho(0)    = rho(1) 
         CP(0)     = CP(1)   
         talpha(0) = talpha(1)   
         sbeta(0)  = sbeta(1)           

c Call to ntflx, put here to allow removal of diagnostic call to vmix
c and to ensure the most recent cp,rho used (consistent with other 
c surface fluxes?)
      call ntflx
c calculate kinematic surface momentum fluxes
      wU(0,1) = -sflux(1,5,0) / rho(0)
      wU(0,2) = -sflux(2,5,0) / rho(0)
      tau     = sqrt( sflux(1,5,0)**2 + sflux(2,5,0)**2 ) +1.e-16
c  1.e-16 added to stop subsequent division by zero if tau=0.0
      ustar   = sqrt( tau / rho(0) )
 
c total turbulent kinematic temperature flux (C m/s)
      wX(0,1)  = - sflux(4,5,0) /rho(0) /CP(0)
c total turbulent kinematic salinity flux (o/oo m/s)
      wX(0,2) =Ssurf*sflux(6,5,0)/rhoh2o+(Ssurf-Sice)*sflux(5,5,0)/rhob

c calculate total kinematic surface buoyancy flux (m**2/s**3)
      B0 = - grav*(talpha(0)*wX(0,1) -sbeta(0)*wX(0,2) )
      wX(0,NSP1) =  - B0
      B0sol = grav * talpha(0) * sflux(3,5,0) /(rho(0) * CP(0))

c     calculate temperature and salt contributions of buoyancy gradients
c               on interfaces for double diffusion
      do n = 1,nz
         alphaDT(n) =0.5 *(talpha(n)+talpha(n+1))*(X(n,1) -X(n+1,1))
         betaDS(n)  =0.5 *(sbeta(n) + sbeta(n+1))*(X(n,2) -X(n+1,2))
      end do

c                      compute buoyancy and shear profiles
      do 115  n = 1,nz
       zref =  epsilon * zm(n) 
c          compute reference buoyancy and velocity
       wz    = AMAX1(zm(1),zref) 
       uref  = U(1,1) * wz / zref
       vref  = U(1,2) * wz / zref
       bref  = buoy(1)* wz / zref
       do 125 kl = 1,nz
         IF(zref.ge.zm(kl)) go to 126
         wz = AMIN1(zm(kl)-zm(kl+1),zm(kl)-zref) 
         del = 0.5 * wz / (zm(kl) - zm(kl+1))
         uref =uref -wz*( U(kl,1)+del*(U(kl+1,1)- U(kl,1)))/zref
         vref =vref -wz*( U(kl,2)+del*(U(kl+1,2)- U(kl,2)))/zref
         bref =bref -wz*(buoy(kl)+del*(buoy(kl+1)-buoy(kl)))/zref
 125   continue
 126   continue

       Ritop(n) = (zref - zm(n)) * (bref - buoy(n))
c     NPK Additions (25/9/2008). Prevent Ritop from going negative.
       IF (Ritop(n) .lt. 0) Ritop(n) = epsln
       dbloc(n) = buoy(n) - buoy(n+1)
       dVsq(n)  = (Uref - U(n,1))**2 + (Vref - U(n,2))**2
       Shsq(n)  = (U(n,1)-U(n+1,1))**2 + (U(n,2)-U(n+1,2))**2
 115  continue

      call kppmix(zm, hm , Shsq, dVsq,
     $            ustar , B0    , B0sol , alphaDT, betaDS ,
     $            dbloc , Ritop , f , jwtype,
     $            difm  , difs  , dift  , ghat , hmixn, kmixn)

! bottom flux turned off
!      if(LNBFLX) then 
!        dlimit = 0.0
!        vlimit = 0.0
!        do n=1,nsclr
!        wxNT(nz,n) = 0.0
!        enddo
!      else
        dlimit = 0.00005
        vlimit = 0.0001
!      endif
      do k=nz,nzp1
        difm(k) = vlimit
        difs(k) = dlimit
        dift(k) = dlimit
      enddo
        ghat(nz) = 0.0


      return
      end
**********************************************************************

