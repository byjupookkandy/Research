      SUBROUTINE  ocnstep(U,X)
!c-----------------------------------------------------------------------
!c Originally ~/KPP/LARGE/ocn.f (from Bill Large)
!c  Modified by SJW to remove those bits which a specific to Bill's 
!c work and just leave the bit we need
!c Started 18/03/02
!c-----------------------------------------------------------------------

!c     Main driver for ocean module.
!c     Integration is performed only on the permanent grid
!c     Written   3 Mar 1991 - WGL
!c     Modified  5 Jun 1992 - jan : implicit scheme
!c              16 Nov      - jan : latest version
!c              16 Nov 1994 - wgl : new KPP codes no temporary grid/

      IMPLICIT NONE
	include 'parameter.inc'
	include 'constants.com'
	include 'times.com'
	include 'timocn.com'
	include 'vert_pgrid.com'
	include 'proc_swit.com'
	include 'ocn_state.com'
	include 'ocn_paras.com'
	include 'ocn_energy.com'
	include 'flx_sfc.com'
	include 'flx_profs.com'
	include 'flx_paras.com'
	include 'kprof_out.com'
	include 'location.com'
c Input/Output
      real U(NZP1,NVEL), X(NZP1,NSCLR)
c
c Local Common Blocks
c
      real hmixd(0:1),     ! storage arrays for extrapolations
     +     Us(NZP1,NVEL ,0:1), ! ..      ..     ..  ..
     +     Xs(NZP1,NSCLR,0:1)  ! ..      ..     ..  ..
      integer old,new ! extrapolation index for Us,Xs,hmixd
      common/ saveUXh / 
     +     old,new,Us,Xs,hmixd
!c
!c Local'
!c
      real Un(NZP1,NVEL),       ! new profiles
     +     Xn(NZP1,NSCLR),      ! ..  ..
     +     hmixe,         ! estimated hmix (integration input )
     +     hmixn,         ! new comp. hmix (    ..      output)
     +     tol                  ! tolerance in hmix iteration
      real Uo(NZP1,NVEL),
     +     Xo(NZP1,NSCLR)
      real Ux(NZP1,NVEL),       ! Additional variables to provide
     +     Xx(NZP1,NSCLR)       ! smoothing in the iteration.
      real lambda                ! Factor to control smoothing
      integer
     +     iter,iconv                ! number of iterations
      integer kmixe,kmixn
!c
!c More Local Variables (to make implicit none)'
!c
      real deltaz,rhonot
      integer k,l,n

      data lambda /0.5/
c initialise and read the constants name list
      spd=86400.                ! secs/day
      dpy=360.                  ! days/year
      twopi=8*atan(1.)          ! 2pi
      onepi=twopi/2.            ! pi
      grav=9.816                ! gravity
      vonk=0.4                  ! Von Karman's constant
      TK0=273.15                ! Kelvin of 0degC
      sbc=5.67e-8               ! Stefan Boltzmann Constant
      epsw=1.0                  ! cor.fac for departure of H2O from B.body
      albocn=0.06               ! albedo for seawater
      sice=4.0                  ! salinity of ice(?)
      EL=2.50e6                 ! Latent heat of evap. at 0C (or constant)
      SL=2512200.               ! Latent heat of evap for ice
      FL=334000.                ! Latent heat of fusion for ice
      FLSN=FL                   ! Latent heat of fusion for snow
      
!      Fix U/V => shear instability due to instant TAU
       U(:,1) = 0.
       U(:,2) = 0. 
c     Copy old profiles (to allow 3D U,X)
      DO k=1,NZP1
         DO l=1,NVEL
            Uo(k,l)=U(k,l)
         ENDDO
         DO l=1,NSCLR
            Xo(k,l)=X(k,l)
         ENDDO
       ! write(*,*)Xs(k,1,0),Xs(k,1,1)
      ENDDO
c Estimate new profiles by  extrapolation
      do 20 k=1,NZP1
         do 22 l=1,NVEL
            Un(k,l)=2.*Us(k,l,new)-Us(k,l,old)
            Ux(k,l)=Un(k,l)
 22      continue
         do 24 l=1,NSCLR
            Xn(k,l)=2.*Xs(k,l,new)-Xs(k,l,old)
            Xx(k,l)=Xn(k,l)
 24      continue
 20   continue      


c Iteration loop for semi-implicit integration
c Reset iteration counter
      iter = 0
      iconv=0
c This loop controls the number of compulsory iterations
c The original value was 2, using (the upper value of the loop +2)
c added by SJW (17 Jan 03) to try an alleviate some non-convergences
      DO iter=0,2
        DO k=1,NZP1
            DO l=1,NVEL
               Un(k,l)=lambda*Ux(k,l)+(1-lambda)*Un(k,l)
               Ux(k,l)=Un(k,l)
            ENDDO
            DO l=1,NSCLR
               Xn(k,l)=lambda*Xx(k,l)+(1-lambda)*Xn(k,l)
               Xx(k,l)=Xn(k,l)
           ENDDO
        ENDDO
        call vmix(Un,Xn,hmixe,kmixe)
        call ocnint(NZ,zm,hm,dm,1,kmixe,Uo,Xo,Un,Xn)
      ENDDO
c The original code can be restored by reseting iter=1 and removing the  
c above loop  
!      iter=1

!      IF (LKPP) THEN
 45      continue
         DO k=1,NZP1
            DO l=1,NVEL
               Un(k,l)=lambda*Ux(k,l)+(1-lambda)*Un(k,l)
               Ux(k,l)=Un(k,l)
            ENDDO
            DO l=1,NSCLR
               Xn(k,l)=lambda*Xx(k,l)+(1-lambda)*Xn(k,l)
               Xx(k,l)=Xn(k,l)
            ENDDO
         ENDDO
         call vmix(Un,Xn,hmixn,kmixn)
         call ocnint(NZ,zm,hm,dm,1,kmixn,Uo,Xo,Un,Xn)
         iter = iter + 1
         
!c     check iteration for convergence
         tol = hmixtolfrac*hm(kmixn)
         if(kmixn.eq.NZP1) tol = hmixtolfrac*hm(NZ)
         if(abs(hmixn-hmixe).gt.tol)  then
!c     Uncommeting the following the lines iconv=0 to IF (iconv ...)
!c     will make the model do two consecutive tests for convergence of the 
!c     hmix (added by SJW 17 Jan 03). This did not work well in testing for
!c     long timestep, high resolution (the model generally failed to satisfy the 
!c     convergence test on two consecutive iterations.
            iconv=0
         ELSE
            iconv=iconv+1
         ENDIF
         IF (iconv .lt. 3) THEN
            if (iter.lt.itermax) then
               hmixe = hmixn
               kmixe = kmixn
               goto 45
            else
!c     use shallower hmix
               if(hmixn.gt.hmixe) then
                  hmixe = hmixn ! comment out for hmix data
                  kmixe = kmixn ! ..      ..  ..  hmix data
                  goto 45       ! ..      ..  ..  hmix data 
               endif  
            endif
         endif
!     ENDIF   ! end of if(LKPP)
      
!c  Output  Results from permanent grid iterations to common.inc
!c Compute diagnostic fluxes for writing to dat file
      do k=1,NZ
        deltaz = 0.5*(hm(k)+hm(k+1))
        do n=1,NSCLR
           wX(k,n)=-difs(k)*((Xn(k,n)-Xn(k+1,n))/deltaz- 
     +          ghat(k)*wX(0,n))
	end do
        wX(k,nsp1)= grav * 
     +       (talpha(k)*wX(k,1) - sbeta(k) * wX(k,2))
        do n=1,NVEL
           wU(k,n)= -difm(k)* (Un(k,n)-Un(k+1,n))/deltaz
	end do 
        end do
 
!c Compute energetics
      rhonot = 1026.
      Eflx = 0.5 * ( (Uo(1,1) + Un(1,1)) * sflux(1,5,0) +
     &               (Uo(1,2) + Un(1,2)) * sflux(2,5,0) )
      Esnk =-0.5*rhonot* ( (Uo(NZ,1) + Un(NZ,1)) * wU(NZ,1) +
     &                        (Uo(NZ,2) + Un(NZ,2)) * wU(NZ,2) )
      Ptke = 0.0
!c     use "amax1" to prevent "underflow" in single precision
      do k=1,NZ-1
      Ptke= Ptke- 0.5*( amax1(wU(k,1),1.E-10)* 
     &        (rhonot   * (Uo(k  ,1) + Un(k  ,1)) -
     &        rhonot   * (Uo(k+1,1) + Un(k+1,1)) ) +
     &        amax1(wU(k,2),1.E-10)*
     &        (rhonot   * (Uo(k  ,2) + Un(k  ,2)) -
     &        rhonot   * (Uo(k+1,2) + Un(k+1,2)) ) )
       end do
      Tmke = 0.0
      do k=1,NZP1
          rmke(k) = 0.5 * rhonot * (Un(k,1)**2 + Un(k,2)**2) * hm(k)
          Tmke  = Tmke + rmke(k)
      end do

!c Set new profiles
      do k=1,NZP1         ! values at NZP1 only change for slab ocean
         do  n=1,NVEL
          U(k,n) = Un(k,n)
	 end do
         do  n=1,NSCLR
          X(k,n) = Xn(k,n)
	 end do
!	write(*,*)X(k,1)
      end do

!DIDI: TEST FIX U/VVEL
!       U(:,1) = 0.
!       U(:,2) = 0. 
!c Set correct surface values, and CP and rho profiles for new profiles
!c Get latest profiles
!c Removed to ensure that hmix,diff,cp,rho diagnostics are consistent with 
!c those used in the final iteration of the timestep 
!c (CP,RHO are not quite correct for updated values, but are the values
!c  by the integration) (SJW 16/01/03)
!c Need to consider improving convergence test!!! (SJW 16/01/03)

!c The final call to vmix is removed to ensure that the diffusion and 
!c boundary layer profiles in the diagnostics are the ones used to calculate
!c the fluxes, as it stands at the moment this means that the CP and rho are
!c also the values used in the timestepping not the values appropriate to the
!c S,T at the new time level.  
      hmix = hmixn
      kmix = kmixn
      uref = U(1,1)
      vref = U(1,2)
      Tref = X(1,1)
!      IF (L_SSref) THEN
!         Ssurf=SSref
!      ELSE
         Ssurf= X(1,2)+Sref
!      ENDIF
       
c Save variables for next timestep
      old = new
      new = 1 - old
      hmixd(new) = hmix
      do k=1,NZP1
         do l=1,NVEL
!          Us(k,l,new)=U(k,l)
          Us(k,l,new)=0.0
         end do  
        do l=1,NSCLR
          Xs(k,l,new)=X(k,l)
        end do
      end do    
      return
      end

***********************************************************************

      SUBROUTINE ocnint(nzi,z,h,d,intri,kmixe,Uo,Xo,Un,Xn)

!c     Integrate the ocn model by backwards Euler(implicit)discretization
!c     On input : Un,Xn are estimated profiles which are used
!c                to estimate diffusivity profiles at new time.
!c              : Updated diffusivities from Un Xn are in common
!c     On output: Un,Xn are new profiles after integration.

!c     Written  19 March 1991 - jan

      IMPLICIT NONE

	include 'parameter.inc'
        include 'constants.com'
        include 'times.com'
        include 'timocn.com'
        include 'vert_pgrid.com'
        include 'proc_swit.com'
        include 'ocn_state.com'
        include 'ocn_paras.com'
        include 'ocn_energy.com'
        include 'flx_sfc.com'
        include 'flx_profs.com'
        include 'kprof_out.com'
        include 'fcorr_in.com'
        include 'location.com'
        include 'dble_diff.com'
        include 'restore.com'
!c Input
      integer intri,            ! index for tri.diag. coeff
     +     nzi
      real  Uo(nzi+1,NVEL), Xo(nzi+1,NSCLR), ! old profiles
     +      z(nzi+1),h(nzi+1),d(0:nzi)      
!c Output
      real  Un(nzi+1,NVEL), Xn(nzi+1,NSCLR)  ! new profiles

!c Common tridiagonal matrix factors (set in "init_ocn")
      real tri(0:NZtmax,0:1,NGRID)    ! dt/dz/dz factors in trid. matrix
      common/ trifac / tri
!c Local
      real cu (NZtmax), ! upper coeff for (k-1) on k line of trid.matrix
     +     cc (NZtmax), ! central ...     (k  ) ..
     +     cl (NZtmax), ! lower .....     (k-1) ..
     +     rhs(NZtmax)  ! right-hand-side terms
      real diff(0:NZtmax),gcap(NZtmax),ntflx(0:NZtmax,NSCLR)
!c More local variables to make implicit none
      integer kmixe,i,npd,imode,n,k
      real ftemp,ghatflux,sturflux
      integer adv_mode
      real adv_mag
      LOGICAL KapFlg
      real Sig0,cp_obs,rho_obs,Sig,cpsw
!c
      real prev_t
      COMMON /save_fcorr_withz/ fcorr_withz , fcorr_saltz
      COMMON /clim_in/ sal_clim, temp_clim
      KapFlg=.False.
!c ********************************************************************
!c U and V solution of tridiagonal matrix
!c                set f = 0 for equatorial application
      ftemp = f  
!c                               set coefficients of tridiagonal matrix

      DO k=0,NZtmax
         diff(k)=difm(k)
      ENDDO
      call tridcof(diff,nzi,intri,cu,cc,cl)
!c                               U right hand side and solution
      rhs(1)= Uo(1,1) + dto*(f*.5*(Uo(1,2)+Un(1,2)) - wU(0,1)/h(1))
      do i=2,nzi-1
        rhs(i)= Uo(i,1) + dto*f*.5*(Uo(i,2)+Un(i,2)) 
      end do
      i=nzi   ! bottom
      rhs(i)= Uo(i,1) + dto*f*.5*(Uo(i,2)+Un(i,2))
     +                + tri(i,1,intri)*difm(i)*Uo(i+1,1)

      call tridmat(cu,cc,cl,rhs,Uo(1,1),nzi,Un(1,1),diff)

!c                                 V rhs and solution
      rhs(1)= Uo(1,2) - dto*(f*.5*(Uo(1,1)+Un(1,1)) + wU(0,2)/h(1))
      do i=2,nzi-1
        rhs(i)= Uo(i,2) - dto*f*.5*(Uo(i,1)+Un(i,1)) 
      end do
      i=nzi
      rhs(i)= Uo(i,2) - dto*f*.5*(Uo(i,1)+Un(i,1))
     +                + tri(i,1,intri)*difm(i)*Uo(i+1,2)

      npd = 1

      call tridmat(cu,cc,cl,rhs,Uo(1,2),nzi,Un(1,2),diff)
      f= ftemp

!c *******************************************************************
!c Scalar solutions of tridiagonal matrix
!c     Temperature (different from other scalars because of ghat-term
!c                  and double diffusion)
!c     ghatflux = wX(0,1) - (1-SWDK(-hmixe,real(time)))
!c    $                     * sflux(3,5,0) / rho(ipt,0) / CP(ipt,0)
!c     ghatflux = wX(0,1) - (1-SWDK(-d(1) ,real(time)))
!c    $                     * sflux(3,5,0) / rho(ipt,0) / CP(ipt,0)
      ghatflux = wX(0,1)
      sturflux = wX(0,1)
      diff(0)=dift(0)
      ntflx(0,1)=wXNT(0,1)
      DO k=1,NZtmax
         diff(k)=dift(k)
         gcap(k)=ghat(k)
         ntflx(k,1)=wXNT(k,1)
      ENDDO
       !write(*,*)'temp diff value',diff(1)
      call tridcof(diff,nzi,intri,cu,cc,cl)
      call tridrhs(npd,h,Xo(1,1),ntflx(0,1),diff,gcap,sturflux,
     >     ghatflux,dto,nzi,intri,rhs)

! restore the temperature back to climatology
      relax_temp=10  !1 W/m2/K 
      do k=1,1  ! temeperature restoring at the surface only
        cp_obs=CPSW(sal_clim(k),temp_clim(k),-zm(i))
        call Sig80(sal_clim(k),temp_clim(k),0.,KapFlg,Sig0,Sig)
        rho_obs=1000.+Sig0
        rhs(k)=rhs(k)+dto*relax_temp*(temp_clim(k)-Xo(k,1))/
     +                         (rho_obs*cp_obs*hm(k))
      end do

! adding flux correction with z
      if (L_FCORR_WITHZ) then
         DO k=1,NZP1
            rhs(k) = rhs(k) + 
     +          dto *fcorr_withz(k)/(rho(k)*cp(k))
         ENDDO
      end if

      call tridmat(cu,cc,cl,rhs,Xo(1,1),nzi,Xn(1,1),diff)

!     Salinity and other scalars
      DO k=0,NZtmax
         diff(k)=difs(k)
      ENDDO   
       !write(*,*)'salt diff value',diff(1)
      call tridcof(diff,nzi,intri,cu,cc,cl)
      do 200 n=2,NSCLR
         DO k=0,NZtmax
            ntflx(k,n)=wXNT(k,n)
         ENDDO
         ghatflux = wX(0,n) 
         sturflux = wX(0,n)
         call tridrhs(npd,h,Xo(1,n),ntflx(0,n),diff,gcap,sturflux,
     >                ghatflux,dto,nzi,intri,rhs)

! modify rhs for advections
!         do imode=1,nmodeadv(2)
!            adv_mode=modeadv(imode,2)
!            adv_mag=advection(imode,2)
!            call rhsmod(2,adv_mode,adv_mag,
!     +           time,dto,dpy,kmixe,d(kmixe),nzi,h,z,rhs)
!         enddo

! modify rhs for salinity restoring
      relax_sal = 1/(50.*360.0*86400)
      do k=1,NZP1
        rhs(k)=rhs(k)+dto*relax_sal*(sal_clim(k)-Xo(k,2))
      end do

      if (L_FCORR_WITHZ) then
         DO k=1,NZP1
            rhs(k) = rhs(k) + dto*fcorr_saltz(k)  !*relax_sal
         ENDDO
        !write(*,*)'fcorr_saltz=',fcorr_saltz(1)
      end if

      call tridmat(cu,cc,cl,rhs,Xo(1,n),nzi,Xn(1,n),diff)
 200  continue
 
      return
      end

************************************************************************

      SUBROUTINE tridcof(diff,nzi,ind,cu,cc,cl)

c     Compute coefficients for tridiagonal matrix (dimension=nzi).
c     Note: cu(1) = 0. and cl(nzi) = 0. are necessary conditions.
c-----     
      IMPLICIT NONE

        include 'parameter.inc'
c Input
      integer nzi,              ! dimension of field
     +     ind                  ! index for tri-coefficients: = kmixo for t-grid,
c                                           =     1 for p-grid.
      real diff(0:nzi) ! diffusivity profile on interfaces
c Output
      real cu(nzi),    ! upper coeff. for (k-1) on k line of trid.matrix
     +     cc(nzi),    ! central ...      (k  ) ..
     +     cl(nzi)     ! lower .....      (k-1) ..
c Common tridiagonal factors (set in "init_ocn<")
      real tri(0:NZtmax,0:1,NGRID)    ! dt/dz/dz factors in trid. matrix
      common/ trifac / tri
c more local variables (to make implicit none)
      integer i
c
c In the surface layer
      cu(1) = 0.
      cc(1) = 1. + tri(1,1,ind)*diff(1)   ! 1.+ dto/h(1)/dzb(1)*diff(1)
      cl(1) =    - tri(1,1,ind)*diff(1)   !   - dto/h(1)/dzb(1)*diff(1)
c Inside the domain
      do i=2,nzi
      cu(i) =    - tri(i,0,ind)*diff(i-1)
      cc(i) = 1. + tri(i,1,ind)*diff(i)   + tri(i,0,ind)*diff(i-1)
      cl(i) =    - tri(i,1,ind)*diff(i)
      end do   
c In the bottom layer
      cl(nzi)= 0.
      return
      end

***********************************************************************
      SUBROUTINE tridrhs(npd,h,yo,ntflux,diff,ghat,sturflux,ghatflux,
     +                   dto,nzi,ind,rhs)

c     Compute right hand side of tridiagonal matrix for scalar fields:
c     =  yo (old field) 
c      + flux-divergence of ghat
c      + flux-divergence of non-turbulant fluxes
c     Note: surface layer needs +dto/h(1) * surfaceflux
c           bottom  ..... ..... +dto/h(nzi)*diff(nzi)/dzb(nzi)*yo(nzi+1)

      IMPLICIT NONE

        include 'parameter.inc'
        include 'times.com'
c Input
      real dto           ! timestep interval (seconds)
      integer nzi,       ! dimension of field
     +        ind        ! index for tri-coefficients:=kmixo for t-grid,
c                              =    1 for p-grid
      real h(nzi+1),     ! layer thickness
     +     yo(nzi+1),    ! old profile
     +     ntflux(0:nzi),! non-turbulent flux = wXNT(0:nzi,1:2)
     +     diff(0:nzi),  ! diffusivity profile on interfaces
     +     ghat(nzi),    ! ghat turbulent flux   
     +     sturflux,     ! surface turbulent (kinematic) flux = wX(0,n)
     +     ghatflux      ! surface flux for ghat: includes solar flux      
      integer npd        ! included in list by sjw for implicit none
c Output
      real rhs(nzi)      ! right hand side
c Common tridiagonal factors (set in "init_ocn")
      real tri(0:NZtmax,0:1,NGRID)    ! dt/dz/dz factors in trid. matrix
      common/ trifac / tri

c more local variables to make implicit none
      integer i
      real divflx
      divflx =  1.0 / float(npd)

!	if(ntime.gt.23683*4.and.ntime.lt.23700*4) then
!	write(*,*)'ntime=',ntime
!        write(*,*)yo(1),dto,h(1),ghatflux,diff(1),ghat(1)
!        write(*,*)sturflux,divflx,ntflux(1),ntflux( 0 )
!	end if 
c In the surface layer (dto/h(1)=tri(0,1,ind)
      rhs(1)= yo(1) + dto/h(1) *
     +        (   ghatflux*diff(1)*ghat(1)
     +       - sturflux*divflx
     +       + ntflux(1) - ntflux( 0 ) )
c  Inside the domain to npd
        if(npd.ge.2) then
      do  i=2,npd
      rhs(i)= yo(i) + dto/h(i) *
     + (                        ghatflux *  diff(i)  *ghat(i)
     +  -                       ghatflux *  diff(i-1)*ghat(i-1)
     +        - sturflux * divflx
     +        + ntflux(i) - ntflux(i-1) )
      end do
        endif

c Inside the rest of the domain
      do i=npd+1,nzi-1
      rhs(i)= yo(i) + dto/h(i) *
     +              ( ghatflux*(diff(i)*ghat(i) - diff(i-1)*ghat(i-1))
     +               +ntflux(i) - ntflux(i-1) )
      end do

c In the bottom layer     
      if(nzi.gt.1) then   ! not for slab ocean
      i=nzi
      rhs(i)= yo(i) + dto/h(i) *
     +              ( ghatflux*(diff(i)*ghat(i) - diff(i-1)*ghat(i-1))
     +               +ntflux(i) - ntflux(i-1) )
     +      + yo(i+1)*tri(i,1,ind)*diff(i)
      endif
	!write(*,*)rhs
      return
      end

c**********************************************************************
      subroutine rhsmod(jsclr,mode,A,time,
     +                  dto,dpy,km,dm,nzi,h,z,rhs)
 
c     Modify rhs to correct scalar, jsclr, 
c     for advection according to mode
c mode = 1 : Steady upper layer horizontal advection
c        2 : Steady mixed layer horizontal advection to km-1 
c        3 : Steady horizontal advection throughout the entire column
c        4 : Steady vertical advection (= deep horizontal) below 100m 
c            to bottom 
c            (Change: start below 100m, instead of at layer 16, and
c            do not advect into bottom layer, 7-1-93)
c        5 : Steady bottom diffusion
c        6 : Seasonal mixed layer horizontal advection to dm
c        7 : Seasonal thermocline horizontal advection to 1.5 dm
      IMPLICIT NONE

        include 'parameter.inc'
        include 'ocn_paras.com'
c Input
      integer nzi,       ! vertical dimension of field
     +         km,       ! index of gridpoint just below h
     +       mode,       ! type of advection
     +      jsclr        ! scalar
      real h(nzi+1),     ! layer thickness
     +     z(nzi+1),     ! z grid levels (added as input on 7-1-93)
     +     rhs(nzi)     ! right hand side from tridrhs
      DOUBLE PRECISION time        ! time in days from jan 1 of any year.
      real dto,          ! ocean time step
     +     dpy,          ! days per year (added as input on 7-1-93)
     +     dm,           ! depth d(km+.5)
     +     A             ! advection of heat(W/m2) or Salt(PSU m/s)      

c Output
c     real rhs(nzi)      ! modified right hand side

      real am,fact,delta,depth,dmax
      integer n,n1,n2,nzend

c Internal
      real f(12)     ! monthly partion of annual advection
      real xsA(21)   ! yearly excess of heat

c     data f/.1,.1,6*0.0,.1,.3,.4,.2/
      data f/.05,.05,5*0.0,.05,.15,.20,.30,.20/
      data xsA/48.26,21.73,29.02,56.59,19.94,15.96,18.28,
     &         40.52,37.06,29.83,29.47,15.77, 1.47,14.55,
     &          4.22,28.19,39.54,19.58,20.27,11.19,21.72/
      if(mode.le.0) return 
c       Am = -12. * f(month) * (xsA(iyr) - 0.0 )       ! Annual
c        Am =  12. * f(month) * A                       ! Seasonal
       Am = A                                          ! Steady

      if(mode.eq.1) then
c                          correct upper layer advection
        if(jsclr.eq.1) fact = dto * Am / (rho(1)*cp(1))
        if(jsclr.eq.2) fact = dto * Am * 0.033
        rhs(1) = rhs(1) + fact / h(1)

      else if(mode.eq.2) then
c                          correct mixed layer advection
        delta = 0.0
        do n=1,km-1
          delta = delta + h(n)
        end do
        do 215 n=1,km-1
        if(jsclr.eq.1) fact = dto * Am / (rho(n)*cp(n)) 
        if(jsclr.eq.2) fact = dto * Am * 0.033
        rhs(n) = rhs(n) + fact  / delta        
 215    continue

      else if (mode.eq.3) then
c                               throughout whole water column
        delta = 0.0
        do n=1,nzi
           delta = delta + h(n)
        end do
        do 315 n=1,nzi
        if(jsclr.eq.1) fact = dto * Am / (rho(n)*cp(n))
        if(jsclr.eq.2) fact = dto * Am * 0.033
         rhs(n) = rhs(n) + fact / delta 
  315   continue

    
      else if (mode.eq.4) then
c                          vertical advection = deep horizontal
        nzend=nzi-1                          ! nzend=nzi (change:7-1-93)
        n1=0                                 ! n1=16     (change:7-1-93)
  401   n1=n1+1
        if(z(n1).ge.-100.) goto 401
        delta = 0.0
        do n=n1,nzend
           delta = delta + h(n)
        end do
        do 415 n=n1,nzend
        if(jsclr.eq.1) fact = dto * Am / (rho(n)*cp(n))
        if(jsclr.eq.2) fact = dto * Am * 0.033
         rhs(n) = rhs(n) + fact / delta 
  415   continue

      else if(mode.eq.5) then
c                          correct bottom layer diffusion
        if(jsclr.eq.1) fact = dto * Am / (rho(nzi)*cp(nzi))
        if(jsclr.eq.2) fact = dto * Am * 0.033
        rhs(nzi) = rhs(nzi) + fact / h(nzi)

      else 

c              seasonal mixed layer or thermocline advection  
c               find ocean year
c       iyr = 1 + idint((time-75.0)/dpy)
c       day = time - dpy * (idint(time/dpy))  ! 365.25))
c       month = 1 + int(12. * day / dpy)    ! 367.)
cdiag
c       if(month.gt.12) then
c          write(nuerr,*) 'STOP rhsmod (ocn.f):'
c          write(nuerr,*) '     rounding error, month gt 12 =',month
c          stop 97
c       endif   
cdiag
c       Am = -12. * f(month) * (xsA(iyr) - 0.0 )       ! Annual
c       Am =  12. * f(month) * A                       ! Seasonal
c       Am = A                                          ! Steady

        if(mode.eq.6) then
c                            mixed layer to dm
          n1 = 1
          depth = h(1)
          dmax  = dm -  0.5 * (h(km) + h(km-1))
          delta = 0.0
          do 605 n =n1,nzi
          n2    = n
          delta = delta + h(n)
          depth = depth + h(n+1)
          if(depth.ge.dmax) go to 606
 605      continue
 606      continue

        else if (mode.eq.7) then
c                                 thermocline to 100m   
          n1 = km - 1
          depth = dm - 0.5 * h(km) 
          dmax = 100.        
          delta = 0.0
          do 705 n=n1,nzi
          n2 = n
          delta = delta + h(n)
          depth = depth + h(n+1)
          if(depth.ge.dmax) go to 706
  705     continue
  706     continue

        else
          write(*,*) 'STOP in rhsmod (ocn.f):'
          write(*,*) '      mode out of range, mode=',mode
        endif

c                          Finish both 6 and 7 here
        do 615 n=n1,n2  
          if(jsclr.eq.1) fact = dto * Am / (rho(n)*cp(n)) 
          if(jsclr.eq.2) fact = dto * Am * 0.033
          rhs(n) = rhs(n) + fact  / delta
 615    continue
      endif
      return
      end

***********************************************************************
      SUBROUTINE tridmat(cu,cc,cl,rhs,yo,nzi,yn,diff)
c
c     Solve tridiagonal matrix for new vector yn, given right hand side
c     vector rhs. Note: yn(nzi+1) = yo(nzi+1).
c-----     
      IMPLICIT NONE

        include 'parameter.inc'
        include 'times.com'
c Input
      integer nzi               ! dimension of matrix
      real cu (nzi),            ! upper coeff. for (k-1) on k line of tridmatrix
     +     cc (nzi),            ! central ...      (k  ) ..
     +     cl (nzi),            ! lower .....      (k-1) ..
     +     rhs(nzi),            ! right hand side
     +     yo(nzi+1),           ! old field
     +     diff(0:nzi)
c Output
      real yn(nzi+1)    ! new field
c Local 
      real gam(NZtmax), ! temporary array for tridiagonal solver
     +     bet          ! ...
c more local for implicit none
      integer i
c Solve tridiagonal matrix.
      !bet   = cc(1)
      bet   = max(cc(1),1.E-12) ! BP added
c     yn(1) = (rhs(1) + tri(0,1,ind)*surflux) / bet    ! surface
      yn(1) =  rhs(1) / bet    ! surface
      do 21 i=2,nzi
         gam(i)= cl(i-1)/bet
         bet   = cc(i) - cu(i)*gam(i)
         if(bet.eq.0.) then
          write(*,*)'* algorithm for solving tridiag matrix fails'
          !write(*,*)'* bet=',bet
          !write(*,*)'*i-1=',i-1,' cc=',cc(i-1),'cl=',cl(i-1)
          !write(*,*)'*i=',i,' cc=',cc(i),' cu=',cu(i),' gam=',gam(i)
          bet=1.E-12
         endif
c        to avoid "Underflow" at single precision on the sun
         yn(i) =      (rhs(i)  - cu(i)  *yn(i-1)  )/bet
 21   continue
      do i=nzi-1,1,-1
         yn(i)  = yn(i) - gam(i+1)*yn(i+1)
      end do
      yn(nzi+1) = yo(nzi+1)

!	if(ntime.gt.23683*4.and.ntime.lt.23700*4) then
!	write(*,*)yn(1),gam(2),yn(2),cc(1)
!	end if
      return
      end

*********************************************************************

      SUBROUTINE init_ocn(U,X)
c
c     Initialize ocean model:
c     Set coefficients for tridiagonal matrix solver.
c     Compute hmix and diffusivity profiles for initial profile.
c     Prepare for first time step.
c     
      IMPLICIT NONE

	include 'parameter.inc'
        include 'constants.com'
        include 'times.com'
        include 'timocn.com'
        include 'vert_pgrid.com'
        include 'flx_profs.com'
        include 'kprof_out.com'
        include 'ocn_paras.com'
        include 'ocn_state.com'

c     Input
      real  U(NZP1,NVEL),X(NZP1,NSCLR)
      real U0(NZP1,NVEL),X0(NZP1,NSCLR)
      
c     Common Blocks
      real tri(0:NZtmax,0:1,NGRID)! dt/dz/dz factors in trid. matrix
      common/ trifac / tri

      real hmixd(0:1),     ! storage arrays for extrapolations
     +     Us(NZP1,NVEL ,0:1),    ! ..      ..     ..  ..
     +     Xs(NZP1,NSCLR,0:1)     ! ..      ..     ..  ..
      integer old,new ! extrapolation index for Us,Xs,hmixd
      common/ saveUXh /
     +     old,new,Us,Xs,hmixd

c Local variables
      real dzb(NZ)              ! diff. between grid-levels below z(j)
c     more local for implicit none
      integer k,kmix0,n,l,i
      real hmix0,deltaz
c     
c     Compute factors for coefficients of tridiagonal matrix elements.
c     tri(0     ,1,.........) : dt/h(1) factor for rhs flux
c     tri(k=1:NZ,0,.........) : dt/h(k)/ {dzb(k-1)=z(k-1)-z(k)=dzabove}
c     tri(k=1:NZ,1,.........) : dt/h(k)/ {dzb(k  )=z(k)-z(k+1)=dzbelow}
c
      do k=1,NZ
         dzb(k)     = zm(k) - zm(k+1)
      end do

      tri(0,1,1)    = dto/hm(1)   !dto=1800 and hm=depth
      tri(1,1,1)    = dto/hm(1)/dzb(1)
      do k=2,NZ
         tri(k,1,1) = dto/hm(k)/dzb(k)
         tri(k,0,1) = dto/hm(k)/dzb(k-1)
      end do
      
c     Determine hmix for initial profile:    
      DO k=1,NZP1
         DO i=1,nvel
            U0(k,i)=U(k,i)
         ENDDO
         DO i=1,nsclr
            X0(k,i)=X(k,i)
         ENDDO
      ENDDO
      call vmix(U0,X0,hmix0,kmix0)
      hmix = hmix0
      kmix = kmix0
      Tref = X(1,1)
c Evaluate initial fluxes (to write to output data file)
      do k=1,NZ
        deltaz = 0.5*(hm(k)+hm(k+1))
        do n=1,NSCLR
           wX(k,n)=-difs(k)*((X(k,n)-X(k+1,n))/deltaz-ghat(k)*wX(0,n))
        end do
        wX(k,nsp1)= grav * (talpha(k)*wX(k,1) - sbeta(k) * wX(k,2))
        do n=1,NVEL
           wU(k,n)= -difm(k)* (U(k,n)-U(k+1,n))/deltaz
        end do
      end do
                     
c     indices for extrapolation
       old = 0
       new = 1
c     initialize array for extrapolating hmixd,Us,Xs
       hmixd(0) = hmix
       hmixd(1) = hmix
       do k=1,NZP1
          do l=1,NVEL
             Us(k,l,0)=U(k,l)
             Us(k,l,1)=U(k,l)
          end do
          do l=1,NSCLR
             Xs(k,l,0)=X(k,l)
             Xs(k,l,1)=X(k,l)
          end do
       end do     
                  
      return
      end

*********************************************************************
      subroutine swfrac(fact, z, jwtype, swdk )
c     compute fraction of solar short-wave flux penetrating to specified
c     depth (times fact) due to exponential decay in  Jerlov water type
c     reference : two band solar absorption model of simpson and 
c     paulson (1977)
      
      IMPLICIT NONE

      integer nwtype
      parameter(nwtype=5) ! max number of different water types 
c
c  model
      integer imt         ! number of horizontal grid points
c  input
      real fact           ! scale  factor to apply to depth array
      real z         ! vertical height ( <0.) for desired sw 
c                           fraction                                 (m)
      integer jwtype ! index for jerlov water type
c  output
      real swdk      !  short wave (radiation) fractional decay
c  local
      real  rfac(nwtype),a1(nwtype),a2(nwtype)
      real rmin,r1,r2
      integer i
      save  rfac,a1,a2,rmin
c
c     jerlov water type :  I       IA      IB      II      III
c                jwtype    1       2       3       4       5
c
      data rfac         /  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /
      data a1           /  0.35 ,  0.6  ,  1.0  ,  1.5  ,  1.4  /
      data a2           / 23.0  , 20.0  , 17.0  , 14.0  ,  7.9  /
      data rmin         / -80. /
c
!        jwtype=3
         r1      = MAX(z*fact/a1(jwtype), rmin)
	 r2      = MAX(z*fact/a2(jwtype), rmin)
         swdk =      rfac(jwtype)  * exp(r1)
     $            + (1.-rfac(jwtype)) * exp(r2)

      return
      end

c *********************************************************


      SUBROUTINE zint(nz,zed,Gout,mz,zm,Gint)
      
      dimension zed(0:nz+1), Gout(0:nz+1)
      dimension  zm(mz+1), Gint(0:mz+1)

      Gint(0) = Gout(0)
      k = 1
      do 35 j=1,mz+1
  45  continue
          if(zed(k).gt.zm(j)) then
            k = k+1
            go to 45
          else
            Gint(j) = Gout(k-1) + (Gout(k)-Gout(k-1)) *
     >                          (zm(j)-zed(k-1))/(zed(k)-zed(k-1))
          endif

 35   continue

      return
      end

************************************************************************
      subroutine swfrac_opt(fact, z, k, jwtype, swdk )
c     compute fraction of solar short-wave flux penetrating to specified
c     depth (times fact) due to exponential decay in  Jerlov water type
c     reference : two band solar absorption model of simpson and 
c     paulson (1977)
      
      IMPLICIT NONE

      integer nwtype
      parameter(nwtype=5) ! max number of different water types 
c
c  model
        include 'parameter.inc'
        include 'times.com'
        include 'vert_pgrid.com'
      
c  input
      real fact           ! scale  factor to apply to depth array
      integer k      ! index of vertical grid point
      real z         ! vertical height ( <0.) for desired sw 
c                           fraction                                 (m)
      integer jwtype ! index for jerlov water type
c  output
      real swdk      !  short wave (radiation) fractional decay
c  local
      real swfrac_save(NZP1)
      real  rfac(nwtype),a1(nwtype),a2(nwtype)
      real rmin,r1,r2
      integer i,l
      save  rfac,a1,a2,rmin
      common /save_swfrac/swfrac_save
c
c     jerlov water type :  I       IA      IB      II      III
c                jwtype    1       2       3       4       5
c
      data rfac         /  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /
      data a1           /  0.35 ,  0.6  ,  1.0  ,  1.5  ,  1.4  /
      data a2           / 23.0  , 20.0  , 17.0  , 14.0  ,  7.9  /
      data rmin         / -80. /
c
!      jwtype=3
        IF ((ntime .eq. 1) .and. (k .eq. 2)) THEN
         DO l=1,NZP1
               r1      = MAX(zm(l)*fact/a1(jwtype), rmin)
               r2      = MAX(zm(l)*fact/a2(jwtype), rmin)
               swfrac_save(l) =      rfac(jwtype)  * exp(r1)
     $              + (1.-rfac(jwtype)) * exp(r2)
               
         ENDDO
      ENDIF
      
         swdk=swfrac_save(k)

      return
      end

!***********************************************************
      SUBROUTINE init_env (dscale,lat)
c     ===================
c     My init_env routine modified from ...
c     Modified  3 March 1991   -  WGL in KPP/input.f

c     the physical environment for the model:
c     vertical grid and geographic location
c
      IMPLICIT NONE

      include 'parameter.inc'
      include 'constants.com'
      include 'vert_pgrid.com'
      include 'location.com'
      include 'ocn_state.com'
c     Inputs
      REAL dscale,lat
c     Local Variables
      REAL sumh,hsum,dfac,sk
      INTEGER i
c define vertical grid fields
	DMAX=1000.0
	dscale=4.0
         sumh = 0.0
         dfac = 1.0 - exp(-dscale)
         do 5 i = 1,(NZ)
            sk = - (float(i)-0.5)/float(NZ)
            hm(i) = DMAX*dfac/float(NZ)/dscale / ( 1.0 + sk*dfac )
            sumh = sumh + hm(i)
    5    continue
c     
c layer thickness h, layer grids zgrid, interface depths d
c
      hsum = 0.0
      dm(0) = 0.0
      do 15 i=1,NZ
            hm(i) = hm(i) * DMAX / sumh
         zm(i) =  0.0 - (hsum + 0.5 * hm(i) )
         hsum = hsum + hm(i)
         dm(i) = hsum
 15   continue
      hm(nzp1) = 1.e-10
      zm(nzp1) = -DMAX

c
c set fraction of land, and initial open water fraction
         focn  = 1.
c     
c compute geographic location
c      
!        rlat=dlat(ipt)*twopi/360.
!        rlon(ipt)=dlon(ipt)*twopi/360.
c
c Coriolis Parameter ! check the necessity of this
c
c
        twopi=8*atan(1.)          ! 2pi    
	f = 2. * (twopi/86164.) * sin(lat*twopi/360.)
      return
      end

