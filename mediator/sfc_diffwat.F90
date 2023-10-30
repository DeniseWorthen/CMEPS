
!> \defgroup GFS_diff_main GFS Surface Layer Module
!> This module calculates surface roughness length.
!> @{
!! This subroutine includes the surface roughness length formulation
!! based on the surface sublayer scheme in
!! Zeng and Dickinson (1998) \cite zeng_and_dickinson_1998.
!> \section arg_table_sfc_diff_run Argument Table
!! \htmlinclude sfc_diff_run.html
!!
!>  \section general_diff GFS Surface Layer Scheme General Algorithm
!! - Calculate the thermal roughness length formulation over the ocean (see eq. (25) and (26)
!!  in Zeng et al. (1998) \cite zeng_et_al_1998).
!! - Calculate Zeng's momentum roughness length formulation over land and sea ice.
!! - Calculate the new vegetation-dependent formulation of thermal roughness length
!! (Zheng et al.(2012) \cite zheng_et_al_2012).
!! Zheng et al. (2012) \cite zheng_et_al_2012 proposed a new formulation on
!! \f$ln(Z_{0m}^,/Z_{0t})\f$ as follows:
!! \f[
!!  ln(Z_{0m}^,/Z_{0t})=(1-GVF)^2C_{zil}k(u*Z_{0g}/\nu)^{0.5}
!! \f]
!! where \f$Z_{0m}^,\f$ is the effective momentum roughness length
!! computed in the following equation for each grid, \f$Z_{0t}\f$
!! is the roughness lenghth for heat, \f$C_{zil}\f$ is a coefficient
!! (taken as 0.8), k is the Von Karman constant (0.4),
!! \f$\nu=1.5\times10^{-5}m^{2}s^{-1}\f$ is the molecular viscosity,
!! \f$u*\f$ is the friction velocity, and \f$Z_{0g}\f$ is the bare
!! soil roughness length for momentum (taken as 0.01).
!! \n In order to consider the convergence of \f$Z_{0m}\f$ between
!! fully vegetated and bare soil, the effective \f$Z_{0m}^,\f$ is
!! computed:
!! \f[
!!  \ln(Z_{0m}^,)=(1-GVF)^{2}\ln(Z_{0g})+\left[1-(1-GVF)^{2}\right]\ln(Z_{0m})
!!\f]
!! - Calculate the exchange coefficients:\f$cm\f$, \f$ch\f$, and \f$stress\f$ as inputs of other \a sfc schemes.
!!
      subroutine sfc_diff_run (im,rvrdm1,eps,epsm1,grav,                &  !intent(in)
     &                    ps,t1,q1,z1,garea,wind,                       &  !intent(in)
     &                    prsl1,prslki,prsik1,prslk1,                   &  !intent(in)
     &                    sigmaf,vegtype,shdmax,ivegsrc,                &  !intent(in)
     &                    z0pert,ztpert,                                &  ! mg, sfc-perts !intent(in)
     &                    flag_iter,redrag,                             &  !intent(in)
     &                    u10m,v10m,sfc_z0_type,                        &  !hafs,z0 type !intent(in)
     &                    wet,dry,icy,                                  &  !intent(in)
     &                    thsfc_loc,                                    &  !intent(in)
     &                    tskin_wat, tskin_lnd, tskin_ice,              &  !intent(in)
     &                    tsurf_wat, tsurf_lnd, tsurf_ice,              &  !intent(in)
     &                     z0rl_wat,  z0rl_lnd,  z0rl_ice,              &  !intent(inout)
     &                     z0rl_wav,                                    &  !intent(inout)
     &                    ustar_wat, ustar_lnd, ustar_ice,              &  !intent(inout)
     &                       cm_wat,    cm_lnd,    cm_ice,              &  !intent(inout)
     &                       ch_wat,    ch_lnd,    ch_ice,              &  !intent(inout)
     &                       rb_wat,    rb_lnd,    rb_ice,              &  !intent(inout)
     &                   stress_wat,stress_lnd,stress_ice,              &  !intent(inout)
     &                       fm_wat,    fm_lnd,    fm_ice,              &  !intent(inout)
     &                       fh_wat,    fh_lnd,    fh_ice,              &  !intent(inout)
     &                     fm10_wat,  fm10_lnd,  fm10_ice,              &  !intent(inout)
     &                      fh2_wat,   fh2_lnd,   fh2_ice,              &  !intent(inout)
     &                    ztmax_wat, ztmax_lnd, ztmax_ice,              &  !intent(inout)
     &                    zvfun,                                        &  !intent(out)
     &                    errmsg, errflg)                                  !intent(out)
!
      implicit none
!
      integer, parameter  :: kp = kind_phys
      integer, intent(in) :: im, ivegsrc
      integer, intent(in) :: sfc_z0_type ! option for calculating surface roughness length over ocean

      integer, dimension(:), intent(in) :: vegtype

      logical, intent(in) :: redrag ! reduced drag coeff. flag for high wind over sea (j.han)
      logical, dimension(:), intent(in) :: flag_iter, dry, icy
      logical, dimension(:), intent(inout) :: wet

      logical, intent(in) :: thsfc_loc ! Flag for reference pressure in theta calculation

      real(kind=kind_phys), dimension(:), intent(in)    :: u10m,v10m
      real(kind=kind_phys), intent(in) :: rvrdm1, eps, epsm1, grav
      real(kind=kind_phys), dimension(:), intent(in)    ::              &
     &                    ps,t1,q1,z1,garea,prsl1,prslki,prsik1,prslk1, &
     &                    wind,sigmaf,shdmax,                           &
     &                    z0pert,ztpert ! mg, sfc-perts
      real(kind=kind_phys), dimension(:), intent(in)    ::              &
     &                    tskin_wat, tskin_lnd, tskin_ice,              &
     &                    tsurf_wat, tsurf_lnd, tsurf_ice

      real(kind=kind_phys), dimension(:), intent(in)    :: z0rl_wav
      real(kind=kind_phys), dimension(:), intent(inout) ::              &
     &                     z0rl_wat,  z0rl_lnd,  z0rl_ice,              &
     &                    ustar_wat, ustar_lnd, ustar_ice,              &
     &                       cm_wat,    cm_lnd,    cm_ice,              &
     &                       ch_wat,    ch_lnd,    ch_ice,              &
     &                       rb_wat,    rb_lnd,    rb_ice,              &
     &                   stress_wat,stress_lnd,stress_ice,              &
     &                       fm_wat,    fm_lnd,    fm_ice,              &
     &                       fh_wat,    fh_lnd,    fh_ice,              &
     &                     fm10_wat,  fm10_lnd,  fm10_ice,              &
     &                      fh2_wat,   fh2_lnd,   fh2_ice,              &
     &                    ztmax_wat, ztmax_lnd, ztmax_ice
      real(kind=kind_phys), dimension(:), intent(out)    :: zvfun
!
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
!
!     locals
!
      integer   i
!
      real(kind=kind_phys) :: rat, tv1, thv1, restar, wind10m,
     &                        czilc, tem1, tem2, virtfac
!

      real(kind=kind_phys) :: tvs, z0, z0max, ztmax, gdx
!
      real(kind=kind_phys), parameter :: z0lo=0.1, z0up=1.0
!
      real(kind=kind_phys), parameter ::
     &        one=1.0_kp, zero=0.0_kp, half=0.5_kp, qmin=1.0e-8_kp
     &,       charnock=.018_kp, z0s_max=.317e-2_kp                      &! a limiting value at high winds over sea
     &,       zmin=1.0e-6_kp                                            &
     &,       vis=1.4e-5_kp, rnu=1.51e-5_kp, visi=one/vis               &
     &,       log01=log(0.01_kp), log05=log(0.05_kp), log07=log(0.07_kp)

!     parameter (charnock=.014,ca=.4)!c ca is the von karman constant
!     parameter (alpha=5.,a0=-3.975,a1=12.32,b1=-7.755,b2=6.041)
!     parameter (a0p=-7.941,a1p=24.75,b1p=-8.705,b2p=7.899,vis=1.4e-5)

!     real(kind=kind_phys) aa1,bb1,bb2,cc,cc1,cc2,arnu
!     parameter (aa1=-1.076,bb1=.7045,cc1=-.05808)
!     parameter (bb2=-.1954,cc2=.009999)
!     parameter (arnu=.135*rnu)
!
!    z0s_max=.196e-2 for u10_crit=25 m/s
!    z0s_max=.317e-2 for u10_crit=30 m/s
!    z0s_max=.479e-2 for u10_crit=35 m/s
!
! mbek -- toga-coare flux algorithm
!     parameter (rnu=1.51e-5,arnu=0.11*rnu)

! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!  initialize variables. all units are supposedly m.k.s. unless specified
!  ps is in pascals, wind is wind speed,
!  surface roughness length is converted to m from cm
!
!       write(0,*)'in sfc_diff, sfc_z0_type=',sfc_z0_type

      do i=1,im
        if(flag_iter(i)) then

          ! Need to initialize ztmax arrays
          ztmax_lnd(i) = 1. ! log(1) = 0
          ztmax_ice(i) = 1. ! log(1) = 0
          ztmax_wat(i) = 1. ! log(1) = 0

          virtfac = one + rvrdm1 * max(q1(i),qmin)

          tv1 = t1(i) * virtfac ! Virtual temperature in middle of lowest layer
          if(thsfc_loc) then ! Use local potential temperature
            thv1  = t1(i) * prslki(i) * virtfac
          else ! Use potential temperature reference to 1000 hPa
            thv1    = t1(i) / prslk1(i) * virtfac
          endif

          zvfun(i) = zero
          gdx = sqrt(garea(i))


! BWG: Everything from here to end of subroutine was after
!      the stuff now put into "stability"

          if (wet(i)) then ! Some open ocean

            zvfun(i) = zero

            if(thsfc_loc) then ! Use local potential temperature
              tvs        = half * (tsurf_wat(i)+tskin_wat(i)) * virtfac
            else
              tvs        = half * (tsurf_wat(i)+tskin_wat(i))/prsik1(i)
     &                          * virtfac
            endif

            z0           = 0.01_kp * z0rl_wat(i)
            z0max        = max(zmin, min(z0,z1(i)))
!           ustar_wat(i) = sqrt(grav * z0 / charnock)
            wind10m      = sqrt(u10m(i)*u10m(i)+v10m(i)*v10m(i))

!**  test xubin's new z0

!           ztmax  = z0max

            restar = max(ustar_wat(i)*z0max*visi, 0.000001_kp)

!           restar = log(restar)
!           restar = min(restar,5.)
!           restar = max(restar,-5.)
!           rat    = aa1 + (bb1 + cc1*restar) * restar
!           rat    = rat    / (1. + (bb2 + cc2*restar) * restar))
!  rat taken from zeng, zhao and dickinson 1997

            rat   = min(7.0_kp, 2.67_kp * sqrt(sqrt(restar)) - 2.57_kp)
            ztmax_wat(i) = max(z0max * exp(-rat), zmin)
!
            if (sfc_z0_type == 6) then
              call znot_t_v6(wind10m, ztmax_wat(i))   ! 10-m wind,m/s, ztmax(m)
            else if (sfc_z0_type == 7) then
              call znot_t_v7(wind10m, ztmax_wat(i))   ! 10-m wind,m/s, ztmax(m)
            else if (sfc_z0_type > 0) then
              write(0,*)'no option for sfc_z0_type=',sfc_z0_type
              errflg = 1
              errmsg = 'ERROR(sfc_diff_run): no option for sfc_z0_type'
              return
            endif
!
            call stability
!  ---  inputs:
     &       (z1(i), zvfun(i), gdx, tv1, thv1, wind(i),
     &        z0max, ztmax_wat(i), tvs, grav, thsfc_loc,
!  ---  outputs:
     &        rb_wat(i), fm_wat(i), fh_wat(i), fm10_wat(i), fh2_wat(i),
     &        cm_wat(i), ch_wat(i), stress_wat(i), ustar_wat(i))
!
!  update z0 over ocean
!
            if (sfc_z0_type >= 0) then
              if (sfc_z0_type == 0) then
!               z0 = (charnock / grav) * ustar_wat(i) * ustar_wat(i)
                tem1 = 0.11 * vis / ustar_wat(i)
                z0 = tem1 + (charnock/grav)*ustar_wat(i)*ustar_wat(i)


! mbek -- toga-coare flux algorithm
!               z0 = (charnock / grav) * ustar(i)*ustar(i) +  arnu/ustar(i)
!  new implementation of z0
!               cc = ustar(i) * z0 / rnu
!               pp = cc / (1. + cc)
!               ff = grav * arnu / (charnock * ustar(i) ** 3)
!               z0 = arnu / (ustar(i) * ff ** pp)

                if (redrag) then
                  z0rl_wat(i) = 100.0_kp * max(min(z0, z0s_max),        &
     &                                                 1.0e-7_kp)
                else
                  z0rl_wat(i) = 100.0_kp * max(min(z0,0.1_kp), 1.e-7_kp)
                endif

              elseif (sfc_z0_type == 6) then   ! wang
                 call znot_m_v6(wind10m, z0)   ! wind, m/s, z0, m
                 z0rl_wat(i) = 100.0_kp * z0   ! cm
              elseif (sfc_z0_type == 7) then   ! wang
                 call znot_m_v7(wind10m, z0)   ! wind, m/s, z0, m
                 z0rl_wat(i) = 100.0_kp * z0   ! cm
              else
                 z0rl_wat(i) = 1.0e-4_kp
              endif

            elseif (z0rl_wav(i) <= 1.0e-7_kp .or.                       &
     &              z0rl_wav(i) > 1.0_kp) then
!             z0 = (charnock / grav) * ustar_wat(i) * ustar_wat(i)
              tem1 = 0.11 * vis / ustar_wat(i)
              z0 = tem1 + (charnock/grav)*ustar_wat(i)*ustar_wat(i)

              if (redrag) then
                z0rl_wat(i) = 100.0_kp * max(min(z0, z0s_max),1.0e-7_kp)
              else
                z0rl_wat(i) = 100.0_kp * max(min(z0,0.1_kp), 1.0e-7_kp)
              endif
            endif

          endif              ! end of if(open ocean)
!
        endif                ! end of if(flagiter) loop
      enddo

      return
      end subroutine sfc_diff_run
