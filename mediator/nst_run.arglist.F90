  call sfc_nst_pre_run(im=GFS_Control%blksz(cdata%blk_no),    &
       wet=GFS_Interstitial(cdata%thrd_no)%wet,               &
       tgice=con_tice,                                        &
       !************************
       tsfco=GFS_Data(cdata%blk_no)%Sfcprop%tsfco,            &
       tsurf_wat=GFS_Interstitial(cdata%thrd_no)%tsurf_water, &
       tseal=GFS_Interstitial(cdata%thrd_no)%tseal,           &
       !************************
       xt=GFS_Data(cdata%blk_no)%Sfcprop%xt,                  &
       xz=GFS_Data(cdata%blk_no)%Sfcprop%xz,                  &
       dt_cool=GFS_Data(cdata%blk_no)%Sfcprop%dt_cool,        &
       z_c=GFS_Data(cdata%blk_no)%Sfcprop%z_c,                &
       tref=GFS_Data(cdata%blk_no)%Sfcprop%tref,              &
       cplflx=GFS_Control%cplflx,                             &
       oceanfrac=GFS_Data(cdata%blk_no)%Sfcprop%oceanfrac,    &
       nthreads=GFS_Control%nthreads,                         &
       errmsg=cdata%errmsg,                                   &
       errflg=cdata%errflg)

  ! pre
  tseal = tsfco  ! what is the point? tseal is over-written by tref + ....
  -> get dtzm_2d
  tref = tsfco + ...dtzm
  tseal = tref + ....
  tsurf_wat = tseal
  ! tref has been calculated, tseal=tsurf_wat=tref +....
  ! tsurf_wat passed out as %tsurf_watER

  ! sfc_nst_run gets tref, tseal (variable tskin) and tsurf_watER (variable tsurf)
  call sfc_nst_run(im =GFS_Control%blksz(cdata%blk_no),                           &
       hvap           =con_hvap,                                                  &
       cp             =con_cp,                                                    &
       hfus           =con_hfus,                                                  &
       jcal           =con_jcal,                                                  &
       eps            =con_eps,                                                   &
       epsm1          =con_epsm1,                                                 &
       rvrdm1         =con_fvirt,                                                 &
       rd             =con_rd,                                                    &
       rhw0           =con_rhw0,                                                  &
       pi             =con_pi,                                                    &
       tgice          =con_tice,                                                  &
       sbc            =con_sbc,                                                   &
!from atm
       ps             =GFS_Data(cdata%blk_no)%Statein%pgr,                        &
       u1             =GFS_Data(cdata%blk_no)%Statein%ugrs(:,1),                  &
       v1             =GFS_Data(cdata%blk_no)%Statein%vgrs(:,1),                  &
       t1             =GFS_Data(cdata%blk_no)%Statein%tgrs(:,1),                  &
       q1             =GFS_Data(cdata%blk_no)%Statein%qgrs(:,1,GFS_Control%ntqv), &
       !************************
       tref           =GFS_Data(cdata%blk_no)%Sfcprop%tref,                       &
       !************************
       cm             =GFS_Interstitial(cdata%thrd_no)%cd_water,                  &
       ch             =GFS_Interstitial(cdata%thrd_no)%cdq_water,                 &
!       lseaspray      =GFS_Control%lseaspray,                                     &
       fm             =GFS_Interstitial(cdata%thrd_no)%ffmm_water,                &
       fm10           =GFS_Interstitial(cdata%thrd_no)%fm10_water,                &

       prsl1          =GFS_Data(cdata%blk_no)%Statein%prsl(:,1),                  &
       prslki         =GFS_Interstitial(cdata%thrd_no)%work3,                     &
       prsik1         =GFS_Data(cdata%blk_no)%Statein%prsik(:,1),                 &
       prslk1         =GFS_Data(cdata%blk_no)%Statein%prslk(:,1),                 &
!       wet            =GFS_Interstitial(cdata%thrd_no)%wet,                       &
!      use_lake_model =GFS_Data(cdata%blk_no)%Sfcprop%use_lake_model,             &
       xlon           =GFS_Data(cdata%blk_no)%Grid%xlon,                          &
       sinlat         =GFS_Data(cdata%blk_no)%Grid%sinlat,                        &
       stress         =GFS_Interstitial(cdata%thrd_no)%stress_water,              &

       sfcemis        =GFS_Data(cdata%blk_no)%Sfcprop%emis_wat,                   &

       dlwflxs        =GFS_Interstitial(cdata%thrd_no)%gabsbdlw_water,            &
       sfcnsw         =GFS_Data(cdata%blk_no)%Intdiag%nswsfci,                    &
       rain           =GFS_Interstitial(cdata%thrd_no)%tprcp_water,               &

       timestep       =GFS_Control%dtf,                                           &
!       kdt            =GFS_Control%kdt,                                           &
       solhr          =GFS_Control%solhr,                                         &
       xcosz          =GFS_Interstitial(cdata%thrd_no)%xcosz,                     &
       wind           =GFS_Interstitial(cdata%thrd_no)%wind,                      &

       flag_iter      =GFS_Interstitial(cdata%thrd_no)%flag_iter,                 &
       flag_guess     =GFS_Interstitial(cdata%thrd_no)%flag_guess,                &

       nstf_name1     =GFS_Control%nstf_name(1),                                  &
       nstf_name4     =GFS_Control%nstf_name(4),                                  &
       nstf_name5     =GFS_Control%nstf_name(5),                                  &
       !
       lprnt          =GFS_Control%lprnt,                                         &
!       ipr            =GFS_Interstitial(cdata%thrd_no)%ipr,                       &
       thsfc_loc      =GFS_Control%thsfc_loc,                                     &
       !************************
       tskin          =GFS_Interstitial(cdata%thrd_no)%tseal,                     &
       tsurf          =GFS_Interstitial(cdata%thrd_no)%tsurf_water,               &
       ! ***********************
       xt             =GFS_Data(cdata%blk_no)%Sfcprop%xt,                         &
       xs             =GFS_Data(cdata%blk_no)%Sfcprop%xs,                         &
       xu             =GFS_Data(cdata%blk_no)%Sfcprop%xu,                         &
       xv             =GFS_Data(cdata%blk_no)%Sfcprop%xv,                         &
       xz             =GFS_Data(cdata%blk_no)%Sfcprop%xz,                         &
       zm             =GFS_Data(cdata%blk_no)%Sfcprop%zm,                         &
       xtts           =GFS_Data(cdata%blk_no)%Sfcprop%xtts,                       &
       xzts           =GFS_Data(cdata%blk_no)%Sfcprop%xzts,                       &
       dt_cool        =GFS_Data(cdata%blk_no)%Sfcprop%dt_cool,                    &
       z_c            =GFS_Data(cdata%blk_no)%Sfcprop%z_c,                        &
!
       c_0            =GFS_Data(cdata%blk_no)%Sfcprop%c_0,                        &
       c_d            =GFS_Data(cdata%blk_no)%Sfcprop%c_d,                        &
       w_0            =GFS_Data(cdata%blk_no)%Sfcprop%w_0,                        &
       w_d            =GFS_Data(cdata%blk_no)%Sfcprop%w_d,                        &
       d_conv         =GFS_Data(cdata%blk_no)%Sfcprop%d_conv,                     &
!       ifd            =GFS_Data(cdata%blk_no)%Sfcprop%ifd,                        &
       qrain          =GFS_Data(cdata%blk_no)%Sfcprop%qrain,                      &
       qsurf          =GFS_Interstitial(cdata%thrd_no)%qss_water,                 &
       gflux          =GFS_Interstitial(cdata%thrd_no)%gflx_water,                &
       cmm            =GFS_Interstitial(cdata%thrd_no)%cmm_water,                 &
       chh            =GFS_Interstitial(cdata%thrd_no)%chh_water,                 &
       eva            =GFS_Interstitial(cdata%thrd_no)%evap_water,                &
       hflx           =GFS_Interstitial(cdata%thrd_no)%hflx_water,                &
       ep             =GFS_Interstitial(cdata%thrd_no)%ep1d_water,                &
!       errmsg         =cdata%errmsg,errflg=cdata%errflg)

 ! temp variables
  tref is used but not changed
  tseal = tsurf, used as variable and updated as tref+dtz
  tskin is updated = tsurf

       call sfc_nst_post_run(im=GFS_Control%blksz(cdata%blk_no),     &
       kdt=GFS_Control%kdt,                                          &
       rlapse=rlapse,                                                &
       tgice=con_tice,                                               &
       wet=GFS_Interstitial(cdata%thrd_no)%wet,                      &
       use_lake_model=GFS_Data(cdata%blk_no)%Sfcprop%use_lake_model, &
       icy=GFS_Interstitial(cdata%thrd_no)%icy,                      &
       oro=GFS_Data(cdata%blk_no)%Sfcprop%oro,                       &
       oro_uf=GFS_Data(cdata%blk_no)%Sfcprop%oro_uf,                 &
       nstf_name1=GFS_Control%nstf_name(1),                          &
       nstf_name4=GFS_Control%nstf_name(4),                          &
       nstf_name5=GFS_Control%nstf_name(5),                          &
       xt=GFS_Data(cdata%blk_no)%Sfcprop%xt,                         &
       xz=GFS_Data(cdata%blk_no)%Sfcprop%xz,                         &
       dt_cool=GFS_Data(cdata%blk_no)%Sfcprop%dt_cool,               &
       z_c=GFS_Data(cdata%blk_no)%Sfcprop%z_c,                       &
       xlon=GFS_Data(cdata%blk_no)%Grid%xlon,                        &
       !************************
       tref=GFS_Data(cdata%blk_no)%Sfcprop%tref,                     &
       tsurf_wat=GFS_Interstitial(cdata%thrd_no)%tsurf_water,        &
       tsfc_wat=GFS_Interstitial(cdata%thrd_no)%tsfc_water,          &
       !************************
       nthreads=GFS_Control%nthreads,                                &
       dtzm=GFS_Interstitial(cdata%thrd_no)%dtzm,                    &
       errmsg=cdata%errmsg,                                          &
       errflg=cdata%errflg)

  !temp variables
  tsfc_wat = tref + ...
  tsurf_wat is unused




  ! pre
  tseal = tsfco  ! what is the point? tseal is over-written by tref + ....
  -> get dtzm_2d
  tref = tsfco + ...dtzm
  tseal = tref + ....
  tsurf_wat = tseal
  ! tref has been calculated, tseal=tsurf_wat=tref +....
  ! tsurf_wat passed out as %tsurf_watER

  ! sfc_nst_run gets tref, tseal (as variable tskin) and tsurf_watER (as variable tsurf)
  ! "tsurf same as tskin but for guess run"
  ! nst run
  tref is used but not changed

  tsurf used to calc evap,hflx "calc latent and sens heat flux over open water with tskin at previous timestep"
  tsurf is used as tsea for LW, density, sen heat due to rain (qrain) ?enthalpy "hrain"
  -> get dtzm_point
  tsurf = ..tref + dtz

  tskin is updated at end w/ variable tsurf
  tskin is used to calculate sens, latent fluxes (evap, hflx)
  tsurf is exported as tsurf_watER


  ! post
  ! gets tsurf_watER (as variable tsurf_wat) but does not use
  ! get tsfc_wat and updates
  ! gets tref and uses but does not change
  calcs dtzm (sends out)
  tsfc_wat = tref + dtzm

  ! tsfc_wat = Interstitial "tsfc_water"
  ! what is this tsfc_wat ? It isn't sent to pre or sfc_nst run
  ! gets passed around as "tskin_wat=GFS_Interstitial(cdata%thrd_no)%tsfc_water or as tsfc_wat"
