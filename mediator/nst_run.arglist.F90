
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

       tref           =GFS_Data(cdata%blk_no)%Sfcprop%tref,                       &
!
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
       tskin          =GFS_Interstitial(cdata%thrd_no)%tseal,                     &
       tsurf          =GFS_Interstitial(cdata%thrd_no)%tsurf_water,               &

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
