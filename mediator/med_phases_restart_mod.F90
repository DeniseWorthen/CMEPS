module med_phases_restart_mod

  !-----------------------------------------------------------------------------
  ! Write/Read mediator restart files
  !-----------------------------------------------------------------------------

  use med_kind_mod            , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_constants_mod       , only : dbug_flag => med_constants_dbug_flag
  use med_utils_mod           , only : chkerr    => med_utils_ChkErr
  use med_internalstate_mod   , only : maintask, logunit, InternalState
  use med_internalstate_mod   , only : ncomps, compname, compocn, complnd, compwav
  use perf_mod                , only : t_startf, t_stopf
  use med_phases_prep_glc_mod , only : FBlndAccum2glc_l, lndAccum2glc_cnt
  use med_phases_prep_glc_mod , only : FBocnAccum2glc_o, ocnAccum2glc_cnt
  use med_phases_prep_rof_mod , only : FBlndAccum2rof_l, lndAccum2rof_cnt
  use pio                     , only : file_desc_t
  implicit none
  private

  public  :: med_phases_restart_read
  public  :: med_phases_restart_write

  private :: med_phases_restart_alarm_init

  logical :: write_restart_at_endofrun = .false.
  logical :: whead(2) = (/.true. , .false./)
  logical :: wdata(2) = (/.false., .true. /)

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine med_phases_restart_alarm_init(gcomp, rc)

    ! --------------------------------------
    ! Initialize mediator restart file alarms (module variables)
    ! --------------------------------------

    use ESMF         , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF         , only : ESMF_Clock, ESMF_ClockGet, ESMF_ClockAdvance, ESMF_ClockSet
    use ESMF         , only : ESMF_Time, ESMF_TimeInterval, ESMF_TimeIntervalGet
    use ESMF         , only : ESMF_Alarm, ESMF_AlarmSet, ESMF_AlarmCreate, ESMF_TimeIntervalSet
    use ESMF         , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR
    use ESMF         , only : ESMF_SUCCESS, ESMF_FAILURE
    use NUOPC        , only : NUOPC_CompAttributeGet
    use NUOPC_Model  , only : NUOPC_ModelGet
    use med_time_mod , only : med_time_AlarmInit
    ! debug
    use ESMF, only : ESMF_ClockPrint, ESMF_AlarmPrint, ESMF_Finalize, ESMF_END_ABORT, ESMF_ClockGetAlarmList
    use ESMF, only : ESMF_ALARMLIST_RINGING, ESMF_ALARMLIST_ALL, ESMF_AlarmList_Flag, ESMF_AlarmIsRinging
    use ESMF, only : ESMF_AlarmGet, ESMF_AlarmRingerOff

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Alarm)        :: alarm
    type(ESMF_Clock)        :: mclock
    type(ESMF_TimeInterval) :: mtimestep
    type(ESMF_Time)         :: mCurrTime
    type(ESMF_TimeInterval) :: AlarmInterval    ! Alarm interval
    type(ESMF_Time)         :: AlarmTime        ! Alarm time for non-interval alarms
    type(ESMF_TimeInterval) :: AlarmTimeStep    ! Clock advance for non-interval alarms
    integer                 :: timestep_length
    character(CL)           :: cvalue          ! attribute string
    character(CL)           :: restart_option  ! freq_option setting (ndays, nsteps, etc)
    integer                 :: restart_n       ! freq_n setting relative to freq_option
    integer                 :: fhcount
    integer                 :: restart_fh(3)
    character(CS)           :: alarm_name
    logical                 :: isPresent
    logical                 :: isSet
    integer                 :: n
    character(len=*), parameter :: subname='(med_phases_restart_alarm_init)'
    ! debug
    character(CL) :: msg
    integer :: alarmcount
    type(ESMF_Alarm), allocatable :: alarms(:)
   !---------------------------------------

    rc = ESMF_SUCCESS

    ! Get model clock
    call NUOPC_ModelGet(gcomp, modelClock=mclock,  rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine restart frequency
    call NUOPC_CompAttributeGet(gcomp, name='restart_option', value=restart_option, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeGet(gcomp, name='restart_n', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) restart_n

    ! Set alarm for instantaneous mediator restart output
    call ESMF_ClockGet(mclock, currTime=mCurrTime,  rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_time_alarmInit(mclock, alarm, option=restart_option, opt_n=restart_n, &
         reftime=mcurrTime, alarmname='alarm_restart', rc=rc)
    call ESMF_AlarmSet(alarm, clock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Advance model clock to trigger alarm then reset model clock back to currtime
    call ESMF_ClockGet(mclock, currTime=mCurrTime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeIntervalGet(mtimestep, s=timestep_length, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockSet(mclock, currTime=mcurrtime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Handle end of run restart
    call NUOPC_CompAttributeGet(gcomp, name="write_restart_at_endofrun", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) .eq. '.true.') write_restart_at_endofrun = .true.
    end if

    !call NUOPC_CompAttributeGet(gcomp, name='restart_fh', isPresent=isPresent, isSet=isSet, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !if (isPresent .and. isSet) then
    !   !call NUOPC_CompAttributeGet(gcomp, name='restart_fh', valueList=restart_fh, itemCount=fhcount, rc=rc)
    !   call NUOPC_CompAttributeGet(gcomp, name='restart_fh:', itemCount=fhcount, rc=rc)
    !   if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !end if
    call ESMF_ClockPrint(mclock, options="currTime", unit=cvalue,rc=rc)
    call ESMF_LogWrite(trim(subname)//": currtime b4 restart_fh alarm set = "//trim(cvalue), ESMF_LOGMSG_INFO)

    restart_fh = (/1,3,5/)
    do n = 1,size(restart_fh)
       write(cvalue,'(i3.3)')restart_fh(n)
       alarm_name = 'alarm_restart'//trim(cvalue)
       ! extra, non-interval alarms at specified forecast hours
       call med_time_alarmInit(mclock, alarm, option='fhours', opt_n=restart_fh(n), &
            reftime=mcurrTime, alarmname=trim(alarm_name), rc=rc)
       call ESMF_AlarmSet(alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Advance model clock to trigger alarm then reset model clock back to currtime
       call ESMF_ClockGet(mclock, currTime=mCurrTime, timeStep=mtimestep, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeIntervalGet(mtimestep, s=timestep_length, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_ClockAdvance(mclock,rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_ClockSet(mclock, currTime=mcurrtime, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do
#ifdef test
    restart_fh = (/1,3,5/)
    do n = 1,size(restart_fh)
       write(cvalue,'(i3.3)')restart_fh(n)
       alarm_name = 'alarm_restart'//trim(cvalue)
       ! set an alarm at a specific time
       call ESMF_TimeIntervalSet(AlarmTimeStep, h=restart_fh(n), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_ClockPrint(mclock, options="currTime", unit=cvalue, rc=rc)
       write(msg,'(a,i6)')"currtime b4 advance  = "//trim(cvalue)//'  '//trim(alarm_name),restart_fh(n)
       call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

       ! advance the clock to get the ring time
       call ESMF_ClockAdvance(mclock, timestep=AlarmTimeStep, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_ClockGet(mclock, currTime=alarmTime, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_ClockPrint(mclock, options="currTime", unit=cvalue, rc=rc)
       write(msg,'(a,i6)')trim(subname)//": currtime at restart_fh alarm  = "//trim(cvalue),restart_fh(n)
       call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

       ! create the alarm and make sure it is off
       alarm = ESMF_AlarmCreate(name=alarm_name, clock=mclock, ringTime=alarmTime, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_AlarmRingerOff( alarm, rc=rc )
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          write(msg,'(a,2i6)')"XXX_alarm_init "//trim(alarm_name)//" is ringing at creation, turn off"
          call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
       end if
       ! set the clock back to the current time and timestep
       call ESMF_ClockSet(mclock, currTime=mcurrTime, timeStep=mtimestep, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_ClockPrint(mclock, options="currTime", unit=cvalue,rc=rc)
       write(msg,'(a,i6)')"currtime after reset  = "//trim(cvalue),restart_fh(n)
       call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

       if (maintask) then
          write(logunit,'(a,i6)') trim(subname) //' creating set time alarm '// trim(alarm_name) &
               //' at hour ',restart_fh(n)
       end if
    end do
    call ESMF_ClockPrint(mclock, options="currTime", unit=cvalue,rc=rc)
    call ESMF_LogWrite(trim(subname)//": currtime af restart_fh alarm set = "//trim(cvalue), ESMF_LOGMSG_INFO)

    ! Write mediator diagnostic output
    if (maintask) then
       write(logunit,*)
       write(logunit,'(a,2x,i8)') trim(subname)//" restart clock timestep = ",timestep_length
       write(logunit,'(a,2x,i8)') trim(subname)//" set restart alarm with option "//&
            trim(restart_option)//" and frequency ",restart_n
       write(logunit,'(a,l7)') trim(subname)//" write_restart_at_endofrun : ", write_restart_at_endofrun
       write(logunit,*)
    end if
#endif
    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_RINGING, alarmCount=alarmcount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(msg,'(a,2i6)')"XXX_alarm_init ringing count = ",alarmcount
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    if (alarmcount > 0) then
       allocate(alarms(alarmcount))
       call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_RINGING, alarmList=alarms, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do n = 1,alarmcount
          call ESMF_AlarmGet(alarms(n), name=alarm_name, rc=rc)
          write(msg,'(a,2i6)')"XXX_alarm_init ringing = "//trim(alarm_name),n,alarmcount
          call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
       end do
       deallocate(alarms)
    end if

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmcount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(msg,'(a,2i6)')"XXX_alarm_init all count = ",alarmcount
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    if (alarmcount > 0) then
       allocate(alarms(alarmcount))
       call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmList=alarms, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do n = 1,alarmcount
          call ESMF_AlarmGet(alarms(n), name=alarm_Name, rc=rc)
          write(msg,'(a,2i6)')"XXX_alarm_init all = "//trim(alarm_name),n,alarmcount
          call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
       end do
       deallocate(alarms)
    end if

  end subroutine med_phases_restart_alarm_init

  !===============================================================================
  subroutine med_phases_restart_write(gcomp, rc)

    ! Write mediator restart

    use ESMF       , only : ESMF_GridComp, ESMF_VM, ESMF_Clock, ESMF_Time, ESMF_Alarm
    use ESMF       , only : ESMF_TimeInterval, ESMF_CalKind_Flag, ESMF_MAXSTR
    use ESMF       , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_FAILURE
    use ESMF       , only : ESMF_LOGMSG_ERROR, operator(==), operator(-)
    use ESMF       , only : ESMF_GridCompGet, ESMF_ClockGet, ESMF_ClockGetNextTime
    use ESMF       , only : ESMF_TimeGet, ESMF_ClockGetAlarm, ESMF_ClockPrint, ESMF_TimeIntervalGet
    use ESMF       , only : ESMF_AlarmIsRinging, ESMF_AlarmRingerOff, ESMF_FieldBundleIsCreated
    use ESMF       , only : ESMF_Calendar, ESMF_AlarmGet, ESMF_AlarmList_Flag, ESMF_ClockGetAlarmList
    use ESMF       , only : ESMF_ALARMLIST_RINGING
    use NUOPC      , only : NUOPC_CompAttributeGet
    use NUOPC_Model, only : NUOPC_ModelGet
    use med_io_mod , only : med_io_define_time, med_io_write_time
    use med_io_mod , only : med_io_write, med_io_wopen, med_io_enddef
    use med_io_mod , only : med_io_close, med_io_date2yyyymmdd, med_io_sec2hms
    use med_phases_history_mod, only : auxcomp
    use med_constants_mod     , only : SecPerDay => med_constants_SecPerDay

    ! Input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(file_desc_t)          :: io_file
    type(ESMF_VM)              :: vm
    type(ESMF_Clock)           :: clock
    type(ESMF_Time)            :: starttime
    type(ESMF_Time)            :: currtime
    type(ESMF_Time)            :: nexttime
    type(ESMF_Time), save      :: lasttimewritten
    type(ESMF_TimeInterval)    :: timediff       ! Used to calculate curr_time
    type(ESMF_Alarm), allocatable :: alarms(:)
    type(ESMF_Alarm)           :: alarm
    type(ESMF_Calendar)        :: calendar
    character(len=CS)          :: alarm_name
    character(len=CS)          :: currtimestr
    character(len=CS)          :: nexttimestr
    type(InternalState)        :: is_local
    integer                    :: m,n,nf,nc      ! counters
    integer                    :: curr_ymd       ! Current date YYYYMMDD
    integer                    :: curr_tod       ! Current time-of-day (s)
    integer                    :: start_ymd      ! Starting date YYYYMMDD
    integer                    :: start_tod      ! Starting time-of-day (s)
    integer                    :: next_ymd       ! Starting date YYYYMMDD
    integer                    :: next_tod       ! Starting time-of-day (s)
    integer                    :: nx,ny          ! global grid size
    integer                    :: yr,mon,day,sec ! time units
    real(R8)                   :: days_since     ! Time interval since start time
    integer                    :: unitn          ! unit number
    character(ESMF_MAXSTR)     :: time_units     ! units of time variable
    character(ESMF_MAXSTR)     :: case_name      ! case name
    character(ESMF_MAXSTR)     :: restart_file   ! Local path to restart filename
    character(ESMF_MAXSTR)     :: restart_pfile  ! Local path to restart pointer filename
    character(ESMF_MAXSTR)     :: cpl_inst_tag   ! instance tag
    character(ESMF_MAXSTR)     :: restart_dir    ! Optional restart directory name
    character(ESMF_MAXSTR)     :: cvalue         ! attribute string
    logical                    :: alarmIsOn      ! generic alarm flag
    real(R8)                   :: tbnds(2)       ! CF1.0 time bounds
    logical                    :: isPresent
    logical                    :: first_time = .true.
    integer                    :: alarmcount
    character(len=*), parameter :: subname='(med_phases_restart_write)'
    !debug
    character(CL) :: msg
    !---------------------------------------

    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeGet(gcomp, name='inst_suffix', isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if(isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='inst_suffix', value=cpl_inst_tag, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       cpl_inst_tag = ""
    endif
    call NUOPC_CompAttributeGet(gcomp, name='restart_dir', isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if(isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='restart_dir', value=restart_dir, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       restart_dir = ""
    endif

    if (first_time) then
       call ESMF_LogWrite("XXX calling alarm_init", ESMF_LOGMSG_INFO)
       call med_phases_restart_alarm_init(gcomp, rc)
       first_time = .false.
       call ESMF_LogWrite("XXX after alarm_init", ESMF_LOGMSG_INFO)
    end if

    !---------------------------------------
    ! --- Get the clock info
    !---------------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!#ifdef test
    ! check if any alarms are ringing and get any ringing alarms
    call ESMF_ClockGetAlarmList(clock, alarmlistflag=ESMF_ALARMLIST_RINGING, alarmCount=alarmcount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(msg,'(a,2i6)')"XXX_restart_mod ringing count = ",alarmcount
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    if (alarmcount > 0) then
       allocate(alarms(alarmcount))
       call ESMF_ClockGetAlarmList(clock, alarmlistflag=ESMF_ALARMLIST_RINGING, alarmList=alarms, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do n = 1,alarmcount
          call ESMF_AlarmGet(alarms(n), name=alarm_Name, rc=rc)
          write(msg,'(a,2i6)')"XXX_restart_mod ringing = "//trim(alarm_name),n,alarmcount
          call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
       end do
       deallocate(alarms)
    end if
       !do n = 1,alarmcount
          !if (ESMF_AlarmIsRinging(alarms(n), rc=rc)) then
             !call ESMF_AlarmGet(alarms(n), name=alarm_Name, rc=rc)
             ! this is a ringing restart alarm
 !             if(trim(alarm_Name(1:13)) == 'alarm_restart') then
!                 alarmIsOn = .true.
!                 call ESMF_AlarmRingerOff( alarms(n), rc=rc )
!                 if (ChkErr(rc,__LINE__,u_FILE_u)) return
!              end if
!           !end if
!           end do
!        end if
! !
!#endif
!#ifdef test
    ! Restart Alarm
    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       alarmIsOn = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       ! Stop Alarm
       call ESMF_ClockGetAlarm(clock, alarmname='alarm_stop', alarm=alarm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (ESMF_AlarmIsRinging(alarm, rc=rc) .and. write_restart_at_endofrun) then
          AlarmIsOn = .true.
       else
          AlarmIsOn = .false.
       endif
    endif
!#endif
    if (alarmIsOn) then
       call ESMF_ClockGet(clock, currtime=currtime, starttime=starttime, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_ClockGetNextTime(clock, nextTime=nexttime, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_ClockGet(clock, calendar=calendar, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(currtime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       write(currtimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
       if (dbug_flag > 1) then
          call ESMF_LogWrite(trim(subname)//": currtime = "//trim(currtimestr), ESMF_LOGMSG_INFO)
       endif
       call ESMF_TimeGet(nexttime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       write(nexttimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
       if (dbug_flag > 1) then
          call ESMF_LogWrite(trim(subname)//": nexttime = "//trim(nexttimestr), ESMF_LOGMSG_INFO)
       endif
       if (maintask) then
          call ESMF_ClockPrint(clock, options="currTime", &
               preString="-------->"//trim(subname)//" mediating for: ", unit=cvalue, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          write(logunit, *) trim(cvalue)
       endif
       timediff = nexttime - starttime
       call ESMF_TimeIntervalGet(timediff, d=day, s=sec, rc=rc)
       days_since = day + sec/real(SecPerDay,R8)

       call ESMF_TimeGet(starttime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ymd2date(yr, mon, day, start_ymd)
       start_tod = sec
       time_units = 'days since '//trim(med_io_date2yyyymmdd(start_ymd))//' '//med_io_sec2hms(start_tod, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_TimeGet(nexttime, yy=yr, mm=mon, dd=day, s=sec)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ymd2date(yr,mon,day,next_ymd)
       next_tod = sec

       call ESMF_TimeGet(currtime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ymd2date(yr,mon,day,curr_ymd)
       curr_tod = sec

       !---------------------------------------
       ! Restart File
       ! Use nexttimestr rather than currtimestr here since that is the time at the end of
       ! the timestep and is preferred for restart file names
       !---------------------------------------

       write(restart_file,"(6a)") trim(restart_dir)//trim(case_name),'.cpl', trim(cpl_inst_tag),'.r.',&
            trim(nexttimestr),'.nc'

       if (maintask) then
          restart_pfile = "rpointer.cpl"//trim(cpl_inst_tag)
          call ESMF_LogWrite(trim(subname)//" write rpointer file = "//trim(restart_pfile), ESMF_LOGMSG_INFO)
          open(newunit=unitn, file=restart_pfile, form='FORMATTED')
          write(unitn,'(a)') trim(restart_file)
          close(unitn)
       endif

       call ESMF_LogWrite(trim(subname)//": write "//trim(restart_file), ESMF_LOGMSG_INFO)
       call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call med_io_wopen(restart_file, io_file, vm, rc, clobber=.true.)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       do m = 1,2
          if (m == 2) then
             call med_io_enddef(io_file)
          end if

          tbnds = days_since
          call ESMF_LogWrite(trim(subname)//": time "//trim(time_units), ESMF_LOGMSG_INFO)
          if (whead(m)) then
             call ESMF_ClockGet(clock, calendar=calendar, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call med_io_define_time(io_file, time_units, calendar, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             call med_io_write_time(io_file, days_since, tbnds=(/days_since,days_since/), nt=1, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

          ! Write out next ymd/tod in place of curr ymd/tod because the
          ! restart represents the time at end of the current timestep
          ! and that is where we want to start the next run.
          call med_io_write(io_file, start_ymd, 'start_ymd', whead(m), wdata(m), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_io_write(io_file, start_tod, 'start_tod', whead(m), wdata(m), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_io_write(io_file, next_ymd , 'curr_ymd' , whead(m), wdata(m), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call med_io_write(io_file, next_tod , 'curr_tod' , whead(m), wdata(m), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          do n = 1,ncomps
             if (is_local%wrap%comp_present(n)) then
                nx = is_local%wrap%nx(n)
                ny = is_local%wrap%ny(n)
                ! Write import field bundles
                if (ESMF_FieldBundleIsCreated(is_local%wrap%FBimp(n,n),rc=rc)) then
                   call med_io_write(io_file, is_local%wrap%FBimp(n,n), whead(m), wdata(m), nx, ny, &
                        nt=1, pre=trim(compname(n))//'Imp', rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                endif
                ! Write export field bundles
                if (ESMF_FieldBundleIsCreated(is_local%wrap%FBexp(n),rc=rc)) then
                   call med_io_write(io_file, is_local%wrap%FBexp(n), whead(m), wdata(m), nx, ny, &
                        nt=1, pre=trim(compname(n))//'Exp', rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                endif
                ! Write fraction field bundles
                if (ESMF_FieldBundleIsCreated(is_local%wrap%FBfrac(n),rc=rc)) then
                   call med_io_write(io_file, is_local%wrap%FBfrac(n), whead(m), wdata(m), nx, ny, &
                        nt=1, pre=trim(compname(n))//'Frac', rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                endif
             end if
          enddo

          ! Write export accumulation to ocn
          if (ESMF_FieldBundleIsCreated(is_local%wrap%FBExpAccumOcn)) then
             nx = is_local%wrap%nx(compocn)
             ny = is_local%wrap%ny(compocn)
             call med_io_write(io_file, is_local%wrap%FBExpAccumOcn, whead(m), wdata(m), nx, ny, &
                  nt=1, pre='ocnExpAccum', rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call med_io_write(io_file, is_local%wrap%ExpAccumOcnCnt, 'ocnExpAccum_cnt', whead(m), wdata(m), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          endif

          ! Write export accumulation to wav
          if (ESMF_FieldBundleIsCreated(is_local%wrap%FBExpAccumWav)) then
             nx = is_local%wrap%nx(compwav)
             ny = is_local%wrap%ny(compwav)
             call med_io_write(io_file, is_local%wrap%FBExpAccumWav, whead(m), wdata(m), nx, ny, &
                  nt=1, pre='wavExpAccum', rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call med_io_write(io_file, is_local%wrap%ExpAccumWavCnt, 'wavExpAccum_cnt', whead(m), wdata(m), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          endif

          ! Write accumulation from lnd to rof if lnd->rof coupling is on
          if (ESMF_FieldBundleIsCreated(FBlndAccum2rof_l)) then
             nx = is_local%wrap%nx(complnd)
             ny = is_local%wrap%ny(complnd)
             call med_io_write(io_file, FBlndAccum2rof_l, whead(m), wdata(m), nx, ny, &
                  nt=1, pre='lndImpAccum2rof', rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call med_io_write(io_file, lndAccum2rof_cnt, 'lndImpAccum2rof_cnt', whead(m), wdata(m), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

          ! Write accumulation from lnd to glc if lnd->glc coupling is on
          if (ESMF_FieldBundleIsCreated(FBlndAccum2glc_l)) then
             nx = is_local%wrap%nx(complnd)
             ny = is_local%wrap%ny(complnd)
             call med_io_write(io_file, FBlndAccum2glc_l, whead(m), wdata(m), nx, ny, &
                  nt=1, pre='lndImpAccum2glc', rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call med_io_write(io_file, lndAccum2glc_cnt, 'lndImpAccum2glc_cnt', whead(m), wdata(m), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

          ! Write accumulation from ocn to glc if ocn->glc coupling is on
          if (ESMF_FieldBundleIsCreated(FBocnAccum2glc_o)) then
             nx = is_local%wrap%nx(compocn)
             ny = is_local%wrap%ny(compocn)
             call med_io_write(io_file, FBocnAccum2glc_o, whead(m), wdata(m), nx, ny, &
                  nt=1, pre='ocnImpAccum2glc_o', rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call med_io_write(io_file, ocnAccum2glc_cnt, 'ocnImpAccum2glc_cnt', whead(m), wdata(m), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

          ! Write ocn albedo field bundle (CESM only)
          if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnalb_o,rc=rc)) then
             nx = is_local%wrap%nx(compocn)
             ny = is_local%wrap%ny(compocn)
             call med_io_write(io_file, is_local%wrap%FBMed_ocnalb_o, whead(m), wdata(m), nx, ny, &
                  nt=1, pre='MedOcnAlb_o', rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

          ! Write auxiliary files accumulation -
          ! For now assume that any time averaged history file has only
          ! one time sample - this will be generalized in the future
          do nc = 2,ncomps
             do nf = 1,auxcomp(nc)%num_auxfiles
                if (auxcomp(nc)%files(nf)%doavg .and. auxcomp(nc)%files(nf)%accumcnt > 0) then
                   nx = is_local%wrap%nx(nc)
                   ny = is_local%wrap%ny(nc)
                   call med_io_write(io_file, auxcomp(nc)%files(nf)%FBaccum, &
                        whead(m), wdata(m), nx, ny, &
                        nt=1, pre=trim(compname(nc))//trim(auxcomp(nc)%files(nf)%auxname), rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   call med_io_write(io_file, auxcomp(nc)%files(nf)%accumcnt, &
                        trim(compname(nc))//trim(auxcomp(nc)%files(nf)%auxname)//'_accumcnt', &
                        whead(m), wdata(m), rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                end if
             end do
          end do

       enddo ! end of whead/wdata loop

       ! Close file
       call med_io_close(io_file, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    !---------------------------------------
    !--- clean up
    !---------------------------------------
    lasttimewritten = currtime
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_phases_restart_write

  !===============================================================================
  subroutine med_phases_restart_read(gcomp, rc)

    ! Read mediator restart

    use ESMF       , only : ESMF_GridComp, ESMF_VM, ESMF_Clock, ESMF_Time, ESMF_MAXSTR
    use ESMF       , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_FAILURE
    use ESMF       , only : ESMF_LOGMSG_ERROR, ESMF_VMBroadCast
    use ESMF       , only : ESMF_GridCompGet, ESMF_ClockGet, ESMF_ClockPrint
    use ESMF       , only : ESMF_FieldBundleIsCreated, ESMF_TimeGet
    use NUOPC      , only : NUOPC_CompAttributeGet
    use med_io_mod , only : med_io_read

    ! Input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_VM)          :: vm
    type(ESMF_Clock)       :: clock
    type(ESMF_Time)        :: currtime
    character(len=CS)      :: currtimestr
    type(InternalState)    :: is_local
    integer                :: n
    integer                :: ierr, unitn
    integer                :: yr,mon,day,sec ! time units
    character(ESMF_MAXSTR) :: case_name      ! case name
    character(ESMF_MAXSTR) :: restart_file   ! Local path to restart filename
    character(ESMF_MAXSTR) :: restart_pfile  ! Local path to restart pointer filename
    character(ESMF_MAXSTR) :: cpl_inst_tag   ! instance tag
    logical                :: isPresent
    character(len=*), parameter :: subname='(med_phases_restart_read)'
    !---------------------------------------
    call t_startf('MED:'//subname)
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get case name and inst suffix
    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeGet(gcomp, name='inst_suffix', isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='inst_suffix', value=cpl_inst_tag, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       cpl_inst_tag = ""
    endif

    ! Get the clock info
    call ESMF_GridCompGet(gcomp, clock=clock)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockGet(clock, currtime=currtime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(currtime,yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(currtimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": currtime = "//trim(currtimestr), ESMF_LOGMSG_INFO)
    endif
    if (maintask) then
       call ESMF_ClockPrint(clock, options="currTime", preString="-------->"//trim(subname)//" mediating for: ", rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    ! Get the restart file name from the pointer file
    restart_pfile = "rpointer.cpl"//trim(cpl_inst_tag)
    if (maintask) then
       call ESMF_LogWrite(trim(subname)//" read rpointer file = "//trim(restart_pfile), ESMF_LOGMSG_INFO)
       open(newunit=unitn, file=restart_pfile, form='FORMATTED', status='old', iostat=ierr)
       if (ierr < 0) then
          call ESMF_LogWrite(trim(subname)//' rpointer file open returns error', ESMF_LOGMSG_INFO)
          rc=ESMF_Failure
          return
       end if
       read (unitn,'(a)', iostat=ierr) restart_file
       if (ierr < 0) then
          call ESMF_LogWrite(trim(subname)//' rpointer file read returns error', ESMF_LOGMSG_INFO)
          rc=ESMF_Failure
          return
       end if
       close(unitn)
       call ESMF_LogWrite(trim(subname)//' restart file from rpointer = '//trim(restart_file), ESMF_LOGMSG_INFO)
    endif
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMBroadCast(vm, restart_file, len(restart_file), 0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(trim(subname)//": read "//trim(restart_file), ESMF_LOGMSG_INFO)

    ! Now read in the restart file
    do n = 1,ncomps
       if (is_local%wrap%comp_present(n)) then
          ! Read import field bundle
          if (ESMF_FieldBundleIsCreated(is_local%wrap%FBimp(n,n),rc=rc)) then
             call med_io_read(restart_file, vm, is_local%wrap%FBimp(n,n), pre=trim(compname(n))//'Imp', rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          endif
          ! Read export field bundle
          if (ESMF_FieldBundleIsCreated(is_local%wrap%FBExp(n),rc=rc)) then
             call med_io_read(restart_file, vm, is_local%wrap%FBexp(n), pre=trim(compname(n))//'Exp', rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          endif
          ! Read fraction field bundles
          if (ESMF_FieldBundleIsCreated(is_local%wrap%FBfrac(n),rc=rc)) then
             call med_io_read(restart_file, vm, is_local%wrap%FBfrac(n), pre=trim(compname(n))//'Frac', rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          endif
       endif
    enddo

    ! Read export field bundle accumulator
    if (ESMF_FieldBundleIsCreated(is_local%wrap%FBExpAccumOcn,rc=rc)) then
       call med_io_read(restart_file, vm, is_local%wrap%FBExpAccumOcn, pre='ocnExpAccum', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call med_io_read(restart_file, vm, is_local%wrap%ExpAccumOcnCnt, 'ocnExpAccum_cnt', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif
    if (ESMF_FieldBundleIsCreated(is_local%wrap%FBExpAccumWav,rc=rc)) then
       call med_io_read(restart_file, vm, is_local%wrap%FBExpAccumWav, pre='wavExpAccum', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call med_io_read(restart_file, vm, is_local%wrap%ExpAccumWavCnt, 'wavExpAccum_cnt', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif
    ! If lnd->rof, read accumulation from lnd to rof (CESM only)
    if (ESMF_FieldBundleIsCreated(FBlndAccum2rof_l)) then
       call med_io_read(restart_file, vm, FBlndAccum2rof_l, pre='lndImpAccum2rof', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call med_io_read(restart_file, vm, lndAccum2rof_cnt, 'lndImpAccum2rof_cnt', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    ! If lnd->glc, read accumulation from lnd to glc (CESM only)
    if (ESMF_FieldBundleIsCreated(FBlndAccum2glc_l)) then
       call med_io_read(restart_file, vm, FBlndAccum2glc_l, pre='lndImpAccum2glc', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call med_io_read(restart_file, vm, lndAccum2glc_cnt, 'lndImpAccum2glc_cnt', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    ! If ocn->glc, read accumulation from ocn to glc (CESM only)
    if (ESMF_FieldBundleIsCreated(FBocnAccum2glc_o)) then
       call med_io_read(restart_file, vm, FBocnAccum2glc_o, pre='ocnImpAccum2glc', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call med_io_read(restart_file, vm, ocnAccum2glc_cnt, 'ocnImpAccum2glc_cnt', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    ! Read ocn albedo field bundle (CESM only)
    if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnalb_o,rc=rc)) then
       call med_io_read(restart_file, vm, is_local%wrap%FBMed_ocnalb_o, pre='MedOcnAlb_o', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    call t_stopf('MED:'//subname)

  end subroutine med_phases_restart_read

  !===============================================================================
  subroutine ymd2date(year,month,day,date)
    ! Converts  year, month, day to coded-date
    ! NOTE: this calendar has a year zero (but no day or month zero)

    integer,intent(in ) :: year,month,day  ! calendar year,month,day
    integer,intent(out) :: date            ! coded (yyyymmdd) calendar date
    !---------------------------------------

    date = abs(year)*10000 + month*100 + day  ! coded calendar date
    if (year < 0) date = -date
  end subroutine ymd2date

end module med_phases_restart_mod
