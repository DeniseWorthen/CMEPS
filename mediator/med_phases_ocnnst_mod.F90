module med_phases_ocnnst_mod

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8, I4=>SHR_KIND_I4
  use med_internalstate_mod , only : InternalState, logunit
  use med_constants_mod     , only : dbug_flag       => med_constants_dbug_flag
  use med_utils_mod         , only : chkerr          => med_utils_chkerr
  use med_methods_mod       , only : FB_diagnose     => med_methods_FB_diagnose
  use med_methods_mod       , only : FB_getFldPtr    => med_methods_FB_getFldPtr
  use med_methods_mod       , only : State_GetScalar => med_methods_State_GetScalar
  use med_internalstate_mod , only : mapconsf, mapnames, compatm, compocn, compice, maintask
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: med_phases_ocnnst_run

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: med_phases_ocnnst_init
  private :: set_ocnnst_pointers

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  type ocnnst_type
     real(r8) , pointer :: lats        (:) => null() ! latitudes  (degrees)
     real(r8) , pointer :: lons        (:) => null() ! longitudes (degrees)
     !real(r8) , pointer :: area        (:) => null() ! area (m2)
     integer  , pointer :: mask        (:) => null() ! ocn domain mask: 0 <=> inactive cell
     !inputs from atm
     real(r8) , pointer :: prsl        (:) => null() ! surface pressure (Pa)
     real(r8) , pointer :: u1          (:) => null() ! zonal component of surface layer wind (m/s)
     real(r8) , pointer :: v1          (:) => null() ! merid component of surface layer wind (m/s)
     real(r8) , pointer :: t1          (:) => null() ! surface layer mean temperature (K)
     real(r8) , pointer :: q1          (:) => null() ! surface layer mean spec humidity (kg/kg)
     real(r8) , pointer :: z1          (:) => null() ! layer 1 height above ground (not MSL) (m)
     !real(r8) , pointer :: u10m        (:) => null() ! zonal component of 10m wind (m/s)
     !real(r8) , pointer :: v10m        (:) => null() ! merid component of 10m wind (m/s)
     real(r8) , pointer :: dlwflx      (:) => null() ! total sky sfc downward lw flux (W/m2)
     real(r8) , pointer :: rain        (:) => null() ! rainfall rate (kg/m2/s)
     real(r8) , pointer :: evap        (:) => null() ! evap rate (!TODO)
     real(r8) , pointer :: sen         (:) => null() ! sensible heat flux (W/m2)
     !real(r8) , pointer :: taux        (:) => null() ! zonal momentum stress (!TODO)
     !real(r8) , pointer :: tauy        (:) => null() ! merid momentum stress (!TODO)
     real(r8) , pointer :: sfcnsw      (:) => null() ! net SW down (calculated from bands)
     real(r8) , pointer :: wind        (:) => null() ! wind magnitude, from Sa_u,Sa_v
     real(r8) , pointer :: stress      (:) => null() ! stress magnitude, from Faxa_taux,Faxa_tauy
     ! input from ice
     real(r8) , pointer :: ifrac       (:) => null() ! sea ice fraction (nd); ofrac=1.0-ifrac
     ! input from ocean
     real(r8) , pointer :: tsfco       (:) => null() ! sea surface temperature (K)

     ! nst variables
     real(r8) , pointer :: tseal       (:) => null() ! ocean surface skin temperature (K)
     real(r8) , pointer :: tsfc_wat    (:) => null() ! surface skin temperature over water (K)
     real(r8) , pointer :: tsurf_wat   (:) => null() ! surface skin temperature after iteration over water (K)
     real(r8) , pointer :: dtzm        (:) => null() ! mean of dT(z)  (z1 to z2) (?)
     real(r8) , pointer :: dtm         (:) => null() ! ?
     ! in sfcf file, need mapping back to atm
     real(r8) , pointer :: c0          (:) => null() ! coefficient1 to calculate d(tz)/d(ts) (nd)
     real(r8) , pointer :: cd          (:) => null() ! coefficient2 to calculate d(tz)/d(ts) (nd)
     real(r8) , pointer :: dconv       (:) => null() ! thickness of free convection layer (m)
     real(r8) , pointer :: dtcool      (:) => null() ! sub-layer cooling amount (K)
     real(r8) , pointer :: qrain       (:) => null() ! sensible heat flux due to rainfall (W)
     real(r8) , pointer :: tref        (:) => null() ! sea surface reference temperature (K)
     real(r8) , pointer :: w0          (:) => null() ! coefficient3 to calculate d(tz)/d(ts) (nd)
     real(r8) , pointer :: wd          (:) => null() ! coefficient4 to calculate d(tz)/d(ts) (nd)

     real(r8) , pointer :: xs          (:) => null() ! salinity  content in diurnal thermocline layer (ppt m)
     real(r8) , pointer :: xt          (:) => null() ! heat content in diurnal thermocline layer (K m)
     real(r8) , pointer :: xtts        (:) => null() ! d(xt)/d(ts) (m)

     real(r8) , pointer :: xu          (:) => null() ! u-current content in diurnal thermocline layer (m2 s-1)
     real(r8) , pointer :: xv          (:) => null() ! v-current  content in diurnal thermocline layer (m2 s-1)
     real(r8) , pointer :: xz          (:) => null() ! diurnal thermocline layer thickness (m)
     real(r8) , pointer :: xzts        (:) => null() ! d(xz)/d(ts) (m K-1)
     real(r8) , pointer :: zc          (:) => null() ! sub-layer cooling thickness (m)
  end type ocnnst_type

  ! used, reused in module
  real(r8) , allocatable, save :: ifd         (:)   ! index to start DTM run or not
  logical  , allocatable, save :: flag_guess  (:)   ! .true.=  guess step to get CD et al
                                                    ! when iter = 1, flag_guess = .true. when wind < 2
                                                    ! when iter = 2, flag_guess = .false. for all grids
  logical  , allocatable, save :: flag_iter   (:)   ! execution or not
                                                    ! when iter = 1, flag_iter = .true. for all grids
                                                    ! when iter = 2, flag_iter = .true. when wind < 2
                                                    ! for both land and ocean (when nstf_name1 > 0)
  real(R8), parameter    :: tgice = 271.20_R8       ! TODO: actually f(sss)
  integer, parameter     :: kp = R8
  real(R8), parameter    :: zero = 0.0_R8, one = 1.0_R8
  character(*),parameter :: u_FILE_u =  __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine med_phases_ocnnst_init(gcomp, ocnnst, rc)

    !-----------------------------------------------------------------------
    ! Initialize pointers to the module variables and then use the module
    ! variables in the med_ocnnst phase
    ! All input field bundles are ASSUMED to be on the ocean grid
    !-----------------------------------------------------------------------

    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_FAILURE
    use ESMF  , only : ESMF_VM, ESMF_VMGet, ESMF_Mesh, ESMF_MeshGet
    use ESMF  , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_MESHLOC_ELEMENT
    use ESMF  , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate, ESMF_FieldBundleGet
    use ESMF  , only : ESMF_TYPEKIND_LOGICAL, ESMF_TYPEKIND_R8, ESMF_TYPEKIND_I4
    use NUOPC , only : NUOPC_CompAttributeGet
    use ESMF  , only : operator(==)

    ! Arguments
    type(ESMF_GridComp)               :: gcomp
    type(ocnnst_type) , intent(inout) :: ocnnst
    integer           , intent(out)   :: rc
    !
    ! Local variables
    type(ESMF_VM)            :: vm
    integer                  :: iam
    type(ESMF_Mesh)          :: lmesh
    integer                  :: n
    integer                  :: lsize
    integer                  :: spatialDim
    integer                  :: numOwnedElements
    integer                  :: fieldcount
    type(InternalState)      :: is_local
    real(R8), pointer        :: ownedElemCoords(:)
    real(r8), pointer        :: dataptr1d(:)
    character(len=CL)        :: tempc1,tempc2
    character(len=CS)        :: cvalue
    logical                  :: isPresent, isSet
    character(CL)            :: msg
    type(ESMF_Field), pointer :: fieldlist(:)
    character(*), parameter  :: subname = '(med_phases_ocnnst_init) '
    !-----------------------------------------------------------------------

    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    ! The following is for debugging
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=iam, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from gcomp
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ocean surface temperature from ocean
    call FB_GetFldPtr(is_local%wrap%FBImp(compocn,compocn), 'So_t', ocnnst%tsfco, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------
    ! Get lat, lon, which are time-invariant
    !----------------------------------

    ! The following assumes that all fields in FBMed_ocnnst_o have the same grid - so
    ! only need to query field 1
    call ESMF_FieldBundleGet(is_local%wrap%FBMed_ocnnst_o, fieldCount=fieldCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldlist(fieldcount))
    call ESMF_FieldBundleGet(is_local%wrap%FBMed_ocnnst_o, fieldlist=fieldlist, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(fieldlist(1), mesh=lmesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    deallocate(fieldlist)
    call ESMF_MeshGet(lmesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lsize = size(ocnnst%tsfco)
    if (numOwnedElements /= lsize) then
       write(tempc1,'(i10)') numOwnedElements
       write(tempc2,'(i10)') lsize
       call ESMF_LogWrite(trim(subname)//": ERROR numOwnedElements "// trim(tempc1) // &
            " not equal to local size "// trim(tempc2), ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if
    allocate(ownedElemCoords(spatialDim*numOwnedElements))
    allocate(ocnnst%lons(numOwnedElements))
    allocate(ocnnst%lats(numOwnedElements))
    call ESMF_MeshGet(lmesh, ownedElemCoords=ownedElemCoords)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do n = 1,lsize
       ocnnst%lons(n) = ownedElemCoords(2*n-1)
       ocnnst%lats(n) = ownedElemCoords(2*n)
    end do

    ! ocean mask
    call FB_GetFldPtr(is_local%wrap%FBImp(compocn,compocn), 'So_omask', dataptr1d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ocnnst%mask(numOwnedElements))
    do n = 1,size(dataptr1d)
       if (dataptr1d(n) == 0._r8) then
          ocnnst%mask(n) = 0
       else
          ocnnst%mask(n) = 1
       end if
    enddo
    ! initialize flags
    allocate(ifd(1:size(dataptr1d)))
    allocate(flag_iter(1:size(dataptr1d)))
    allocate(flag_guess(1:size(dataptr1d)))
    ifd = 0.0_r8
    flag_iter(:) = .true.
    flag_guess(:) = .false.

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_phases_ocnnst_init

  !===============================================================================

  subroutine med_phases_ocnnst_run(gcomp, rc)

    !-----------------------------------------------------------------------
    ! Compute ocean NST (on the ocean grid)
    !-----------------------------------------------------------------------

    use NUOPC_Mediator, only : NUOPC_MediatorGet
    use ESMF          , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_TimeInterval
    use ESMF          , only : ESMF_Clock, ESMF_ClockGet, ESMF_Time, ESMF_TimeGet
    use ESMF          , only : ESMF_ClockIsCreated, ESMF_ClockGetNextTime
    use ESMF          , only : ESMF_VM, ESMF_VMGet
    use ESMF          , only : ESMF_LogWrite, ESMF_LogFoundError
    use ESMF          , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_LOGMSG_INFO
    use ESMF          , only : ESMF_Field, ESMF_FieldGet
    use ESMF          , only : ESMF_FieldBundleGet, ESMF_FieldBundleIsCreated
    use ESMF          , only : operator(+)
    use NUOPC         , only : NUOPC_CompAttributeGet
    use med_phases_history_mod, only : med_phases_history_write_med
    use module_nst_water_prop , only : get_dtzm_point

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    ! local variables
    type(ocnnst_type), save :: ocnnst
    type(ESMF_VM)           :: vm
    integer                 :: iam
    integer                 :: iter, day, jday
    type(InternalState)     :: is_local
    type(ESMF_Clock)        :: clock
    type(ESMF_Clock)        :: dclock
    type(ESMF_Time)         :: currTime
    type(ESMF_Time)         :: nextTime
    type(ESMF_TimeInterval) :: timeStep
    character(CL)           :: cvalue
    character(CS)           :: starttype        ! config start type
    character(CL)           :: runtype          ! initial, continue, hybrid, branch
    real(R8)                :: nextsw_cday      ! calendar day of next atm shortwave
    real(R8)                :: solhr            ! fcst hour at the end of prev time step (currTime)
    real(R8)                :: z_c_0, zsea1, zsea2
    real(R8)                :: tem2
    integer                 :: lsize            ! local size
    integer                 :: i                ! indices
    real(R8), parameter     :: omz1 = 2.0_R8
    character(CL)           :: msg
    logical         , save  :: ocnnst_created
    logical         , save  :: first_call = .true.
    character(len=*)  , parameter :: subname='(med_phases_ocnnst_run)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! Determine main task
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=iam, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! TODO: revisit this initialization
    ! Determine if ocnnst data type will be initialized - and if not return
    if (first_call) then
       if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnnst_o, rc=rc)) then

          ! Initialize ocean NST calculation
          call med_phases_ocnnst_init(gcomp, ocnnst, rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ocnnst_created = .true.
       else
          ocnnst_created = .false.
       end if

       ! Now set first_call to .false.
       first_call = .false.
    end if

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    call t_startf('MED:'//subname)

    if (ocnnst_created) then
       ! get clock, current time and timestep
       call ESMF_GridCompGet(gcomp, clock=clock)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_ClockGet( clock, currTime=currTime, timeStep=timeStep, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ! Clock is not advanced until the end of ModelAdvance
       call ESMF_TimeGet( currTime, dd=day, h_r8=solhr, rc=rc )
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet( currTime, dayOfYear=jday, rc=rc )
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write(msg,*)trim(subname)//' jday = ',jday,' solhr = ',solhr
       call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
       !
       ! Calculate ocean NST on the ocean grid
       !
       lsize = size(ocnnst%mask)
       ! ice fraction
       call FB_GetFldPtr(is_local%wrap%FBImp(compice,compice), 'Si_ifrac', ocnnst%ifrac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! ocean surface temperature from ocean
       call FB_GetFldPtr(is_local%wrap%FBImp(compocn,compocn),  'So_t', ocnnst%tsfco, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! set pointers
       call set_ocnnst_pointers(is_local%wrap%FBImp(compatm,compocn), is_local%wrap%FBImp(compice,compocn), &
            is_local%wrap%FBMed_ocnalb_o, is_local%wrap%FBMed_ocnnst_o, lsize, ocnnst, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBMed_ocnnst_o, string=trim(subname)//' b4 FBMed_ocnnst_o', rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if

       ! NST
       do iter = 1,2
          write(msg,*)trim(subname)//' julian day ',jday,' solhr = ',solhr,' iter = ',iter
          call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

          ! loop_control_part1
          do i = 1,lsize
             if (ocnnst%mask(i) == 1) then
                if (iter == 1 .and. ocnnst%wind(i) < 2.0d0) then
                   flag_guess(i) = .true.
                end if
             end if
          end do

          ! nst_pre
          z_c_0 = zero
          do i = 1,lsize
             if (ocnnst%mask(i) == 1) then
                call get_dtzm_point(ocnnst%xt(i), ocnnst%xz(i), ocnnst%dtcool(i), z_c_0, zero, omz1, ocnnst%dtzm(i))
                ocnnst%tref(i) = max(tgice, ocnnst%tsfco(i) - ocnnst%dtzm(i))
                if (abs(ocnnst%xz(i)) > zero) then
                   tem2 = one / ocnnst%xz(i)
                else
                   tem2 = zero
                endif
                ocnnst%tseal(i)     = ocnnst%tref(i) + (ocnnst%xt(i)+ocnnst%xt(i)) * tem2 - ocnnst%dtcool(i)
                ocnnst%tsurf_wat(i) = ocnnst%tseal(i)
             endif
          enddo

          ! nst_run
          !zsea1 = 0.001_kp*real(nstf_name4)
          !zsea2 = 0.001_kp*real(nstf_name5)
          ! nstf_name4,5 are both 0 right now
          zsea1 = 0.0_r8
          zsea2 = 0.0_r8
          call sfc_nst_run(lsize, mask=ocnnst%mask, flag_iter=flag_iter, flag_guess=flag_guess, ifd=ifd,             &
               zsea1=zsea1, zsea2=zsea2, xlon=ocnnst%lons, solhr=solhr, z1=ocnnst%z1, t1=ocnnst%t1, ps=ocnnst%prsl,  &
               u1=ocnnst%u1, v1=ocnnst%v1, q1=ocnnst%q1, hflx=ocnnst%sen, evap=ocnnst%evap, rain=ocnnst%rain,        &
               dlwflx=ocnnst%dlwflx, sfcnsw=ocnnst%sfcnsw, stress=ocnnst%stress, wind=ocnnst%wind, tref=ocnnst%tref, &
               tskin=ocnnst%tseal, tsurf=ocnnst%tsurf_wat, xt=ocnnst%xt, xs=ocnnst%xs, xu=ocnnst%xu,                 &
               xv=ocnnst%xv, xz=ocnnst%xz, xtts=ocnnst%xtts, xzts=ocnnst%xzts, dt_cool=ocnnst%dtcool,                &
               z_c=ocnnst%zc, c_0=ocnnst%c0, c_d=ocnnst%cd, w_0=ocnnst%w0, w_d=ocnnst%wd, d_conv=ocnnst%dconv)
          ! nst_post
          do i = 1,lsize
             if (ocnnst%mask(i) == 1) then
                call get_dtzm_point(ocnnst%xt(i), ocnnst%xz(i), ocnnst%dtcool(i), ocnnst%zc(i), zsea1, zsea2, &
                     ocnnst%dtzm(i))
                ocnnst%tsfc_wat(i) = max(tgice, ocnnst%tref(i) + ocnnst%dtzm(i))
             end if
          end do

          ! loop_control_part2
          do i = 1, lsize
             if (ocnnst%mask(i) == 0) then
                flag_iter(i)  = .false.
                flag_guess(i) = .false.
             else
                if (iter == 1 .and. ocnnst%wind(i) < 2.0d0) then
                   flag_iter(i) = .true.
                endif
             endif
          end do
       end do
       write(msg,*)trim(subname)//' done iter loop '
       call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

       if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnnst_o, rc=rc)) then
          call NUOPC_MediatorGet(gcomp, driverClock=dClock, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (ESMF_ClockIsCreated(dclock)) then
             call med_phases_history_write_med(gcomp, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if

       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBMed_ocnnst_o, string=trim(subname)//' FBMed_ocnnst_o', rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
       call t_stopf('MED:'//subname)
    end if

  end subroutine med_phases_ocnnst_run

!===============================================================================

!===============================================================================

  subroutine set_ocnnst_pointers(fldbun_a, fldbun_o, fldbun_alb, fldbun_m, lsize, ocnnst, rc)

    use ESMF  , only : ESMF_FieldBundle, ESMF_SUCCESS

    ! Set pointers for ocnnst

    use med_methods_mod , only : FB_fldchk    => med_methods_FB_FldChk

    ! input/output variables
    type(ESMF_FieldBundle)     , intent(inout) :: fldbun_a
    type(ESMF_FieldBundle)     , intent(inout) :: fldbun_o
    type(ESMF_FieldBundle)     , intent(inout) :: fldbun_alb
    type(ESMF_FieldBundle)     , intent(inout) :: fldbun_m
    integer                    , intent(in)    :: lsize
    type(ocnnst_type)          , intent(inout) :: ocnnst
    integer                    , intent(out)   :: rc

    ! local variables
    real(R8), pointer   :: Faxa_taux(:)
    real(R8), pointer   :: Faxa_tauy(:)
    real(R8), pointer   :: avsdr(:)
    real(R8), pointer   :: avsdf(:)
    real(R8), pointer   :: anidr(:)
    real(R8), pointer   :: anidf(:)
    real(R8), pointer   :: Faxa_swvdf(:)
    real(R8), pointer   :: Faxa_swndf(:)
    real(R8), pointer   :: Faxa_swvdr(:)
    real(R8), pointer   :: Faxa_swndr(:)
    real(R8), pointer   :: Fioi_swpen_vdr(:)
    real(R8), pointer   :: Fioi_swpen_vdf(:)
    real(R8), pointer   :: Fioi_swpen_idr(:)
    real(R8), pointer   :: Fioi_swpen_idf(:)

    real(R8) :: fswabsv, fswabsi, fswpen
    integer  :: n

    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !----------------------------------
    ! Set pointers to fields needed for NST calculations
    !----------------------------------

    ! Import from atm
    call FB_GetFldPtr(fldbun_a, 'Sa_pbot' ,   ocnnst%prsl,   rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Sa_z' ,      ocnnst%z1,     rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Sa_shum' ,   ocnnst%q1,     rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Sa_tbot' ,   ocnnst%t1,     rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Sa_u' ,      ocnnst%u1,     rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Sa_v' ,      ocnnst%v1,     rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Faxa_evap' , ocnnst%evap,   rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Faxa_rain' , ocnnst%rain,   rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Faxa_sen' ,  ocnnst%sen,    rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Faxa_rain' , ocnnst%rain,   rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Faxa_lwdn' , ocnnst%dlwflx, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Compute |stress| and |wind| for ocean
    !---------------------------------------

    call FB_GetFldPtr(fldbun_m, 'Snst_wind' , ocnnst%wind, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ocnnst%wind = 0.0_R8
    ocnnst%wind = sqrt(ocnnst%u1*ocnnst%u1 + ocnnst%v1*ocnnst%v1)

    call FB_GetFldPtr(fldbun_a, 'Faxa_taux' , Faxa_taux, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Faxa_tauy' , Faxa_tauy, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call FB_GetFldPtr(fldbun_m, 'Snst_stress' , ocnnst%stress, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ocnnst%stress = 0.0_R8
    ocnnst%stress = sqrt(Faxa_taux*Faxa_taux + Faxa_tauy*Faxa_tauy)

    !---------------------------------------
    ! Compute netsw for ocean
    !---------------------------------------
    ! netsw_for_ocn = downsw_from_atm * (1-ocn_albedo) * (1-ice_fraction) + pensw_from_ice * (ice_fraction)

    ! Input from atm
    call FB_GetFldPtr(fldbun_a, 'Faxa_swvdr', Faxa_swvdr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Faxa_swndr', Faxa_swndr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Faxa_swvdf', Faxa_swvdf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Faxa_swndf', Faxa_swndf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! Input from mediator, ocean albedos
    call FB_GetFldPtr(fldbun_alb, 'So_avsdr' , avsdr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_alb, 'So_anidr' , anidr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_alb, 'So_avsdf' , avsdf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_alb, 'So_anidf' , anidf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! Input from ice, sw pen through ice
    call FB_GetFldPtr(fldbun_o, 'Fioi_swpen_vdr', Fioi_swpen_vdr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_o, 'Fioi_swpen_vdf', Fioi_swpen_vdf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_o, 'Fioi_swpen_idr', Fioi_swpen_idr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_o, 'Fioi_swpen_idf', Fioi_swpen_idf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call FB_GetFldPtr(fldbun_m, 'Snst_sfcnsw' , ocnnst%sfcnsw, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ocnnst%sfcnsw = 0.0_R8
    do n = 1,lsize
       ! Compute total swnet to ocean
       fswpen   = Fioi_swpen_vdr(n) + Fioi_swpen_vdf(n) + Fioi_swpen_idr(n) + Fioi_swpen_idf(n)
       fswabsv  = Faxa_swvdr(n) * (1.0_R8 - avsdr(n)) + Faxa_swvdf(n) * (1.0_R8 - avsdf(n))
       fswabsi  = Faxa_swndr(n) * (1.0_R8 - anidr(n)) + Faxa_swndf(n) * (1.0_R8 - anidf(n))
       ocnnst%sfcnsw(n) = (fswabsv + fswabsi)*(1.0_R8 - ocnnst%ifrac(n)) + fswpen*ocnnst%ifrac(n)
    end do

    ! Ocean NST fields (mapped back to ATM)
    call FB_GetFldPtr(fldbun_m, 'Snst_tref' ,        ocnnst%tref,   rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_dconv' ,       ocnnst%dconv,  rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_dtcool' ,      ocnnst%dtcool, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_qrain' ,       ocnnst%qrain,  rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_xtts' ,        ocnnst%xtts,   rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_xzts' ,        ocnnst%xzts,   rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_c0' ,          ocnnst%c0,     rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_cd' ,          ocnnst%cd,     rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_w0' ,          ocnnst%w0,     rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_wd' ,          ocnnst%wd,     rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_xs' ,          ocnnst%xs,     rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_xt' ,          ocnnst%xt,     rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_xu' ,          ocnnst%xu,     rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_xv' ,          ocnnst%xv,     rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_xz' ,          ocnnst%xz,     rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_zc' ,          ocnnst%zc,     rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call FB_GetFldPtr(fldbun_m, 'Snst_tseal' ,       ocnnst%tseal,     rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_tsfc_water' ,  ocnnst%tsfc_wat,  rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_tsurf_water' , ocnnst%tsurf_wat, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_dtzm' ,        ocnnst%dtzm,      rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_dtm' ,         ocnnst%dtm,       rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine set_ocnnst_pointers
  ! call sfc_nst_run(lsize, mask=ocnnst%mask, flag_iter=flag_iter, flag_guess=flag_guess,              &
  !      zsea1=zsea1, zsea2=zsea2, xlon=ocnnst%lons, solhr=solhr, t1=ocnnst%t1, prsl1=ocnnst%prsl,     &
  !      u1=ocnnst%u1, v1=ocnnst%v1, q1=ocnnst%q1, hflx=ocnnst%sen, evap=ocnnst%evap, sfcnsw=ocnnst%sfcnsw, &
  !      stress=ocnnst%stress, wind=ocnnst%wind, &
  !      tskin=ocnnst%tseal, tsurf=ocnnst%tsurf_wat, xt=ocnnst%xt, xs=ocnnst%xs, xu=ocnnst%xu,         &
  !      xv=ocnnst%xv, xz=ocnnst%xz, xtts=ocnnst%xtts, xzts=ocnnst%xzts, dt_cool=ocnnst%dt_cool,       &
  !      z_c=ocnnst%zc, c_0=ocnnst%c_0, c_d=ocnnst%c_d, w_0=ocnnst%w_0, w_d=ocnnst%w_d, d_conv=ocnnst%d_conv)

  !TODO: add SSS
  subroutine sfc_nst_run(im, mask, flag_iter, flag_guess, ifd, zsea1, zsea2, xlon, solhr, z1, t1, ps,  &
       u1, v1, q1, evap, rain, hflx, dlwflx, sfcnsw, stress, wind, tref, tskin, tsurf, xt, xs, xu, xv, &
       xz, xtts, xzts, dt_cool, z_c, c_0, c_d, w_0, w_d, d_conv)

    use funcphys              , only : fpvs
    use module_nst_water_prop , only : get_dtzm_point, density, rhocoef, sw_ps_9b, sw_ps_9b_aw, grv
    use module_nst_parameters , only : eps, epsm1, t0k,cp_w,omg_m,omg_sh, sigma_r,solar_time_6am,ri_c, &
         z_w_max,delz,wd_max, rad2deg,const_rot,tau_min,tw_max,sst_max,rvrdm1,rd,cp_a, hvap
    use nst_module, only : cool_skin,dtm_1p,cal_w,cal_ttop, convdepth,dtm_1p_fca,dtm_1p_tla,           &
         dtm_1p_mwa,dtm_1p_mda,dtm_1p_mta, dtl_reset

    integer, intent(in) :: im
    integer, intent(in) :: mask(:)
    logical, intent(in) :: flag_iter(:)
    logical, intent(in) :: flag_guess(:)

    real (R8), intent(in) :: zsea1, zsea2, solhr
    real (R8), intent(in) :: xlon(:)
    real (R8), intent(in) :: z1(:)
    real (R8), intent(in) :: t1(:)
    real (R8), intent(in) :: u1(:)
    real (R8), intent(in) :: v1(:)
    real (R8), intent(in) :: q1(:)
    real (R8), intent(in) :: ps(:)
    real (R8), intent(in) :: hflx(:)
    real (R8), intent(in) :: evap(:)
    real (R8), intent(in) :: rain(:)
    real (R8), intent(in) :: dlwflx(:)
    real (R8), intent(in) :: sfcnsw(:)
    real (R8), intent(in) :: stress(:)
    real (R8), intent(in) :: wind(:)
    real (R8), intent(in) :: tref(:)

    real (R8), intent(inout) :: ifd(:)
    real (R8), intent(inout) :: tskin(:)
    real (R8), intent(inout) :: tsurf(:)
    real (R8), intent(inout) :: xt(:)  , xs(:)  , xu(:)     , xv(:)    , xz(:)
    real (R8), intent(inout) :: xtts(:), xzts(:), dt_cool(:), d_conv(:)
    real (R8), intent(inout) :: z_c(:) , c_0(:) , c_d(:)
    real (R8), intent(inout) :: w_0(:) , w_d(:)

    ! local variables
    integer :: i, kdt
    ! TODO
    real (R8) :: timestep
    logical   :: flag(im)
    real (R8) :: tsea, t12, alon, soltim, sss, le, dwat, dtmp, wetc, alfac, tem
    real (R8) :: f_nsol, sep, ustar_a, rnl_ts, hs_ts, sbc, rho_w, cp, grav, alpha, beta
    real (R8) :: fw, fc, taux, tauy, q_ts, sina, cosa, dz, q_warm, ttop0, t0, dta, dtz
    real (R8) :: cpinv, hvapi, elocp, hl_ts, rf_ts, rich, ttop, sstc
    real (R8), dimension(im) :: xt_old, xs_old, xu_old, xv_old, xz_old, xtts_old
    real (R8), dimension(im) :: xzts_old, ifd_old, tskin_old, dt_cool_old, z_c_old
    real (R8), dimension(im) :: wndmag, ulwflx, nswsfc, q0, tv1, rho_a, sfcemis, sinlat
    real (R8), dimension(im) :: qss, rch

    !TODO: do what with this?
    real (R8), dimension(im) :: qrain

    kdt = 1 ! only used for prints
    timestep = 1.0_R8 ! only used for prints
    cp = cp_a
    sss = 34.0_R8
    sbc = sigma_r
    sfcemis = 0.97_R8
    sinlat = sin(xlon)
    cpinv = one/cp
    hvapi = one/hvap
    elocp = hvap/cp
    q0 = zero
    tv1 = zero
    rho_a = zero
    rch = zero
    qss = zero
!
! flag for open water and where the iteration is on
!
      flag = .false.
      do i = 1,im
         flag(i) = (mask(i) == 1 .and. flag_iter(i))
      end do
!
!  save nst-related prognostic fields for guess run
!
      do i=1, im
         if (mask(i) == 1 .and. flag_guess(i)) then
          xt_old(i)      = xt(i)
          xs_old(i)      = xs(i)
          xu_old(i)      = xu(i)
          xv_old(i)      = xv(i)
          xz_old(i)      = xz(i)
          !zm_old(i)      = zm(i)
          xtts_old(i)    = xtts(i)
          xzts_old(i)    = xzts(i)
          ifd_old(i)     = ifd(i)
          tskin_old(i)   = tskin(i)
          dt_cool_old(i) = dt_cool(i)
          z_c_old(i)     = z_c(i)
       endif
     enddo
     !  --- ...  initialize variables. all units are m.k.s. unless specified.
     !           ps is in pascals, wind is wind speed, theta1 is surface air
     !           estimated from level 1 temperature, rho_a is air density and
     !           qss is saturation specific humidity at the water surface
     do i = 1, im
        if ( flag(i) ) then
           nswsfc(i) = sfcnsw(i) ! net solar radiation at the air-sea surface (positive downward)
           wndmag(i) = wind(i)
           q0(i)     = max(q1(i), 1.0e-8_kp)
           tv1(i)    = t1(i) * (one + rvrdm1*q0(i))
           print *,'XXX ',i,tv1(i),q0(i),ps(i),rd*tv1(i)
           rho_a(i)  = ps(i) / (rd*tv1(i))
           qss(i)    = fpvs(tsurf(i))                          ! pa
           qss(i)    = eps*qss(i) / (ps(i) + epsm1*qss(i))     ! pa
           rch(i)    = evap(i)/(elocp * (qss(i) - q0(i)))      ! iffy
           tem       = 0.00001_kp/z1(i)                        ! from sfc_diff
           rch(i)    = max(rch(i), tem)
        end if
     end do
     !
     ! run nst model: dtm + slm
     !
     !> - Call module_nst_water_prop::density() to compute sea water density.
     !> - Call module_nst_water_prop::rhocoef() to compute thermal expansion
     !! coefficient (\a alpha) and saline contraction coefficient (\a beta).
     do i = 1, im
        if ( flag(i) ) then
           tsea      = tsurf(i)
           t12       = tsea*tsea
           ulwflx(i) = sfcemis(i) * sbc * t12 * t12
           alon      = xlon(i)*rad2deg
           grav      = grv(sinlat(i))
           soltim  = mod(alon/15.0_kp + solhr, 24.0_kp)*3600.0_kp
           call density(tsea,sss,rho_w)                     ! sea water density
           call rhocoef(tsea,sss,rho_w,alpha,beta)          ! alpha & beta
           !
           !> - Calculate sensible heat flux (\a qrain) due to rainfall.
           !
           le       = (2.501_kp-0.00237_kp*tsea)*1e6_kp
           dwat     = 2.11e-5_kp*(t1(i)/t0k)**1.94_kp               ! water vapor diffusivity
           dtmp     = (one+3.309e-3_kp*(t1(i)-t0k)-1.44e-6_kp*(t1(i)-t0k) &
                * (t1(i)-t0k))*0.02411_kp/(rho_a(i)*cp)             ! heat diffusivity
           wetc     = 622.0_kp*le*qss(i)/(rd*t1(i)*t1(i))
           alfac    = one / (one + (wetc*le*dwat)/(cp*dtmp))        ! wet bulb factor
           tem      = (1.0e3_kp * rain(i) / rho_w) * alfac * cp_w
           qrain(i) =  tem * (tsea-t1(i)+1.0e3_kp*(qss(i)-q0(i))*le/cp)

           !> - Calculate input non solar heat flux as upward = positive to models here

           f_nsol   = hflx(i) + evap(i) + ulwflx(i) - dlwflx(i) + omg_sh*qrain(i)

           sep      = sss*(evap(i)/le-rain(i))/rho_w
           ustar_a  = sqrt(stress(i)/rho_a(i))          ! air friction velocity
           !
           !  sensitivities of heat flux components to ts
           !
           rnl_ts = 4.0_kp*sfcemis(i)*sbc*tsea*tsea*tsea     ! d(rnl)/d(ts)
           hs_ts  = rch(i)
           hl_ts  = rch(i)*elocp*eps*hvap*qss(i)/(rd*t12)
           rf_ts  = tem * (one+rch(i)*hl_ts)
           q_ts   = rnl_ts + hs_ts + hl_ts + omg_sh*rf_ts
           !
           !> - Call cool_skin(), which is the sub-layer cooling parameterization
           !! (Fairfall et al. (1996) \cite fairall_et_al_1996).
           ! & calculate c_0, c_d
           !
           call cool_skin(ustar_a,f_nsol,nswsfc(i),evap(i),sss,alpha,beta, &
                rho_w,rho_a(i),tsea,q_ts,hl_ts,grav,le,                    &
                dt_cool(i),z_c(i),c_0(i),c_d(i))

           tem  = one / wndmag(i)
           cosa = u1(i)*tem
           sina = v1(i)*tem
           taux = max(stress(i),tau_min)*cosa
           tauy = max(stress(i),tau_min)*sina
           fc   = const_rot*sinlat(i)
           !
           !  Run DTM-1p system.
           !
           if ( (soltim > solar_time_6am .and. ifd(i) == zero) ) then
           else
              ifd(i) = one
              !
              !     calculate fcl thickness with current forcing and previous time's profile
              !

              !> - Call convdepth() to calculate depth for convective adjustments.
              if ( f_nsol > zero .and. xt(i) > zero ) then
                 call convdepth(kdt,timestep,nswsfc(i),f_nsol,sss,sep,rho_w, &
                      alpha,beta,xt(i),xs(i),xz(i),d_conv(i))
              else
                 d_conv(i) = zero
              endif
              !
              !    determine rich: wind speed dependent (right now)
              !
              !           if ( wind(i) < 1.0 ) then
              !             rich = 0.25 + 0.03*wind(i)
              !           elseif ( wind(i) >= 1.0 .and. wind(i) < 1.5 ) then
              !             rich = 0.25 + 0.1*wind(i)
              !           elseif ( wind(i) >= 1.5 .and. wind(i) < 6.0 ) then
              !             rich = 0.25 + 0.6*wind(i)
              !           elseif ( wind(i) >= 6.0 ) then
              !             rich = 0.25 + min(0.8*wind(i),0.50)
              !           endif

              rich = ri_c

              !> - Call the diurnal thermocline layer model dtm_1p().
              call dtm_1p(kdt,timestep,rich,taux,tauy,nswsfc(i),    &
                   f_nsol,sss,sep,q_ts,hl_ts,rho_w,alpha,beta,alon, &
                   sinlat(i),soltim,grav,le,d_conv(i),              &
                   xt(i),xs(i),xu(i),xv(i),xz(i),xzts(i),xtts(i))


              !  apply mda
              if ( xt(i) > zero ) then
                 !>  - If \a dtl heat content \a xt > 0.0, call dtm_1p_mda() to apply
                 !!  minimum depth adjustment (mda).
                 call dtm_1p_mda(xt(i),xtts(i),xz(i),xzts(i))
                 if ( xz(i) >= z_w_max ) then
                    !>   - If \a dtl thickness >= module_nst_parameters::z_w_max, call dtl_reset()
                    !! to reset xt/xs/x/xv to zero, and xz to module_nst_parameters::z_w_max.
                    call dtl_reset(xt(i),xs(i),xu(i),xv(i),xz(i),xtts(i), xzts(i))
                 endif

                 !  apply fca
                 if ( d_conv(i) > zero ) then
                    !>  - If thickness of free convection layer > 0.0, call dtm_1p_fca()
                    !! to apply free convection adjustment.
                    !>   - If \a dtl thickness >= module_nst_parameters::z_w_max(), call dtl_reset()
                    !! to reset xt/xs/x/xv to zero, and xz to module_nst_parameters::z_w_max().
                    call dtm_1p_fca(d_conv(i),xt(i),xtts(i),xz(i),xzts(i))
                    if ( xz(i) >= z_w_max ) then
                       call dtl_reset (xt(i),xs(i),xu(i),xv(i),xz(i),xzts(i),xtts(i))
                    endif
                 endif

                 !  apply tla
                 dz = min(xz(i),max(d_conv(i),delz))
                 !
                 !>  - Call sw_ps_9b() to compute the fraction of the solar radiation
                 !! absorbed by the depth \a delz (Paulson and Simpson (1981) \cite paulson_and_simpson_1981).
                 !! And calculate the total heat absorbed in warm layer.
                 call sw_ps_9b(delz,fw)
                 q_warm = fw*nswsfc(i)-f_nsol    !total heat absorbed in warm layer

                 !>  - Call cal_ttop() to calculate the diurnal warming amount at the top layer with
                 !! thickness of \a dz.
                 if ( q_warm > zero ) then
                    call cal_ttop(kdt,timestep,q_warm,rho_w,dz, xt(i),xz(i),ttop0)
                    ttop = ((xt(i)+xt(i))/xz(i))*(one-dz/((xz(i)+xz(i))))

                    !>  - Call dtm_1p_tla() to apply top layer adjustment.
                    if ( ttop > ttop0 ) then
                       call dtm_1p_tla(dz,ttop0,xt(i),xtts(i),xz(i),xzts(i))
                       if ( xz(i) >= z_w_max ) then
                          call dtl_reset (xt(i),xs(i),xu(i),xv(i),xz(i),xzts(i),xtts(i))
                       endif
                    endif
                 endif           ! if ( q_warm > 0.0 ) then

                 !  apply mwa
                 !>  - Call dt_1p_mwa() to apply maximum warming adjustment.
                 t0 = (xt(i)+xt(i))/xz(i)
                 if ( t0 > tw_max ) then
                    call dtm_1p_mwa(xt(i),xtts(i),xz(i),xzts(i))
                    if ( xz(i) >= z_w_max ) then
                       call dtl_reset (xt(i),xs(i),xu(i),xv(i),xz(i),xzts(i),xtts(i))
                    endif
                 endif

                 !  apply mta
                 !>  - Call dtm_1p_mta() to apply maximum temperature adjustment.
                 sstc = tref(i) + (xt(i)+xt(i))/xz(i) - dt_cool(i)

                 if ( sstc > sst_max ) then
                    dta = sstc - sst_max
                    call  dtm_1p_mta(dta,xt(i),xtts(i),xz(i),xzts(i))
                    if ( xz(i) >= z_w_max ) then
                       call dtl_reset (xt(i),xs(i),xu(i),xv(i),xz(i),xzts(i),xtts(i))
                    endif
                 endif
                 !
              endif             ! if ( xt(i) > 0.0 ) then
              !           reset dtl at midnight and when solar zenith angle > 89.994 degree
              if ( abs(soltim) < 2.0_kp*timestep ) then
                 call dtl_reset (xt(i),xs(i),xu(i),xv(i),xz(i),xzts(i),xtts(i))
              endif

           endif             ! if (solar_time > solar_time_6am .and. ifd(i) == 0.0 ) then: too late to start the first day


           !     update tsurf  (when flag(i) .eqv. .true. )
           !>  - Call get_dtzm_point() to computes \a dtz and \a tsurf.
           call get_dtzm_point(xt(i),xz(i),dt_cool(i),z_c(i), zsea1,zsea2,dtz)
           tsurf(i) = max(tgice, tref(i) + dtz )

           !>  - Call cal_w() to calculate \a w_0 and \a w_d.
           if ( xt(i) > zero ) then
              call cal_w(kdt,xz(i),xt(i),xzts(i),xtts(i),w_0(i),w_d(i))
           else
              w_0(i) = zero
              w_d(i) = zero
           endif

           !         if ( xt(i) > 0.0 ) then
           !           rig(i) = grav*xz(i)*xz(i)*(alpha*xt(i)-beta*xs(i))
           !    &             /(2.0*(xu(i)*xu(i)+xv(i)*xv(i)))
           !         else
           !           rig(i) = 0.25
           !         endif

           !         qrain(i) = rig(i)
           !zm(i) = wind(i)

        endif
     enddo

     ! restore nst-related prognostic fields for guess run
     do i=1, im
        if (mask(i) == 1) then
           if (flag_guess(i)) then    ! when it is guess of
              xt(i)      = xt_old(i)
              xs(i)      = xs_old(i)
              xu(i)      = xu_old(i)
              xv(i)      = xv_old(i)
              xz(i)      = xz_old(i)
              !zm(i)      = zm_old(i)
              xtts(i)    = xtts_old(i)
              xzts(i)    = xzts_old(i)
              ifd(i)     = ifd_old(i)
              tskin(i)   = tskin_old(i)
              dt_cool(i) = dt_cool_old(i)
              z_c(i)     = z_c_old(i)
           else
              !
              !         update tskin when coupled and not guess run
              !         (all other NSST variables have been updated in this case)
              !
              tskin(i) = tsurf(i)
           endif                 ! if flag_guess(i) then
        endif                   ! if mask==1
     enddo

  end subroutine sfc_nst_run
end module med_phases_ocnnst_mod
