module med_phases_ocnnst_mod

  use ESMF                  , only : ESMF_FieldBundle

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8, I4=>SHR_KIND_I4
  use med_internalstate_mod , only : InternalState, logunit
  use med_constants_mod     , only : dbug_flag       => med_constants_dbug_flag
  use med_utils_mod         , only : chkerr          => med_utils_chkerr
  use med_methods_mod       , only : FB_diagnose     => med_methods_FB_diagnose
  use med_methods_mod       , only : FB_getFldPtr    => med_methods_FB_getFldPtr
  use med_methods_mod       , only : State_GetScalar => med_methods_State_GetScalar
  use med_internalstate_mod , only : mapconsf, mapnames, compatm, compocn, maintask
  use perf_mod              , only : t_startf, t_stopf
  use module_nst_water_prop , only : get_dtzm_point

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
     real(r8) , pointer :: area        (:) => null() ! area (m2)
     integer  , pointer :: mask        (:) => null() ! ocn domain mask: 0 <=> inactive cell
     !inputs
     real(r8) , pointer :: ifrac       (:) => null() ! sea ice fraction (nd); ofrac=1.0-ifrac
     real(r8) , pointer :: tsfco       (:) => null() ! sea surface temperature (K)
     real(r8) , pointer :: ps          (:) => null() ! surface pressure (Pa)
     real(r8) , pointer :: u1          (:) => null() ! zonal component of surface layer wind (m/s)
     real(r8) , pointer :: v1          (:) => null() ! merid component of surface layer wind (m/s)
     real(r8) , pointer :: t1          (:) => null() ! surface layer mean temperature (K)
     real(r8) , pointer :: q1          (:) => null() ! surface layer mean spec humidity (kg/kg)
     real(r8) , pointer :: z1          (:) => null() ! layer 1 height above ground (not MSL) (m)
     real(r8) , pointer :: u10m        (:) => null() ! zonal component of 10m wind (m/s)
     real(r8) , pointer :: v10m        (:) => null() ! merid component of 10m wind (m/s)
     real(r8) , pointer :: dlwflx      (:) => null() ! total sky sfc downward lw flux (W/m2)
     real(r8) , pointer :: rain        (:) => null() ! rainfall rate (kg/m2/s)
     ! ?
     real(r8) , pointer :: tseal       (:) => null() ! ocean surface skin temperature (K)
     real(r8) , pointer :: tsfc_wat    (:) => null() ! surface skin temperature over water (K)
     real(r8) , pointer :: tsurf_wat   (:) => null() ! surface skin temperature after iteration over water (K)
     real(r8) , pointer :: dtzm        (:) => null() ! mean of dT(z)  (z1 to z2) (?)
     real(r8) , pointer :: dtm         (:) => null() ! ?
     ! in sfcf file
     real(r8) , pointer :: c_0         (:) => null() ! coefficient1 to calculate d(tz)/d(ts) (nd)
     real(r8) , pointer :: c_d         (:) => null() ! coefficient2 to calculate d(tz)/d(ts) (nd)
     real(r8) , pointer :: d_conv      (:) => null() ! thickness of free convection layer (m)
     real(r8) , pointer :: dt_cool     (:) => null() ! sub-layer cooling amount (K)
     real(r8) , pointer :: qrain       (:) => null() ! sensible heat flux due to rainfall (W)
     real(r8) , pointer :: tref        (:) => null() ! sea surface reference temperature (K)
     real(r8) , pointer :: w_0         (:) => null() ! coefficient3 to calculate d(tz)/d(ts) (nd)
     real(r8) , pointer :: w_d         (:) => null() ! coefficient4 to calculate d(tz)/d(ts) (nd)

     real(r8) , pointer :: xs          (:) => null() ! salinity  content in diurnal thermocline layer (ppt m)
     real(r8) , pointer :: xt          (:) => null() ! heat content in diurnal thermocline layer (K m)
     real(r8) , pointer :: xtts        (:) => null() ! d(xt)/d(ts) (m)

     real(r8) , pointer :: xu          (:) => null() ! u-current content in diurnal thermocline layer (m2 s-1)
     real(r8) , pointer :: xv          (:) => null() ! v-current  content in diurnal thermocline layer (m2 s-1)
     real(r8) , pointer :: xz          (:) => null() ! diurnal thermocline layer thickness (m)
     real(r8) , pointer :: xzts        (:) => null() ! d(xz)/d(ts) (m K-1)
     real(r8) , pointer :: zc          (:) => null() ! sub-layer cooling thickness (m)
     logical            :: created   ! has memory been allocated here
  end type ocnnst_type

  ! local
  type(ESMF_FieldBundle) :: FBnst
  real(r8) , pointer     :: wind        (:) => null() ! wind speed (m/s)
  real(r8) , pointer     :: sfcnsw      (:) => null() ! total sfc netsw flx into ocean (W/m2)
  ! logical  , pointer     :: flag_guess  (:) => null() ! .true.=  guess step to get CD et al
  !                                                     ! when iter = 1, flag_guess = .true. when wind < 2
  !                                                     ! when iter = 2, flag_guess = .false. for all grids
  ! logical  , pointer     :: flag_iter   (:) => null() ! execution or not
  !                                                     ! when iter = 1, flag_iter = .true. for all grids
  !                                                     ! when iter = 2, flag_iter = .true. when wind < 2
  !                                                     ! for both land and ocean (when nstf_name1 > 0)
  integer(i4) , pointer :: flag_guess (:) => null() ! 1 = true, 0 = false
  integer(i4) , pointer :: flag_iter  (:) => null() ! 1 = true, 0 = false

  character(*),parameter :: u_FILE_u =  __FILE__

  !--- potential temperature definition in surface layer physics
  logical              :: thsfc_loc      = .true.     !< flag for local vs. standard potential temperature

  ! used, reused in module
  !logical  :: use_nextswcday  ! use the scalar field for next time (otherwise, will be set using clock)

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
    use ESMF  , only : ESMF_FieldBundleCreate, ESMF_FieldBundleAdd, ESMF_FieldBundleGet
    use ESMF  , only : ESMF_FieldRegridGetArea
    use ESMF  , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate
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
    type(ESMF_Field)         :: lfield
    type(ESMF_Mesh)          :: lmesh
    integer                  :: n
    integer                  :: lsize
    integer                  :: spatialDim
    integer                  :: numOwnedElements
    type(InternalState)      :: is_local
    real(R8), pointer        :: ownedElemCoords(:)
    real(r8), pointer        :: dataptr1d(:)
    integer(i4), pointer     :: intptr1d(:)
    character(len=CL)        :: tempc1,tempc2
    character(len=CS)        :: cvalue
    logical                  :: isPresent, isSet
    integer                  :: fieldCount
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

    ! ocean grid ara
    allocate(ocnnst%area(numOwnedElements))
    call ESMF_FieldBundleGet(is_local%wrap%FBImp(compocn,compocn), fieldname='So_omask', field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridGetArea(lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=dataptr1d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ocnnst%area(:) = dataptr1d(:)

    ! Allow setting of NST timestep using the clock instead of the atm's next timestep
    ! use_nextswcday = .true.
    ! call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxNextSwCday", isPresent=isPresent, rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! if (.not. isPresent ) then
    !    use_nextswcday = .false.
    ! endif
    ! write(msg,'(A,l)') trim(subname)//': use_nextswcday setting is ',use_nextswcday
    ! call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    !----------------------------------
    ! local FB needed for NST calculation
    !----------------------------------

    FBnst = ESMF_FieldBundleCreate(name='FBnst', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lfield = ESMF_FieldCreate(lmesh, ESMF_TYPEKIND_R8, name='wind', meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBnst, (/lfield/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=dataptr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr1d = 0.0

    lfield = ESMF_FieldCreate(lmesh, ESMF_TYPEKIND_R8, name='sfcnsw', meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBnst, (/lfield/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=dataptr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr1d = 0.0

    lfield = ESMF_FieldCreate(lmesh, ESMF_TYPEKIND_I4, name='flag_iter', meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBnst, (/lfield/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=intptr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    intptr1d = 1
    where (ocnnst%mask == 0) intptr1d = 0

    lfield = ESMF_FieldCreate(lmesh, ESMF_TYPEKIND_I4, name='flag_guess', meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(FBnst, (/lfield/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=intptr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    intptr1d = 0

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
    use med_constants_mod , only : shr_const_pi
    use med_phases_history_mod, only : med_phases_history_write_med

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    ! local variables
    type(ocnnst_type), save :: ocnnst
    type(ESMF_VM)           :: vm
    type(ESMF_Field)        :: lfield
    integer                 :: iam
    integer                 :: iter
    logical                 :: update_nst
    type(InternalState)     :: is_local
    type(ESMF_Clock)        :: clock
    type(ESMF_Clock)        :: dclock
    type(ESMF_Time)         :: currTime
    type(ESMF_Time)         :: nextTime
    type(ESMF_TimeInterval) :: timeStep
    real(R8), pointer       :: swnet_vdr(:)
    real(R8), pointer       :: swnet_vdf(:)
    real(R8), pointer       :: swnet_idr(:)
    real(R8), pointer       :: swnet_idf(:)
    character(CL)           :: cvalue
    character(CS)           :: starttype        ! config start type
    character(CL)           :: runtype          ! initial, continue, hybrid, branch
    real(R8)                :: nextsw_cday      ! calendar day of next atm shortwave
    real(R8)                :: solhr            ! fcst hour at the end of prev time step (currTime)
    real(R8)                :: z_c_0, zsea1, zsea2
    real(R8)                :: tem2
    integer                 :: lsize            ! local size
    integer                 :: i                ! indices
    real(R8)                :: rlat             ! gridcell latitude in radians
    real(R8)                :: rlon             ! gridcell longitude in radians
    real(R8), parameter     :: const_deg2rad = shr_const_pi/180.0_R8  ! deg to rads
    real(R8), parameter     :: tgice = 271.20_R8  ! TODO? actually f(s)
    real(R8), parameter     :: zero = 0.0_R8, one = 1.0_R8, omz1 = 2.0_R8
    character(CL)           :: msg
    logical                 :: first_call = .true.
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

    ! Determine if ocnnst data type will be initialized - and if not return
    if (first_call) then
       if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnnst_o, rc=rc)) then
          ocnnst%created = .true.
       else
          ocnnst%created = .false.
       end if
    end if
    if (.not. ocnnst%created) then
       return
    end if

    ! Note that in the mct version the atm was initialized first so
    ! that nextsw_cday could be passed to the other components - this
    ! assumed that atmosphere component was ALWAYS initialized first.
    ! In the nuopc version it will be easier to assume that on startup
    ! - nextsw_cday is just what cam was setting it as the current calendar day

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    call t_startf('MED:'//subname)

    ! get clock
    call ESMF_GridCompGet(gcomp, clock=clock)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (first_call) then

       ! Initialize ocean NST calculation
       call med_phases_ocnnst_init(gcomp, ocnnst, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) starttype

       if (trim(starttype) == trim('startup')) then
          runtype = "initial"
       else if (trim(starttype) == trim('continue') ) then
          runtype = "continue"
       else if (trim(starttype) == trim('branch')) then
          runtype = "continue"
       else
          call ESMF_LogWrite( subname//' ERROR: unknown starttype', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if

       call ESMF_ClockGet( clock, currTime=currTime, timeStep=timeStep, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

    !    if (trim(runtype) == 'initial') then
    !       call ESMF_TimeGet( currTime, dayOfYear_r8=nextsw_cday, rc=rc )
    !       if (chkerr(rc,__LINE__,u_FILE_u)) return
    !    else
    !       ! obtain nextsw_cday from atm if it is in the import state
    !       if (use_nextswcday) then
    !          call State_GetScalar(&
    !               state=is_local%wrap%NstateImp(compatm), &
    !               flds_scalar_name=is_local%wrap%flds_scalar_name, &
    !               flds_scalar_num=is_local%wrap%flds_scalar_num, &
    !               scalar_id=is_local%wrap%flds_scalar_index_nextsw_cday, &
    !               scalar_value=nextsw_cday, rc=rc)
    !          if (chkerr(rc,__LINE__,u_FILE_u)) return
    !       else
    !          call ESMF_TimeGet( currTime, dayOfYear_r8=nextsw_cday, rc=rc )
    !          if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !       end if
    !    end if
    !    first_call = .false.
    ! else
    !    ! Note that med_methods_State_GetScalar includes a broadcast to all other pets
    !    if (use_nextswcday) then
    !       call State_GetScalar(&
    !            state=is_local%wrap%NstateImp(compatm), &
    !            flds_scalar_name=is_local%wrap%flds_scalar_name, &
    !            flds_scalar_num=is_local%wrap%flds_scalar_num, &
    !            scalar_id=is_local%wrap%flds_scalar_index_nextsw_cday, &
    !            scalar_value=nextsw_cday, rc=rc)
    !       if (chkerr(rc,__LINE__,u_FILE_u)) return
    !    else
    !       call ESMF_ClockGetNextTime(clock, nextTime, rc=rc)
    !       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !       call ESMF_TimeGet(nextTime, dayOfYear_r8=nextsw_cday, rc=rc)
    !       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !    end if
    ! end if

       if (trim(runtype) == 'initial') then
          call ESMF_TimeGet( currTime, dayOfYear_r8=nextsw_cday, rc=rc )
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          first_call = .false.
       else
          call ESMF_ClockGetNextTime(clock, nextTime, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_TimeGet(nextTime, dayOfYear_r8=nextsw_cday, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    ! Clock is not advanced until the end of ModelAdvance
    call ESMF_TimeGet( currTime, h_r8=solhr, rc=rc )
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !write(msg,*)trim(subname)//' nextsw_cday = ',nextsw_cday
    !call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
    !
    ! Calculate ocean NST on the ocean grid
    !
    update_nst = .false.
    lsize = size(ocnnst%mask)

    ! ice fraction on ocean
    call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ifrac', ocnnst%ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! ocean surface temperature from ocean
    call FB_GetFldPtr(is_local%wrap%FBImp(compocn,compocn), 'So_t', ocnnst%tsfco, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! ocean netsw from prep_ocn_custom (nst must run after ocn_accum, which calls prep_ocn_custom)
    call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdr',  swnet_vdr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdf',  swnet_vdf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_idr',  swnet_idr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_idf',  swnet_idf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! set pointers
    call set_ocnnst_pointers(is_local%wrap%FBImp(compatm,compocn), is_local%wrap%FBMed_ocnnst_o, &
         ocnnst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call FB_GetFldPtr(FBnst, 'wind', wind, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    wind = sqrt(ocnnst%u1*ocnnst%u1 + ocnnst%v1*ocnnst%v1)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call FB_GetFldPtr(FBnst, 'sfcnsw', sfcnsw, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    sfcnsw = swnet_vdr + swnet_vdf + swnet_idr + swnet_idf

    call ESMF_FieldBundleGet(FBnst, fieldname='flag_iter', field=lfield, rc=rc)
    call ESMF_FieldGet(lfield, farrayPtr=flag_iter, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleGet(FBnst, fieldname='flag_guess', field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=flag_guess, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! NST
    do iter = 1,2

       ! loop_control_part1
       do i = 1,lsize
          if (ocnnst%mask(i) == 1) then
             if (iter == 1 .and. wind(i) < 2.0d0) then
                flag_guess(i) = 1
             end if
          end if
       end do

       ! nst_pre
       z_c_0 = zero
       do i = 1,lsize
          if (ocnnst%mask(i) == 1) then
             call get_dtzm_point(ocnnst%xt(i), ocnnst%xz(i), ocnnst%dt_cool(i), z_c_0, zero, omz1, ocnnst%dtzm(i))
             ocnnst%tref(i) = max(tgice, ocnnst%tsfco(i) - ocnnst%dtzm(i))
             if (abs(ocnnst%xz(i)) > zero) then
                tem2 = one / ocnnst%xz(i)
             else
                tem2 = zero
             endif
             ocnnst%tseal(i)     = ocnnst%tref(i) + (ocnnst%xt(i)+ocnnst%xt(i)) * tem2 - ocnnst%dt_cool(i)
             ocnnst%tsurf_wat(i) = ocnnst%tseal(i)
          endif
       enddo

       ! nst_run
       ! Compute NST
       !do n = 1,lsize
       !   if (ocnnst%mask(n) == 1) then
       !      rlat = const_deg2rad * ocnnst%lats(n)
       !      rlon = const_deg2rad * ocnnst%lons(n)
       !   end if
       !end do
       !update_nst = .true.

       ! nst_post
       !zsea1 = 0.001_kp*real(nstf_name4)
       !zsea2 = 0.001_kp*real(nstf_name5)
       ! nstf_name4,5 are both 0 right now
       zsea1 = 0.0_r8
       zsea2 = 0.0_r8
       do i = 1,lsize
          if (ocnnst%mask(i) == 1) then
             call get_dtzm_point(ocnnst%xt(i), ocnnst%xz(i), ocnnst%dt_cool(i), ocnnst%zc(i), zsea1, zsea2, ocnnst%dtzm(i))
             ocnnst%tsfc_wat(i) = max(tgice, ocnnst%tref(i) + ocnnst%dtzm(i))
          end if
       end do

       ! loop_control_part2
       do i = 1, lsize
          if (ocnnst%mask(i) == 0) then
             flag_iter(i)  = 0
             flag_guess(i) = 0
          else
             if (iter == 1 .and. wind(i) < 2.0d0) then
                flag_iter(i) = 1
             endif
          endif
       end do
    end do

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

  end subroutine med_phases_ocnnst_run

!===============================================================================

!===============================================================================

  subroutine set_ocnnst_pointers(fldbun_a, fldbun_m, ocnnst, rc)

    use ESMF  , only : ESMF_FieldBundle, ESMF_SUCCESS

    ! Set pointers for ocnnst

    use med_methods_mod , only : FB_fldchk    => med_methods_FB_FldChk

    ! input/output variables
    type(ESMF_FieldBundle)     , intent(inout) :: fldbun_a
    type(ESMF_FieldBundle)     , intent(inout) :: fldbun_m
    type(ocnnst_type)          , intent(inout) :: ocnnst
    integer                    , intent(out)   :: rc
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !----------------------------------
    ! Set pointers to fields needed for NST calculations
    !----------------------------------

    ! Import from atm, 'inst, height_lowest fields'
    call FB_GetFldPtr(fldbun_a, 'Sa_pbot' ,           ocnnst%ps,         rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Sa_u' ,              ocnnst%u1,         rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Sa_v' ,              ocnnst%v1,         rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Sa_shum' ,           ocnnst%q1,         rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Sa_tbot' ,           ocnnst%t1,         rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Sa_z' ,              ocnnst%z1,         rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Sa_u10m' ,           ocnnst%u10m,       rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Sa_v10m' ,           ocnnst%v10m,       rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Faxa_rain' ,         ocnnst%rain,       rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_a, 'Faxa_lwdn' ,         ocnnst%dlwflx,     rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Ocean NST fields (mapped back to ATM)
    call FB_GetFldPtr(fldbun_m, 'Snst_tref' ,         ocnnst%tref,       rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_dconv' ,        ocnnst%d_conv,     rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_dtcool' ,       ocnnst%dt_cool,    rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_qrain' ,        ocnnst%qrain,      rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_xtts' ,         ocnnst%xtts,       rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_xzts' ,         ocnnst%xzts,       rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_c0' ,           ocnnst%c_0,        rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_cd' ,           ocnnst%c_d,        rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_w0' ,           ocnnst%w_0,        rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_wd' ,           ocnnst%w_d,        rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_xs' ,           ocnnst%xs,         rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_xt' ,           ocnnst%xt,         rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_xu' ,           ocnnst%xu,         rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_xv' ,           ocnnst%xv,         rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_xz' ,           ocnnst%xz,         rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_zc' ,           ocnnst%zc,         rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call FB_GetFldPtr(fldbun_m, 'Snst_tseal' ,       ocnnst%tseal,       rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_tsfc_water' ,  ocnnst%tsfc_wat,    rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_tsurf_water' , ocnnst%tsurf_wat,   rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_dtzm' ,        ocnnst%dtzm,        rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(fldbun_m, 'Snst_dtm' ,         ocnnst%dtm,         rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine set_ocnnst_pointers

  subroutine stability                                                &
 !  ---  inputs:
       ( z1, zvfun, gdx, tv1, thv1, wind, z0max, ztmax, tvs, grav,    &
       thsfc_loc,                                                     &
 !  ---  outputs:
       rb, fm, fh, fm10, fh2, cm, ch, stress, ustar)

    integer, parameter :: kind_phys = R8
    integer, parameter :: kp = kind_phys
    real (kind=kind_phys), parameter :: ca=0.4_kind_phys  ! ca - von karman constant

    !  ---  inputs:
    real(kind=kind_phys), intent(in) ::                               &
         z1, zvfun, gdx, tv1, thv1, wind, z0max, ztmax, tvs, grav
    logical,              intent(in) :: thsfc_loc

    !  ---  outputs:
    real(kind=kind_phys), intent(out) ::                              &
         rb, fm, fh, fm10, fh2, cm, ch, stress, ustar

    !  ---  locals:
    real(kind=kind_phys), parameter :: alpha=5.0_kp, a0=-3.975_kp,    &
         a1=12.32_kp, alpha4=4.0_kp*alpha,                            &
         b1=-7.755_kp, b2=6.041_kp,                                   &
         xkrefsqr=0.3_kp, xkmin=0.05_kp,                              &
         xkgdx=3000.0_kp,                                             &
         a0p=-7.941_kp, a1p=24.75_kp, b1p=-8.705_kp, b2p=7.899_kp,    &
         zolmin=-10.0_kp, zero=0.0_kp, one=1.0_kp

    real(kind=kind_phys) :: aa,     aa0,    bb,     bb0, dtv,   adtv, &
         hl1,    hl12,   pm,     ph,  pm10,  ph2,                     &
         z1i,                                                         &
         fms,    fhs,    hl0,    hl0inf, hlinf,                       &
         hl110,  hlt,    hltinf, olinf,                               &
         tem1,   tem2,   zolmax

    real(kind=kind_phys) :: xkzo

    z1i = one / z1

    !
    !  set background diffusivities with one for gdx >= xkgdx and
    !   as a function of horizontal grid size for gdx < xkgdx
    !   (i.e., gdx/xkgdx for gdx < xkgdx)
    !
    if(gdx >= xkgdx) then
       xkzo = one
    else
       xkzo = gdx / xkgdx
    endif

    tem1 = tv1 - tvs
    if(tem1 > zero) then
       tem2 = xkzo * zvfun
       xkzo = min(max(tem2, xkmin), xkzo)
    endif

    zolmax = xkrefsqr / sqrt(xkzo)

    !  compute stability indices (rb and hlinf)

    dtv     = thv1 - tvs
    adtv    = max(abs(dtv),0.001_kp)
    dtv     = sign(1.0_kp,dtv) * adtv

    if(thsfc_loc) then ! Use local potential temperature
       rb      = max(-5000.0_kp, (grav+grav) * dtv * z1 / ((thv1 + tvs) * wind * wind))
    else ! Use potential temperature referenced to 1000 hPa
       rb      = max(-5000.0_kp, grav * dtv * z1 / (tv1 * wind * wind))
    endif

    tem1    = one / z0max
    tem2    = one / ztmax
    fm      = log((z0max+z1)  * tem1)
    fh      = log((ztmax+z1)  * tem2)
    fm10    = log((z0max+10.0_kp) * tem1)
    fh2     = log((ztmax+2.0_kp)  * tem2)
    hlinf   = rb * fm * fm / fh
    hlinf   = min(max(hlinf,zolmin),zolmax)
    !
    !  stable case
    !
    if (dtv >= zero) then
       hl1 = hlinf
       if(hlinf > 0.25_kp) then
          tem1   = hlinf * z1i
          hl0inf = z0max * tem1
          hltinf = ztmax * tem1
          aa     = sqrt(one + alpha4 * hlinf)
          aa0    = sqrt(one + alpha4 * hl0inf)
          bb     = aa
          bb0    = sqrt(one + alpha4 * hltinf)
          pm     = aa0 - aa + log( (aa + one)/(aa0 + one) )
          ph     = bb0 - bb + log( (bb + one)/(bb0 + one) )
          fms    = fm - pm
          fhs    = fh - ph
          hl1    = fms * fms * rb / fhs
          hl1    = min(hl1, zolmax)
       endif
       !
       !  second iteration
       !
       tem1  = hl1 * z1i
       hl0   = z0max * tem1
       hlt   = ztmax * tem1
       aa    = sqrt(one + alpha4 * hl1)
       aa0   = sqrt(one + alpha4 * hl0)
       bb    = aa
       bb0   = sqrt(one + alpha4 * hlt)
       pm    = aa0 - aa + log( (one+aa)/(one+aa0) )
       ph    = bb0 - bb + log( (one+bb)/(one+bb0) )
       hl110 = hl1 * 10.0_kp * z1i
       aa    = sqrt(one + alpha4 * hl110)
       pm10  = aa0 - aa + log( (one+aa)/(one+aa0) )
       hl12  = (hl1+hl1) * z1i
       !           aa    = sqrt(one + alpha4 * hl12)
       bb    = sqrt(one + alpha4 * hl12)
       ph2   = bb0 - bb + log( (one+bb)/(one+bb0) )
       !
       !  unstable case - check for unphysical obukhov length
       !
    else                          ! dtv < 0 case
       olinf = z1 / hlinf
       tem1  = 50.0_kp * z0max
       if(abs(olinf) <= tem1) then
          hlinf = -z1 / tem1
          hlinf = max(hlinf, zolmin)
       endif
       !
       !  get pm and ph
       !
       if (hlinf >= -0.5_kp) then
          hl1   = hlinf
          pm    = (a0  + a1*hl1)  * hl1   / (one+ (b1+b2*hl1)  *hl1)
          ph    = (a0p + a1p*hl1) * hl1   / (one+ (b1p+b2p*hl1)*hl1)
          hl110 = hl1 * 10.0_kp * z1i
          pm10  = (a0 + a1*hl110) * hl110/(one+(b1+b2*hl110)*hl110)
          hl12  = (hl1+hl1) * z1i
          ph2   = (a0p + a1p*hl12) * hl12/(one+(b1p+b2p*hl12)*hl12)
       else                       ! hlinf < 0.05
          hl1   = -hlinf
          tem1  = one / sqrt(hl1)
          pm    = log(hl1) + 2.0_kp * sqrt(tem1) - .8776_kp
          ph    = log(hl1) + 0.5_kp * tem1 + 1.386_kp
          !             pm    = log(hl1) + 2.0 * hl1 ** (-.25) - .8776
          !             ph    = log(hl1) + 0.5 * hl1 ** (-.5) + 1.386
          hl110 = hl1 * 10.0_kp * z1i
          pm10  = log(hl110) + 2.0_kp/sqrt(sqrt(hl110)) - 0.8776_kp
          !             pm10  = log(hl110) + 2. * hl110 ** (-.25) - .8776
          hl12  = (hl1+hl1) * z1i
          ph2   = log(hl12) + 0.5_kp / sqrt(hl12) + 1.386_kp
          !             ph2   = log(hl12) + .5 * hl12 ** (-.5) + 1.386
       endif

    endif          ! end of if (dtv >= 0 ) then loop
    !
    !  finish the exchange coefficient computation to provide fm and fh
    !
    fm        = fm - pm
    fh        = fh - ph
    fm10      = fm10 - pm10
    fh2       = fh2 - ph2
    cm        = ca * ca / (fm * fm)
    ch        = ca * ca / (fm * fh)
    tem1      = 0.00001_kp/z1
    cm        = max(cm, tem1)
    ch        = max(ch, tem1)
    stress    = cm * wind * wind
    ustar     = sqrt(stress)

  end subroutine stability

end module med_phases_ocnnst_mod
