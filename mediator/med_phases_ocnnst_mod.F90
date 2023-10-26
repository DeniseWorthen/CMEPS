module med_phases_ocnnst_mod

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_internalstate_mod , only : InternalState, logunit
  use med_constants_mod     , only : dbug_flag       => med_constants_dbug_flag
  use med_utils_mod         , only : chkerr          => med_utils_chkerr
  use med_methods_mod       , only : FB_diagnose     => med_methods_FB_diagnose
  use med_methods_mod       , only : FB_getFldPtr    => med_methods_FB_getFldPtr
  use med_methods_mod       , only : State_GetScalar => med_methods_State_GetScalar
  use med_internalstate_mod , only : mapconsf, mapnames, compatm, compocn, maintask
  use perf_mod              , only : t_startf, t_stopf
  use shr_orb_mod           , only : shr_orb_cosz, shr_orb_decl
  use shr_orb_mod           , only : shr_orb_params, SHR_ORB_UNDEF_INT, SHR_ORB_UNDEF_REAL
  use shr_log_mod           , only : shr_log_unit

  implicit none
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public med_phases_ocnnst_run

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private med_phases_ocnnst_init
  private med_phases_ocnnst_orbital_update
  private med_phases_ocnnst_orbital_init

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  type ocnnst_type
     real(r8) , pointer :: lats        (:) => null() ! latitudes  (degrees)
     real(r8) , pointer :: lons        (:) => null() ! longitudes (degrees)
     integer  , pointer :: mask        (:) => null() ! ocn domain mask: 0 <=> inactive cell
     !inputs
     real(r8) , pointer :: tsfco       (:) => null() ! sea surface temperature (K)
     real(r8) , pointer :: ps          (:) => null() ! surface pressure (Pa)
     real(r8) , pointer :: u1          (:) => null() ! zonal component of surface layer wind (m/s)
     real(r8) , pointer :: v1          (:) => null() ! merid component of surface layer wind (m/s)
     real(r8) , pointer :: t1          (:) => null() ! surface layer mean temperature (K)
     real(r8) , pointer :: q1          (:) => null() ! surface layer mean spec humidity (kg/kg)
     real(r8) , pointer :: dlwflx      (:) => null() ! total sky sfc downward lw flux (W/m2)
     real(r8) , pointer :: sfcnsw      (:) => null() ! total sfc netsw flx into ocean (W/m2)
     real(r8) , pointer :: rain        (:) => null() ! rainfall rate (kg/m2/s)
     real(r8) , pointer :: wind        (:) => null() ! wind speed (m/s)
     ! outputs
     real(r8) , pointer :: tseal       (:) => null() ! ocean surface skin temperature (K)
     real(r8) , pointer :: tsfc_water  (:) => null() ! surface skin temperature over water (K)
     real(r8) , pointer :: tsurf_water (:) => null() ! surface skin temperature after iteration over water (K)
     real(r8) , pointer :: dtzm        (:) => null() ! mean of dT(z)  (z1 to z2) (?)
     ! in sfcf file
     real(r8) , pointer :: c0          (:) => null() ! coefficient1 to calculate d(tz)/d(ts) (nd)
     real(r8) , pointer :: c2          (:) => null() ! coefficient2 to calculate d(tz)/d(ts) (nd)
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

     logical            :: created   ! has memory been allocated here
  end type ocnnst_type

  character(*),parameter :: u_FILE_u = &
       __FILE__
  character(len=CL)      :: orb_mode        ! attribute - orbital mode
  integer                :: orb_iyear       ! attribute - orbital year
  integer                :: orb_iyear_align ! attribute - associated with model year
  real(R8)               :: orb_obliq       ! attribute - obliquity in degrees
  real(R8)               :: orb_mvelp       ! attribute - moving vernal equinox longitude
  real(R8)               :: orb_eccen       ! attribute and update-  orbital eccentricity

  character(len=*) , parameter :: orb_fixed_year       = 'fixed_year'
  character(len=*) , parameter :: orb_variable_year    = 'variable_year'
  character(len=*) , parameter :: orb_fixed_parameters = 'fixed_parameters'
  ! used, reused in module
  logical  :: use_nextswcday  ! use the scalar field for next time (otherwise, will be set using clock)

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
    use ESMF  , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF  , only : ESMF_FieldBundleGet, ESMF_Field, ESMF_FieldGet
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

    !----------------------------------
    ! Set pointers to fields needed for NST calculations
    !----------------------------------

    ! These must must be on the ocean grid since the ocean NST computation is on the ocean grid
    ! The following sets pointers to the module arrays

    call FB_GetFldPtr(is_local%wrap%,FBMed_ocnnst_o, 'So_xt', ocnnst%xt, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%,FBMed_ocnnst_o, 'So_xz', ocnnst%xz, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%,FBMed_ocnnst_o, 'So_dtcool', ocnnst%dtcool, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%,FBMed_ocnnst_o, 'So_dtzm', ocnnst%dtzm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call FB_GetFldPtr(is_local%wrap%,FBMed_ocnnst_o, 'So_tref', ocnnst%tref, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%,FBMed_ocnnst_o, 'So_tseal', ocnnst%tseal, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%,FBMed_ocnnst_o, 'So_tsurf_water', ocnnst%tsurf_water, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !call FB_GetFldPtr(is_local%wrap%FBMed_ocnnst_o, 'So_t', ocnnst%tsfco, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compocn,compocn), 'So_t', ocnnst%tsfco, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compocn,compocn), 'So_omask', dataptr1d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    lsize = size(ocnnst%tsfco)
    do n = 1,lsize
       if (dataptr1d(n) == 0._r8) then
          ocnnst%mask(n) = 0
       else
          ocnnst%mask(n) = 1
       end if
    enddo

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

    ! Initialize orbital values
    call  med_phases_ocnnst_orbital_init(gcomp, logunit, iam==0, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Allow setting of NST timestep using the clock instead of the atm's next timestep
    use_nextswcday = .true.
    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxNextSwCday", isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (.not. isPresent ) then
       use_nextswcday = .false.
    endif
    write(msg,'(A,l)') trim(subname)//': use_nextswcday setting is ',use_nextswcday
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

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
    logical                 :: update_nst
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
    real(R8)                :: z_c_0
    real(R8)                :: tem2
    !real(R8), pointer       :: ofrac(:)
    !real(R8), pointer       :: ofrad(:)
    !real(R8), pointer       :: ifrac(:)
    !real(R8), pointer       :: ifrad(:)
    integer                 :: lsize            ! local size
    integer                 :: n                ! indices
    real(R8)                :: rlat             ! gridcell latitude in radians
    real(R8)                :: rlon             ! gridcell longitude in radians
    real(R8)                :: cosz             ! Cosine of solar zenith angle
    real(R8)                :: eccen            ! Earth orbit eccentricity
    real(R8)                :: mvelpp           ! Earth orbit
    real(R8)                :: lambm0           ! Earth orbit
    real(R8)                :: obliqr           ! Earth orbit
    real(R8)                :: delta            ! Solar declination angle  in radians
    real(R8)                :: eccf             ! Earth orbit eccentricity factor
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

       if (trim(runtype) == 'initial') then
          call ESMF_TimeGet( currTime, dayOfYear_r8=nextsw_cday, rc=rc )
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else
          ! obtain nextsw_cday from atm if it is in the import state
          if (use_nextswcday) then
             call State_GetScalar(&
                  state=is_local%wrap%NstateImp(compatm), &
                  flds_scalar_name=is_local%wrap%flds_scalar_name, &
                  flds_scalar_num=is_local%wrap%flds_scalar_num, &
                  scalar_id=is_local%wrap%flds_scalar_index_nextsw_cday, &
                  scalar_value=nextsw_cday, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
          else
             call ESMF_TimeGet( currTime, dayOfYear_r8=nextsw_cday, rc=rc )
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if
       first_call = .false.
    else
       ! Note that med_methods_State_GetScalar includes a broadcast to all other pets
       if (use_nextswcday) then
          call State_GetScalar(&
               state=is_local%wrap%NstateImp(compatm), &
               flds_scalar_name=is_local%wrap%flds_scalar_name, &
               flds_scalar_num=is_local%wrap%flds_scalar_num, &
               scalar_id=is_local%wrap%flds_scalar_index_nextsw_cday, &
               scalar_value=nextsw_cday, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else
          call ESMF_ClockGetNextTime(clock, nextTime, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_TimeGet(nextTime, dayOfYear_r8=nextsw_cday, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    ! Get orbital values
    call med_phases_ocnnst_orbital_update(clock, logunit, iam==0, eccen, obliqr, lambm0, mvelpp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Calculate ocean NST on the ocean grid
    update_nst = .false.
    lsize = size(ocnnst%tsfco)

    ! Clock is not advanced until the end of ModelAdvance
    call ESMF_TimeGet( currTime, hr_r8=solhr, rc=rc )
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Solar declination
    ! Will only do albedo calculation if nextsw_cday is not -1.
    write(msg,*)trim(subname)//' nextsw_cday = ',nextsw_cday
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    ! nst_pre
    z_c_0 = zero
    do n = 1,lsize
       if (ocnnst%mask(n) > 0) then
          call get_dtzm_point(ocnnst%xt(n), ocnnst%xz(n), ocnnst%dtcool(n), z_c_0, zero, omz1, ocnnst%dtm(n))
          ocnnst%tref(n) = max(tgice, ocnnst%tsfco(n) - ocnnst%dtmz(n))
          if (abs(ocnnst%xz(n)) > zero) then
             tem2 = one / ocnnst%xz(n)
          else
             tem2 = zero
          endif
          ocnnst%tseal(n)     = ocnnst%tref(n) + (ocnnst%xt(n)+ocnnst%xt(n)) * tem2 - ocnnst%dtcool(n)
          ocnnst%tsurf_wat(n) = ocnnst%tseal(n)
       endif
    enddo



    if (nextsw_cday >= -0.5_r8) then
       call shr_orb_decl(nextsw_cday, eccen, mvelpp,lambm0, obliqr, delta, eccf)

       ! Compute NST
       do n = 1,lsize
          if (ocnnst%mask(n) > 0.0) then
             rlat = const_deg2rad * ocnnst%lats(n)
             rlon = const_deg2rad * ocnnst%lons(n)
             cosz = shr_orb_cosz( nextsw_cday, rlat, rlon, delta )
             !if (cosz  >  0.0_r8) then !--- sun hit --
             !end if
          end if
       end do
       update_nst = .true.

    endif    ! nextsw_cday

    ! ! Update current ifrad/ofrad values if albedo was updated in field bundle
    ! if (update_nst) then
    !    call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compocn), fieldname='ifrac', field=lfield, rc=rc)
    !    if (chkerr(rc,__LINE__,u_FILE_u)) return
    !    call ESMF_FieldGet(lfield, farrayptr=ifrac, rc=rc)
    !    if (chkerr(rc,__LINE__,u_FILE_u)) return
    !    call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compocn), fieldname='ifrad', field=lfield, rc=rc)
    !    if (chkerr(rc,__LINE__,u_FILE_u)) return
    !    call ESMF_FieldGet(lfield, farrayptr=ifrad, rc=rc)
    !    if (chkerr(rc,__LINE__,u_FILE_u)) return
    !    call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compocn), fieldname='ofrac', field=lfield, rc=rc)
    !    if (chkerr(rc,__LINE__,u_FILE_u)) return
    !    call ESMF_FieldGet(lfield, farrayptr=ofrac, rc=rc)
    !    if (chkerr(rc,__LINE__,u_FILE_u)) return
    !    call ESMF_FieldBundleGet(is_local%wrap%FBFrac(compocn), fieldname='ofrad', field=lfield, rc=rc)
    !    if (chkerr(rc,__LINE__,u_FILE_u)) return
    !    call ESMF_FieldGet(lfield, farrayptr=ofrad, rc=rc)
    !    if (chkerr(rc,__LINE__,u_FILE_u)) return
    !    ifrad(:) = ifrac(:)
    !    ofrad(:) = ofrac(:)
    ! endif

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

  subroutine med_phases_ocnnst_orbital_init(gcomp, logunit, maintask, rc)

    !----------------------------------------------------------
    ! Obtain orbital related values
    !----------------------------------------------------------

    use ESMF  , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF  , only : ESMF_LogWrite, ESMF_LogFoundError, ESMF_LogSetError
    use ESMF  , only : ESMf_SUCCESS, ESMF_FAILURE, ESMF_LOGMSG_INFO, ESMF_RC_NOT_VALID
    use NUOPC , only : NUOPC_CompAttributeGet

    ! input/output variables
    type(ESMF_GridComp)                 :: gcomp
    integer             , intent(in)    :: logunit         ! output logunit
    logical             , intent(in)    :: maintask
    integer             , intent(out)   :: rc              ! output error

    ! local variables

    character(len=CL) :: msgstr          ! temporary
    character(len=CL) :: cvalue          ! temporary
    character(len=*) , parameter :: subname = "(med_phases_ocnnst_orbital_init)"
    !-------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine orbital attributes from input
    call NUOPC_CompAttributeGet(gcomp, name="orb_mode", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_mode

    call NUOPC_CompAttributeGet(gcomp, name="orb_iyear", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_iyear

    call NUOPC_CompAttributeGet(gcomp, name="orb_iyear_align", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_iyear_align

    call NUOPC_CompAttributeGet(gcomp, name="orb_obliq", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_obliq

    call NUOPC_CompAttributeGet(gcomp, name="orb_eccen", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_eccen

    call NUOPC_CompAttributeGet(gcomp, name="orb_mvelp", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_mvelp

    ! Error checks
    if (trim(orb_mode) == trim(orb_fixed_year)) then
       orb_obliq = SHR_ORB_UNDEF_REAL
       orb_eccen = SHR_ORB_UNDEF_REAL
       orb_mvelp = SHR_ORB_UNDEF_REAL
       if (orb_iyear == SHR_ORB_UNDEF_INT) then
          write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
          write(logunit,*) trim(subname),' ERROR: fixed_year settings = ',orb_iyear
          write (msgstr, *) ' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       endif
    elseif (trim(orb_mode) == trim(orb_variable_year)) then
       orb_obliq = SHR_ORB_UNDEF_REAL
       orb_eccen = SHR_ORB_UNDEF_REAL
       orb_mvelp = SHR_ORB_UNDEF_REAL
       if (orb_iyear == SHR_ORB_UNDEF_INT .or. orb_iyear_align == SHR_ORB_UNDEF_INT) then
          write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
          write(logunit,*) trim(subname),' ERROR: variable_year settings = ',orb_iyear, orb_iyear_align
          write (msgstr, *) subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       endif
    elseif (trim(orb_mode) == trim(orb_fixed_parameters)) then
       !-- force orb_iyear to undef to make sure shr_orb_params works properly
       orb_iyear = SHR_ORB_UNDEF_INT
       orb_iyear_align = SHR_ORB_UNDEF_INT
       if (orb_eccen == SHR_ORB_UNDEF_REAL .or. &
           orb_obliq == SHR_ORB_UNDEF_REAL .or. &
           orb_mvelp == SHR_ORB_UNDEF_REAL) then
          write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
          write(logunit,*) trim(subname),' ERROR: orb_eccen = ',orb_eccen
          write(logunit,*) trim(subname),' ERROR: orb_obliq = ',orb_obliq
          write(logunit,*) trim(subname),' ERROR: orb_mvelp = ',orb_mvelp
          write (msgstr, *) subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       endif
    else
       write (msgstr, *) subname//' ERROR: invalid orb_mode '//trim(orb_mode)
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       rc = ESMF_FAILURE
       return  ! bail out
    endif
  end subroutine med_phases_ocnnst_orbital_init

  !===============================================================================

  subroutine med_phases_ocnnst_orbital_update(clock, logunit,  maintask, eccen, obliqr, lambm0, mvelpp, rc)

    !----------------------------------------------------------
    ! Update orbital settings
    !----------------------------------------------------------

    use ESMF, only : ESMF_Clock, ESMF_ClockGet, ESMF_Time, ESMF_TimeGet
    use ESMF, only : ESMF_LogSetError, ESMF_RC_NOT_VALID, ESMF_SUCCESS

    ! input/output variables
    type(ESMF_Clock) , intent(in)    :: clock
    integer          , intent(in)    :: logunit
    logical          , intent(in)    :: maintask
    real(R8)         , intent(inout) :: eccen  ! orbital eccentricity
    real(R8)         , intent(inout) :: obliqr ! Earths obliquity in rad
    real(R8)         , intent(inout) :: lambm0 ! Mean long of perihelion at vernal equinox (radians)
    real(R8)         , intent(inout) :: mvelpp ! moving vernal equinox long of perihelion plus pi (rad)
    integer          , intent(out)   :: rc     ! output error

    ! local variables
    type(ESMF_Time)   :: CurrTime ! current time
    integer           :: year     ! model year at current time
    integer           :: orb_year ! orbital year for current orbital computation
    character(len=CL) :: msgstr   ! temporary
    logical           :: lprint
    logical           :: first_time = .true.
    character(len=*) , parameter :: subname = "(med_phases_ocnnst_orbital_update)"
    !-------------------------------------------

    rc = ESMF_SUCCESS
    lprint = .false.
    if (trim(orb_mode) == trim(orb_variable_year)) then
       call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(CurrTime, yy=year, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       orb_year = orb_iyear + (year - orb_iyear_align)
       lprint = maintask
    else
       orb_year = orb_iyear
       if (first_time) then
          lprint = maintask
          first_time = .false.
       else
          lprint = .false.
       end if
    end if

    eccen = orb_eccen
    shr_log_unit = logunit
    call shr_orb_params(orb_year, eccen, orb_obliq, orb_mvelp, obliqr, lambm0, mvelpp, lprint)

    if ( eccen  == SHR_ORB_UNDEF_REAL .or. obliqr == SHR_ORB_UNDEF_REAL .or. &
         mvelpp == SHR_ORB_UNDEF_REAL .or. lambm0 == SHR_ORB_UNDEF_REAL) then
       write (msgstr, *) subname//' ERROR: orb params incorrect'
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       return  ! bail out
    endif

  end subroutine med_phases_ocnnst_orbital_update

!===============================================================================

end module med_phases_ocnnst_mod
