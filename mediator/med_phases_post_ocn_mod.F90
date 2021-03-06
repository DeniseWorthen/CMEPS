module med_phases_post_ocn_mod

  !-----------------------------------------------------------------------------
  ! Mediator post ocn phase - maps ocn->ice, accumulate glc input from ocn
  !-----------------------------------------------------------------------------

  implicit none
  private

  public  :: med_phases_post_ocn

  logical :: ocn2glc_coupling

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_post_ocn(gcomp, rc)

    use med_kind_mod            , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
    use ESMF                    , only : ESMF_GridComp
    use ESMF                    , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use med_utils_mod           , only : chkerr      => med_utils_ChkErr
    use med_constants_mod       , only : dbug_flag   => med_constants_dbug_flag
    use med_map_mod             , only : med_map_field_packed
    use med_internalstate_mod   , only : InternalState, logunit, mastertask
    use med_phases_prep_glc_mod , only : med_phases_prep_glc_accum_ocn
    use esmFlds                 , only : compice, compglc, compocn, num_icesheets
    use perf_mod                , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: ns
    logical             :: first_call = .true.
    character(len=*),parameter :: subname='(med_phases_post_ocn)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Map ocn->ice
    if (is_local%wrap%med_coupling_active(compocn,compice)) then
       call t_startf('MED:'//trim(subname)//' map_ocn2ice')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compocn,compocn), &
            FBDst=is_local%wrap%FBImp(compocn,compice), &
            FBFracSrc=is_local%wrap%FBFrac(compocn), &
            field_normOne=is_local%wrap%field_normOne(compocn,compice,:), &
            packed_data=is_local%wrap%packed_data(compocn,compice,:), &
            routehandles=is_local%wrap%RH(compocn,compice,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_ocn2ice')
    end if

    ! Accumulate ocn input for glc if there is ocn->glc coupling
    if (first_call) then
       do ns = 1,num_icesheets
          if (is_local%wrap%med_coupling_active(compocn,compglc(ns))) then
             ocn2glc_coupling = .true.
             exit
          end if
       end do
       first_call = .false.
    end if
    if (ocn2glc_coupling) then
       call med_phases_prep_glc_accum_ocn(gcomp, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_post_ocn

end module med_phases_post_ocn_mod
