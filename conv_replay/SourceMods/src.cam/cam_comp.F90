module cam_comp
!-----------------------------------------------------------------------
!
! Purpose:  The CAM Community Atmosphere Model component. Interfaces with 
!           a merged surface field that is provided outside of this module. 
!           This is the atmosphere only component. It can interface with a 
!           host of surface components.
!
!-----------------------------------------------------------------------
   use shr_kind_mod,      only: r8 => SHR_KIND_R8, cl=>SHR_KIND_CL, cs=>SHR_KIND_CS
   use pmgrid,            only: plat, plev
   use spmd_utils,        only: masterproc
   use abortutils,        only: endrun
   use camsrfexch,        only: cam_out_t, cam_in_t     
   use shr_sys_mod,       only: shr_sys_flush
   use infnan,            only: nan
   use physics_types,     only: physics_state, physics_tend
   use cam_control_mod,   only: nsrest, print_step_cost, obliqr, lambm0, mvelpp, eccen
   use dyn_comp,          only: dyn_import_t, dyn_export_t
   use ppgrid,            only: begchunk, endchunk, pver, pcols
   use perf_mod
   use cam_logfile,       only: iulog
     use physics_buffer, only: physics_buffer_desc, pbuf_read_restart, pbuf_init_restart, pbuf_deallocate, pbuf_get_index, pbuf_get_field, pbuf_get_chunk, pbuf_old_tim_idx
   implicit none
   private
   save
   !
   ! Public access methods
   !
   public cam_init      ! First phase of CAM initialization
   public cam_run1      ! CAM run method phase 1
   public cam_run2      ! CAM run method phase 2
   public cam_run3      ! CAM run method phase 3
   public cam_run4      ! CAM run method phase 4
   public cam_final     ! CAM Finalization
   !
   ! Private module data
   !
#if ( defined SPMD )
   real(r8) :: cam_time_beg              ! Cam init begin timestamp
   real(r8) :: cam_time_end              ! Cam finalize end timestamp
   real(r8) :: stepon_time_beg = -1.0_r8 ! Stepon (run1) begin timestamp
   real(r8) :: stepon_time_end = -1.0_r8 ! Stepon (run4) end timestamp
   integer  :: nstep_beg = -1            ! nstep at beginning of run
#else
   integer  :: mpicom = 0
#endif

  real(r8) :: gw(plat)           ! Gaussian weights
  real(r8) :: etamid(plev) = nan ! vertical coords at midpoints
  real(r8) :: dtime              ! Time step for either physics or dynamics (set in dynamics init)

  type(dyn_import_t) :: dyn_in   ! Dynamics import container
  type(dyn_export_t) :: dyn_out  ! Dynamics export container

  type(physics_state), pointer :: phys_state(:)
  type(physics_tend ), pointer :: phys_tend(:)
  type(physics_buffer_desc), pointer :: pbuf2d(:,:)
  type(physics_buffer_desc), pointer :: pbuf(:)

  real(r8) :: wcstart, wcend     ! wallclock timestamp at start, end of timestep
  real(r8) :: usrstart, usrend   ! user timestamp at start, end of timestep
  real(r8) :: sysstart, sysend   ! sys timestamp at start, end of timestep
  integer :: nstep

!-----------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
!
subroutine cam_init( cam_out, cam_in, mpicom_atm, &
	             start_ymd, start_tod, ref_ymd, ref_tod, stop_ymd, stop_tod, &
                    perpetual_run, perpetual_ymd, calendar)

   !-----------------------------------------------------------------------
   !
   ! Purpose:  CAM initialization.
   !
   !-----------------------------------------------------------------------

   use history_defaults, only: bldfld
   use cam_initfiles,    only: cam_initfiles_open
   use inital,           only: cam_initial
   use cam_restart,      only: cam_read_restart
   use stepon,           only: stepon_init
   use physpkg,          only: phys_init, phys_register
   
   use dycore,           only: dycore_is
#if (defined BFB_CAM_SCAM_IOP)
   use history_defaults, only: initialize_iop_history
#endif
!   use shr_orb_mod,      only: shr_orb_params
   use camsrfexch,       only: hub2atm_alloc, atm2hub_alloc
   use cam_history,      only: addfld, add_default, phys_decomp, intht, init_masterlinkedlist
   use history_scam,     only: scm_intht
   use scamMod,          only: single_column
   use cam_pio_utils,    only: init_pio_subsystem
   use cam_instance,     only: inst_suffix

#if ( defined SPMD )   
   real(r8) :: mpi_wtime  ! External
#endif
   !-----------------------------------------------------------------------
   !
   ! Arguments
   !
   type(cam_out_t), pointer        :: cam_out(:)       ! Output from CAM to surface
   type(cam_in_t) , pointer        :: cam_in(:)        ! Merged input state to CAM
!   type(cam_in_t)        :: cam_in(:)        ! Merged input state to
   integer            , intent(in) :: mpicom_atm       ! CAM MPI communicator
   integer            , intent(in) :: start_ymd        ! Start date (YYYYMMDD)
   integer            , intent(in) :: start_tod        ! Start time of day (sec)
   integer            , intent(in) :: ref_ymd          ! Reference date (YYYYMMDD)
   integer            , intent(in) :: ref_tod          ! Reference time of day (sec)
   integer            , intent(in) :: stop_ymd         ! Stop date (YYYYMMDD)
   integer            , intent(in) :: stop_tod         ! Stop time of day (sec)
   logical            , intent(in) :: perpetual_run    ! If in perpetual mode or not
   integer            , intent(in) :: perpetual_ymd    ! Perpetual date (YYYYMMDD)
   character(len=cs)  , intent(in) :: calendar         ! Calendar type
   !
   ! Local variables
   !
   integer :: dtime_cam        ! Time-step
   logical :: log_print        ! Flag to print out log information or not
   character(len=cs) :: filein ! Input namelist filename
   !-----------------------------------------------------------------------
   !
   ! Initialize CAM MPI communicator, number of processors, master processors 
   !	
#if ( defined SPMD )
   cam_time_beg = mpi_wtime()
#endif
   !
   ! Initialization needed for cam_history
   ! 
   call init_masterlinkedlist()
   !
   ! Set up spectral arrays
   !
   call trunc()
   !
   ! Initialize index values for advected and non-advected tracers
   !
   call phys_register ()
   !
   ! Determine input namelist filename
   !
   filein = "atm_in" // trim(inst_suffix)
   !
   ! Do appropriate dynamics and history initialization depending on whether initial, restart, or 
   ! branch.  On restart run intht need not be called because all the info is on restart dataset.
   !
   call init_pio_subsystem(filein)

   if ( nsrest == 0 )then

      call cam_initfiles_open()
      call cam_initial(dyn_in, dyn_out, NLFileName=filein)

      ! Allocate and setup surface exchange data
      call atm2hub_alloc(cam_out)
      call hub2atm_alloc(cam_in)

   else

      call cam_read_restart ( cam_out, dyn_in, dyn_out, pbuf2d, stop_ymd, stop_tod, NLFileName=filein )

      call hub2atm_alloc( cam_in )
#if (defined BFB_CAM_SCAM_IOP)
      call initialize_iop_history()
#endif
   end if


   call phys_init( phys_state, phys_tend, pbuf2d,  cam_out )

   call bldfld ()       ! master field list (if branch, only does hash tables)

   !
   ! Setup the characteristics of the orbit
   ! (Based on the namelist parameters)
   !
   if (masterproc) then
      log_print = .true.
   else
      log_print = .false.
   end if

   call stepon_init( gw, etamid, dyn_in, dyn_out ) ! dyn_out necessary?

   if (single_column) call scm_intht()
   call intht()


end subroutine cam_init

!
!-----------------------------------------------------------------------
!
subroutine cam_run1(cam_in, cam_out)
!-----------------------------------------------------------------------
!
! Purpose:   First phase of atmosphere model run method.
!            Runs first phase of dynamics and first phase of
!            physics (before surface model updates).
!
!-----------------------------------------------------------------------
   
   use physpkg,          only: phys_run1
   use stepon,           only: stepon_run1
#if ( defined SPMD )
   use mpishorthand,     only: mpicom
#endif
   use time_manager,     only: get_nstep

   type(cam_in_t)  :: cam_in(begchunk:endchunk)
   type(cam_out_t) :: cam_out(begchunk:endchunk)

#if ( defined SPMD )
   real(r8) :: mpi_wtime
#endif
!-----------------------------------------------------------------------

#if ( defined SPMD )
   if (stepon_time_beg == -1.0_r8) stepon_time_beg = mpi_wtime()
   if (nstep_beg == -1) nstep_beg = get_nstep()
#endif
   if (masterproc .and. print_step_cost) then
      call t_stampf (wcstart, usrstart, sysstart)
   end if
   !----------------------------------------------------------
   ! First phase of dynamics (at least couple from dynamics to physics)
   ! Return time-step for physics from dynamics.
   !----------------------------------------------------------
   call t_barrierf ('sync_stepon_run1', mpicom)
   call t_startf ('stepon_run1')
   call stepon_run1( dtime, phys_state, phys_tend, pbuf2d, dyn_in, dyn_out )
   call t_stopf  ('stepon_run1')

   !
   !----------------------------------------------------------
   ! PHYS_RUN Call the Physics package
   !----------------------------------------------------------
   !
   call t_barrierf ('sync_phys_run1', mpicom)
   call t_startf ('phys_run1')
   call phys_run1(phys_state, dtime, phys_tend, pbuf2d,  cam_in, cam_out)
   call t_stopf  ('phys_run1')

end subroutine cam_run1

!
!-----------------------------------------------------------------------
!

subroutine cam_run2( cam_out, cam_in )
!-----------------------------------------------------------------------
!
! Purpose:   Second phase of atmosphere model run method.
!            Run the second phase physics, run methods that
!            require the surface model updates.  And run the
!            second phase of dynamics that at least couples
!            between physics to dynamics.
!
!-----------------------------------------------------------------------
   
   use physpkg,          only: phys_run2
   use stepon,           only: stepon_run2
   use time_manager,     only: is_first_step, is_first_restart_step
#if ( defined SPMD )
   use mpishorthand,     only: mpicom
#endif

   type(cam_out_t), intent(inout) :: cam_out(begchunk:endchunk)
   type(cam_in_t),  intent(inout) :: cam_in(begchunk:endchunk)

   !
   ! Second phase of physics (after surface model update)
   !
   call t_barrierf ('sync_phys_run2', mpicom)
   call t_startf ('phys_run2')
   call phys_run2(phys_state, dtime, phys_tend, pbuf2d,  cam_out, cam_in )
   call t_stopf  ('phys_run2')

   !
   ! Second phase of dynamics (at least couple from physics to dynamics)
   !
   call t_barrierf ('sync_stepon_run2', mpicom)
   call t_startf ('stepon_run2')
   call stepon_run2( phys_state, phys_tend, dyn_in, dyn_out )

   call t_stopf  ('stepon_run2')

   if (is_first_step() .or. is_first_restart_step()) then
      call t_startf ('cam_run2_memusage')
      call t_stopf  ('cam_run2_memusage')
   end if
end subroutine cam_run2

!
!-----------------------------------------------------------------------
!

subroutine cam_run3( cam_out )
!-----------------------------------------------------------------------
!
! Purpose:  Third phase of atmosphere model run method. This consists
!           of the third phase of the dynamics. For some dycores
!           this will be the actual dynamics run, for others the
!           dynamics happens before physics in phase 1.
!
!-----------------------------------------------------------------------
   use stepon,           only: stepon_run3
   use time_manager,     only: is_first_step, is_first_restart_step
#if ( defined SPMD )
   use mpishorthand,     only: mpicom
#endif

   type(cam_out_t), intent(inout) :: cam_out(begchunk:endchunk)
!-----------------------------------------------------------------------
   !
   ! Third phase of dynamics
   !
   call t_barrierf ('sync_stepon_run3', mpicom)
   call t_startf ('stepon_run3')
   call stepon_run3( dtime, etamid, cam_out, phys_state, dyn_in, dyn_out )

   call t_stopf  ('stepon_run3')

   if (is_first_step() .or. is_first_restart_step()) then
      call t_startf ('cam_run3_memusage')
      call t_stopf  ('cam_run3_memusage')
   end if
end subroutine cam_run3

!
!-----------------------------------------------------------------------
!

subroutine cam_run4( cam_out, cam_in, rstwr, nlend, &
                     yr_spec, mon_spec, day_spec, sec_spec )

!-----------------------------------------------------------------------
!
! Purpose:  Final phase of atmosphere model run method. This consists
!           of all the restart output, history writes, and other
!           file output.
!
!-----------------------------------------------------------------------
   use cam_history,      only: wshist, wrapup
   use cam_restart,      only: cam_write_restart, cam_read_restart
   use dycore,           only: dycore_is
       use constituents,    only: pcnst
 use phys_grid        , only: get_ncols_p
use pio              , only: file_desc_t, pio_closefile
  use filenames        , only: interpret_filename_spec
      use cam_pio_utils,    only: cam_pio_openfile
    use restart_dynamics, only: read_restart_dynamics
    use restart_physics, only: read_restart_physics
use camsrfexch,       only: hub2atm_alloc,hub2atm_deallocate


#if ( defined SPMD )
   use mpishorthand,     only: mpicom
#endif

!   type(cam_out_t), intent(inout)        :: cam_out(begchunk:endchunk)
   type(cam_out_t), pointer        :: cam_out(:)
   type(cam_in_t), pointer          :: cam_in(:) 
!   type(cam_in_t) , intent(inout)        :: cam_in(begchunk:endchunk)
   logical            , intent(in)           :: rstwr           ! true => write restart file
   logical            , intent(in)           :: nlend           ! true => this is final timestep
   integer            , intent(in), optional :: yr_spec         ! Simulation year
   integer            , intent(in), optional :: mon_spec        ! Simulation month
   integer            , intent(in), optional :: day_spec        ! Simulation day
   integer            , intent(in), optional :: sec_spec        ! Seconds into current simulation day
   logical            , save                 :: do_restart=.FALSE.
   logical            , save                 :: first_time=.true.
    integer :: c, p, q, ncol, lchnk, i, k
    character(len=50) :: locfn
    character(len=50) :: filein 
   type(file_desc_t) :: File
  character(len=50) :: fname_pbuf_cam    
   character(len=50) :: rfilename_spec_cam = '%c.cam.r.%y-%m-%d-%s.nc'
   logical, save :: is_first=.TRUE.
   integer :: itim

real(r8), pointer, dimension(:,:) :: t_ttend   
real(r8), pointer, dimension(:,:) :: TEOUT     
real(r8), pointer, dimension(:,:) :: DTCORE    
real(r8), pointer, dimension(:,:) :: QCWAT     
real(r8), pointer, dimension(:,:) :: LCWAT     
real(r8), pointer, dimension(:,:) :: TCWAT     
real(r8), pointer, dimension(:,:) :: CLD       
real(r8), pointer, dimension(:,:) :: AST       
real(r8), pointer, dimension(:,:) :: CONCLD    
real(r8), pointer, dimension(:,:) :: DP_FLXPRC 
real(r8), pointer, dimension(:,:) :: DP_FLXSNW 
real(r8), pointer, dimension(:,:) :: DP_CLDLIQ 
real(r8), pointer, dimension(:,:) :: DP_CLDICE 
real(r8), pointer, dimension(:) :: cush      
real(r8), pointer, dimension(:,:) :: QRS       
real(r8), pointer, dimension(:,:) :: QRL       
real(r8), pointer, dimension(:,:) :: pblh      
real(r8), pointer, dimension(:,:) :: tke       
real(r8), pointer, dimension(:,:) :: kvh       
real(r8), pointer, dimension(:,:) :: kvm       
real(r8), pointer, dimension(:,:) :: Turbtype  
real(r8), pointer, dimension(:,:) :: smaw      
real(r8), pointer, dimension(:) :: Tauresx   
real(r8), pointer, dimension(:) :: Tauresy   
real(r8), pointer, dimension(:) :: tpert     
real(r8), pointer, dimension(:,:) :: qpert


real(r8), pointer, dimension(:,:) :: t_ttend_old  
real(r8), pointer, dimension(:,:) :: TEOUT_old    
real(r8), pointer, dimension(:,:) :: DTCORE_old   
real(r8), pointer, dimension(:,:) :: QCWAT_old    
real(r8), pointer, dimension(:,:) :: LCWAT_old    
real(r8), pointer, dimension(:,:) :: TCWAT_old    
real(r8), pointer, dimension(:,:) :: CLD_old      
real(r8), pointer, dimension(:,:) :: AST_old      
real(r8), pointer, dimension(:,:) :: CONCLD_old   
real(r8), pointer, dimension(:,:) :: DP_FLXPRC_old
real(r8), pointer, dimension(:,:) :: DP_FLXSNW_old
real(r8), pointer, dimension(:,:) :: DP_CLDLIQ_old
real(r8), pointer, dimension(:,:) :: DP_CLDICE_old
real(r8), pointer, dimension(:) :: cush_old     
real(r8), pointer, dimension(:,:) :: QRS_old      
real(r8), pointer, dimension(:,:) :: QRL_old      
real(r8), pointer, dimension(:,:) :: pblh_old     
real(r8), pointer, dimension(:,:) :: tke_old      
real(r8), pointer, dimension(:,:) :: kvh_old      
real(r8), pointer, dimension(:,:) :: kvm_old      
real(r8), pointer, dimension(:,:) :: Turbtype_old 
real(r8), pointer, dimension(:,:) :: smaw_old     
real(r8), pointer, dimension(:) :: Tauresx_old  
real(r8), pointer, dimension(:) :: Tauresy_old  
real(r8), pointer, dimension(:) :: tpert_old     
real(r8), pointer, dimension(:,:) :: qpert_old


    integer :: t_ttend_idx       = 0
    integer ::  TEOUT_idx      = 0 
    integer ::  DTCORE_idx     = 0
    integer ::  QCWAT_idx     = 0
    integer ::  LCWAT_idx     = 0
    integer ::  TCWAT_idx     = 0
    integer ::  CLD_idx      = 0
    integer ::  AST_idx      = 0
    integer ::  CONCLD_idx     = 0
    integer ::  DP_FLXPRC_idx  = 0 
    integer ::  DP_FLXSNW_idx     = 0
    integer ::  DP_CLDLIQ_idx     = 0
    integer ::  DP_CLDICE_idx     = 0
    integer ::  cush_idx     = 0
    integer ::  QRS_idx     =  0
    integer ::  QRL_idx      =  0
    integer ::  pblh_idx     = 0
    integer ::  tke_idx      = 0
    integer ::  kvh_idx      = 0
    integer ::  kvm_idx       = 0
    integer ::  Turbtype_idx     = 0
    integer ::  smaw_idx     = 0
    integer ::  Tauresx_idx     = 0
    integer ::  Tauresy_idx     = 0
    integer ::  tpert_idx     = 0
    integer ::  qpert_idx     = 0


  integer :: t_ttend_oldid       = 0
    integer ::  TEOUT_oldid      = 0 
    integer ::  DTCORE_oldid     = 0
    integer ::  QCWAT_oldid     = 0
    integer ::  LCWAT_oldid     = 0
    integer ::  TCWAT_oldid     = 0
    integer ::  CLD_oldid      = 0
    integer ::  AST_oldid      = 0
    integer ::  CONCLD_oldid     = 0
    integer ::  DP_FLXPRC_oldid  = 0 
    integer ::  DP_FLXSNW_oldid     = 0
    integer ::  DP_CLDLIQ_oldid     = 0
    integer ::  DP_CLDICE_oldid     = 0
    integer ::  cush_oldid     = 0
    integer ::  QRS_oldid     =  0
    integer ::  QRL_oldid      =  0
    integer ::  pblh_oldid     = 0
    integer ::  tke_oldid      = 0
    integer ::  kvh_oldid      = 0
    integer ::  kvm_oldid       = 0
    integer ::  Turbtype_oldid     = 0
    integer ::  smaw_oldid     = 0
    integer ::  Tauresx_oldid     = 0
    integer ::  Tauresy_oldid     = 0
    integer ::  tpert_oldid     = 0
    integer ::  qpert_oldid     = 0


#if ( defined SPMD )
   real(r8) :: mpi_wtime
#endif
!-----------------------------------------------------------------------
! print_step_cost

   !
   !----------------------------------------------------------
   ! History and restart logic: Write and/or dispose history tapes if required
   !----------------------------------------------------------
   !
   call t_barrierf ('sync_wshist', mpicom)
   call t_startf ('wshist')
   call wshist ()
   call t_stopf  ('wshist')

#if ( defined SPMD )
   stepon_time_end = mpi_wtime()
#endif
   !
   ! Write restart files
   !
   if (rstwr) then
      call t_startf ('cam_write_restart')
      if (present(yr_spec).and.present(mon_spec).and.present(day_spec).and.present(sec_spec)) then
         call cam_write_restart( cam_out, dyn_out, pbuf2d, &
              yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
      else
         call cam_write_restart( cam_out, dyn_out, pbuf2d )
      end if
      call t_stopf  ('cam_write_restart') 
   end if

!when it's time to do restart, find the right file name and read it, zeroing out
!some more tendencies

if ( mod(sec_spec,21600)==10800 .AND. do_restart ) then


        do_restart=.FALSE.

        filein = "atm_in"
        fname_pbuf_cam = interpret_filename_spec( rfilename_spec_cam, yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec=sec_spec-10800 )
        locfn= './' // trim(fname_pbuf_cam)
do lchnk=begchunk,endchunk
    pbuf=> pbuf_get_chunk(pbuf2d,lchnk)

t_ttend_idx   = pbuf_get_index('T_TTEND') 
TEOUT_idx     = pbuf_get_index('TEOUT') 
DTCORE_idx    = pbuf_get_index('DTCORE') 
QCWAT_idx     = pbuf_get_index('QCWAT') 
LCWAT_idx     = pbuf_get_index('LCWAT') 
TCWAT_idx     = pbuf_get_index('TCWAT') 
CLD_idx       = pbuf_get_index('CLD') 
AST_idx       = pbuf_get_index('AST') 
CONCLD_idx    = pbuf_get_index('CONCLD') 
DP_FLXPRC_idx = pbuf_get_index('DP_FLXPRC') 
DP_FLXSNW_idx = pbuf_get_index('DP_FLXSNW') 
DP_CLDLIQ_idx = pbuf_get_index('DP_CLDLIQ') 
DP_CLDICE_idx = pbuf_get_index('DP_CLDICE') 
cush_idx      = pbuf_get_index('cush') 
QRS_idx       = pbuf_get_index('QRS') 
QRL_idx       = pbuf_get_index('QRL') 
pblh_idx      = pbuf_get_index('pblh') 
tke_idx       = pbuf_get_index('tke') 
kvh_idx       = pbuf_get_index('kvh') 
kvm_idx       = pbuf_get_index('kvm') 
Turbtype_idx  = pbuf_get_index('turbtype') 
smaw_idx      = pbuf_get_index('smaw') 
Tauresx_idx   = pbuf_get_index('tauresx') 
Tauresy_idx   = pbuf_get_index('tauresy') 
tpert_idx     = pbuf_get_index('tpert') 
qpert_idx     = pbuf_get_index('qpert')

t_ttend_oldid   = pbuf_get_index('T_TTEND_OLD') 
TEOUT_oldid     = pbuf_get_index('TEOUT_OLD')
DTCORE_oldid    = pbuf_get_index('DTCORE_OLD') 
QCWAT_oldid     = pbuf_get_index('QCWAT_OLD') 
LCWAT_oldid     = pbuf_get_index('LCWAT_OLD') 
TCWAT_oldid     = pbuf_get_index('TCWAT_OLD') 
CLD_oldid       = pbuf_get_index('CLD_OLD') 
AST_oldid       = pbuf_get_index('AST_OLD') 
CONCLD_oldid    = pbuf_get_index('CONCLD_OLD') 
DP_FLXPRC_oldid = pbuf_get_index('DP_FLXPRC_OLD') 
DP_FLXSNW_oldid = pbuf_get_index('DP_FLXSNW_OLD') 
DP_CLDLIQ_oldid = pbuf_get_index('DP_CLDLIQ_OLD') 
DP_CLDICE_oldid = pbuf_get_index('DP_CLDICE_OLD') 
cush_oldid      = pbuf_get_index('cush_OLD') 
QRS_oldid       = pbuf_get_index('QRS_OLD') 
QRL_oldid       = pbuf_get_index('QRL_OLD') 
pblh_oldid      = pbuf_get_index('pblh_OLD') 
tke_oldid       = pbuf_get_index('tke_OLD') 
kvh_oldid       = pbuf_get_index('kvh_OLD') 
kvm_oldid       = pbuf_get_index('kvm_OLD') 
Turbtype_oldid  = pbuf_get_index('turbtype_OLD') 
smaw_oldid      = pbuf_get_index('smaw_OLD') 
Tauresx_oldid   = pbuf_get_index('tauresx_OLD') 
Tauresy_oldid   = pbuf_get_index('tauresy_OLD') 
tpert_oldid     = pbuf_get_index('tpert_OLD') 
qpert_oldid     = pbuf_get_index('qpert_OLD')

call pbuf_get_field(pbuf,t_ttend_oldid  , t_ttend_old  ) 
call pbuf_get_field(pbuf,TEOUT_oldid    , TEOUT_old    )
call pbuf_get_field(pbuf,DTCORE_oldid   , DTCORE_old   )
call pbuf_get_field(pbuf,QCWAT_oldid    , QCWAT_old    )
call pbuf_get_field(pbuf,LCWAT_oldid    , LCWAT_old    )
call pbuf_get_field(pbuf,TCWAT_oldid    , TCWAT_old    )
call pbuf_get_field(pbuf,CLD_oldid      , CLD_old      )
call pbuf_get_field(pbuf,AST_oldid      , AST_old      )
call pbuf_get_field(pbuf,CONCLD_oldid   , CONCLD_old   )
call pbuf_get_field(pbuf,DP_FLXPRC_oldid, DP_FLXPRC_old)
call pbuf_get_field(pbuf,DP_FLXSNW_oldid, DP_FLXSNW_old)
call pbuf_get_field(pbuf,DP_CLDLIQ_oldid, DP_CLDLIQ_old)
call pbuf_get_field(pbuf,DP_CLDICE_oldid, DP_CLDICE_old)
call pbuf_get_field(pbuf,cush_oldid     , cush_old     )
call pbuf_get_field(pbuf,QRS_oldid      , QRS_old      )
call pbuf_get_field(pbuf,QRL_oldid      , QRL_old      )
call pbuf_get_field(pbuf,pblh_oldid     , pblh_old     )
call pbuf_get_field(pbuf,tke_oldid      , tke_old      )
call pbuf_get_field(pbuf,kvh_oldid      , kvh_old      )
call pbuf_get_field(pbuf,kvm_oldid      , kvm_old      )
call pbuf_get_field(pbuf,turbtype_oldid , Turbtype_old )
call pbuf_get_field(pbuf,smaw_oldid     , smaw_old     )
call pbuf_get_field(pbuf,tauresx_oldid  , Tauresx_old  )
call pbuf_get_field(pbuf,tauresy_oldid  , Tauresy_old  )
call pbuf_get_field(pbuf,tpert_oldid    , tpert_old    )
call pbuf_get_field(pbuf,qpert_oldid    , qpert_old    )

call pbuf_get_field(pbuf,t_ttend_idx  , t_ttend   ) 
call pbuf_get_field(pbuf,TEOUT_idx    , TEOUT     )
call pbuf_get_field(pbuf,DTCORE_idx   , DTCORE    )
call pbuf_get_field(pbuf,QCWAT_idx    , QCWAT     )
call pbuf_get_field(pbuf,LCWAT_idx    , LCWAT     )
call pbuf_get_field(pbuf,TCWAT_idx    , TCWAT     )
call pbuf_get_field(pbuf,CLD_idx      , CLD       )
call pbuf_get_field(pbuf,AST_idx      , AST       )
call pbuf_get_field(pbuf,CONCLD_idx   , CONCLD    )
call pbuf_get_field(pbuf,DP_FLXPRC_idx, DP_FLXPRC )
call pbuf_get_field(pbuf,DP_FLXSNW_idx, DP_FLXSNW )
call pbuf_get_field(pbuf,DP_CLDLIQ_idx, DP_CLDLIQ )
call pbuf_get_field(pbuf,DP_CLDICE_idx, DP_CLDICE )
call pbuf_get_field(pbuf,cush_idx     , cush      )
call pbuf_get_field(pbuf,QRS_idx      , QRS       )
call pbuf_get_field(pbuf,QRL_idx      , QRL       )
call pbuf_get_field(pbuf,pblh_idx     , pblh      )
call pbuf_get_field(pbuf,tke_idx      , tke       )
call pbuf_get_field(pbuf,kvh_idx      , kvh       )
call pbuf_get_field(pbuf,kvm_idx      , kvm       )
call pbuf_get_field(pbuf,turbtype_idx , Turbtype  )
call pbuf_get_field(pbuf,smaw_idx     , smaw      )
call pbuf_get_field(pbuf,tauresx_idx  , Tauresx   )
call pbuf_get_field(pbuf,tauresy_idx  , Tauresy   )
call pbuf_get_field(pbuf,tpert_idx    , tpert     )
call pbuf_get_field(pbuf,qpert_idx    , qpert     )

ncol = get_ncols_p(lchnk)
do i = 1, ncol 
 do k = 1, pver
t_ttend(i,k)  = t_ttend_old(i,k)   
TEOUT(i,k)    = TEOUT_old(i,k)     
DTCORE(i,k)   = DTCORE_old(i,k)    
QCWAT(i,k)    = QCWAT_old(i,k)     
LCWAT(i,k)    = LCWAT_old(i,k)     
TCWAT(i,k)    = TCWAT_old(i,k)     
CLD(i,k)      = CLD_old(i,k)       
AST(i,k)      = AST_old(i,k)       
CONCLD(i,k)   = CONCLD_old(i,k)    
DP_FLXPRC(i,k)= DP_FLXPRC_old(i,k) 
DP_FLXSNW(i,k)= DP_FLXSNW_old(i,k) 
DP_CLDLIQ(i,k)= DP_CLDLIQ_old(i,k) 
DP_CLDICE(i,k)= DP_CLDICE_old(i,k) 
QRS(i,k)      = QRS_old(i,k)       
QRL(i,k)      = QRL_old(i,k)       
pblh(i,k)     = pblh_old(i,k)      
tke(i,k)      = tke_old(i,k)       
kvh(i,k)      = kvh_old(i,k)       
kvm(i,k)      = kvm_old(i,k)       
turbtype(i,k) = turbtype_old(i,k)  
smaw(i,k)     = smaw_old(i,k)      
end do
tauresx(i)  = tauresx_old(i)   
tauresy(i)  = tauresy_old(i)   
tpert(i)    = tpert_old(i)     
cush(i)     = cush_old(i)      
qpert(i,1)  = qpert_old(i,1)
end do
end do

!        call cam_read_restart( cam_out, dyn_in, dyn_out, pbuf2d, -1, 0, NLFileName=filein, passed_loc=locfn)

        do lchnk = begchunk,endchunk
             ncol = get_ncols_p(lchnk)
             do k = 1, pver
                do i = 1, ncol
                    phys_tend(lchnk)%dudt (i,k)=0
                    phys_tend(lchnk)%dvdt (i,k)=0
                end do
             end do
        end do
        do lchnk = begchunk,endchunk
             ncol = get_ncols_p(lchnk)
             do i = 1, ncol
                    phys_tend(lchnk)%flx_net (i)=0
                    phys_tend(lchnk)%te_tnd (i)=0
                    phys_tend(lchnk)%tw_tnd (i)=0
             end do
        end do
end if

if ( mod(sec_spec,21600)==0  .AND. .NOT. do_restart ) then 
        do_restart=.TRUE. 

do lchnk=begchunk,endchunk 
    pbuf=> pbuf_get_chunk(pbuf2d,lchnk)

      itim = pbuf_old_tim_idx()

t_ttend_idx   = pbuf_get_index('T_TTEND') 
TEOUT_idx     = pbuf_get_index('TEOUT') 
DTCORE_idx    = pbuf_get_index('DTCORE') 
QCWAT_idx     = pbuf_get_index('QCWAT') 
LCWAT_idx     = pbuf_get_index('LCWAT') 
TCWAT_idx     = pbuf_get_index('TCWAT') 
CLD_idx       = pbuf_get_index('CLD') 
AST_idx       = pbuf_get_index('AST') 
CONCLD_idx    = pbuf_get_index('CONCLD') 
DP_FLXPRC_idx = pbuf_get_index('DP_FLXPRC') 
DP_FLXSNW_idx = pbuf_get_index('DP_FLXSNW') 
DP_CLDLIQ_idx = pbuf_get_index('DP_CLDLIQ') 
DP_CLDICE_idx = pbuf_get_index('DP_CLDICE') 
cush_idx      = pbuf_get_index('cush') 
QRS_idx       = pbuf_get_index('QRS') 
QRL_idx       = pbuf_get_index('QRL') 
pblh_idx      = pbuf_get_index('pblh') 
tke_idx       = pbuf_get_index('tke') 
kvh_idx       = pbuf_get_index('kvh') 
kvm_idx       = pbuf_get_index('kvm') 
Turbtype_idx  = pbuf_get_index('turbtype') 
smaw_idx      = pbuf_get_index('smaw') 
Tauresx_idx   = pbuf_get_index('tauresx') 
Tauresy_idx   = pbuf_get_index('tauresy') 
tpert_idx     = pbuf_get_index('tpert') 
qpert_idx     = pbuf_get_index('qpert')

t_ttend_oldid   = pbuf_get_index('T_TTEND_OLD') 
TEOUT_oldid     = pbuf_get_index('TEOUT_OLD')
DTCORE_oldid    = pbuf_get_index('DTCORE_OLD') 
QCWAT_oldid     = pbuf_get_index('QCWAT_OLD') 
LCWAT_oldid     = pbuf_get_index('LCWAT_OLD') 
TCWAT_oldid     = pbuf_get_index('TCWAT_OLD') 
CLD_oldid       = pbuf_get_index('CLD_OLD') 
AST_oldid       = pbuf_get_index('AST_OLD') 
CONCLD_oldid    = pbuf_get_index('CONCLD_OLD') 
DP_FLXPRC_oldid = pbuf_get_index('DP_FLXPRC_OLD') 
DP_FLXSNW_oldid = pbuf_get_index('DP_FLXSNW_OLD') 
DP_CLDLIQ_oldid = pbuf_get_index('DP_CLDLIQ_OLD') 
DP_CLDICE_oldid = pbuf_get_index('DP_CLDICE_OLD') 
cush_oldid      = pbuf_get_index('cush_OLD') 
QRS_oldid       = pbuf_get_index('QRS_OLD') 
QRL_oldid       = pbuf_get_index('QRL_OLD') 
pblh_oldid      = pbuf_get_index('pblh_OLD') 
tke_oldid       = pbuf_get_index('tke_OLD') 
kvh_oldid       = pbuf_get_index('kvh_OLD') 
kvm_oldid       = pbuf_get_index('kvm_OLD') 
Turbtype_oldid  = pbuf_get_index('turbtype_OLD') 
smaw_oldid      = pbuf_get_index('smaw_OLD') 
Tauresx_oldid   = pbuf_get_index('tauresx_OLD') 
Tauresy_oldid   = pbuf_get_index('tauresy_OLD') 
tpert_oldid     = pbuf_get_index('tpert_OLD') 
qpert_oldid     = pbuf_get_index('qpert_OLD')

call pbuf_get_field(pbuf,t_ttend_oldid  , t_ttend_old  ) 
call pbuf_get_field(pbuf,TEOUT_oldid    , TEOUT_old    )
call pbuf_get_field(pbuf,DTCORE_oldid   , DTCORE_old   )
call pbuf_get_field(pbuf,QCWAT_oldid    , QCWAT_old    )
call pbuf_get_field(pbuf,LCWAT_oldid    , LCWAT_old    )
call pbuf_get_field(pbuf,TCWAT_oldid    , TCWAT_old    )
call pbuf_get_field(pbuf,CLD_oldid      , CLD_old      )
call pbuf_get_field(pbuf,AST_oldid      , AST_old      )
call pbuf_get_field(pbuf,CONCLD_oldid   , CONCLD_old   )
call pbuf_get_field(pbuf,DP_FLXPRC_oldid, DP_FLXPRC_old)
call pbuf_get_field(pbuf,DP_FLXSNW_oldid, DP_FLXSNW_old)
call pbuf_get_field(pbuf,DP_CLDLIQ_oldid, DP_CLDLIQ_old)
call pbuf_get_field(pbuf,DP_CLDICE_oldid, DP_CLDICE_old)
call pbuf_get_field(pbuf,cush_oldid     , cush_old     )
call pbuf_get_field(pbuf,QRS_oldid      , QRS_old      )
call pbuf_get_field(pbuf,QRL_oldid      , QRL_old      )
call pbuf_get_field(pbuf,pblh_oldid     , pblh_old     )
call pbuf_get_field(pbuf,tke_oldid      , tke_old      )
call pbuf_get_field(pbuf,kvh_oldid      , kvh_old      )
call pbuf_get_field(pbuf,kvm_oldid      , kvm_old      )
call pbuf_get_field(pbuf,turbtype_oldid , Turbtype_old )
call pbuf_get_field(pbuf,smaw_oldid     , smaw_old     )
call pbuf_get_field(pbuf,tauresx_oldid  , Tauresx_old  )
call pbuf_get_field(pbuf,tauresy_oldid  , Tauresy_old  )
call pbuf_get_field(pbuf,tpert_oldid    , tpert_old    )
call pbuf_get_field(pbuf,qpert_oldid    , qpert_old    )

call pbuf_get_field(pbuf,t_ttend_idx  , t_ttend   ) 
call pbuf_get_field(pbuf,TEOUT_idx    , TEOUT     )
call pbuf_get_field(pbuf,DTCORE_idx   , DTCORE    )
call pbuf_get_field(pbuf,QCWAT_idx    , QCWAT     )
call pbuf_get_field(pbuf,LCWAT_idx    , LCWAT     )
call pbuf_get_field(pbuf,TCWAT_idx    , TCWAT     )
call pbuf_get_field(pbuf,CLD_idx      , CLD       )
call pbuf_get_field(pbuf,AST_idx      , AST       )
call pbuf_get_field(pbuf,CONCLD_idx   , CONCLD    )
call pbuf_get_field(pbuf,DP_FLXPRC_idx, DP_FLXPRC )
call pbuf_get_field(pbuf,DP_FLXSNW_idx, DP_FLXSNW )
call pbuf_get_field(pbuf,DP_CLDLIQ_idx, DP_CLDLIQ )
call pbuf_get_field(pbuf,DP_CLDICE_idx, DP_CLDICE )
call pbuf_get_field(pbuf,cush_idx     , cush      )
call pbuf_get_field(pbuf,QRS_idx      , QRS       )
call pbuf_get_field(pbuf,QRL_idx      , QRL       )
call pbuf_get_field(pbuf,pblh_idx     , pblh      )
call pbuf_get_field(pbuf,tke_idx      , tke       )
call pbuf_get_field(pbuf,kvh_idx      , kvh       )
call pbuf_get_field(pbuf,kvm_idx      , kvm       )
call pbuf_get_field(pbuf,turbtype_idx , Turbtype  )
call pbuf_get_field(pbuf,smaw_idx     , smaw      )
call pbuf_get_field(pbuf,tauresx_idx  , Tauresx   )
call pbuf_get_field(pbuf,tauresy_idx  , Tauresy   )
call pbuf_get_field(pbuf,tpert_idx    , tpert     )
call pbuf_get_field(pbuf,qpert_idx    , qpert     )

ncol = get_ncols_p(lchnk)
do i = 1, ncol 
 do k = 1, pver
t_ttend_old(i,k)  = t_ttend(i,k)   
TEOUT_old(i,k)    = TEOUT(i,k)     
DTCORE_old(i,k)   = DTCORE(i,k)    
QCWAT_old(i,k)    = QCWAT(i,k)     
LCWAT_old(i,k)    = LCWAT(i,k)     
TCWAT_old(i,k)    = TCWAT(i,k)     
CLD_old(i,k)      = CLD(i,k)       
AST_old(i,k)      = AST(i,k)       
CONCLD_old(i,k)   = CONCLD(i,k)    
DP_FLXPRC_old(i,k)= DP_FLXPRC(i,k) 
DP_FLXSNW_old(i,k)= DP_FLXSNW(i,k) 
DP_CLDLIQ_old(i,k)= DP_CLDLIQ(i,k) 
DP_CLDICE_old(i,k)= DP_CLDICE(i,k) 
QRS_old(i,k)      = QRS(i,k)       
QRL_old(i,k)      = QRL(i,k)       
pblh_old(i,k)     = pblh(i,k)      
tke_old(i,k)      = tke(i,k)       
kvh_old(i,k)      = kvh(i,k)       
kvm_old(i,k)      = kvm(i,k)       
turbtype_old(i,k) = turbtype(i,k)  
smaw_old(i,k)     = smaw(i,k)      
end do
tauresx_old(i)  = tauresx(i)   
tauresy_old(i)  = tauresy(i)   
tpert_old(i)    = tpert(i)     
cush_old(i)     = cush(i)      
qpert_old(i,1)  = qpert(i,1)
end do
end do
end if

   call t_startf ('cam_run4_wrapup')
   call wrapup(rstwr, nlend)
   call t_stopf  ('cam_run4_wrapup')

   if (masterproc .and. print_step_cost) then
      call t_startf ('cam_run4_print')
      call t_stampf (wcend, usrend, sysend)
      write(iulog,'(a,3f8.3,a)')'Prv timestep wallclock, usr, sys=', &
                            wcend-wcstart, usrend-usrstart, sysend-sysstart, &
                            ' seconds'
      call t_stopf  ('cam_run4_print')
   end if

#ifndef UNICOSMP
   call t_startf ('cam_run4_flush')
   call shr_sys_flush(iulog)
   call t_stopf  ('cam_run4_flush')
#endif

end subroutine cam_run4

!
!-----------------------------------------------------------------------
!

subroutine cam_final( cam_out, cam_in )
!-----------------------------------------------------------------------
!
! Purpose:  CAM finalization.
!
!-----------------------------------------------------------------------
   use units,            only: getunit
   use time_manager,     only: get_nstep, get_step_size
#if ( defined SPMD )
   use mpishorthand,     only: mpicom, mpiint, &
                               nsend, nrecv, nwsend, nwrecv
   use spmd_utils,       only: iam, npes
#endif
   use stepon,           only: stepon_final
   use physpkg,          only: phys_final
   use cam_initfiles,    only: cam_initfiles_close
   use camsrfexch,       only: atm2hub_deallocate, hub2atm_deallocate
   !
   ! Arguments
   !
   type(cam_out_t), pointer :: cam_out(:) ! Output from CAM to surface
   type(cam_in_t), pointer :: cam_in(:)   ! Input from merged surface to CAM
!-----------------------------------------------------------------------
   !
   ! Local variables
   !
   integer :: nstep           ! Current timestep number.

#if ( defined SPMD )   
   integer :: iu              ! SPMD Statistics output unit number
   character*24 :: filenam    ! SPMD Stats output filename
   integer :: signal          ! MPI message buffer
!------------------------------Externals--------------------------------
   real(r8) :: mpi_wtime
#endif

   call phys_final( phys_state, phys_tend , pbuf2d)
   call stepon_final(dyn_in, dyn_out)

   if(nsrest==0) then
      call cam_initfiles_close()
   end if

   call hub2atm_deallocate(cam_in)
   call atm2hub_deallocate(cam_out)

#if ( defined SPMD )
   if (.false.) then
      write(iulog,*)'The following stats are exclusive of initialization/boundary datasets'
      write(iulog,*)'Number of messages sent by proc ',iam,' is ',nsend
      write(iulog,*)'Number of messages recv by proc ',iam,' is ',nrecv
   end if
#endif

   ! This flush attempts to ensure that asynchronous diagnostic prints from all 
   ! processes do not get mixed up with the "END OF MODEL RUN" message printed 
   ! by masterproc below.  The test-model script searches for this message in the 
   ! output log to figure out if CAM completed successfully.  This problem has 
   ! only been observed with the Linux Lahey compiler (lf95) which does not 
   ! support line-buffered output.  
#ifndef UNICOSMP
   call shr_sys_flush( 0 )   ! Flush all output to standard error
   call shr_sys_flush( iulog )   ! Flush all output to standard output
#endif

   if (masterproc) then
#if ( defined SPMD )
      cam_time_end = mpi_wtime()
#endif
      nstep = get_nstep()
      write(iulog,9300) nstep-1,nstep
9300  format (//'Number of completed timesteps:',i6,/,'Time step ',i6, &
                ' partially done to provide convectively adjusted and ', &
                'time filtered values for history tape.')
      write(iulog,*)'------------------------------------------------------------'
#if ( defined SPMD )
      write(iulog,*)
      write(iulog,*)' Total run time (sec) : ', cam_time_end-cam_time_beg
      write(iulog,*)' Time Step Loop run time(sec) : ', stepon_time_end-stepon_time_beg
      if (((nstep-1)-nstep_beg) > 0) then
         write(iulog,*)' SYPD : ',  &
         236.55_r8/((86400._r8/(dtime*((nstep-1)-nstep_beg)))*(stepon_time_end-stepon_time_beg))
      endif
      write(iulog,*)
#endif
      write(iulog,*)'******* END OF MODEL RUN *******'
   end if


#if ( defined SPMDSTATS )
   if (t_single_filef()) then
      write(filenam,'(a17)') 'spmdstats_cam.all'
      iu = getunit ()
      if (iam .eq. 0) then
         open (unit=iu, file=filenam, form='formatted', status='replace')
         signal = 1
      else
         call mpirecv(signal, 1, mpiint, iam-1, iam, mpicom) 
         open (unit=iu, file=filenam, form='formatted', status='old', position='append')
      endif
      write (iu,*)'************ PROCESS ',iam,' ************'
      write (iu,*)'iam ',iam,' msgs  sent =',nsend
      write (iu,*)'iam ',iam,' msgs  recvd=',nrecv
      write (iu,*)'iam ',iam,' words sent =',nwsend
      write (iu,*)'iam ',iam,' words recvd=',nwrecv
      write (iu,*)
      close(iu)
      if (iam+1 < npes) then
         call mpisend(signal, 1, mpiint, iam+1, iam+1, mpicom)
      endif
   else
      iu = getunit ()
      write(filenam,'(a14,i5.5)') 'spmdstats_cam.', iam
      open (unit=iu, file=filenam, form='formatted', status='replace')
      write (iu,*)'************ PROCESS ',iam,' ************'
      write (iu,*)'iam ',iam,' msgs  sent =',nsend
      write (iu,*)'iam ',iam,' msgs  recvd=',nrecv
      write (iu,*)'iam ',iam,' words sent =',nwsend
      write (iu,*)'iam ',iam,' words recvd=',nwrecv
      close(iu)
   endif
#endif

end subroutine cam_final

!
!-----------------------------------------------------------------------
!

end module cam_comp
