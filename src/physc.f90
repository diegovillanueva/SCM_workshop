#ifdef __xlC__
@PROCESS HOT
@PROCESS SPILLSIZE(5000)
#endif
!OCL NOALIAS

!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE physc(krow)
  !
  ! Description:
  !
  ! Controls the calls to the various physical subroutines.
  !
  !
  !  *physc* is called from *gpc*.
  !
  !  Externals:
  !
  !  *geopot*    computes full level geopotentials.
  !  *pres*      computes half level pressure.
  !  *presf*     computes full level pressure.
  !  *radiation* controls radiation computations.
  !  *radheat*   computes radiation tendencies.
  !  *vdiff*     computes vertical exchange by turbulence.
  !  *ssodrag*   computes gravity wave drag.
  !  *cucall*    controls mass-flux scheme.
  !  *cloud*     computes large scale water phase changes and cloud cover.
  !  *ml_ocean*  computes mixed layer ocean
  !
  !  Authors:
  !
  !  M. Jarraud, ECMWF, January 1982, original source
  !
  !     Modifications.
  !     --------------
  !  M.A. Giorgetta, MPI-Hamburg, May 2000, modified for ECHAM5
  !  A.Tompkins,     MPI-Hamburg, June  2000, cloud cover scheme
  !  U. Schlese,     MPI-Hamburg, July  2000, calling sequence changed
  !  I. Kirchner,    MPI, July 2000, tendency diagnostics revision
  !  L. Kornblueh,   MPI, August 2001, changes for different advection
  !                       schemes, adapt polefilter call 
  !  U. Schulzweida, MPI, May 2002, blocking (nproma)
  !  I. Kirchner,    MPI, August 2002, nmi revision/extension
  !  U. Schlese, M. Esch, MPI, September 2002, mixed layer ocean
  !  L. Kornblueh,   MPI, August 2004, new tropopause calculation
  !  S.J. Lorenz,    MPI, November 2007, volcanic forcing
  !  M. Esch,        MPI, November 2007, switch for double radiation
  !  J. Kazil        MPI, October 2008,
  !                  Modified M7 call, export of PBL top level,
  !                  tropopause level conversion integer -> double
  !  D.O'Donnell     MPI-Met, February 2009, added routine 'call_chem0' to
  !                  allow submodels to do some work before radiation is called
  !                  modified call_chem1 and call_chem2 calls, see tag '>>dod' 
  !  L. Kornblueh    MPI, Februrary 2012, remove spitfire
  !
  USE mo_kind,              ONLY: dp
  USE mo_memory_g1a,        ONLY: xlm1, xim1, tm1, qm1, alpsm1, xtm1
  USE mo_memory_g2a,        ONLY: vm1, um1
  USE mo_memory_g3a,        ONLY: aprfluxm, geospm, ustrm, ustrgwm,            &
       seaicem, acdncm
  USE mo_memory_g3b        
USE mo_activ,             ONLY: icnc_instantan
  USE mo_control,           ONLY: ltdiag, lcouple, lmidatm, lhd,               &
       nlev, nlevp1, ldiagamip, ltimer,             &
       lmeltpond, lcouple_co2, lfractional_mask
  USE mo_submodel,          ONLY: lanysubmodel
  USE mo_vphysc,            ONLY: vphysc, set_vphysc_var !++s.stadtler
  USE mo_hyb,               ONLY: delb, nlevm1
  USE mo_param_switches,    ONLY: lcond, lconv, lgwdrag, lrad, &
                                  lcdnc_progn, nic_cirrus
  USE mo_math_constants,    ONLY: pi
  USE mo_physical_constants,ONLY: grav, rd, vtmpc1, tmelt, cpd, vtmpc2
  USE mo_radiation_parameters,  ONLY: iaero, flx_ratio_cur, solc
  USE mo_radiation,         ONLY: radiation
  USE mo_cover,             ONLY: cover
  USE mo_cloud,             ONLY: cloud
  USE mo_scan_buffer,       ONLY: vol, vom, qte, xlte, xite, tte, xtte,      &
       alnpr, alpha, alpste, vervel
  USE mo_physc2,            ONLY: cqsncr, cwlmax
  USE mo_tracdef,           ONLY: ntrac
  !ham_ps: interface via call_submodels: xtdriver1, xtdriver2, xtdiagn
  !
  USE mo_midatm,            ONLY: gwspectrum
  USE mo_ssortns,           ONLY: ssodrag
  !
  USE mo_cumastr,           ONLY: cucall
  !
  USE mo_decomposition,     ONLY: ldc => local_decomposition
  USE mo_geoloc,            ONLY: amu0_x, rdayl_x, sqcst_2d, &
       philat_2d, set_geo_loop
  !
  USE mo_time_control,      ONLY: lstart, lresume, delta_time, l_trigrad,&
       time_step_len, l_getocean
  USE mo_advection

  USE mo_timer,             ONLY: timer_start,   timer_stop,    timer_radiation, &
       timer_cloud,   timer_cover,   timer_vdiff,     &
       timer_radheat,   &
       timer_gwdrag,  timer_cucall
  !
  USE mo_nmi,               ONLY: nmi_phase, NMI_ACCU, NMI_USE_AVG,              &
       dh_t, dh_m, dh_l, buf_t, buf_m, buf_l,         &
       lnmi_run, lnmi_cloud
  USE mo_diag_amip2,        ONLY: collect_amip2_diag 
  USE mo_tropopause,        ONLY: WMO_tropopause
  USE mo_cosp_simulator,    ONLY: locosp, call_cospsimulator
  USE mo_cosp_offline,      ONLY: locospoffl, cospoffl_geom1, cospoffl_geohm1, &
       cospoffl_p, cospoffl_ph
  ! JSBACH interface add-ons
  USE mo_surface_memory,    ONLY: jrsfl, jrsfc, jssfl, jssfc, ztrfli,            &
       zsofli,                                        &
       jsswnir, jsswdifnir, jsswvis, jsswdifvis,      &
       jsswpar, jsswdifpar, jsswniracc, jsswvisacc,   &
       jsswparacc, jsswdifniracc, jsswdifvisacc,      &
       jsswdifparacc
  USE mo_co2,               ONLY: co2m1, co2_flux_ocean, co2atmos, co2flux_cpl,  &
       co2_flux_atmosphere_ocean
  USE mo_surface_ice,       ONLY: meltpond, meltpond_ice
  USE mo_diag_tendency_new, ONLY: tdiag_vars, set_tendency
  !USE mo_test

  USE mo_submodel_interface,ONLY: physc_subm_1, &
       physc_subm_2, &
       physc_subm_3, &
       physc_subm_4

  USE mo_hydrology,      ONLY: hydrology_collect_lake
  USE mo_station_diag,   ONLY : lostation, collect_station_diag
  USE mo_memory_cfdiag,  ONLY : locfdiag, calc_cfdiag

!>>SF
  USE mo_cloud_micro_2m, ONLY: cloud_micro_interface
!<<SF

  IMPLICIT NONE
  !
#ifdef _EUROHACK_2015
  LOGICAL, SAVE :: ldycore_only = .TRUE. 
#endif
  !
  INTEGER :: krow
  ! 
  !  Local array bounds
  INTEGER :: nproma, nbdim
  !
  !  Local scalars:
  REAL(dp) :: zcst, zrcst, ztwodt, zprat,  zcdnc,  zn1,  zn2                    &
       , zepsec, zsigfac, zsigh, zsn_mm, zcons1
  INTEGER :: jlev, jl, kfdia, kidia, ktdia, jk, nexp
  LOGICAL :: loconv, locond
  ! 
  !  Local arrays:
  REAL(dp) ::  zdadc(ldc%nproma), zdpsdt(ldc%nproma), zgeo(ldc%nproma)          &
       ,ztvm1(ldc%nproma,nlev), ztslnew(ldc%nproma)                          &
       ,zcair(ldc%nproma,nlev),zi0(ldc%nproma)
  REAL(dp) ::  zqtec(ldc%nproma,nlev)
  REAL(dp) ::  zxtecl(ldc%nproma,nlev),zxteci(ldc%nproma,nlev)
  REAL(dp) ::  zqtold(ldc%nproma) !!$, zflux_corr_acc(ldc%nproma)
  !
  LOGICAL :: lonorth(ldc%nproma)
  !
  !  Surface fluxes over land/water/ice
  !
  REAL(dp) :: zhfsl(ldc%nproma),  zhfsw(ldc%nproma),  zhfsi(ldc%nproma),         &
       zhflw(ldc%nproma),  zhfli(ldc%nproma),   zevapw(ldc%nproma),           &
       zevapi(ldc%nproma), ztrflw(ldc%nproma),  zsofll(ldc%nproma),           &
       zsoflw(ldc%nproma)

  REAL(dp) ::  zfrl(ldc%nproma),  zfrw(ldc%nproma),   zfri(ldc%nproma)           &
       ,zcvsc(ldc%nproma),  zwlmx(ldc%nproma)                                 &
       ,zcvs(ldc%nproma),  zcvw(ldc%nproma)
  REAL(dp) :: zteffl4(ldc%nproma),ztsnew(ldc%nproma), zradtemp_old(ldc%nproma)
  REAL(dp) :: zrain(ldc%nproma), zsnow(ldc%nproma)
  ! 
  INTEGER :: ilab(ldc%nproma,nlev), itype(ldc%nproma), ictop(ldc%nproma)
  INTEGER :: invb(ldc%nproma)
  !
  INTEGER :: itrpwmo(ldc%nproma), itrpwmop1(ldc%nproma)
  !
  !
  !    Arrays internal to physics
  !
  REAL(dp) :: qhfla(ldc%nproma), evapot(ldc%nproma), zprecip(ldc%nproma)
  !
  REAL(dp) :: zdtime
  !
!---Included for prognostic CDNC/IC scheme (Ulrike Lohmann, 11/02/2007)---
REAL(dp):: zicnc(ldc%nproma,nlev)
!--- End included -----------------------------------------------------------
  ! Local arrays
  !
  LOGICAL  :: loland(ldc%nproma), loglac(ldc%nproma), lolake(ldc%nproma)
  REAL(dp) :: geom1(ldc%nproma,nlev)
  REAL(dp) :: aphm1(ldc%nproma,nlevp1), apm1(ldc%nproma,nlev) 
  REAL(dp) :: aphp1(ldc%nproma,nlevp1), app1(ldc%nproma,nlev) 
  REAL(dp) :: rsfc(ldc%nproma), ssfc(ldc%nproma)
  REAL(dp) :: rsfl(ldc%nproma), ssfl(ldc%nproma)
  REAL(dp) :: zpbl(ldc%nproma)
  REAL(dp) :: alpha0(ldc%nproma,nlev)
  REAL(dp) :: geohm1(ldc%nproma,nlevp1)
  REAL(dp) :: zxlp1(ldc%nproma,nlev)
  !
  ! Local Array for ssodrag
  REAL(dp)::    zlat(ldc%nproma)      ! latitude (radian)
  ! 
  !  External subroutines
  EXTERNAL :: geopot, pres, presf, vdiff, radheat, collect
  ! 
#ifdef _EUROHACK_2015
  IF (.NOT. ldycore_only) THEN
#endif
  !  Local array bounds

  nbdim = ldc% nproma

  IF ( krow == ldc% ngpblks ) THEN
    nproma = ldc% npromz
  ELSE
    nproma = ldc% nproma
  END IF
  !
  !*    COMPUTATIONAL CONSTANTS.
  !     ------------- ----------
  !
  zdtime = delta_time
  zepsec=1.E-12_dp
  zsigfac=0.15_dp
  zcons1 = cpd*vtmpc2
  !
  !     ----------------------------------------------------------------
  !
  !*        2.    ALLOCATE STORAGE 
  !               -------- ------- 
  !
  !     ------------------------------------------------------------
  !
  !*        3.    COMPUTE SOME FIELDS NEEDED BY THE PHYSICAL ROUTINES.
  !               ------- ---- ------ ------ -- --- -------- ---------
  !               set parameters depending on the location

  CALL set_geo_loop (krow)

  !
  !*        3.0   Set start conditions for JSBACH interface variables
  !
  ! INITIALIZE SURFACE PARAMETERS WHICH COME FROM MO_SURFACE IN THE
  !  NEXT TIME STEP
  IF (lstart) THEN
    jrsfl(:,krow) = 0.0_dp  ! convectiv, large scale
    jrsfc(:,krow) = 0.0_dp  ! snow and rain fall
    jssfl(:,krow) = 0.0_dp  ! with units in
    jssfc(:,krow) = 0.0_dp  ! kg/m2/s, set for initial values
  END IF
  zrain(:) = 0._dp
  zsnow(:) = 0._dp
  zrain(1:nproma) = jrsfl(1:nproma,krow) + jrsfc(1:nproma,krow)
  zsnow(1:nproma) = jssfl(1:nproma,krow) + jssfc(1:nproma,krow)
  !
  !*        3.1   COMPUTE VIRTUAL TEMPERATURE AT T-DT AND SET *ZGEO* TO 0.
  !
  ztvm1(1:nproma,:) = tm1(1:nproma,:,krow)*(1._dp+vtmpc1*qm1(1:nproma,:,krow)  &
       -(xlm1(1:nproma,:,krow)+xim1(1:nproma,:,krow)))
  !
  zgeo(1:nproma)=0._dp
  !
  !*        3.2   COMPUTE (PHI-PHIS) AT T-DT USING LN(P) AT T.
  !
  CALL geopot(geom1,ztvm1,alnpr(:,:,krow),alpha(:,:,krow),zgeo,nbdim,nproma)
  IF (lanysubmodel) THEN
    ! store geopotential height in vphysc for use in bounday condition scheme
    ! and compute geopotential at interfaces (e.g. for vertical interpolation)
    vphysc%geom1(1:nproma,:,krow) = geom1(1:nproma,:)
    alpha0 = 0._dp
    CALL geopot(geohm1(:,2:nlevp1),ztvm1,alnpr(:,:,krow),alpha0,zgeo,nbdim,nproma) ! half level geopotential
    geohm1(:,1) = 1.e30_dp                              ! set upper boundary to a huge value
    vphysc%geohm1(1:nproma,:,krow) = geohm1(1:nproma,:)
  END IF
  !
  ! 3.2b Specific heat of moist air
  !
  DO jk=1,nlev
    DO jl=1,nproma
      zcair(jl,jk)=cpd+zcons1*MAX(qm1(jl,jk,krow),0.0_dp)
    END DO
  END DO
  !
  !*        3.3    COMPUTE PRESSURE AT FULL AND HALF LEVELS AT T-DT.
  !
  aphm1(1:nproma,nlevp1)=EXP(alpsm1(1:nproma,krow))
  !
  CALL pres(aphm1,nbdim,aphm1(1,nlevp1),nproma)
  IF (lanysubmodel) THEN
     vphysc%aphm1(1:nproma,:,krow) = aphm1(1:nproma,:)
  END IF

  !
  CALL presf(apm1,nbdim,aphm1,nproma)
  !
  !*        3.4   COMPUTE REAL WINDS AND WIND TENDENCIES.
  !
  DO jlev = 1, nlev
!DIR$ CONCURRENT
!IBM* novector
    DO jl = 1, nproma
      zrcst=1._dp/sqcst_2d(jl,krow)
      um1(jl,jlev,krow) = zrcst*um1(jl,jlev,krow)
      vm1(jl,jlev,krow) = zrcst*vm1(jl,jlev,krow)
      vol(jl,jlev,krow) = zrcst*vol(jl,jlev,krow)
      vom(jl,jlev,krow) = zrcst*vom(jl,jlev,krow)
    END DO
  END DO
  !
  ! Potential temperature
  !
  DO jlev=1,nlev
    DO jl=1,nproma
      tpot(jl,jlev,krow) = tm1(jl,jlev,krow)*(1.E05_dp/apm1(jl,jlev))**(rd/cpd)
    END DO
  END DO
  !
  locond=lcond
  loconv=lconv

  IF (lnmi_run) THEN
    ! control cloud parameterisation in nmi initialization mode
    locond=lnmi_cloud
    loconv=lnmi_cloud
    SELECT CASE(nmi_phase)
    CASE(NMI_ACCU)     ! store tendencies
!DIR$ CONCURRENT
      buf_t(:,:) = tte(:,:,krow)
!DIR$ CONCURRENT
      buf_m(:,:) = vom(:,:,krow)
!DIR$ CONCURRENT
      buf_l(:,:) = vol(:,:,krow)
    CASE(NMI_USE_AVG)  ! store tendencies
!DIR$ CONCURRENT
      buf_t(:,:) = tte(:,:,krow)
!DIR$ CONCURRENT
      buf_m(:,:) = vom(:,:,krow)
!DIR$ CONCURRENT
      buf_l(:,:) = vol(:,:,krow)
    CASE default
    END SELECT
  END IF
  !
  ! ------------------------------------------------------------------
  !
  !*        3.5   ESTIMATE ADIABATIC CONVERSION OF POTENTIAL ENERGY.
  !
  zdadc(1:nproma)=0._dp
  !
  DO 351 jl=1,nproma
    zdpsdt(jl)=aphm1(jl,nlevp1)*alpste(jl,krow)
    zdadc(jl)=zdadc(jl)+geospm(jl,krow)*zdpsdt(jl)
351 END DO
  !
  DO 353 jlev=1,nlev
    DO 352 jl=1,nproma
      zdadc(jl)=zdadc(jl)+(1._dp+vtmpc2*qm1(jl,jlev,krow))*cpd*              &
           (tte(jl,jlev,krow)*(aphm1(jl,jlev+1)-aphm1(jl,jlev))           &
           +tm1(jl,jlev,krow)*delb(jlev)*zdpsdt(jl))
352 END DO
353 END DO
  !
  !
  !*        3.6   COMPUTE LOGICAL MASK FOR LAND AND GLACIER.
  !
  ! special handling for fractional/non-fractional masks
  !
  IF (lfractional_mask) THEN
     DO jl=1,nproma
       loland(jl)=slf(jl,krow).GT.0._dp
       lolake(jl)=alake(jl,krow).GT.0._dp
       zfrl(jl)=slf(jl,krow)
     END DO
  ELSE
     DO jl=1,nproma
       loland(jl)=slm(jl,krow).GT.0._dp
       lolake(jl)=alake(jl,krow).GE.0.5_dp
       zfrl(jl)=slm(jl,krow)
     END DO
  ENDIF
  !
  DO 365 jl=1,nproma
    loglac(jl)=loland(jl).AND.glac(jl,krow).GT.0._dp
    lonorth(jl)=philat_2d(jl,krow).GT.0._dp ! true in northern hemisphere
365 END DO
  !
  !       3.7 Weighting factors for fractional surface coverage
  !           Accumulate ice portion for diagnostics
  !
!DIR$ CONCURRENT
  DO jl=1,nproma
!    zfrl(jl)=slf(jl,krow) ! see special handling !
    zfrw(jl)=(1._dp-zfrl(jl))*(1._dp-seaice(jl,krow))
    zfri(jl)=1._dp-zfrl(jl)-zfrw(jl)
!!$      friac(jl,krow)=friac(jl,krow)+zdtime*zfri(jl)
    !! friac now directly imported from mo_surface into g3b stream
  END DO
  IF(lcouple) THEN
    DO jl=1,nproma
      IF(slf(jl,krow).GT.1.0_dp-zepsec) THEN
        tsi(jl,krow)=tmelt
        tsw(jl,krow)=tmelt
      END IF
    END DO
  ENDIF
  !
  !      3.8  Skin reservoir, wet skin fraction and snow cover
  !           (bare land, canopy, lake ice)
  !
!IBM* novector
  DO jl=1,nproma
    IF (.NOT.loglac(jl)) THEN
      zwlmx(jl)=cwlmax*(1._dp+vlt(jl,krow))
      zcvw(jl)=MIN(wl(jl,krow)/zwlmx(jl),1.0_dp)
      zsn_mm=1000._dp*sn(jl,krow)
      zsigh=SQRT(zsn_mm/(zsn_mm+zepsec+zsigfac*orostd(jl,krow)))
      zcvs(jl)=cqsncr*TANH(zsn_mm/10._dp)*zsigh
      zcvsc(jl)=MIN(1._dp,snc(jl,krow)/(zwlmx(jl)-cwlmax+EPSILON(1._dp)))
      IF (zcvs(jl).LT.EPSILON(1._dp) .AND. zcvsc(jl).GE.EPSILON(1._dp)) THEN
        zcvs(jl)=zcvsc(jl)
      END IF
    ELSE
      zwlmx(jl)=0._dp
      zcvw(jl)=0._dp
      zcvs(jl)=1._dp
      zcvsc(jl)=0._dp
    END IF
  END DO
  !
  !
  !*         3.8    SET LOOP  VALUES FOR PHYSICAL PARAMETERIZATIONS
  !
  kidia=1
  kfdia=nproma
  ktdia=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! following is a submodel interface !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                                            !!
  !! call to submodels that need to update before radiation                     !!
  !!                                                                            !!
  !! a) set diagnostic tracers for SOA and calculate gas/aerosol equilibrium    !!
  !!                                                                            !!
  !!                                                                            !!
  !IF (ltimer) call timer_start(timer_physci_1)                                 !!
  !!                                                                            !!
  IF (lanysubmodel) THEN
    CALL physc_subm_1 (      &
         nproma, nbdim, nlev,    &
         nlevp1, ntrac,          &
         krow,                   &
         apm1(:,:),              &
         aphm1(:,:),             &
         tm1(:,:,krow),          &
         tte(:,:,krow),          &
         xtm1(:,:,:,krow),       &
         xtte(:,:,:,krow),       &
         qm1(:,:,krow),          &
         qte(:,:,krow)           )
  END IF
  !!                                                                            !!
  !IF (ltimer) call timer_stop(timer_physci_1)                                  !!
  !ENDIF                                                                        !!
  !!                                                                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !
  !         3.9   DETERMINE TROPOPAUSE HEIGHT AND MASS BUDGETS
  !
  !
  CALL WMO_tropopause (nproma, nbdim, nlev,                                    &
       tm1(:,:,krow),  apm1, tropo(:,krow),                    &
       itrpwmo, itrpwmop1)
  !
  !
  !*    3.12 INITIALISATION OF CLOUD DROPLET NUMBER CONCENTRATION 
  !          (1/M**3) USED IN RADLSW AND CLOUD 
  !
  IF (lstart) THEN
    DO 4103 jk=ktdia,nlev        
!DIR$ CONCURRENT
      DO 4102 jl=kidia,kfdia       
        nexp=2
        zprat=(MIN(8._dp,80000._dp/apm1(jl,jk)))**nexp
        IF (loland(jl).AND.(.NOT.loglac(jl)).OR.lolake(jl)) THEN
          zn1= 20._dp
          zn2=180._dp
        ELSE 
          zn1= 20._dp
          zn2= 80._dp
        ENDIF
        IF (apm1(jl,jk).LT.80000._dp) THEN
          zcdnc=1.e6_dp*(zn1+(zn2-zn1)*(EXP(1._dp-zprat)))
        ELSE
          zcdnc=zn2*1.e6_dp
        ENDIF
        acdnc(jl,jk,krow)=zcdnc
        acdncm(jl,jk,krow)=zcdnc
4102  END DO
4103 END DO
  ENDIF

!*    3.12a Initialisation of ice nuclei number concentration (1/M**3) acc. to Pruppacher/Klett (1997)
!           See Seinfeld/Pandis p.808

  IF (nic_cirrus>0) THEN
    IF (lstart) THEN
      zicnc(1:nproma,1:nlev)=1.e3_dp*exp(0.6_dp*(253._dp-tm1(1:nproma,1:nlev,krow)))
    ELSE
      zicnc(1:nproma,1:nlev)=icnc_instantan(1:nproma,1:nlev,krow)
    END IF
  ELSE
    zicnc(1:nproma,1:nlev)=0._dp
  END IF

  !
  DO 4104 jl=kidia,kfdia       
    itype(jl)=NINT(rtype(jl,krow))
4104 END DO
  !
  DO 4106 jk=ktdia,nlev
    DO 4105 jl=kidia,kfdia
      zqtec(jl,jk)=0.0_dp
4105 END DO
4106 END DO
  !
  !*        3.13   DIAGNOSE CURRENT CLOUD COVER
  !
  IF(locond) THEN

    IF (ltimer) CALL timer_start(timer_cover)

    CALL cover( nproma, nbdim, ktdia, nlev, nlevp1                            &
         , itype,             zfrw,             zfri                     &
         , aphm1,             apm1,             geom1                    &
         , tm1(:,:,krow),     qm1(:,:,krow),    xim1(:,:,krow)           &
         , aclc(:,:,krow)                                                &
         , invb,              rintop(:,krow)                             &
         )

    IF (ltimer) CALL timer_stop(timer_cover)

  ENDIF
  !
  !*        4.    RADIATION PARAMETERISATION.
  !               --------- -----------------
  !
  ! If radiation is active (lrad == True) then call of radiation computation on radiation timestep (l_trigrad == True)
  ! Shortwave fluxes are normalized by the solar constant and zenith angle at the time at which radiation is calculated.
  ! These are then renormalized by the value of the solar constant and zenith angle each time the heating rates are
  ! calculated
  !
  IF (lrad) THEN
    IF (l_trigrad) THEN
      IF (ltimer) CALL timer_start(timer_radiation)
      CALL radiation( &
           &  nproma                 ,nbdim                   ,nlev                   ,nlevp1                   &
           & ,krow                   ,ntrac                   ,itype                  ,loland                   &
           & ,loglac                 ,albedo_vis_dir(:,krow)  ,albedo_nir_dir(:,krow) ,albedo_vis_dif(:,krow)   &
           & ,albedo_nir_dif(:,krow) ,radtemp(:,krow)         ,aphm1                  ,apm1                     &
           & ,tm1(:,:,krow)          ,qm1(:,:,krow)           ,xlm1(:,:,krow)         ,xim1(:,:,krow)           &
!>>SF
           & ,geom1                  ,co2m1(:,:,krow)         ,acdnc(:,:,krow) , zicnc(:,:)       ,aclc(:,:,krow)           &
!           & ,geom1                  ,co2m1(:,:,krow)         ,acdnc(:,:,krow)        ,aclc(:,:,krow)           &
!<<SF
           & ,xtm1(:,:,:,krow)       ,aclcv(:,krow)           ,sswvis_frc(:,krow)     ,sswpar(:,krow)           &
           & ,sswdifnir(:,krow)      ,sswdifvis(:,krow)       ,sswdifpar(:,krow)      ,emtef(:,:,krow)          &
           & ,trsof(:,:,krow)        ,emtef0(:,:,krow)        ,trsof0(:,:,krow)       ,emter(:,:,krow)          &
           & ,trsol(:,:,krow)        ,ao3(:,:,krow)                                                             )
      IF (ltimer) CALL timer_stop(timer_radiation)
    END IF
    !
    DO  jl=1,nproma
      zi0(jl)=flx_ratio_cur*solc*amu0_x(jl,krow)*rdayl_x(jl,krow)
    END DO
  ELSE
    !
    ! No radiative fluxes, which implies no radiative effect of the atmosphere or the surface
    !
    emter(:,:,krow)  = 0._dp
    trsol(:,:,krow)  = 1._dp
    emtef(:,:,krow)  = 0._dp
    trsof(:,:,krow)  = 1._dp
    emtef0(:,:,krow) = 0._dp
    trsof0(:,:,krow) = 1._dp

    DO  jl=1,nproma
      zi0(jl)=0._dp
    END DO
  END IF

  IF ( locfdiag ) THEN
    CALL calc_cfdiag( nproma, nbdim, zi0, krow )
  END IF
  !
  !       Compute diffuse, direct, NIR and visible fluxes for JSBACH
  !
  DO  jl=1,nproma
    jsswvis(jl,krow) = zi0(jl) * trsol(jl,nlevp1,krow) * sswvis_frc(jl,krow) 
    jsswnir(jl,krow) = zi0(jl) * trsol(jl,nlevp1,krow) - jsswvis(jl,krow)
    jsswpar(jl,krow) = zi0(jl) * sswpar(jl,krow)

    jsswdifnir(jl,krow) = sswdifnir(jl,krow)
    jsswdifvis(jl,krow) = sswdifvis(jl,krow)
    jsswdifpar(jl,krow) = sswdifpar(jl,krow)

    jsswniracc(jl,krow)    = jsswniracc(jl,krow) + jsswnir(jl,krow) * zdtime
    jsswvisacc(jl,krow)    = jsswvisacc(jl,krow) + jsswvis(jl,krow) * zdtime
    jsswparacc(jl,krow)    = jsswparacc(jl,krow) + jsswpar(jl,krow) * zdtime
    jsswdifparacc(jl,krow) = jsswdifparacc(jl,krow) + jsswdifpar(jl,krow)  &
         * jsswpar(jl,krow) * zdtime
    jsswdifniracc(jl,krow) = jsswdifniracc(jl,krow) + jsswdifnir(jl,krow)  &
         * jsswnir(jl,krow) * zdtime
    jsswdifvisacc(jl,krow) = jsswdifvisacc(jl,krow) + jsswdifvis(jl,krow)  &
         * jsswvis(jl,krow) * zdtime
  END DO

  !
  !     ------------------------------------------------------------
  !
  !*              VERTICAL EXCHANGE OF U,V,T,Q BY TURBULENCE.
  !               -------- -------- -- - - - - -- -----------
  !
  !
  !      COMPUTE PRESSURE AT FULL AND HALF LEVELS AT T+DT.
  !
  ztwodt=time_step_len
  DO 522 jl=1,nproma
    aphp1(jl,nlevp1)=EXP(alpsm1(jl,krow)+ztwodt*alpste(jl,krow))
522 END DO
  !
  CALL pres(aphp1,nbdim,aphp1(1,nlevp1),nproma)
  !
  CALL presf(app1,nbdim,aphp1,nproma)
  !

  IF(lcouple_co2) CALL co2_flux_atmosphere_ocean(krow, nproma, nlev)

  IF (ltimer) CALL timer_start(timer_vdiff)

  IF(ltdiag) THEN
    IF (ASSOCIATED(tdiag_vars%dudt_vdiff)) THEN
      CALL set_tendency(tdiag_vars%dudt_vdiff(:,:,krow)  ,vom(:,:,krow)   ,nproma &
           ,nbdim                            ,nlev            ,'sub'  )
    END IF
    IF (ASSOCIATED(tdiag_vars%dvdt_vdiff)) THEN
      CALL set_tendency(tdiag_vars%dvdt_vdiff(:,:,krow)  ,vol(:,:,krow)   ,nproma &
           ,nbdim                            ,nlev            ,'sub'  )
    END IF
    IF (ASSOCIATED(tdiag_vars%dtdt_vdiff)) THEN
      CALL set_tendency(tdiag_vars%dtdt_vdiff(:,:,krow)  ,tte(:,:,krow)   ,nproma &
           ,nbdim                            ,nlev            ,'sub'  )
    END IF
    IF (ASSOCIATED(tdiag_vars%dqdt_vdiff)) THEN
      CALL set_tendency(tdiag_vars%dqdt_vdiff(:,:,krow)  ,qte(:,:,krow)   ,nproma &
           ,nbdim                            ,nlev            ,'sub'  )
    END IF
    IF (ASSOCIATED(tdiag_vars%dxldt_vdiff)) THEN
      CALL set_tendency(tdiag_vars%dxldt_vdiff(:,:,krow) ,xlte(:,:,krow)  ,nproma &
           ,nbdim                            ,nlev            ,'sub'  )
    END IF
    IF (ASSOCIATED(tdiag_vars%dxidt_vdiff)) THEN
      CALL set_tendency(tdiag_vars%dxidt_vdiff(:,:,krow) ,xite(:,:,krow)  ,nproma &
           ,nbdim                            ,nlev            ,'sub'  )
    END IF
  END IF

  CALL vdiff (nproma, nbdim, ktdia, nlev, nlevm1, nlevp1, ntrac                &
       , krow                                                             &
       , xtm1(:,:,:,krow)                                                 &
       ! contains also CO2 concentration
       , qm1(:,:,krow),        tm1(:,:,krow),        um1(:,:,krow)        &
       , vm1(:,:,krow),        xlm1(:,:,krow),       xim1(:,:,krow)       &
       , ahfl(:,krow),         ahfs(:,krow),         az0(:,krow)          &
       , dew2(:,krow),         evap(:,krow),         forest(:,krow)       &
       , temp2(:,krow),        t2max(:,krow)                              &
       , t2min(:,krow),        wind10w(:,krow),      vdis(:,krow)         &
       , u10(:,krow),          v10(:,krow),          ustr(:,krow)         &
       , vstr(:,krow),         wimax(:,krow),        wind10(:,krow)       &
       , vgrat(:,krow),        vlt(:,krow),          zcvw(:)              &
       , tsw(:,krow),          tsi(:,krow)                                &
       , ocu(:,krow),          ocv(:,krow)                                &
       , az0l(:,krow),         az0w(:,krow),         az0i(:,krow)         &
       , zhfsl,                zhfsw,                zhfsi                &
       , zhflw,                zevapw,               zevapi               &
       , ahfslac(:,krow),      ahfswac(:,krow),      ahfsiac(:,krow)      &
       , ahfllac(:,krow),      ahflwac(:,krow),      ahfliac(:,krow)      &
       , evaplac(:,krow),      evapwac(:,krow),      evapiac(:,krow)      &
       , ustrl(:,krow),        ustrw(:,krow),        ustri(:,krow)        &
       , vstrl(:,krow),        vstrw(:,krow),        vstri(:,krow)        &
       , albedo(:,krow),       albedo_vis(:,krow)                         &
       , albedo_nir(:,krow),   alsol(:,krow)                              &
       , albedo_vis_dir(:,krow),     albedo_nir_dir(:,krow)               &
       , albedo_vis_dif(:,krow),     albedo_nir_dif(:,krow)               &
       , zpbl(:)                                                          &
       , tke(:,:,krow),        tkem1(:,:,krow),      tkem(:,:,krow)       &
       , aclc(:,:,krow),       emter(:,:,krow)                            &
       , thvvar(:,:,krow),     thvsig(:,krow)                             &
       , sh_vdiff(:,krow),     ev_vdiff(:,krow)                           &
       , aphm1,                apm1,                 geom1                &
       , ztvm1                                                            &
       , qhfla,                evapot                                     &
       , ztslnew,              zfrl,                 loland               &
       , lonorth,              zfrw,                 zfri                 &
       , alake(:,krow),        xtte(:,:,:,krow)                           &
       , vol(:,:,krow),        vom(:,:,krow),        qte(:,:,krow)        &
       , tte(:,:,krow),        xlte(:,:,krow),       xite(:,:,krow)       &
       , siced(:,krow)                                                    &
       ! Rain and snow over last time step
       , zrain(:),             zsnow(:)                                   &
       , jsswnir(:,krow),      jsswdifnir(:,krow)                         &
       , jsswvis(:,krow),      jsswdifvis(:,krow)                         &
       , jsswpar(:,krow),      jsswdifpar(:,krow)                         &
       , zi0,                  trsol(:,:,krow)                            &
       , ztrflw,               ztrfli(:,krow),       zhfli                &
       , zsofll,               zsoflw,               zsofli(:,krow)       &
       , trfllac(:,krow),      trflwac(:,krow),      trfliac(:,krow)      &
       , sofllac(:,krow),      soflwac(:,krow),      sofliac(:,krow)      &
       , alsoi(:,krow),        alsow(:,krow),        radtemp(:,krow)      &
       , alsobs(:,krow),       alsom(:,krow),        taus(:,krow)         &
       , fage(:,krow),         snifrac(:,krow),      barefrac(:,krow)     &
       , ameltdepth(:,krow),   ameltfrac(:,krow)                          &
       , amlcorr(:,krow),      amlcorac(:,krow),     amlheatac(:,krow)    &
       , sni(:,krow),          ahfice(:,krow),       fluxres(:,krow)      &
       , qres(:,krow),         ahfcon(:,krow),       ahfres(:,krow)       &
       , zteffl4,              ztsnew,               tsurf(:,krow)        &
       , zradtemp_old                                                     &
       , aphp1                                                            &
       , seaice(:,krow)      )
  ! 
  IF (ltimer) call timer_stop(timer_vdiff)
  IF(ltdiag) THEN
    IF (ASSOCIATED(tdiag_vars%dudt_vdiff)) THEN
      CALL set_tendency(tdiag_vars%dudt_vdiff(:,:,krow)  ,vom(:,:,krow)   ,nproma &
           ,nbdim                            ,nlev            ,'add'  )
    END IF
    IF (ASSOCIATED(tdiag_vars%dvdt_vdiff)) THEN
      CALL set_tendency(tdiag_vars%dvdt_vdiff(:,:,krow)  ,vol(:,:,krow)   ,nproma &
           ,nbdim                            ,nlev            ,'add'  )
    END IF
    IF (ASSOCIATED(tdiag_vars%dtdt_vdiff)) THEN
      CALL set_tendency(tdiag_vars%dtdt_vdiff(:,:,krow)  ,tte(:,:,krow)   ,nproma &
           ,nbdim                            ,nlev            ,'add'  )
    END IF
    IF (ASSOCIATED(tdiag_vars%dqdt_vdiff)) THEN
      CALL set_tendency(tdiag_vars%dqdt_vdiff(:,:,krow)  ,qte(:,:,krow)   ,nproma &
           ,nbdim                            ,nlev            ,'add'  )
    END IF
    IF (ASSOCIATED(tdiag_vars%dxldt_vdiff)) THEN
      CALL set_tendency(tdiag_vars%dxldt_vdiff(:,:,krow) ,xlte(:,:,krow)  ,nproma &
           ,nbdim                            ,nlev            ,'add'  )
    END IF
    IF (ASSOCIATED(tdiag_vars%dxidt_vdiff)) THEN
      CALL set_tendency(tdiag_vars%dxidt_vdiff(:,:,krow) ,xite(:,:,krow)  ,nproma &
           ,nbdim                            ,nlev            ,'add'  )
    END IF
  END IF

  !
  !------------------------------------------------------------------------------
  !
  !*            ADD RADIATION TENDENCIES EVERY TIME STEP.
  !
  IF (lrad) THEN
    IF (ltimer) call timer_start(timer_radheat)
    CALL radheat (nproma, nbdim, nlev, nlevp1,                                   &
         krow,                                                                &
         zi0,                                                                 &
         tm1(:,:,krow)  ,   qm1(:,:,krow),                                    &
         trsof(:,:,krow),   trsol(:,:,krow),                                  &
         emtef(:,:,krow),   emter(:,:,krow),                                  &
         emtef0(:,:,krow),  trsof0(:,:,krow),                                 &
         srad0(:,krow),     srads(:,krow),                                    &
         sradl(:,krow),     srafl(:,krow),                                    &
         srad0u(:,krow),    sradsu(:,krow),                                   &
         sraf0(:,krow),     srafs(:,krow),                                    &
         srad0d(:,krow),                                                      &
         trad0(:,krow),     trads(:,krow),                                    &
         tradl(:,krow),     trafl(:,krow),                                    &
         traf0(:,krow),     trafs(:,krow),                                    &
         tradsu(:,krow),                                                      &
         albedo(:,krow),                                                      &
         aphm1,             apm1,                                             &
         tte(:,:,krow),                                                       &
         zradtemp_old,      ztsnew  )
    IF (ltimer) call timer_stop(timer_radheat)

  ELSE
    ! lrad=.FALSE.
    ! --> no radiative effect at the surface
    ztrflw(:)=0._dp
    ztrfli(:,krow)=0._dp
    zsofll(:)=0._dp
    zsoflw(:)=0._dp
    zsofli(:,krow)=0._dp
    ! --> lake, ml_ocean, licetemp and sicetemp get zero fluxes 
  END IF


  !
  !     ------------------------------------------------------------
  !
  !        ***  GRAVITY WAVE DRAG PARAMETERISATION  ***
  !
  IF (lmidatm) THEN
    IF (lresume) THEN
      aprflux(1:nproma,krow) = 0._dp
      aprfluxm(1:nproma,krow) = 0._dp
    ENDIF

    IF(ltdiag) THEN
      IF (ASSOCIATED(tdiag_vars%dudt_hines)) THEN
        CALL set_tendency(tdiag_vars%dudt_hines(:,:,krow)  ,vom(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'sub'  )
      END IF
      IF (ASSOCIATED(tdiag_vars%dvdt_hines)) THEN
        CALL set_tendency(tdiag_vars%dvdt_hines(:,:,krow)  ,vol(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'sub'  )
      END IF
      IF (ASSOCIATED(tdiag_vars%dtdt_hines)) THEN
        CALL set_tendency(tdiag_vars%dtdt_hines(:,:,krow)  ,tte(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'sub'  )
      END IF
    END IF
    CALL gwspectrum ( krow, nproma,  nbdim,  nlev,                             &
         aphm1,            apm1,                                      &
         tm1(:,:,krow),    um1(:,:,krow),   vm1(:,:,krow),            &
         aprflux(:,krow),                                             &
         tte(:,:,krow),    vol(:,:,krow),   vom(:,:,krow) )
    IF(ltdiag) THEN
      IF (ASSOCIATED(tdiag_vars%dudt_hines)) THEN
        CALL set_tendency(tdiag_vars%dudt_hines(:,:,krow)  ,vom(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'add'  )
      END IF
      IF (ASSOCIATED(tdiag_vars%dvdt_hines)) THEN
        CALL set_tendency(tdiag_vars%dvdt_hines(:,:,krow)  ,vol(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'add'  )
      END IF
      IF (ASSOCIATED(tdiag_vars%dtdt_hines)) THEN
        CALL set_tendency(tdiag_vars%dtdt_hines(:,:,krow)  ,tte(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'add'  )
      END IF
    END IF
    !  
  END IF


  IF (lgwdrag) THEN
    IF (ltimer) call timer_start(timer_gwdrag)

    zlat(:) = philat_2d(:,krow)/180._dp*pi

    IF(ltdiag) THEN
      IF (ASSOCIATED(tdiag_vars%dudt_sso)) THEN
        CALL set_tendency(tdiag_vars%dudt_sso(:,:,krow)    ,vom(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'sub'  )
      END IF
      IF (ASSOCIATED(tdiag_vars%dvdt_sso)) THEN
        CALL set_tendency(tdiag_vars%dvdt_sso(:,:,krow)    ,vol(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'sub'  )
      END IF
      IF (ASSOCIATED(tdiag_vars%dtdt_sso)) THEN
        CALL set_tendency(tdiag_vars%dtdt_sso(:,:,krow)    ,tte(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'sub'  )
      END IF
    END IF
    CALL ssodrag (nproma,  nbdim,  nlev,  zlat(:),   time_step_len,            &
         aphm1,            apm1,            geom1,                    &
         tm1(:,:,krow),    um1(:,:,krow),   vm1(:,:,krow),            &
         oromea(:,krow),   orostd(:,krow),  orosig(:,krow),           &
         orogam(:,krow),                                              &
         orothe(:,krow),   oropic(:,krow),  oroval(:,krow),           &
         ustrgw(:,krow),   vstrgw(:,krow),  vdisgw(:,krow),           &
         tte(:,:,krow),    vol(:,:,krow),   vom(:,:,krow) )

    IF(ltdiag) THEN
      IF (ASSOCIATED(tdiag_vars%dudt_sso)) THEN
        CALL set_tendency(tdiag_vars%dudt_sso(:,:,krow)    ,vom(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'add'  )
      END IF
      IF (ASSOCIATED(tdiag_vars%dvdt_sso)) THEN
        CALL set_tendency(tdiag_vars%dvdt_sso(:,:,krow)    ,vol(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'add'  )
      END IF
      IF (ASSOCIATED(tdiag_vars%dtdt_sso)) THEN
        CALL set_tendency(tdiag_vars%dtdt_sso(:,:,krow)    ,tte(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'add'  )
      END IF
    END IF
    IF (ltimer) call timer_stop(timer_gwdrag)
  END IF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! following is a submodel interface !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (lanysubmodel) THEN
    CALL physc_subm_2( &
         nproma, nbdim, nlev,     &
         nlevp1, ntrac, krow,     &
         itrpwmo,                 &
         itrpwmop1,               &
         aphm1(:,:),              &
         apm1(:,:),               &
         aphp1(:,:),              &
         app1(:,:),               &
         tm1(:,:,krow),           &
         tte(:,:,krow),           &
         tsurf(:,krow),           &
         qm1(:,:,krow),           &
         qte(:,:,krow),           &
         xlm1(:,:,krow),          &
         xlte(:,:,krow),          &
         xim1(:,:,krow),          &
         xite(:,:,krow),          &
         xtm1(:,:,:,krow),        &
         xtte(:,:,:,krow),        &
         aclc(:,:,krow),          &
         zpbl(:),                 &
         loland(:),               &
         loglac(:)              )
  END IF
  !
  !     ------------------------------------------------------------------
  !
  !*        6.    CONVECTION PARAMETERISATION.
  !               ---------- -----------------
  !
  !*        6.3    COMPUTE *T* AND *Q* TENDENCIES BY MOIST CONVECTION.
  !*                AND SHALLOW CONVECTION.
  !
  itype(1:nproma)=0  !!!!!
  !
  !
  !
  !*         6.3.1   INITIALIZE ARRAYS FOR CONVECTIVE PRECIPITATION
  !*                 AND COPY ARRAYS FOR CONVECTIVE CLOUD PARAMETERS
  !*                 -----------------------------------------------
  !
  zxtecl(1:nproma,:) = 0._dp
  zxteci(1:nproma,:) = 0._dp

  !
  DO 632 jl=1,nproma
    rsfc(jl)=0._dp
    ssfc(jl)=0._dp
632 END DO
  !
  !
  !*         6.3.3   CALL SUBROUTINE CUCALL FOR CUMULUS PARAMETERIZATION
  !                  ---------------------------------------------------
  !


  IF (loconv) THEN
    IF (ltimer) call timer_start(timer_cucall)
    !
    IF(ltdiag) THEN
      IF (ASSOCIATED(tdiag_vars%dudt_cucall)) THEN
        CALL set_tendency(tdiag_vars%dudt_cucall(:,:,krow) ,vom(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'sub'  )
      END IF
      IF (ASSOCIATED(tdiag_vars%dvdt_cucall)) THEN
        CALL set_tendency(tdiag_vars%dvdt_cucall(:,:,krow) ,vol(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'sub'  )
      END IF
      IF (ASSOCIATED(tdiag_vars%dtdt_cucall)) THEN
        CALL set_tendency(tdiag_vars%dtdt_cucall(:,:,krow) ,tte(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'sub'  )
      END IF
      IF (ASSOCIATED(tdiag_vars%dqdt_cucall)) THEN
        CALL set_tendency(tdiag_vars%dqdt_cucall(:,:,krow) ,qte(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'sub'  )
      END IF
    END IF

    CALL cucall(nproma, nbdim, nlev, nlevp1, nlevm1,                    &!in
         ntrac,                                                         &!in
         time_step_len,                                                 &!in
         loland,                                                        &!in
         tm1(:,:,krow),    um1(:,:,krow),    vm1(:,:,krow),             &!in
         qm1(:,:,krow),    xlm1(:,:,krow),   xim1(:,:,krow),            &!in
         xtm1(:,:,:,krow),                                              &!in
         qte(:,:,krow),    xlte(:,:,krow),   xite(:,:,krow),            &!in
         vervel(:,:,krow), qhfla,            geom1,                     &!in
         apm1,             aphm1,            thvsig(:,krow),            &!in
         tte(:,:,krow),    vom(:,:,krow),    vol(:,:,krow),             &!in
         xtte(:,:,:,krow),                                              &!in
         zqtec,                                                         &!inout
         ch_concloud(:,krow), cw_concloud(:,krow),                      &!inout
         rsfc,             ssfc,                                        &!out
         zxtecl,           zxteci,                                      &!out(?)
         itype,            ictop,            ilab,                      &!out
         krow,                                                          &!in
         aprc(:,krow),     aprs(:,krow),                                &!out(?)
         topmax(:,krow)                                                ) !inout
    !
      IF(ltdiag) THEN
      IF (ASSOCIATED(tdiag_vars%dudt_cucall)) THEN
        CALL set_tendency(tdiag_vars%dudt_cucall(:,:,krow) ,vom(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'add'  )
      END IF
      IF (ASSOCIATED(tdiag_vars%dvdt_cucall)) THEN
        CALL set_tendency(tdiag_vars%dvdt_cucall(:,:,krow) ,vol(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'add'  )
      END IF
      IF (ASSOCIATED(tdiag_vars%dtdt_cucall)) THEN
        CALL set_tendency(tdiag_vars%dtdt_cucall(:,:,krow) ,tte(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'add'  )
      END IF
      IF (ASSOCIATED(tdiag_vars%dqdt_cucall)) THEN
        CALL set_tendency(tdiag_vars%dqdt_cucall(:,:,krow) ,qte(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'add'  )
      END IF
    END IF

    IF (ltimer) call timer_stop(timer_cucall)
  ELSE
    !       NECESSARY COMPUTATIONS IF MASSFLUX IS BY-PASSED
    !
    ilab(1:nproma,1:nlev)=0
    ictop(1:nproma)      =nlev-1
    !
  ENDIF

  !
  !     ------------------------------------------------------------
  !
  !*       7.    LARGE SCALE CONDENSATION.
  !              ----- ----- -------------
  !
  IF(locond) THEN
    !
    IF (ltimer) CALL timer_start(timer_cloud)

    IF(ltdiag) THEN
      IF (ASSOCIATED(tdiag_vars%dtdt_cloud)) THEN
        CALL set_tendency(tdiag_vars%dtdt_cloud(:,:,krow)  ,tte(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'sub'  )
      END IF
      IF (ASSOCIATED(tdiag_vars%dqdt_cloud)) THEN
        CALL set_tendency(tdiag_vars%dqdt_cloud(:,:,krow)  ,qte(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'sub'  )
      END IF
      IF (ASSOCIATED(tdiag_vars%dxldt_cloud)) THEN
        CALL set_tendency(tdiag_vars%dxldt_cloud(:,:,krow) ,xlte(:,:,krow)  ,nproma &
             ,nbdim                            ,nlev            ,'sub'  )
      END IF
      IF (ASSOCIATED(tdiag_vars%dxidt_cloud)) THEN
        CALL set_tendency(tdiag_vars%dxidt_cloud(:,:,krow) ,xite(:,:,krow)  ,nproma &
             ,nbdim                            ,nlev            ,'sub'  )
      END IF
    END IF

    IF (.NOT. lcdnc_progn) THEN

       CALL cloud( nproma, nbdim, ktdia, nlev, nlevp1                       & ! in
            , delta_time,        time_step_len                              & ! in, time steps
            , ntrac,             krow                                       & ! in
            , invb,              ictop                                      & ! in
            , aphm1,             vervel(:,:,krow)                           & ! in
            , apm1,              app1,             acdnc(:,:,krow)          & ! in
            , qm1(:,:,krow),     tm1(:,:,krow),    ztvm1                    & ! in
            , xlm1(:,:,krow),    xim1(:,:,krow)                             & ! in
            , zcair(:,:),        geom1,            aphp1                    & ! in
            , xtm1(:,:,:,krow)                                              & ! in
            , aclcov(:,krow),    aprl(:,krow),     qvi(:,krow)              & ! inout
            , xlvi(:,krow),      xivi(:,krow)                               & ! inout
            , aprs(:,krow),      itype                                      & ! inout
            , ch_concloud(:,krow), cw_concloud(:,krow)                      & ! inout
            , zxtecl,            zxteci,           zqtec                    & ! inout
            , qte(:,:,krow),     tte(:,:,krow)                              & ! inout
            , xlte(:,:,krow),    xite(:,:,krow)                             & ! inout
            , xtte(:,:,:,krow)                                              & ! inout
            , aclc(:,:,krow),    aclcac(:,:,krow)                           & ! inout
            , ssfl,              rsfl                                       & ! out
            , relhum(:,:,krow)   )                                            ! out

    ELSE 

       CALL cloud_micro_interface(nproma, nbdim, nlev, nlevp1, ntrac, ktdia, krow,       & ! in
                                  invb,                                                  & ! in
                                  aphm1,           apm1,              app1,              & ! in
                                  qm1(:,:,krow),   tm1(:,:,krow),     ztvm1,             & ! in
                                  xlm1(:,:,krow),    xim1(:,:,krow),                     & ! in
                                  vervel(:,:,krow),  geom1,           xtm1(:,:,:,krow),  & ! in
                                  aphp1,           tkem1(:,:,krow),                      & ! in
                                  zqtec,             aclc(:,:,krow),                     & ! inout
                                  aclcac(:,:,krow),  aclcov(:,krow),  aprl(:,krow),      & ! inout
                                  qvi(:,krow),       xlvi(:,krow),    xivi(:,krow),      & ! inout
                                  qte(:,:,krow),     tte(:,:,krow),   xlte(:,:,krow),    & ! inout
                                  xite(:,:,krow),    xtte(:,:,:,krow),aprs(:,krow),      & ! inout
                                  acdnc(:,:,krow),                                       & ! inout
                                  zicnc(:,:),        relhum(:,:,krow),                   & ! out
                                  ssfl(:), rsfl(:))                                        ! out
        
       IF (nic_cirrus > 0) THEN
          icnc_instantan(1:nproma,1:nlev,krow)=zicnc(1:nproma,1:nlev)
       END IF


    END IF

    ! >>s.stadtler
    !--- save cloud condensation nuclei number concentration for
    !--- heterogeneous chemistry
    zxlp1(1:nproma,:) = xlm1(1:nproma,:,krow) + time_step_len * xlte(1:nproma,:,krow)
    CALL set_vphysc_var(nproma, nlev, krow, pcdnc=acdnc(1:nproma,:,krow))
    CALL set_vphysc_var(nproma, nlev, krow, pclw=zxlp1)
    CALL set_vphysc_var(nproma, nlev, krow, paclc=aclc(1:nproma,:,krow))
    ! <<s.stadtler

    DO jl=kidia,kfdia       
      rtype(jl,krow)=REAL(itype(jl),dp)
    END DO
    
    IF(ltdiag) THEN
      IF (ASSOCIATED(tdiag_vars%dtdt_cloud)) THEN
        CALL set_tendency(tdiag_vars%dtdt_cloud(:,:,krow)  ,tte(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'add'  )
      END IF
      IF (ASSOCIATED(tdiag_vars%dqdt_cloud)) THEN
        CALL set_tendency(tdiag_vars%dqdt_cloud(:,:,krow)  ,qte(:,:,krow)   ,nproma &
             ,nbdim                            ,nlev            ,'add'  )
      END IF
      IF (ASSOCIATED(tdiag_vars%dxldt_cloud)) THEN
        CALL set_tendency(tdiag_vars%dxldt_cloud(:,:,krow) ,xlte(:,:,krow)  ,nproma &
             ,nbdim                            ,nlev            ,'add'  )
      END IF
      IF (ASSOCIATED(tdiag_vars%dxidt_cloud)) THEN
        CALL set_tendency(tdiag_vars%dxidt_cloud(:,:,krow) ,xite(:,:,krow)  ,nproma &
             ,nbdim                            ,nlev            ,'add'  )
      END IF
    END IF

    IF (ltimer) CALL timer_stop(timer_cloud)

  ELSE
    !
    !              NECESSARY COMPUTATIONS IF *CLOUD* IS BY-PASSED.
    !
    ssfl(1:nproma) = 0._dp
    rsfl(1:nproma) = 0._dp
    !
    aclc(1:nproma,:,krow) = 0._dp
    !
  ENDIF
  !

  IF (lmidatm) THEN
    DO jl=1,nproma
      aprflux(jl,krow)=rsfl(jl)+ssfl(jl)+rsfc(jl)+ssfc(jl)
    END DO
  END IF
  DO jl=1,nproma
    zprecip(jl)=rsfl(jl)+ssfl(jl)+rsfc(jl)+ssfc(jl)
  END DO
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! following is a submodel interface !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (lanysubmodel) THEN
    CALL physc_subm_3( &
         nproma, nbdim, nlev,     &
         nlevp1, ntrac,           &
         krow,                    &
         aphm1(:,:),              &
         apm1(:,:),               &
         aphp1(:,:),              &
         app1(:,:),               &
         tm1(:,:,krow),           &
         tte(:,:,krow),           &
         tsurf(:,krow),           &
         qm1(:,:,krow),           &
         qte(:,:,krow),           &
         xlm1(:,:,krow),          &
         xlte(:,:,krow),          &
         xim1(:,:,krow),          &
         xite(:,:,krow),          &
         xtm1(:,:,:,krow),        &
         xtte(:,:,:,krow),        &
         geom1(:,:),              &
         geohm1(:,:),             &
         aclc(:,:,krow),          &
         zpbl(:),                 &
         vervel(:,:,krow),        &
!>>>mgs
         zfrl(:),                 &
         zfrw(:),                 &
         seaice(:, krow),         &       ! use seaice here, not zfri!
         glac(:, krow)            )
!<<<
  END IF


  !!                                                                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !*        8.    Computation of new surface values over land points
  !

  IF (lnmi_run) THEN
    SELECT CASE(nmi_phase)
    CASE(NMI_ACCU)    ! accumulate tendencies
!DIR$ CONCURRENT
      dh_t(:,:,krow) = dh_t(:,:,krow) + tte(:,:,krow) - buf_t(:,:)
!DIR$ CONCURRENT
      dh_m(:,:,krow) = dh_m(:,:,krow) + vom(:,:,krow) - buf_m(:,:)
!DIR$ CONCURRENT
      dh_l(:,:,krow) = dh_l(:,:,krow) + vol(:,:,krow) - buf_l(:,:)

    CASE(NMI_USE_AVG) ! prepare/reset tendencies
!DIR$ CONCURRENT
      tte(:,:,krow) = buf_t(:,:) + dh_t(:,:,krow)
!DIR$ CONCURRENT
      vom(:,:,krow) = buf_m(:,:) + dh_m(:,:,krow)
!DIR$ CONCURRENT
      vol(:,:,krow) = buf_l(:,:) + dh_l(:,:,krow)

    END SELECT
  END IF
  !
  !      p-e budget correction for coupling
  !
  DO jl=1,nproma
    zqtold(jl)=qtnew(jl,krow)
    qtnew(jl,krow)=0._dp
  END DO
  DO jk=ktdia,nlev
!DIR$ CONCURRENT
    DO jl=1,nproma
      qtnew(jl,krow)=qtnew(jl,krow)+(qm1(jl,jk,krow)                       &
           +xlm1(jl,jk,krow)+xim1(jl,jk,krow))                   &
           *(aphm1(jl,jk+1)-aphm1(jl,jk))/grav
    END DO
  END DO
  !     
  !     Accumulate p-e correction for standard diagnostics
  !
!DIR$ CONCURRENT
  DO jl=1,nproma
    apmeb(jl,krow)=apmeb(jl,krow)-(qtnew(jl,krow)-zqtold(jl))              &
         -(zprecip(jl)+qhfla(jl))*zdtime
  END DO
  !
  !     Vertical integral of anthropogenic sulfur burden 
  !
  IF(iaero == 4) THEN
    DO jk=ktdia,nlev
      DO jl=1,nproma
        abso4(jl,krow)=abso4(jl,krow)                                      &
             +(so4all(jl,jk,krow)-so4nat(jl,jk,krow))                   &
             *(aphm1(jl,jk+1)-aphm1(jl,jk))/grav*zdtime
      END DO
    END DO
  ENDIF
  !
  !
!!$ HILFESTELLUNG ZUR KOPPLUNG: kalle, 020904.
!!$      zhflw(jl) = ocean%latent_heat_flux_inst(1:nproma,krow)
!!$      zhfsw(jl) = ocean%sensible_heat_flux_inst(1:nproma,krow)
!!$      ahfice(jl,krow) = ice%ahfice(1:nproma,krow)
!!$      ztrflw(jl) = ocean%trflw(1:nproma,krow)
!!$      zsoflw(jl)= ocean%soflw(1:nproma,krow)
!!$      qres(jl,krow) = ice%qres(1:nproma,krow)
!!$      zevapw(jl) = ocean%evaporation_inst(1:nproma, krow)
!!$      zevapi(jl)= ice%evaporation_inst(1:nproma, krow)
!!$      ustrw(jl,krow) = ocean%u_stress(1:nproma,krow)
!!$      vstrw(jl,krow) = ocean%v_stress(1:nproma,krow)
!!$
  !
  IF(lmeltpond) THEN
    CALL meltpond_ice ( nproma , ameltdepth(:,krow)                           &
         , sicepdw(:,krow),  sicepdi(:,krow), tsicepdi(:,krow)                &
         , sicepres(:,krow), siced(:,krow),   zevapi                          &
         , zhflw,            zhfsw,           ztrflw,         zsoflw          &
         , zhfli,            zhfsi,           ztrfli(:,krow), zsofli(:,krow)  &
         , alsow(:,krow),   alsom(:,krow)                                     &
         , alsoi(:,krow) )
    CALL meltpond ( nproma                                                    &
         , sicepdw(:,krow),   sicepdi(:,krow),   tsicepdi(:,krow)             &
         , siced(:,krow),     seaice(:,krow)                                  &
         , sicemin(:,krow),   ameltdepth(:,krow),ameltfrac(:,krow)            &
         , qres(:,krow),      sni(:,krow) )
  END IF
  !
  IF (lcouple) THEN
    !
    IF (l_getocean) THEN
      apmebco(:,krow) = 0._dp ! set to zero after coupling
      rain(:,krow) = 0._dp    !       "         "
    ENDIF
    !
    !     Accumulate p-e correction for coupling
!DIR$ CONCURRENT
    DO jl=1,nproma
      apmebco(jl,krow)=apmebco(jl,krow)-(qtnew(jl,krow)-zqtold(jl))          &
           -(zprecip(jl)+qhfla(jl))*zdtime
      rain(jl,krow)=rain(jl,krow)+(rsfl(jl)+rsfc(jl))*zdtime
    END DO
    !
    !       collect data needed as input for the ocean model
    !
    CALL collect ( nproma                                                       &
         , zhflw,          zhfsw,         ahfice(:,krow)                        &
         , ztrflw,         zsoflw                                               &
         , qres(:,krow),   zevapw,        zevapi                                &
         , ustrw(:,krow),  vstrw(:,krow), ustri(:,krow), vstri(:,krow)          &
         , alake(:,krow),  slf(:,krow),   seaice(:,krow)                        &
         , wind10w(:,krow), co2m1(:,nlev,krow), co2_flux_ocean(:,krow)          &
         , awhea(:,krow),  awsol(:,krow), awust(:,krow)                         &
         , awvst(:,krow),  aicon(:,krow), aiqre(:,krow), aifre(:,krow)          &
         , aiust(:,krow),  aivst(:,krow), awsta(:,krow)                         &
         , co2atmos(:,krow), co2flux_cpl(:,krow)                                &
         , rsfc,           ssfc,          rsfl,          ssfl         )
  END IF
  !
  !    collect data needed (and modified) within the hydrological discharge model,
  !     that is also needed for the coupling
  !
  IF (lhd) THEN
    CALL hydrology_collect_lake ( nproma, krow, slf(:,krow), seaice(:,krow), alake(:,krow)  &
         , rsfc, ssfc, rsfl, ssfl, zevapw)
  END IF

  IF (lanysubmodel) THEN
    CALL physc_subm_4(nproma,                nbdim,                   nlev,             &
         nlevp1,                ntrac,                   krow,             &
         aphm1,                 zfrl,                    zfrw,             &
         zfri,                  loland,                  xtm1(:,:,:,krow), &
         xtte(:,:,:,krow))   
  ENDIF
  !
  ! cosp simulator diagnostics
  !
  IF (locosp ) THEN
    CALL call_cospsimulator  &
         (nproma,                   nlev,                           krow,                     &
         app1(1:nproma,:),         aphp1(1:nproma,:),              geom1(1:nproma,:),        &
         geospm(1:nproma,krow),    slf(1:nproma,krow),             tm1(1:nproma,:,krow),     &
         qm1(1:nproma,:,krow),                                     xlm1(1:nproma,:,krow),    &
         xim1(1:nproma,:,krow),    zfrl(1:nproma),                 zfrw(1:nproma),           &
         zfri(1:nproma),           tslm1(1:nproma,krow),           tsw(1:nproma,krow),       &
         tsi(1:nproma,krow)                                                                  )
  END IF
  !
  ! collect data for high frequency output at individual locations (cfSites)
  !
  IF (  lostation ) THEN
    CALL  collect_station_diag(nbdim,nlev,krow,geom1 )
    !!, vervel(:,:,krow) ) !!, geom1    )  

  END IF
  IF ( locospoffl ) THEN
    IF ( .NOT. lanysubmodel) THEN
      alpha0 = 0._dp
      CALL geopot(geohm1(:,2:nlevp1),ztvm1,alnpr(:,:,krow),alpha0,zgeo,nbdim,nproma) ! half level geopotential
      geohm1(:,1) = 1.e30_dp                              ! set upper boundary to a huge value
    END IF
    cospoffl_geom1(1:nproma,:,krow) = geom1(1:nproma,:)
    cospoffl_geohm1(1:nproma,:,krow) = geohm1(1:nproma,:)
    cospoffl_p(1:nproma,:,krow) = app1(1:nproma,:)
    cospoffl_ph(1:nproma,:,krow) = aphp1(1:nproma,:)
  END IF

  !
  ! daily block statistics - AMIP2 global diagnostics
  !  
  IF(ldiagamip) THEN
    CALL collect_amip2_diag(nproma,nbdim,nlev,nlevp1,krow                      &
         ,tm1(:,:,krow),um1(:,:,krow),vm1(:,:,krow),aphm1(:,:),geospm(:,krow)    &
         ,ustr(:,krow),ustrgw(:,krow),ustrm(:,krow),ustrgwm(:,krow)              &
         ,tslm1(:,krow),seaicem(:,krow),loland)
  END IF

  !
  ! Save large-scale/convective rain/snow rate for next time step
  !
  jrsfl(1:nproma,krow) = rsfl(1:nproma)
  jrsfc(1:nproma,krow) = rsfc(1:nproma)
  jssfl(1:nproma,krow) = ssfl(1:nproma)
  jssfc(1:nproma,krow) = ssfc(1:nproma)
  !
  !      
  !*        9.    RESTORE WINDS AND WIND TENDENCIES.
  !
  !*        9.2   RESTORE WINDS AND WIND TENDENCIES.
  !
  DO jlev = 1, nlev
!DIR$ CONCURRENT
    DO jl = 1, nproma
      zcst = sqcst_2d(jl,krow)
      um1(jl,jlev,krow) = zcst*um1(jl,jlev,krow)
      vm1(jl,jlev,krow) = zcst*vm1(jl,jlev,krow)
      vol(jl,jlev,krow) = zcst*vol(jl,jlev,krow)
      vom(jl,jlev,krow) = zcst*vom(jl,jlev,krow)
    END DO
  END DO
  !     ------------------------------------------------------------
  !
  !*       10.    RELEASE SPACE.
  !               ------- ------
  !
  !     ------------------------------------------------------------
  !
#ifdef _EUROHACK_2015
  ENDIF
#endif
END SUBROUTINE physc
