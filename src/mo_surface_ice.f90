#ifdef __xlC__
@PROCESS STRICT
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_surface_ice

  USE mo_echam_convect_tables, ONLY: lookup_ua_list_spline
  USE mo_exception,      ONLY: finish  
  USE mo_kind,           ONLY: wp

  IMPLICIT NONE

  REAL(wp), SAVE :: calbmni     ! minimum (bare sea ice)
  REAL(wp), SAVE :: calbmxi     ! maximum (bare sea ice)
  REAL(wp), SAVE :: calbmns     ! minimum (snow on ice)
  REAL(wp), SAVE :: calbmxs     ! maximum (snow on ice)
  REAL(wp), PARAMETER :: cb      = 2.0_wp  ! constant for zenith angle dependence (snow on ice)
  REAL(wp), PARAMETER :: cfresh  = 1.0_wp  ! mm snow water equivalent (constant for snow aging)
  REAL(wp), PARAMETER :: cs1alb  = 0.95_wp ! maximum albedo of fresh snow on ice (vis)
  REAL(wp), PARAMETER :: cs2alb  = 0.65_wp ! maximum albedo of fresh snow on ice (nir)
  REAL(wp), PARAMETER :: cs1age  = 0.10_wp ! constant for impact of snow aging on albedo  (vis)
  REAL(wp), PARAMETER :: cs2age  = 0.20_wp ! constant for impact of snow aging on albedo  (nir)
  REAL(wp), PARAMETER :: csola   = 0.13_wp ! constant for bare-ice albedo  (vis)
  REAL(wp), PARAMETER :: cnira   = 0.05_wp ! constant for bare-ice albedo  (nir)
  REAL(wp), PARAMETER :: csolb   = 0.20_wp ! constant for bare-ice albedo  (vis)
  REAL(wp), PARAMETER :: cnirb   = 0.20_wp ! constant for bare-ice albedo  (nir)

  REAL(wp) :: cseep         ! seepage rate in m/s
  REAL(wp) :: cagetime      ! inverse time scale of snow aging (s**-1)
  REAL(wp) :: cmpfac        ! tuning parameter for meltpond fraction
  
CONTAINS
  !
  !-------------------------------------------------------------------------------------------------

  SUBROUTINE init_surface_ice(nn)

    INTEGER, INTENT(in) :: nn      ! spectral truncation

    SELECT CASE (nn)

    CASE (31)
      
      cseep    = 0.016_wp/86400._wp
      cagetime = 4.E-7_wp
      cmpfac   = 0.25_wp
      
    CASE (63)
      
      cseep    = 0.016_wp/86400._wp
      cagetime = 4.E-7_wp
      cmpfac   = 0.50_wp
      
    CASE (127)
      
      cseep    = 0.008_wp/86400._wp
      cagetime = 1.0E-6_wp
      cmpfac   = 1.00_wp
      
    CASE (255)
      
      cseep    = 0.008_wp/86400._wp
      cagetime = 1.0E-6_wp
      cmpfac   = 1.00_wp

    CASE DEFAULT
      
      CALL finish ('mo_surface_ice::init_surface_ice', 'Truncation not supported.')

    END SELECT

  END SUBROUTINE init_surface_ice
    
  SUBROUTINE init_albedo_ice

    USE mo_control,   ONLY: nn
    USE mo_exception, ONLY: finish

    ! Define albedo parameters depending on resolution
    calbmxi = 0.75_wp
    IF (nn == 31) THEN
       calbmns = 0.65_wp
       calbmxs = 0.8_wp
       calbmni = 0.55_wp
    ELSE IF (nn == 63 .OR. nn == 127 .OR. nn==255) THEN
       calbmns = 0.70_wp
       calbmxs = 0.85_wp
       calbmni = 0.60_wp
    ELSE
       CALL finish ('init_albedo_ice', 'Truncation not supported.')
    ENDIF

  END SUBROUTINE init_albedo_ice
  !
  !---------------------------------------------------------------------------------------------------------
  SUBROUTINE update_albedo_ice(mask, surf_temp, snow, p_albedo)
    
    USE mo_physical_constants, ONLY:  &
         tmelt                 ! Melting temperature of ice/snow
    
    LOGICAL, INTENT(in)  :: mask(:)
    REAL(wp),    INTENT(in)  :: surf_temp(:), snow(:)
    REAL(wp),    INTENT(out) :: p_albedo(:)
    
    ! Local variables
    REAL(wp), DIMENSION(SIZE(mask)) :: alb_min, alb_max              & ! Minimum and maximum snow albedo
         , temp_upper_limit                ! Upper temperature limit for cold snow albedo

    temp_upper_limit = tmelt-1.0_wp

    p_albedo = MERGE(0.55_wp,0.0_wp,mask)

    WHERE (mask)
       
       ! Minimum and maximum albedo
       WHERE (snow > 0.01_wp)               ! sea ice covered with snow
          alb_min = calbmns
          alb_max = calbmxs
       ELSEWHERE                         ! bare sea ice
          alb_min = calbmni
          alb_max = calbmxi
       END WHERE
       
       ! Temperature-dependent snow albedo
       WHERE (surf_temp >= tmelt)
          p_albedo = alb_min
       ELSEWHERE (surf_temp < temp_upper_limit)
          p_albedo = alb_max
       ELSEWHERE
          p_albedo = alb_min + &
               (((alb_max-alb_min) / (tmelt-temp_upper_limit)) *    & ! Change of snow albedo per deg C
               (tmelt-surf_temp))
       END WHERE
       
    END WHERE
    
  END SUBROUTINE update_albedo_ice
  !
  !---------------------------------------------------------------------------------------------------------
  SUBROUTINE update_albedo_ice_meltpond(kdim, mask                            &
                             , ptsi, psnow                                    &
                             , pcvsi, psiced, pmu0                            &
                             , pameltdepth, pameltfrac                        &
                             , ptaus, pfage, psnifrac, pbarefrac              &
                             , pswvis, pswdifvis, pswnir, pswdifnir           &
                             , palbedo_vis_dir, palbedo_vis_dif               &
                             , palbedo_nir_dir, palbedo_nir_dif               &
                             , pali1dir, pali2dir, pali1dif, pali2dif         &
                             , pali1, pali2, palsom, palsobs, palsoi)
!
    USE mo_time_control,   ONLY: delta_time
    USE mo_physical_constants,      ONLY: tmelt, dice
!        
    INTEGER,  INTENT(in) :: kdim              ! number of longitudes
    LOGICAL,  INTENT(in) :: mask(kdim)        ! ice mask
    REAL(wp), INTENT(in) :: ptsi(kdim)        ! ice temperature
    REAL(wp), INTENT(in) :: psnow(kdim)       ! total snowfall
    REAL(wp), INTENT(in) :: psiced(kdim)      ! sea ice thickness
    REAL(wp), INTENT(in) :: pcvsi(kdim)       ! snow cover (ice)
    REAL(wp), INTENT(in) :: pameltdepth(kdim) ! melt pond depth
    REAL(wp), INTENT(in) :: pameltfrac(kdim)  ! melt pond fraction
    REAL(wp), INTENT(in) :: pswvis(kdim)      ! net surface visible
    REAL(wp), INTENT(in) :: pswdifvis(kdim)   ! fraction of diffuse visible
    REAL(wp), INTENT(in) :: pswnir(kdim)      ! net surface near infrared
    REAL(wp), INTENT(in) :: pswdifnir(kdim)   ! fraction of diffuse near infrared
    REAL(wp), INTENT(in) :: palbedo_vis_dir(kdim)   ! grid-mean visible albedo (dir)
    REAL(wp), INTENT(in) :: palbedo_vis_dif(kdim)   ! grid-mean visible albedo (dif)
    REAL(wp), INTENT(in) :: palbedo_nir_dir(kdim)   ! grid-mean near infrared ! albedo (dir)
    REAL(wp), INTENT(in) :: palbedo_nir_dif(kdim)   ! grid-mean near infrared albedo (dif)
    REAL(wp), INTENT(in) :: pmu0(kdim)        ! cos of solar zenith angle
    REAL(wp), INTENT(inout) :: ptaus(kdim)    !  non-dimensional age of snow on ice

    ! OUTPUT
    ! ------

    REAL(wp), INTENT(inout) :: pfage(kdim)      ! aging factor of snow on ice
    REAL(wp), INTENT(inout) :: psnifrac(kdim)   ! fraction of ice covered with snow
    REAL(wp), INTENT(inout) :: pbarefrac(kdim)  ! bare ice fraction
    REAL(wp), INTENT(inout) :: pali1dir(kdim)   ! ice albedo (vis, dir)
    REAL(wp), INTENT(inout) :: pali2dir(kdim)   ! ice albedo (nir, dir)
    REAL(wp), INTENT(inout) :: pali1dif(kdim)   ! ice albedo (vis, dif)
    REAL(wp), INTENT(inout) :: pali2dif(kdim)   ! ice albedo (nir, dif)
    REAL(wp), INTENT(inout) :: palsom(kdim)     ! albedo of melt ponds
    REAL(wp), INTENT(inout) :: palsoi(kdim)     ! ice albedo including ponds
    REAL(wp), INTENT(inout) :: palsobs(kdim)    ! albedo of bare ice and snow without ponds
    REAL(wp), INTENT(inout) :: pali1(kdim)      ! albedo of ice, mean of 3 surfaces (vis, dir+dif)
    REAL(wp), INTENT(inout) :: pali2(kdim)      ! albedo of ice, mean of 3 surfaces (nir, dir+dif)

    ! LOCAL
    ! -----

    REAL(wp) :: zalsom(kdim)     ! albedo of melt ponds
    REAL(wp) :: zalsoi(kdim)     ! ice albedo including ponds
    REAL(wp) :: zalsobs(kdim)    ! albedo of bare ice and snow without ponds

    REAL(wp) :: zeps        ! safety
    REAL(wp) :: zsnifrac    ! fraction of ice covered with snow
    REAL(wp) :: zbarefrac   ! fraction of bare ice
    REAL(wp) :: zalmp1dir   ! albedo of melt ponds (vis, dir)
    REAL(wp) :: zalmp2dir   ! albedo of melt ponds (nir, dir)
    REAL(wp) :: zalmp1dif   ! albedo of melt ponds (vis, dif)
    REAL(wp) :: zalmp2dif   ! albedo of melt ponds (nir, dif)
    REAL(wp) :: zalmp1      ! albedo of melt ponds (vis)
    REAL(wp) :: zalmp2      ! albedo of melt ponds (nir)
    REAL(wp) :: zalbare1dir ! albedo of bare ice (vis, dir)
    REAL(wp) :: zalbare2dir ! albedo of bare ice (nir, dir)
    REAL(wp) :: zalbare1dif ! albedo of bare ice (vis, dif)
    REAL(wp) :: zalbare2dif ! albedo of bare ice (nir, dif)
    REAL(wp) :: zd          ! area-mean pond depth / pond fraction
    REAL(wp) :: ztheta      ! function of solar zenith angle (for snow on ice)
    REAL(wp) :: zalsni1     ! albedo of snow on ice (vis, dir+dif)
    REAL(wp) :: zalsni2     ! albedo of snow on ice (nir, dir+dif)
    REAL(wp) :: zalsob      ! albedo of bare ice (vis+nir, dir+dif)
    REAL(wp) :: zalsos      ! albedo of snow on ice (vis+nir, dir+dif)
    REAL(wp) :: zvisp_down  ! downward visible direct (parallel) radiation at surface
    REAL(wp) :: zvisd_down  ! downward visible diffuse radiation at surface
    REAL(wp) :: zvis_down   ! downward visible radiation at surface
    REAL(wp) :: znirp_down  ! downward near-infrared direct (parallel) radiation at surfac
    REAL(wp) :: znird_down  ! downward near-infrared diffuse radiation at surface
    REAL(wp) :: znir_down   ! downward near-infrared radiation at surface
    REAL(wp) :: zrad_down   ! downward visible + near-infrared radiation at surface
    REAL(wp) :: zsiced, zalbare1, zalbare2, zr1, zr2, zr3, zdtaus, zsnow, zdtime

    INTEGER :: jl    ! loop index

    zeps        = EPSILON(1.0_wp)
    zdtime      = delta_time
    pfage       = 0.0_wp
    psnifrac    = (1.0_wp-pameltfrac)*pcvsi
    pbarefrac   = 1.0_wp-psnifrac-pameltfrac
    zalsoi      = palsoi
    zalsobs     = palsobs
    zalsom      = palsom

    DO jl = 1, kdim

!    IF (psiced(jl).GE.dice) THEN
    IF(mask(jl)) THEN

!   Visible and near infrared downwelling solar radiation at the surface
!   (separately for parallel and diffuse radiation)

    zvisp_down    = pswvis(jl)/(1.0_wp-palbedo_vis_dir(jl)) * (1.0_wp-pswdifvis(jl))
    zvisd_down    = pswvis(jl)/(1.0_wp-palbedo_vis_dif(jl)) * pswdifvis(jl)
    znirp_down    = pswnir(jl)/(1.0_wp-palbedo_nir_dir(jl)) * (1.0_wp-pswdifnir(jl))
    znird_down    = pswnir(jl)/(1.0_wp-palbedo_nir_dif(jl)) * pswdifnir(jl)
    zvis_down     = zvisp_down + zvisd_down
    znir_down     = znirp_down + znird_down 
    zrad_down     = zvis_down  + znir_down          
    
    ! Area fractions of bare ice and snow covered ice
    ! -----------------------------------------------
     
    zbarefrac     = pbarefrac(jl)
    zsnifrac      = psnifrac(jl) 
    
    !  Snow on ice
    !  -----------
    
    IF (zsnifrac > 0.0_wp) THEN
       
       !   Snow aging
       !   ----------
       
       zr1        = MIN(1.0_wp,EXP(5000.0_wp*(1.0_wp/tmelt-1.0_wp/ptsi(jl))))
       zr2        = zr1**10.0_wp
       zr3        = 0.01_wp
       zdtaus     = (zr1+zr2+zr3)*zdtime*cagetime
       zsnow      = psnow(jl)*zdtime/cfresh
       ptaus(jl)  = MAX(0.0_wp,(ptaus(jl)+zdtaus)*(1.0_wp-zsnow))
       pfage(jl)  = ptaus(jl)/(1.0_wp+ptaus(jl))
       IF (pmu0(jl) < 0.5_wp) THEN
          ztheta  = ((1.0_wp+cb)/(1.0_wp+2.0_wp*cb*pmu0(jl))-1.0_wp)/cb
       ELSE
          ztheta  = 0.0_wp
       END IF
       zalsni1    = cs1alb*(1.0_wp-cs1age*pfage(jl))
       zalsni2    = cs2alb*(1.0_wp-cs2age*pfage(jl))
       zalsni1    = zalsni1+0.4*ztheta*(1.0_wp-zalsni1)
       zalsni2    = zalsni2+0.4*ztheta*(1.0_wp-zalsni2)
    ELSE
       ptaus(jl)  = 0.0_wp
       pfage(jl)  = 0.0_wp
       zalsni1    = 0.0_wp
       zalsni2    = 0.0_wp
    END IF  !  snow on ice
    
    IF (zrad_down .GT. 0.0_wp) THEN
       
       !   Albedo of snow on ice (1/vis + 2/nir)
       !   -------------------------------------
       
       zalsos       = (zalsni1*zvis_down+zalsni2*znir_down)/(zrad_down+zeps)
       
       !   bare ice
       !   --------
       
       IF (zbarefrac > 0.0_wp) THEN
          zsiced         = MAX(dice,psiced(jl))*100.0_wp
          IF (zsiced < 100.0_wp) THEN
             zalbare1dir = csola*LOG(zsiced)+csolb
             zalbare2dir = cnira*LOG(zsiced)+cnirb
          ELSE
             zalbare1dir = 0.80_wp
             zalbare2dir = 0.43_wp
          END IF
          zalbare1dif    = zalbare1dir
          zalbare2dif    = zalbare2dir
!!$         IF (zsiced < 100.0_wp) THEN
!!$            zalbare2dir = 0.047_wp*LOG(zsiced)+0.074_wp
!!$            zalbare2dif = 0.049_wp*LOG(zsiced)+0.085_wp
!!$         ELSE
!!$            zalbare2dir = 0.29_wp
!!$            zalbare2dif = 0.31_wp
!!$          END IF
          
          ! Albedo of bare ice
          ! ------------------
          !   dir/p + dif/d

          zalbare1     = (zalbare1dir*zvisp_down+zalbare1dif*zvisd_down)/(zvis_down+zeps)
          zalbare2     = (zalbare2dir*znirp_down+zalbare2dif*znird_down)/(znir_down+zeps)

          !   1/vis + 2/nir
          !   -------------

          zalsob       = (zalbare1*zvis_down+zalbare2*znir_down)/(zrad_down+zeps)
       ELSE
          zalbare1dir  = 0.0_wp
          zalbare2dir  = 0.0_wp
          zalbare1dif  = 0.0_wp
          zalbare2dif  = 0.0_wp
          zalsob       = 0.0_wp
       END IF     ! bare ice
       
       ! Albedo of bare ice + snow covered ice
       ! -------------------------------------
              
       palsobs(jl)     = (zbarefrac*zalsob+zsnifrac*zalsos)/(zbarefrac+zsnifrac+zeps)
       
       !   Albedo of melt ponds
       !   --------------------
       
       IF (pameltfrac(jl) > 0.0_wp) THEN
          zd           = MIN(psiced(jl), pameltdepth(jl)/pameltfrac(jl))
          zalmp1dir    = 0.336_wp+EXP( -9.457_wp*zd-1.061_wp)
          zalmp2dir    = 0.017_wp+EXP(-18.904_wp*zd-0.909_wp)
          zalmp1dif    = 0.413_wp+EXP(-24.014_wp*zd-1.086_wp)
          zalmp2dif    = 0.061_wp+EXP(-17.449_wp*zd-1.075_wp)
          
          ! dir/p + dif/d
          ! -------------

          zalmp1       = (zalmp1dir*zvisp_down+zalmp1dif*zvisd_down)/(zvis_down+zeps)
          zalmp2       = (zalmp2dir*znirp_down+zalmp2dif*znird_down)/(znir_down+zeps) 
 
          ! 1/vis + 2/nir
          ! -------------

          palsom(jl)   = (zalmp1*zvis_down+zalmp2*znir_down)/(zrad_down+zeps)
       ELSE
          zalmp1dir    = 0.0_wp
          zalmp2dir    = 0.0_wp
          zalmp1dif    = 0.0_wp
          zalmp2dif    = 0.0_wp
          palsom(jl)   = 0.0_wp
       END IF   ! melt ponds
       
       ! Albedo of bare ice + snow covered ice + melt ponds
       ! --------------------------------------------------
       
       pali1dir(jl) = zbarefrac*zalbare1dir+zsnifrac*zalsni1+pameltfrac(jl)*zalmp1dir
       pali2dir(jl) = zbarefrac*zalbare2dir+zsnifrac*zalsni2+pameltfrac(jl)*zalmp2dir
       pali1dif(jl) = zbarefrac*zalbare1dif+zsnifrac*zalsni1+pameltfrac(jl)*zalmp1dif
       pali2dif(jl) = zbarefrac*zalbare2dif+zsnifrac*zalsni2+pameltfrac(jl)*zalmp2dif

       ! dir/p + dif/d
       ! -------------

       pali1(jl)    = (pali1dir(jl)*zvisp_down+pali1dif(jl)*zvisd_down)/(zvis_down+zeps)
       pali2(jl)    = (pali2dir(jl)*znirp_down+pali2dif(jl)*znird_down)/(znir_down+zeps)
       
       !  Mean ice albedo (all surfaces)
       !  ------------------------------
       
       palsoi(jl)   = zbarefrac*zalsob+zsnifrac*zalsos+pameltfrac(jl)*palsom(jl)

    ELSE  ! night: albedos from last timestep
     
       pali1dir(jl) = zalsoi(jl)
       pali2dir(jl) = zalsoi(jl)
       pali1dif(jl) = zalsoi(jl)
       pali2dif(jl) = zalsoi(jl)
       pali1(jl)    = zalsoi(jl)
       pali2(jl)    = zalsoi(jl)
       palsoi(jl)   = zalsoi(jl)
       palsobs(jl)  = zalsobs(jl)
       palsom(jl)   = zalsom(jl)

    END IF  ! zrad_down > 0
   END IF   ! mask
  END DO
    
  END SUBROUTINE update_albedo_ice_meltpond
  !
  !---------------------------------------------------------------------------------------------------------
  SUBROUTINE update_z0_ice(kdim, mask, paz0i)

    USE mo_physc2,         ONLY: cz0ice
    
    INTEGER, INTENT(in)     :: kdim
    LOGICAL, INTENT(in)     :: mask(kdim)
    REAL(wp), INTENT(out)       :: paz0i(kdim)

    paz0i = cz0ice

    WHERE(mask)
       paz0i(:) = cz0ice
    END WHERE

  END SUBROUTINE UPDATE_Z0_ICE
  !---------------------------------------------------------------------------------------------------------
  !
  SUBROUTINE update_stress_ice(     &
       kdim          , mask         &
     , zcfmi       , zudif        &
     , zvdif       , pustri       &
     , pvstri                     &
     )

    USE mo_time_control, ONLY: time_step_len
    USE mo_physical_constants,    ONLY: grav
    
    INTEGER, INTENT(in)    :: kdim
    LOGICAL, INTENT(in)    :: mask(kdim)
    REAL(wp), INTENT(in)       :: zcfmi(kdim)
    REAL(wp), INTENT(in)       :: zudif(kdim), zvdif(kdim)
    REAL(wp), INTENT(out)      :: pustri(kdim), pvstri(kdim)

    REAL(wp) :: ztmst, zcons15

    pustri = 0._wp
    pvstri = 0._wp

    ztmst   = time_step_len
    zcons15 = 1._wp / (grav * ztmst)

    WHERE(mask)
       pustri(:) = zcons15 * zcfmi(:) * zudif(:)
       pvstri(:) = zcons15 * zcfmi(:) * zvdif(:)
    END WHERE

  END SUBROUTINE update_stress_ice

  !---------------------------------------------------------------------------------------------------------
  !
  SUBROUTINE precalc_ice( &
       kdim , mask &
       , ptsi, paphm1 &
       , pqm1, zx &
       , ptm1, zqss &
       , zteta1, ztvir1 &
       , zfaxe, paclc &
       , zlteta1, pgeom1 &
       , pum1, pvm1 &
       , pocu, pocv &
       , paz0i, zghabl &
       , zqsi, zcpti &
       , zrii, zcfmi &
       , zchi, zcfhi &
       , zbni, zbmi &
       , zbhi, zustari &
       , ztkevi &
!---------------------------------------------------------------------------
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
       , zcfnci, ztvi, zcdni   & 
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
!---------------------------------------------------------------------------
       , zusti &
       )

    USE mo_physical_constants, ONLY: grav, rd, vtmpc1, cpd, vtmpc2
    USE mo_physc2,       ONLY: ckap, cc, cb, cvdifts
    USE mo_time_control, ONLY: time_step_len

    INTEGER, INTENT(in)    :: kdim
    LOGICAL, INTENT(in)    :: mask(kdim)
    REAL(wp),    INTENT(IN)    :: ptsi(kdim), paphm1(kdim), pqm1(kdim), zx(kdim)
    REAL(wp),    INTENT(in)    :: ptm1(kdim), zqss(kdim), zteta1(kdim)
    REAL(wp),    INTENT(IN)    :: ztvir1(kdim), zfaxe(kdim), paclc(kdim)
    REAL(wp),    INTENT(in)    :: zlteta1(kdim), pum1(kdim), pvm1(kdim), paz0i(kdim)
    REAL(wp),    INTENT(in)    :: pocu(kdim), pocv(kdim)
    REAL(wp),    INTENT(IN)    :: pgeom1(kdim), zghabl(kdim)
    
    ! ACHTUNG klevp1:  paphm1
    ! ACHTUNG klev:    pqm1, ptm1, zqss, zx, ztvir1, zteta1,zface, paclc, zlteta1, pgeom1, 
    REAL(wp),    INTENT(out)   :: zqsi(kdim), zcpti(kdim), zrii(kdim)
    REAL(wp),    INTENT(out)   :: zcfmi(kdim), zchi(kdim), zcfhi(kdim)
    REAL(wp),    INTENT(out)   :: zbni(kdim), zbmi(kdim), zbhi(kdim)
    REAL(wp),    INTENT(out)   :: zustari(kdim), ztkevi(kdim), zusti(kdim)
    
!---------------------------------------------------------------------------
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
    REAL(wp),    INTENT(out)   :: zcfnci(kdim) 
    REAL(wp),    INTENT(out)   :: zcdni(kdim) 
    REAL(wp),    INTENT(out)   :: ztvi(kdim) 
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
!---------------------------------------------------------------------------

    !  SAVE VARIABLE FOR DIAGNOSE AFTER VERTICAL COLUMN IST CALCULATED WITHIN VDIFF
    
    
    ! Local Variables
    INTEGER  :: loidx(kdim), is, nl, jl
    
    REAL(wp)     :: zes(kdim), zqmitte(kdim), zqtmit(kdim), ztmit(kdim), zqsmit(kdim), zvirmitte(kdim)
    REAL(wp)     :: ztemitte(kdim), zfux(kdim), zfox(kdim), zmult1(kdim), zqddif(kdim), zbuoy(kdim)
    REAL(wp)     :: zmult2(kdim), zmult3(kdim), zmult4(kdim), zmult5(kdim), zdus1(kdim), zdus2(kdim)
    REAL(wp)     :: zteldif(kdim), zalo(kdim), zcons(kdim), zdthv(kdim) 
    REAL(wp)     :: zucfi(kdim), zscfi(kdim)
    REAL(wp)     :: zcdn2m(kdim), zcfm2m(kdim)
    REAL(wp)     :: ztesi(kdim), zdu2oc(kdim), zcdnr(kdim), zconvs(kdim), zwsti(kdim)
    REAL(wp)     :: zmonob(kdim), zstabf(kdim)
    REAL(wp)     :: ua(kdim)
  
    REAL(wp)     :: zkappa, zkap, zrvrd, zcons9, zcons8, zcons11, zcons12, zepsec, zcons17, zepdu2
    REAL(wp)     :: zrdrv, zepz0o, zepsr, zcons6, zustf
    REAL(wp)     :: zsmn, zshn, zm1, zm2, zm4, zh1, zh2, zwstf
    
    ! Constants
    zkappa          = rd / cpd
    zkap            = ckap
    zrvrd           = vtmpc1 + 1._wp
    zepdu2          = 1.0_wp
    zcons8          = 2._wp * cb
    zcons9          = 3._wp * cb
    zepsec          = 1.e-2_wp
    zrdrv           = 1._wp / zrvrd
    zepz0o          = 2._wp
    zepsr           = 1.e-10_wp
    
    zcons11         = 3._wp * cb * cc
    zcons12         = cvdifts * time_step_len * grav / rd
    zcons17         = 1._wp / zkap**2
    zcons6          = 1._wp / 3._wp
    zh1             = 2.22_wp
    zh2             = 0.22_wp
    zm1             = 1.24_wp
    zm2             = 2.37_wp
    zm4             = 3.69_wp
    zshn            = zh1 * zh2 * SQRT(2._wp)
    zsmn            = zshn * zm1 * zm2 / zm4
    zustf           = 1._wp / zsmn**2
    zwstf           = 0.2_wp

    zqsi = 0._wp
    zcpti = 0._wp
    zrii = 0._wp
    zcfmi = 0._wp
    zchi = 0._wp
    zcfhi = 0._wp
    zbni = 0._wp
    zbmi = 0._wp
    zbhi = 0._wp
    zustari = 0._wp
    ztkevi = 0._wp
    zusti = 0._wp
    loidx = 0

!grazia  #78
    ztvi(:)   = 0._wp
    zcdni(:)  = 0._wp
    zcfnci(:) = 0._wp
!grazia #78

    !------------------------------------------------------
    !      2.2   surface humidity and virtual temperature
    is = 0
    DO jl=1,kdim
      IF(mask(jl)) THEN
        is = is + 1
        loidx(is) = jl
      ENDIF
    ENDDO
   
      CALL lookup_ua_list_spline('precalc_ice',kdim,is,loidx(1), ptsi(1), ua(1))

    DO nl=1,is
      jl = loidx(nl)
      zes(jl)       = ua(nl) / paphm1(jl)
      zqsi(jl)      = zes(jl) / (1._wp- vtmpc1 * zes(jl))
      zcpti(jl)     = ptsi(jl) * cpd * (1._wp+ vtmpc2 * zqsi(jl))
      ztesi(jl)     = ptsi(jl) * (1.e5_wp / paphm1(jl))**zkappa
      ztvi(jl)      = ztesi(jl) * (1._wp + vtmpc1 * zqsi(jl))
    ENDDO
    
    !     ------------------------------------------------------------------
    !        3.     COMPUTATION OF THE EXCHANGE COEFFICIENTS.
    !        3.1       COMPUTATION OF BASIC QUANTITIES: WIND SHEAR,
    !                  RICHARDSON NUMBER,SQUARED MIXING LENGTHS, UNSTABLE
    !                  AND STABLE CASE COMMON FACTORS AND NEUTRAL CASE
    !                  COMMON PART OF THE DRAG COEFFICIENTS.
    !
    WHERE(mask)
       zdu2oc(:)       = MAX(zepdu2 , pum1(:)**2 + pvm1(:)**2)
    !
    ! correction for water and ice points
    !
       zdu2oc(:)    = MAX(zepdu2,(pum1(:)-pocu(:))**2                  &
                                +(pvm1(:)-pocv(:))**2)
       zcons(:)      = zcons12 * paphm1(:) / (ptm1(:) * (1._wp + vtmpc1 * pqm1(:) - zx(:)))       
       zqmitte(:)    = (pqm1(:) + zqsi(:)) / 2._wp
       zqtmit(:)     = zx(:) * 0.5_wp + zqmitte(:)
       ztmit(:)      = (ptm1(:) + ptsi(:)) / 2._wp
       zqsmit(:)     = (zqss(:) + zqsi(:)) / 2._wp
       ztemitte(:)   = (zteta1(:) + ztesi(:)) / 2._wp
       zvirmitte(:)  = (ztvir1(:) + ztvi(:)) / 2._wp
       zfux(:)       = zfaxe(:) / (cpd * ztmit(:))
       zfox(:)       = zfaxe(:) / (rd * ztmit(:))
       zmult1(:)     = 1._wp + vtmpc1 * zqtmit(:)
       zmult2(:)     = zfux(:) * zmult1(:) - zrvrd
       zmult3(:)     = zrdrv * zfox(:) * zqsmit(:) / (1._wp + zrdrv * zfox(:) * zfux(:) * zqsmit(:))
       zmult5(:)     = zmult1(:) - zmult2(:) * zmult3(:)
       zmult4(:)     = zfux(:) * zmult5(:) - 1._wp
       zdus1(:)      = paclc(:) * zmult5(:) + (1._wp - paclc(:)) * zmult1(:)
       zdus2(:)      = paclc(:) * zmult4(:) + (1._wp - paclc(:)) * vtmpc1
       zteldif(:)    = zlteta1(:) - ztesi(:)
       zqddif(:)     = pqm1(:) + zx(:) - zqsi(:)
       zbuoy(:)      = zdus1(:) * zteldif(:) + zdus2(:) * ztemitte(:) * zqddif(:)
       zrii(:)       = pgeom1(:) * zbuoy(:) / (zvirmitte(:) * zdu2oc(:))
    END WHERE
    WHERE(mask)
       zalo(:)       = LOG(1._wp+ pgeom1(:) / (grav * paz0i(:)))
    END WHERE
    WHERE(mask)
       zcdni(:)      = (ckap / zalo(:))**2
    END WHERE
    WHERE(mask)
       zucfi(:)      = 1._wp / (1._wp + zcons11 * zcdni(:) * SQRT(ABS(zrii(:)) &
            * (1._wp + pgeom1(:) / (grav * paz0i(:)))))
    END WHERE
    WHERE(mask)
       zscfi(:)      = SQRT(1._wp+ ABS(zrii(:)))
       zcfnci(:)     = zcons(:) * SQRT(zdu2oc(:)) * zcdni(:)
    END WHERE
    WHERE(mask)
       zdthv(:)      = MAX(0._wp,(ztvi(:) - ztvir1(:)))
       zwsti(:)      = zdthv(:) * SQRT(zdu2oc(:)) / zvirmitte(:)
    ENDWHERE
    !----------------------------------------------------------------------------------------------
    !     3.2  DIMENSIONLESS HEAT TRANSFER COEFFICIENTS MULTIPLIED
    !          BY PRESSURE THICKNESSES FOR MOMENTUM AND HEAT EXCHANGE
    !
    WHERE(mask) 
       WHERE(zrii(:).GT.0._wp)
          zcfmi(:)  = zcfnci(:) / (1._wp + zcons8 * zrii(:) / zscfi(:))
          zcfhi(:)  = zcfnci(:) / (1._wp + zcons8 * zrii(:) * zscfi(:))
          zchi(:)   = zcfhi(:) / zcfnci(:) * zcdni(:)
       ELSEWHERE
          zcfmi(:)  = zcfnci(:) * (1._wp - zcons8 * zrii(:) * zucfi(:))
          zcfhi(:)  = zcfnci(:) * (1._wp - zcons9 * zrii(:) * zucfi(:))
          zchi(:)   = zcfhi(:) / zcfnci(:) * zcdni(:)
       ENDWHERE
    ENDWHERE
    !-----------------------------------------------------------------------------------------------
    !     interpolation functions for diagnostics
    !
    WHERE(mask)
       zbni(:)      = zkap / SQRT(zcdni(:))
       zbmi(:)      = MAX(zepsec, SQRT(zcfmi(:) * zcdni(:) * zcons17 / zcfnci(:)))
       zbhi(:)      = MAX(zepsec, zchi(:) / zbmi(:) * zcons17)
       zbmi(:)      = 1._wp / zbmi(:)
       zbhi(:)      = 1._wp / zbhi(:)
    ENDWHERE
    
    !-----------------------------------------------------------------------------------------------
    !*       3.4       COMPUTATION OF THE PBL EXTENSION.
    !
    WHERE(mask)
       WHERE(paz0i(:).GT.zepz0o)
          zcdn2m(:)  = (zkap / LOG(1._wp+ pgeom1(:) / (grav * zepz0o)))**2
       ELSEWHERE
          zcdn2m(:)   = zcdni(:)
       ENDWHERE
       zcdnr(:)       = zcdn2m(:) / zcdni(:)
       WHERE( paz0i(:).GT.zepz0o.AND.zrii(:).LT.0._wp)
          zcfm2m(:)   = zcfnci(:) * zcdnr(:) * (1._wp - zcons8 *zrii(:) &
               / (1._wp + zcons11 * zcdn2m(:) * SQRT(ABS(zrii(:)) &
               * (1._wp + pgeom1(:) / (grav * zepz0o)))))
       ELSEWHERE
          zcfm2m(:)   = zcfmi(:) * zcdnr(:)
       ENDWHERE
       zusti(:)       = zcfm2m(:) * SQRT(zdu2oc(:))
       zustari(:)     = SQRT(zusti(:) * ptm1(:) &
            * (1._wp + vtmpc1 * pqm1(:) - zx(:)) &
            / (zcons12 * paphm1(:)))
    ENDWHERE
    
    !----------------------------------------------------------------------------------------------
    !      CONVECTIVE VELOCITY SCALE, MONIN-OBUKHOV LENGTH AND
    !      TKE BOUNDARY CONDITION (MAILHOT/BENOIT, 1982)
    WHERE(mask)
       WHERE(zwsti(:) .GT. zepsr)
          zconvs(:)    = (zwsti(:) * zchi(:) * zghabl(:))**zcons6
          zmonob(:)    = (zustari(:)**3) / (zkap * grav * zwsti(:) * zchi(:))
          zstabf(:)    = (pgeom1(:) / (grav * zmonob(:)))**(zcons6 * 2._wp)
          zstabf(:)    = MIN(zustf * 3._wp, zstabf(:))
       ELSEWHERE
          zconvs=0._wp
          zstabf=0._wp
       ENDWHERE
       ztkevi(:)       = (zustf + zstabf(:)) * (zustari(:)**2) + zwstf * (zconvs(:)**2)
       
    ENDWHERE
  END SUBROUTINE PRECALC_ICE
  !------------------------------------------------------------------------------------------------------
  !
  SUBROUTINE richtmyer_ice(kdim, mask, & 
       klev, klevp1, klevm1, & 
       paphm1, zcfh, &
       zebsh, zqdif, &
       ztdif, zcfhi, &
       zetni, zftni, &
       zeqni, zfqni &
       )

    USE mo_physc2,       ONLY: cvdifts

    INTEGER,   INTENT(IN)  :: kdim, klev, klevp1, klevm1
    LOGICAL,   INTENT(IN)  :: mask(kdim)
    REAL(wp),      INTENT(IN)  :: zebsh(kdim, klev), zcfhi(kdim)
    REAL(wp),      INTENT(IN)  :: paphm1(kdim,klevp1), zcfh(kdim,klev), ztdif(kdim,klev)
    REAL(wp),      INTENT(IN)  :: zqdif(kdim,klev)
    REAL(wp),      INTENT(OUT) :: zetni(kdim), zftni(kdim), zeqni(kdim), zfqni(kdim)
    
    REAL(wp)   :: zdisci(kdim), zdisqi(kdim)
    REAL(wp)   :: zqdp(kdim), zfac(kdim)
    REAL(wp)   :: ztpfac1
    
    !------------------------
    ! CONSTANTS

    ztpfac1   = cvdifts
    zetni = 0._wp
    zftni = 0._wp
    zeqni = 0._wp
    zfqni = 0._wp

    WHERE(mask)
       zqdp(:)       = 1._wp / (paphm1(:,klevp1) - paphm1(:,klev))
       zfac(:)       = zcfh(:,klevm1) * zqdp(:)
       !-------------------------------------------------------------------
       !*  CALCULATION OF THE EN AND FN COEFFICIENTS OF THE RICHTMYER-
       !*  MORTON-SCHEME CONCERNING THE EQUATION:
       !
       !*  XN = EN * XS + FN
       !
       !*  WITH XN = S_ATM  OR  XN = QATM : ATM. VALUE OF S OR Q
       !*  AND  XS = SSURF  OR  XS = QSAT : SURFACE VALUE OF S OR SAT. SPEC.
       !*                                   HUM. AT THE SURFACE
       !
       zdisci(:)     = 1._wp / (1._wp + zfac(:) * (1._wp - zebsh(:,klevm1)) + zcfhi(:) * zqdp(:))
       zdisqi(:)     = zdisci(:)
       zetni(:)      = zdisci(:) * zcfhi(:) * zqdp(:)
       zftni(:)      = zdisci(:) * (ztdif(:,klev) + zfac(:) * ztdif(:,klevm1)) * ztpfac1
       zeqni(:)      = zdisqi(:) * zcfhi(:) * zqdp(:)
       zfqni(:)      = zdisqi(:) * (zqdif(:,klev) + zfac(:) * zqdif(:,klevm1)) * ztpfac1
    ENDWHERE    
  END SUBROUTINE richtmyer_ice
  !-----------------------------------------------------------------------------------
  !  
  
  SUBROUTINE update_ice (kdim, mask                     &
       , zetni, zftni &
       , zeqni, zfqni                     &
       , zcpti, zqsi                                    &
       , ztklevi, zqklevi &
       )

    INTEGER, INTENT(IN) :: kdim
    LOGICAL, INTENT(IN) :: mask(kdim)
    REAL(wp), INTENT(IN)    :: zetni(kdim), zftni(kdim), zeqni(kdim), zfqni(kdim)
    REAL(wp), INTENT(IN)    :: zcpti(kdim), zqsi(kdim)
    REAL(wp), INTENT(OUT)   :: ztklevi(kdim), zqklevi(kdim)
    
    !
    !*  CALCULATION OF SKLEV AND QKLEV USING THE NEW SURFACE VALUES
    !*  ZSNEW AND ZQSNEW WHICH WERE CALCULATED IN SUBROUTINE SURFTEMP
    !
    ztklevi = 0._wp
    zqklevi = 0._wp

    WHERE(MASK)
       ztklevi(:)   = zetni(:) * zcpti(:) + zftni(:)
       zqklevi(:)   = zeqni(:) * zqsi(:) + zfqni(:)
    ENDWHERE
  END SUBROUTINE update_ice
  !--------------------------------------------------------------------------------------
  !
  SUBROUTINE postproc_ice(kdim, klevp1, mask     &
       , zcfhi, zqsi &
       , zqdif, pqm1 &
       , pgeom1, ztdif &
       , zcptgz, zcpti &
       , ptsi , zbni &
       , zbhi, zrii &
       , ptm1, papm1 &
       , zx, pum1 &
       , pvm1, pevapi &
       , pahfsi, pahfli &
       , zspeedi, zt2i &
       , paphm1, zbmi &
       , zdew2i, zu10i &
       , zv10i &
       )

    USE mo_physical_constants,    ONLY: grav, vtmpc1, als, tmelt, rd,                   &
                                        c3les, c4les, c3ies, c4ies, c2es, vtmpc2, cpd
    USE mo_physc2,       ONLY: cvdifts
    USE mo_time_control, ONLY: time_step_len

    INTEGER, INTENT(in)  :: kdim, klevp1
    LOGICAL, INTENT(in)  :: mask(kdim)
    REAL(wp), INTENT(in)     :: zcfhi(kdim), zqsi(kdim)
    REAL(wp), INTENT(in)     :: zqdif(kdim), pqm1(kdim), pgeom1(kdim) 
    REAL(wp), INTENT(in)     :: ztdif(kdim), zcptgz(kdim)
    REAL(wp), INTENT(in)     :: zcpti(kdim), ptsi(kdim), zbni(kdim), zbhi(kdim), zrii(kdim)
    REAL(wp), INTENT(in)     :: ptm1(kdim), papm1(kdim), zx(kdim)
    REAL(wp), INTENT(in)     :: pum1(kdim), pvm1(kdim), paphm1(kdim,klevp1)
    REAL(wp), INTENT(in)     :: zbmi(kdim)

    REAL(wp), INTENT(out)    :: pevapi(kdim), pahfsi(kdim), pahfli(kdim), zspeedi(kdim), zt2i(kdim)
    REAL(wp), INTENT(out)    :: zdew2i(kdim), zu10i(kdim), zv10i(kdim)


    INTEGER  :: loidx(kdim), is, nl, jl
    REAL(wp)     :: ztmst, zcons15, ztpfac1, ztpfac2, zcons16, zhuv, zhtq, zephum
    REAL(wp)     :: zcoefi(kdim), zqhfli(kdim), zthfli(kdim)
    REAL(wp)     :: zrat(kdim), zcbn(kdim), zcbs(kdim), zcbu(kdim), zmerge(kdim), zred(kdim)
    REAL(wp)     :: zh2m(kdim), zqs1(kdim), zrh2m(kdim), zcvm3(kdim), zcvm4(kdim)
    REAL(wp)     :: zaph2m(kdim), zqs2(kdim), zq2m(kdim), zfrac(kdim)
    REAL(wp)     :: zmerge1(kdim)
    REAL(wp)     :: ua(kdim)

    !CONSTANTS
    ztmst         = time_step_len
    ztpfac1       = cvdifts
    ztpfac2       = 1._wp / ztpfac1
!    ztpfac3       = 1._wp - ztpfac2
    zcons15       = 1._wp / (grav * ztmst)
    zcons16       = cpd * vtmpc2
    zhuv          = 10._wp * grav
    zhtq          = 2._wp * grav
    zephum        = 5.e-2_wp

    pevapi = 0._wp
    pahfsi = 0._wp
    pahfli = 0._wp
    zspeedi = 0._wp
    zt2i = 0._wp
    zdew2i = 0._wp
    zu10i =  0._wp
    zv10i = 0._wp
    loidx = 0


       !*         5.8     Surface fluxes of heat and moisture
       !*    Moisture fluxes
       !
    WHERE(mask)
       zcoefi(:)      = zcons15 * zcfhi(:) * ztpfac2
       zqhfli(:)      = zcoefi(:) * (zqdif(:) - zqsi(:))
       !*    Sensible heat fluxes
       zthfli(:)      = zcoefi(:) * (ztdif(:) - zcpti(:))
       !    Accumulated sensible heat flux and evaporation
       pevapi(:)      = zqhfli(:)
       pahfsi(:)      = zthfli(:)
       !     Latent heat fluxes
       pahfli(:)      = als * zqhfli(:)
    ENDWHERE
       !
       !     Compute new t2m, t2m_max t2m_min
       !
    WHERE(mask)
       zrat(:)        = zhtq / pgeom1(:)
       zcbn(:)        = LOG(1._wp + (EXP (zbni(:)) - 1._wp) * zrat(:))
       zcbs(:)        = -(zbni(:) - zbhi(:)) * zrat(:)
       zcbu(:)        = -LOG(1._wp + (EXP (zbni(:) - zbhi(:)) - 1._wp) * zrat(:))
       WHERE(zrii(:).GT.0._wp)
          zmerge(:)   = zcbs(:)
       ELSEWHERE
          zmerge(:)   = zcbu(:)
       ENDWHERE
       zred(:)        = (zcbn(:) + zmerge(:)) / zbhi(:)
       zh2m(:)        = zcpti(:) + zred(:) * (zcptgz(:) - zcpti(:))
       zt2i(:)        = (zh2m(:) - zhtq) / (cpd * (1._wp + vtmpc2 * pqm1(:)))  
    ENDWHERE
       !
       !           5.96   2M DEW POINT
       !
    is = 0
    DO jl=1,kdim
      IF(mask(jl)) THEN
        is = is + 1
        loidx(is) = jl
      ENDIF
    ENDDO
   
    CALL lookup_ua_list_spline('postproc_ice(1)',kdim,is,loidx(1), ptm1(1), ua(1))

    DO nl=1,is
      jl = loidx(nl)
      zqs1(jl)      = ua(nl) / papm1(jl)
      zqs1(jl)      = zqs1(jl) / (1._wp- vtmpc1 * zqs1(jl))
      zrh2m(jl)     = MAX(zephum, pqm1(jl) / zqs1(jl))
    ENDDO

    WHERE(mask)
       WHERE(zt2i(:) .GT. tmelt)
          zcvm3(:)    = c3les
          zcvm4(:)    = c4les
       ELSEWHERE
          zcvm3(:)    = c3ies
          zcvm4(:)    = c4ies
       ENDWHERE
       zaph2m(:)      = paphm1(:,klevp1) * &
            (1._wp - zhtq / (rd * zt2i(:) * (1._wp + vtmpc1 * pqm1(:) - zx(:))))
    ENDWHERE

    CALL lookup_ua_list_spline('postproc_ice(2)',kdim,is,loidx(1), zt2i(1), ua(1))

    DO nl=1,is
      jl = loidx(nl)
      zqs2(jl)      = ua(nl) / zaph2m(jl)
      zqs2(jl)      = zqs2(jl) / (1._wp- vtmpc1 * zqs2(jl))
      zq2m(jl)      = zrh2m(jl) * zqs2(jl)
      zfrac(jl)     = LOG(zaph2m(jl) * zq2m(jl) / (c2es * (1._wp + vtmpc1 * zq2m(jl)))) / zcvm3(jl)
      zdew2i(jl)    = MIN(zt2i(jl), (tmelt - zfrac(jl) * zcvm4(jl)) / (1._wp - zfrac(jl)))
    ENDDO
       !
       !*          5.97   10M WIND COMPONENTS, MAX 10M WINDSPEED
       !
    WHERE(mask)
       zrat(:)        = zhuv / pgeom1(:)
       zcbn(:)        = LOG(1._wp + (EXP (zbni(:)) - 1._wp) * zrat(:))
       zcbs(:)        = -(zbni(:) - zbmi(:)) * zrat(:)
       zcbu(:)        = -LOG(1._wp + (EXP (zbni(:) - zbmi(:)) - 1._wp) * zrat(:))
       WHERE(zrii(:).GT.0._wp)
          zmerge1(:)  = zcbs(:)
       ELSEWHERE
          zmerge1(:)  = zcbu(:)
       ENDWHERE
       zred(:)        = (zcbn(:) + zmerge1(:)) / zbmi(:)
       zu10i(:)       = zred(:) * pum1(:)
       zv10i(:)       = zred(:) * pvm1(:)
       zspeedi(:)     = SQRT(zu10i(:)**2 + zv10i(:)**2)
      
    ENDWHERE

  END SUBROUTINE postproc_ice

  SUBROUTINE s_licetemp(            &
       kdim        , psiced         &
     , psni        , palake         &
     , ptsi        , ptrfli         &
     , psofli      , pahfice        &
     , pqres       , pahfcon        &
     , pahfres     , pevapi         &
     , psnow       , pahfsi         &
     , pahfli      , pcvsi          &
     , pfri                         &
     , palsoi      , palsobs        &
     )

    USE mo_param_switches, ONLY: lice
    USE mo_physical_constants, ONLY: tmelt, rhoh2o, alice, cpsno, snicecond, dice,   &
                                 hcapice, rhowlf 
    USE mo_time_control,   ONLY: delta_time
    USE mo_control,        ONLY: lfractional_mask

    INTEGER, INTENT(in)     ::  kdim
    REAL(wp), INTENT(in)    ::  psiced(kdim),  palake(kdim)
    REAL(wp), INTENT(inout) ::  psni(kdim),    ptsi(kdim)          
    REAL(wp), INTENT(in)    ::  ptrfli(kdim),  psofli(kdim)          
    REAL(wp), INTENT(inout) ::  pahfice(kdim), pqres(kdim)
    REAL(wp), INTENT(inout) ::  pahfcon(kdim), pahfres(kdim)
    REAL(wp), INTENT(in)    ::  palsoi(kdim),  palsobs(kdim)
    REAL(wp), INTENT(in)    ::  pahfsi(kdim),  pahfli(kdim),   pcvsi(kdim)           
    REAL(wp), INTENT(in)    ::  pfri(kdim),    psnow(kdim),    pevapi(kdim)          
    REAL(wp) :: zdtime, zcpdt, zcprosn, zcpdte, zdtlfro, zlfrodt, zdtro
    REAL(wp) :: zevsnd, zicefl, zsflx, zmelfac, zsmelt, zsmelres, zsniced, zsnowd
    INTEGER  :: jl

    ! CONSTANTS

    zdtime   = delta_time
    zcprosn  = rhoh2o*cpsno/zdtime
    zcpdt    = hcapice/zdtime
    zdtlfro  = zdtime/rhowlf
    zlfrodt  = rhowlf/zdtime
    zdtro    = zdtime/rhoh2o
    !
    IF (lice) THEN
    !
       DO jl=1,kdim
          IF ((lfractional_mask .AND. palake(jl).GT.0.0_wp) .OR.            &
              (.NOT. lfractional_mask .AND. palake(jl).GE.0.5_wp)) THEN
             IF (psiced(jl).GE.dice) THEN                               ! ice
                zsnowd=psnow(jl)*zdtro
                zevsnd=pcvsi(jl)*pevapi(jl)*zdtro
                psni(jl)=MAX(psni(jl)+zsnowd+zevsnd,0._wp)
                zsniced=psiced(jl)+snicecond*psni(jl)
                zicefl=alice*tmelt/zsniced
                zsflx=ptrfli(jl)+psofli(jl)+pahfsi(jl)+pahfli(jl)
                zcpdte=zcpdt+zcprosn*psni(jl)
                ptsi(jl)=(zcpdte*ptsi(jl)+zsflx+zicefl)/(zcpdte+alice/zsniced)
                IF (ptsi(jl).GT.tmelt .AND. psni(jl).GT.0.0_wp) THEN
                   zmelfac=(zcpdte+alice/zsniced)*zdtlfro               ! m/K
                   zsmelt=MIN(zmelfac*(ptsi(jl)-tmelt),psni(jl))        ! m
                   psni(jl)=psni(jl)-zsmelt                             ! m
                   zsniced=psiced(jl)+snicecond*psni(jl)                ! m
                   ptsi(jl)=ptsi(jl)-zsmelt/zmelfac                     ! K
                   zsmelres=zsmelt*zlfrodt                              ! W/(m**2)
                ELSE
                   zsmelres=0.0_wp
                END IF
                IF (ptsi(jl).GT.tmelt) THEN
                   zcpdte=zcpdt+zcprosn*psni(jl)
                   pqres(jl)=(zcpdte+alice/zsniced)*(ptsi(jl)-tmelt)    ! W/(m**2)
                   ptsi(jl)=tmelt
                ELSE
                   pqres(jl)=0.0_wp
                END IF
                pahfres(jl)=pahfres(jl)+zdtime*pfri(jl)*(pqres(jl)+zsmelres)
                pahfice(jl)=alice*(ptsi(jl)-tmelt)/zsniced
                pahfcon(jl)=pahfcon(jl)+zdtime*pfri(jl)*pahfice(jl)
             ELSE                                                        ! water
                pqres(jl)=0.0_wp
                pahfice(jl)=0.0_wp
                ptsi(jl)=tmelt
                psni(jl)=0.0_wp
             END IF
          END IF
       END DO
       !
       !        Necessary computations if subroutine is bypassed
       !
    ELSE
       DO jl = 1, kdim
          ptsi(jl)=tmelt
          psni(jl)=0.0_wp
       END DO
    END IF

  END SUBROUTINE s_licetemp

  SUBROUTINE s_sicetemp(          &
       kdim        , psiced       &
     , psni        , palake       &
     , pslf                       &
     , ptsi        , ptrfli       &
     , psofli      , pahfice      &
     , pqres       , pevapi       &
     , pahfcon     , pahfres      &
     , pahfsi      , pahfli       &
     , psnow       , pcvsi        &
     , pfri                       &
     , palsoi      , palsobs      &
     )

    USE mo_param_switches, ONLY: lice
    USE mo_physical_constants, ONLY: tmelt, rhoh2o, alice, cpsno, snicecond, dice,    &
                                 hcapice, rhowlf
    USE mo_time_control,   ONLY: delta_time
    USE mo_physc2,         ONLY: ctfreez
    USE mo_control,        ONLY: lcouple, lfractional_mask

    INTEGER,  INTENT(in )    :: kdim
    REAL(wp), INTENT(in)     :: psiced(kdim),     palake(kdim)          
    REAL(wp), INTENT(in)     :: pslf(kdim)
    REAL(wp), INTENT(inout)  :: psni(kdim)
    REAL(wp), INTENT(inout)  :: ptsi(kdim),       pqres(kdim)  
    REAL(wp), INTENT(in)     :: ptrfli(kdim),     psofli(kdim)          
    REAL(wp), INTENT(inout)  :: pahfice(kdim)         
    REAL(wp), INTENT(inout)  :: pahfcon(kdim),    pahfres(kdim)   
    REAL(wp), INTENT(in)     :: pahfsi(kdim),     pahfli(kdim)
    REAL(wp), INTENT(in)     :: psnow(kdim),      pcvsi(kdim)        
    REAL(wp), INTENT(in)     :: palsoi(kdim),     palsobs(kdim)        
    REAL(wp), INTENT(in)     :: pfri(kdim),       pevapi(kdim)

    REAL(wp) :: zdtime, zcpdt, zsniced, zicefl, zsflx
    REAL(wp) :: zsnowd, zevsnd, zmelfac, zsmelt, zsmelres
    REAL(wp) :: zcprosn, zcpdte, zdtlfro, zlfrodt, zdtro
    INTEGER :: jl

    ! CONSTANTS

    zdtime   = delta_time
    zcprosn  = rhoh2o*cpsno/zdtime
    zcpdt    = hcapice/zdtime
    zdtlfro  = zdtime/rhowlf
    zlfrodt  = rhowlf/zdtime
    zdtro    = zdtime/rhoh2o

!-- 2. Compute new skin-temperature
    
    IF (lice) THEN
       !
       IF (.NOT.lcouple) THEN     ! uncoupled model or coupling to slab ocean
          DO jl=1,kdim
             IF ((lfractional_mask .AND. palake(jl).EQ.0.0_wp) .OR.            &
                 (.NOT. lfractional_mask .AND. palake(jl).LT.0.5_wp)) THEN
                IF (psiced(jl).GE.dice) THEN
                   zsnowd=psnow(jl)*zdtro
                   zevsnd=pcvsi(jl)*pevapi(jl)*zdtro
                   psni(jl)=MAX(psni(jl)+zsnowd+zevsnd,0.0_wp)
                   psni(jl)=MIN(psni(jl),0.5_wp)
                   zsniced=psiced(jl)+snicecond*psni(jl)
                   zicefl=alice*ctfreez/zsniced
                   zsflx=ptrfli(jl)+psofli(jl)+pahfsi(jl)+pahfli(jl)
                   zcpdte=zcpdt+zcprosn*psni(jl)
                   ptsi(jl)=(zcpdte*ptsi(jl)+zsflx+zicefl)/(zcpdte+alice/zsniced)
                   IF (ptsi(jl).GT.tmelt .AND. psni(jl).GT.0.0_wp) THEN
                      zmelfac=(zcpdte+alice/zsniced)*zdtlfro                ! m/K
                      zsmelt=MIN(zmelfac*(ptsi(jl)-tmelt),psni(jl))         ! m
                      psni(jl)=psni(jl)-zsmelt                              ! m
                      zsniced=psiced(jl)+snicecond*psni(jl)                 ! m
                      ptsi(jl)=ptsi(jl)-zsmelt/zmelfac                      ! K
                      zsmelres=zsmelt*zlfrodt                               ! W/(m**2)
                   ELSE
                      zsmelres=0.0_wp
                   END IF
                   IF (ptsi(jl).GT.tmelt) THEN
                      zcpdte=zcpdt+zcprosn*psni(jl)
                      pqres(jl)=(zcpdte+alice/zsniced)*(ptsi(jl)-tmelt)     ! W/(m**2)
                      ptsi(jl)=tmelt
                   ELSE
                      pqres(jl)=0.0_wp
                   END IF
                   pahfice(jl)=alice*(ptsi(jl)-ctfreez)/zsniced
                   pahfcon(jl)=pahfcon(jl)+zdtime*pfri(jl)*pahfice(jl)
                   pahfres(jl)=pahfres(jl)+zdtime*pfri(jl)*(pqres(jl)+zsmelres)
                ELSE   ! psiced < 0.1
                   pqres(jl)=0.0_wp
                   pahfice(jl)=0.0_wp
                   ptsi(jl)=tmelt
                   psni(jl)=0.0_wp
                END IF  ! psiced
             END IF     ! palake
          END DO
       ELSE                        ! coupled model
          DO jl=1,kdim
             IF ((lfractional_mask .AND. palake(jl).EQ.0.0_wp) .OR.            &
                 (.NOT. lfractional_mask .AND. palake(jl).LT.0.5_wp)) THEN
                IF (pslf(jl).LT.1.0_wp) THEN   ! not completely land
                   zsniced=MAX(psiced(jl),dice)+snicecond*psni(jl)
                   zicefl=alice*ctfreez/zsniced
                   zsflx=ptrfli(jl)+psofli(jl)+pahfsi(jl)+pahfli(jl)
                   zcpdte=zcpdt+zcprosn*psni(jl)
                   ptsi(jl)=(zcpdte*ptsi(jl)+zsflx+zicefl)/(zcpdte+alice/zsniced)
                   IF (ptsi(jl).GT.tmelt) THEN
                      pqres(jl)=(alice/zsniced+zcpdte)*(ptsi(jl)-tmelt)
                      ptsi(jl)=tmelt
                   ELSE
                      pqres(jl)=0.0_wp
                   END IF
                   pahfice(jl)=alice*(ptsi(jl)-ctfreez)/zsniced
                   pahfcon(jl)=pahfcon(jl)+zdtime*pfri(jl)*pahfice(jl)
                   pahfres(jl)=pahfres(jl)+zdtime*pfri(jl)*pqres(jl)
                ELSE
                   pqres(jl)=0.0_wp
                   pahfice(jl)=0.0_wp
                   ptsi(jl)=tmelt
                END IF  ! pslf
             END IF  ! palake
          END DO
       END IF   ! cases: uncoupled, coupled
!
!       Necessary computations if subroutine is bypassed
!
    ELSE
       DO jl = 1, kdim
          ptsi(jl)=tmelt
          psni(jl)=0.0_wp
       END DO
    END IF

  END SUBROUTINE s_sicetemp

  SUBROUTINE ice_rad(               &
       kdim          , mask         &
      ,ztrdown      , ptsi         &
      ,pi0          , ptrsol       &
      ,palsoi       , palbedo      &
      ,psofli       , ptrfli       &
      )

    USE mo_radiation_parameters,         ONLY: cemiss
    USE mo_physical_constants,           ONLY: stbo

    INTEGER, INTENT(in)    :: kdim
    LOGICAL, INTENT(in)    :: mask(kdim)
    REAL(wp),    INTENT(in)    :: ptsi(kdim), pi0(kdim), ptrsol(kdim)
    REAL(wp),    INTENT(in)    :: palsoi(kdim), ztrdown(kdim), palbedo(kdim)
    REAL(wp),    INTENT(out)   :: psofli(kdim), ptrfli(kdim)

    REAL(wp) :: zflxs(kdim)

    psofli = 0._wp
    ptrfli = 0._wp

    WHERE(mask)    
       ptrfli(:)  = ztrdown(:) - cemiss * stbo * ptsi(:)**4
       zflxs(:)   = pi0(:) * ptrsol(:)     
       psofli(:)  = (1._wp - palsoi(:)) * zflxs(:) / (1._wp - palbedo(:))
    END WHERE

  END SUBROUTINE ice_rad
  !
  !---------------------------------------------------------------------------------------------------------
SUBROUTINE meltpond ( kproma                                         & 
                    , psicepdw,    psicepdi,   ptsicepdi             & 
                    , psiced,      pseaice,    psicemin              & 
                    , pameltdepth, pameltfrac, pqres,     psni       & 
                    ) 
! 
!  Description 
! 
!  Calculation of melt ponds on sea-ice 
! 
!  Method:  
! 
!  *meltpond* called from physc
!                               
!  *physc* called gpc 
! 
!  Authors:   
!  E. Roeckner, MPI, July 2007 
!  ----------------------------------------------------------------- 
  USE mo_kind,           ONLY: wp 
  USE mo_physical_constants, ONLY: tmelt, rhoiw, rhoilf, rhowlf, dice, dicepond,      &
                               dicelim, dpondmin
  USE mo_time_control,   ONLY: delta_time
! 
  IMPLICIT NONE 
! 
  INTEGER, INTENT (IN):: kproma 
! 
! Arguments 
! 
  REAL(wp) :: psicepdw(kproma),    psicepdi(kproma),    ptsicepdi(kproma)    & 
            , psiced(kproma),      pseaice(kproma)                           & 
            , psicemin(kproma),    pameltdepth(kproma), pameltfrac(kproma)   & 
            , psni(kproma),        pqres(kproma)                      
! 
! psicepdw     melt pond depth on sea-ice   
! psicepdi     ice thickness on melt pond  <= rhoh2o/rhoice * psicepdw
! ptsicepdi    ice temperature of frozen melt pond 
! psiced       sea-ice thickness
! pseaice      sea-ice fraction oftotal water fraction  
! psicemin     minimum annual ice thickness (used in the definition of MYI/FYI)  
! pameltdepth  total melt pond depth 
! pameltfrac   fractional area of melt ponds on sea-ice (parameterized) 
! psni         snow depth on sea-ice 
! pqres        residual heat flux on ice = melt rate of ice
! 
! Local variables 
! 
  REAL(wp) :: zdtime, zdtrwlf, zrilfdt, zeps, zaiqre, zpondnew, zgrowth    & 
            , zd, z1, z2, z3, z4, z5, z6, z7, z8, z9  
! 
  INTEGER  :: jl    
!
  LOGICAL :: lopondopen(kproma)      !  .true. if melt pond open        
!      
! Executable statements 
! 
! 1. Set up constants 
! 
  zdtime   = delta_time               
  zdtrwlf  = zdtime/rhowlf 
  zrilfdt  = rhoilf/zdtime 
  zeps     = EPSILON(1.0_wp) 
  zgrowth  = dice/4.7304E7_wp  ! growth rate of psicemin (m/s = �dice� per 1.5 years)
! 
! 2. Formation of a new melt pond or further growth or closing by freezing 
! 
  DO jl=1, kproma 
!    
    IF (psiced(jl) < dice) THEN
!      First-year ice (valid in the next melting period)
       psicemin(jl) = 0.0_wp
    ELSE
!      Increase of minimum ice thickness with time: IF the FYI threshold is not
!      reached in the next summer, pcicemin > zdseaice after about 18 months 
!      ==> Multi-year ice in the following summer. 
       psicemin(jl) = MIN(1.0_wp,psicemin(jl)+zdtime*zgrowth)
    END IF
!   
    IF (pseaice(jl) > zeps .AND. psiced(jl) >= dice) THEN
! 
!      a) New formation of a melt pond (no snow on sea-ice and melt rate > 0)
!
       IF (pameltdepth(jl) < zeps .AND. psni(jl) < zeps .AND. pqres(jl) > zeps) THEN
!
          pameltdepth(jl) = MIN(psiced(jl),zdtrwlf*pqres(jl))
          psicepdw(jl)    = pameltdepth(jl) 
          psicepdi(jl)    = 0.0_wp
          ptsicepdi(jl)   = tmelt 
          lopondopen(jl)  = .true.
          IF (pameltdepth(jl) < zeps) THEN
             pameltdepth(jl) = 0.0_wp
             psicepdw(jl)    = 0.0_wp
             lopondopen(jl)  = .false. 
          END IF
!
!      b) A melt pond already exists and may or may not grow further depending only on
!         the melt rate, but independent of whether there is snow on top of the sea-ice.
! 
       ELSE IF (pameltdepth(jl) >= zeps .AND. psicepdi(jl) <= dicelim) THEN
!
          zaiqre = pqres(jl)-psicepdi(jl)*zrilfdt
!
          IF (zaiqre > zeps) THEN                   ! further growth of melt pond
             pameltdepth(jl) = pameltdepth(jl)+zdtrwlf*zaiqre 
             pameltdepth(jl) = MIN(psiced(jl), pameltdepth(jl)) 
             psicepdw(jl)    = pameltdepth(jl)      ! Input for *meltpond_ice* 
             psicepdi(jl)    = 0.0_wp 
             ptsicepdi(jl)   = tmelt 
             lopondopen(jl)  = .true.
          ELSE IF (psicepdi(jl) < dicepond) THEN    ! pond remains open
             lopondopen(jl)  = .true.  
          ELSE                                      ! pond temporarily closed 
             lopondopen(jl)  = .false. 
          END IF
!
!       c) New melt pond formation is not possible because the melt rate is zero 
!          or there is snow on top of the sea-ice. The pond is closed when the pond-ice 
!          exceeds the 'dicelim' threshold
! 
       ELSE                                     
          pameltdepth(jl) = 0.0_wp 
          psicepdw(jl)    = 0.0_wp 
          psicepdi(jl)    = 0.0_wp 
          ptsicepdi(jl)   = tmelt 
          lopondopen(jl)   = .false.    
       END IF
!
! 3. Seepage
!
       zpondnew = psicepdw(jl) - cseep*zdtime
       IF (zpondnew > 0.0_wp) THEN
          psicepdw(jl)    = zpondnew
          pameltdepth(jl) = zpondnew + rhoiw*psicepdi(jl)
       ELSE
          pameltdepth(jl) = 0.0_wp 
          psicepdw(jl)    = 0.0_wp 
          psicepdi(jl)    = 0.0_wp 
          ptsicepdi(jl)   = tmelt 
          lopondopen(jl)   = .false. 
       END IF     
! 
! 4. Parameterization of melt-pond fraction depending on open pond depth 
!    and ice type (multi-year / first-year ice) 
!     
       IF (lopondopen(jl) .AND. pameltdepth(jl) >= dpondmin) THEN  
!
          zd = pameltdepth(jl)
          IF (psicemin(jl) > dice) THEN                   ! multi-year ice 
             z1 = -0.0072463_wp * zd**8 
             z2 =  0.1443800_wp * zd**7 
             z3 = -1.1914000_wp * zd**6 
             z4 =  5.2599500_wp * zd**5 
             z5 = -13.371010_wp * zd**4 
             z6 =  19.530300_wp * zd**3 
             z7 = -15.270190_wp * zd**2 
             z8 =  5.2667400_wp * zd 
             z9 = -0.1254900_wp 
             pameltfrac(jl) = MAX(zeps,z1+z2+z3+z4+z5+z6+z7+z8+z9)  
          ELSE                                            ! first-year ice 
             z1 = 30.0_wp 
             z2 =  2.5_wp 
             pameltfrac(jl) = MAX(zeps,0.5_wp*(1.0_wp+TANH(z1*zd-z2))) 
          END IF 
          pameltfrac(jl) = pameltfrac(jl)*cmpfac
       ELSE 
          pameltfrac(jl) = 0.0_wp     
       END IF
    ELSE 
! 
! 5. Land, lake, or open ocean 
! 
       pameltdepth(jl) = 0.0_wp 
       pameltfrac (jl) = 0.0_wp 
       psicepdw(jl)    = 0.0_wp 
       psicepdi(jl)    = 0.0_wp 
       ptsicepdi(jl)   = tmelt 
       lopondopen(jl)  = .false. 
    END IF 
  END DO 
! 
!  ----------------------------------------------------------------------- 
! 
     RETURN 
  END SUBROUTINE meltpond
  !
  !---------------------------------------------------------------------------------------------------------
SUBROUTINE meltpond_ice ( kproma                                      &
                        , pameltdepth, psicepdw,  psicepdi            &
                        , ptsicepdi,   psicepres, psiced,    pevapi   &
                        , pahflw,      pahfsw,    ptrflw,    psoflw   &
                        , pahfli,      pahfsi,    ptrfli,    psofli   &
                        , palsow,      palsom,    palsoi              &
                        )
!
!  Description
!
!  Prognostic calculation of ice formation/melting on melt ponds
!
!  Method: 
!
!  *meltpond_ice* called from mo_surface
!  *mo_surface* called from vdiff
!
!  Authors:  
!  E. Roeckner, MPI, July 2007
!  ---------------------------------------------------------------------
!
  USE mo_kind,           ONLY: wp
  USE mo_physical_constants, ONLY: tmelt, alice, rhoice, rhoiw, rhoilf, &
                                   rhoh2o, cpliq, dice, dicepond,       &
                                   dpondmin, albpondi, hcapicep 
  USE mo_time_control,   ONLY: delta_time
!
  IMPLICIT NONE
!
  INTEGER, INTENT (IN):: kproma
!
! Arguments
!
  REAL(wp) :: pameltdepth(kproma), psicepdw(kproma), psicepdi(kproma) &
            , ptsicepdi(kproma),   psicepres(kproma)                  &
            , psiced(kproma),      pevapi(kproma)                     &
            , pahflw(kproma),      pahfsw(kproma)                     &
            , ptrflw(kproma),      psoflw(kproma)                     &
            , pahfli(kproma),      pahfsi(kproma)                     &
            , ptrfli(kproma),      psofli(kproma)                     &
            , palsow(kproma),      palsom(kproma), palsoi(kproma)
!
! pameltdepth  total pond depth = psicepdw + zrhoiw*psicepdi
! psicepdw     melt pond depth on sea ice  
! psicepdi     ice thickness for frozen melt pond
! ptsicepdi    ice temperature of frozen melt pond
! psicepres    residual heat flux 
! psiced       sea ice thickness
! pevapi       sublimation rate
! pahflw       surface latent heat flux (water)
! pahfsw       surface sensible heat flux (water)                 
! ptrflw       surface net longwave radiation (water)
! psoflw       surface net shortwave radiation (water)
! pahfli       surface latent heat flux (sea ice)
! pahfsi       surface sensible heat flux (sea ice)                 
! ptrfli       surface net longwave radiation (sea ice)
! psofli       surface net shortwave radiation (sea ice)
! palsow       albedo of water
! palsom       albedo of melt ponds
! palsoi       mean sea-ice albedo
!
! Local variables
!
  REAL(wp) :: zdtime, zdtrilf, zcpdt, zeps                                   & 
           ,  zfreez, zpondcap, zmcapdt, zmcaprilf, zfluxi, zsflx, zfluxw    &
           ,  zts, zfres, zconhflx, zhi, zsoflm, zsofli, zsubice
!
  INTEGER  :: jl   
!     
! Executable statements
!
! 1. Set up constants
!
  zdtime   = delta_time                
  zdtrilf  = zdtime/rhoilf
  zcpdt    = hcapicep/zdtime
  zeps     = EPSILON(1.0_wp)
  zfreez   = -dicepond/zdtrilf
!
! 2. Ice formation and skin temperature of ice on melt pond
!
  DO jl=1, kproma
!
   IF (psiced(jl) >= dice .AND. pameltdepth(jl) >= dpondmin                   &
                              .AND. psicepdw(jl) > zeps) THEN
!
      IF (psicepdi(jl) < dicepond) THEN         ! no ice on melt pond
         zpondcap  = rhoh2o*cpliq*pameltdepth(jl)    
         zmcapdt   = zdtime/zpondcap
         zmcaprilf = zpondcap/rhoilf
!     net solar radiation over meltponds
         zsoflm    = psoflw(jl)*(1.0_wp-palsom(jl))/(1.0_wp-palsow(jl))
         zsoflm    = MIN(zsoflm, psoflw(jl))
         zfluxw = pahflw(jl)+pahfsw(jl)+ptrflw(jl)+zsoflm
         zts    = tmelt+zmcapdt*(zfluxw+psicepres(jl))
         psicepres(jl) = 0.0_wp  
         ptsicepdi(jl) = tmelt                ! initial value
         zfres = (zts-tmelt)/zmcapdt        
!
         IF (zfres < 0.0_wp) THEN             ! check ice formation
            IF (zfres <= zfreez) THEN         ! ice formation (>= zdpice)
               psicepdi(jl)  = MIN(pameltdepth(jl)/rhoiw,zmcaprilf*(tmelt-zts))
               psicepdw(jl)  = psicepdw(jl)-psicepdi(jl)*rhoiw
               psicepdw(jl)  = MAX(0.0_wp, psicepdw(jl))
            ELSE                              ! ice formation not realized
               psicepres(jl) = zfres          ! residual < 0 for next time step
            END IF
         END IF
!
!  -----------------------------------------------------------------------------------
      ELSE                                    ! psicepdi >= dicepond
!
!        Skin temperature of pond-ice (decides on freezing/melting of pond-ice)
!
         zsofli        = psofli(jl)*(1.0_wp-albpondi)/(1.0_wp-palsoi(jl))
         zsflx         = ptrfli(jl)+zsofli+pahfsi(jl)+pahfli(jl)
         zsflx         = zsflx+psicepres(jl)+alice*tmelt/psicepdi(jl)
         psicepres(jl) = 0.0_wp
         ptsicepdi(jl) = (zcpdt*ptsicepdi(jl)+zsflx)/(zcpdt+alice/psicepdi(jl))
         IF (ptsicepdi(jl) > tmelt) THEN                    ! melt rate �psicepres� > 0
            psicepres(jl) = (zcpdt+alice/psicepdi(jl))*(ptsicepdi(jl)-tmelt)
            ptsicepdi(jl) = tmelt
         END IF
!        
!        New ice thickness (freezing if zconhflx < 0 or melting if psicepres > 0)
!
         zconhflx      = alice*(ptsicepdi(jl)-tmelt)/psicepdi(jl) 
         zsubice       = pevapi(jl)*zdtime/rhoice
         zfluxi        = zdtrilf*(zconhflx+psicepres(jl))+zsubice     ! thickness change
         zhi           = psicepdi(jl)-zfluxi    ! preliminary thickness
!
         IF (zhi >= dicepond) THEN              ! ice growth or melting
            psicepdi(jl)  = MIN(pameltdepth(jl)/rhoiw, zhi)
            psicepdw(jl)  = MAX(0.0_wp, psicepdw(jl)+rhoiw*zfluxi)
            psicepres(jl) = 0.0_wp
         ELSE IF (zhi <= 0.0_wp) THEN           ! complete melting
            psicepdw(jl)  = psicepdw(jl)+rhoiw*psicepdi(jl)
            psicepdi(jl)  = 0.0_wp
            psicepres(jl) = -zhi/zdtrilf        ! excess heat if zhi < 0
         ELSE                                   ! incomplete melting
            psicepdw(jl)  = psicepdw(jl)+(psicepdi(jl)-dicepond)*rhoiw
            psicepdi(jl)  = dicepond 
            psicepres(jl) = (dicepond-zhi)/zdtrilf  !         > 0.0 (melting next timestep)
         END IF
      END IF     ! ice thickness on ponds
    END IF       ! palake
  END DO
!  
!  ----------------------------------------------------------------------------------
!
     RETURN
  END SUBROUTINE meltpond_ice

END MODULE mo_surface_ice
