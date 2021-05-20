!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_tr_gp_ffsl

  ! Joerg Behrens, DKRZ, February 2010, initial version
  ! Luis Kornblueh, MPIM, February 2010, update to conform to programming rules
  ! Joerg Behrens, DKRZ, November 2010, debug case nproma>=nglon*nglat
  ! Luis Kornblueh, MPIM, February 2012, fixed INTENT bug 
  !
  !
  ! Experimental code for the transposition between gp- and ffsl-decomposition
  ! this code is not intended for production purposes but to explore the scaling
  ! of some aspects of the MPI communication and the combination with OpenMP
  !
  ! usage (gridpoint to ffsl):
  ! ==========================
  ! CALL pack_gp(tracer1_gp_field)
  ! CALL pack_gp(tracer2_gp_field)
  ! pack some more tracers
  !
  ! CALL gp2ffsl
  !
  ! CALL unpack_ffsl(tracer1_ffsl_field)
  ! CALL unpack_ffsl(tracer2_ffsl_field)
  ! unpack some more tracers
  !
  !
  ! usage (ffsl to gridpoint):
  ! ==========================
  !
  ! CALL pack_ffsl(tracer1_ffsl_field)
  ! CALL pack_ffsl(tracer2_ffsl_field)
  ! pack some more tracers
  !
  ! CALL ffsl2gp
  !
  ! CALL unpack_gp(tracer1_gp_field)
  ! CALL unpack_gp(tracer2_gp_field)
  ! unpack some more tracers
  !
  ! restrictions:
  ! =============
  !    - no dynamic number of tracers yet (tracer_max set once at init)
  !    - 3d and 2d transformations cannot be mixed yet
  !
  !





  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: finish
  USE mo_decomposition, ONLY: ldc => local_decomposition, &
       &                      gdc => global_decomposition
  USE mo_control,       ONLY: ltimer
  USE mo_mpi,           ONLY: p_nprocs, p_all_comm, p_real_dp, &
                              p_isend, p_recv, p_irecv, p_send, p_wait
  USE mo_real_timer
  USE mo_utils
  USE mpi

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pack_gp, unpack_ffsl, gp2ffsl
  PUBLIC :: unpack_gp, pack_ffsl, ffsl2gp
  PUBLIC :: init_tr_gp_ffsl
  
  INTERFACE pack_gp
    MODULE PROCEDURE pack_gp3
    MODULE PROCEDURE pack_gp2
  END INTERFACE

  INTERFACE unpack_gp
    MODULE PROCEDURE unpack_gp3
    MODULE PROCEDURE unpack_gp4
    ! slow reference implementation for the 4d (multi-tracer) case
    ! using unpack_gp3:
    !!MODULE PROCEDURE unpack_gp4_slow 
  END INTERFACE

  INTERFACE pack_ffsl
    MODULE PROCEDURE pack_ffsl3
    MODULE PROCEDURE pack_ffsl4
  END INTERFACE

  INTERFACE unpack_ffsl
    MODULE PROCEDURE unpack_ffsl3
    MODULE PROCEDURE unpack_ffsl2
  END INTERFACE

  !number of possible tracers
  INTEGER :: tracer_max

  !message tags:
  INTEGER, PARAMETER :: gp2ffsl_tag = 1001, ffsl2gp_tag = 1002

  TYPE msg_buf_2d_t
    REAL(dp),ALLOCATABLE:: buf(:,:) ! (nglon*nglat/2,sendlevels*tracers)
    INTEGER :: offset
  END TYPE msg_buf_2d_t

  ! ffsl message logic and buffer; wording (src,dest) fits the gp->ffsl 
  ! transposition:
  INTEGER, ALLOCATABLE :: kstack(:,:)
  INTEGER, ALLOCATABLE :: kstack_n(:)
  INTEGER :: kstack_max = 0
  INTEGER :: dest_n, src_n
  TYPE(msg_buf_2d_t), ALLOCATABLE :: gp_msg(:) ! (nglon*nglat/2,sendlevels*tracers,idest)
  INTEGER :: gp_unpack_p0, gp_pack_p0          ! read-, write-offset for the next tracer field in gp_msg
  INTEGER :: ffsl_unpack_p0, ffsl_pack_p0      ! read-,write-offset for the next tracer field in ffsl_msg
  TYPE(msg_buf_2d_t), ALLOCATABLE :: ffsl_msg(:) !(nglon*nglat/2,recvlevels*tracers,isrc)
  INTEGER :: tr_dim=0                          ! do we use a 2d or 3d transposition

  ! aux info for splitting blocked gp data into two regions 
  INTEGER :: bsplit_n1, bsplit_nb1, bsplit_nz1, bsplit_ncross
  INTEGER :: bsplit_n2, bsplit_nb2, bsplit_nz2
  ! type for indirect access to meta data about destination (or source) related data
  TYPE ix_type
    INTEGER :: xpe ! mpi-rank  of communication-partner within dcomm
    INTEGER :: idx ! index winthin global decomposition table gdc
    INTEGER :: ir  ! gp-region_flag: 1=exchange 1. region; 2=exchange 2. region; 3=exchange both
    INTEGER :: nr  ! number of transferred gp-regions
    INTEGER :: jr  ! ffsl-region flag: 1: 1. substripe, 2: 2. substripe (gp-space-orientation), 3: both (restriction: 1=>1, 2=>2)
    INTEGER :: hsize ! one-region horizontal data size (nglon*nglat/2), i.e., size of first dim of msg buffer
  END TYPE ix_type
  TYPE(ix_type), ALLOCATABLE :: dest(:), src(:) ! dest contains local region size, src contains remote region sizes

  ! timer:
  INTEGER :: timer_gp2ffsl, timer_ffsl2gp
  INTEGER :: timer_pack_gp3, timer_pack_gp2, timer_unpack_ffsl3, timer_unpack_ffsl2
  INTEGER :: timer_pack_ffsl3, timer_unpack_gp3, timer_unpack_gp4

  ! do we want to use collectives?
  LOGICAL, PARAMETER :: use_mpi_collectives = .FALSE.

  ! debug control:
  ! we are generous with conditional checks on the fortran level and rely on the
  ! compiler to remove dead code if we set a parameter-switch to .false. 
  LOGICAL, PARAMETER :: ldebug  = .FALSE. 
  LOGICAL, PARAMETER :: lcheck  = .FALSE. 
  LOGICAL, PARAMETER :: lassert = .FALSE. 

  ! module state control:
  LOGICAL, SAVE :: require_init = .TRUE.
  
  INTEGER, SAVE :: d_nprocs, nproca, nprocb, dcomm_size
  INTEGER, SAVE :: my_idx, my_pe, pe_shift
  INTEGER, SAVE, ALLOCATABLE :: pe2idx(:) ! maps dcomm rank to gdc index

  ! domain communicator to separate non-debug demain from debug domain,
  INTEGER, SAVE :: dcomm = 0

CONTAINS

  !------------------------------------------------------------------------------

  SUBROUTINE init_tr_gp_ffsl(default_buffer_size)    
    INTEGER, INTENT(in), OPTIONAL :: default_buffer_size

    CHARACTER(len=*), PARAMETER :: context = 'mo_tr_gp_ffsl:init_tr_gp_ffsl'
    CHARACTER(len=mpi_max_error_string) :: err_str
    INTEGER :: err_len
    INTEGER :: color, key
    INTEGER :: ip1, ip2, idx, pe, iastat, ierror, error
    INTEGER :: pe_shift_vec(0:ldc%d_nprocs-1)

    IF (.NOT. require_init) RETURN
    require_init = .FALSE.

    timer_pack_gp2      = new_timer('tr_gp_ffsl:pack_gp2')
    timer_pack_gp3      = new_timer('tr_gp_ffsl:pack_gp3')
    timer_unpack_ffsl2  = new_timer('tr_gp_ffsl:unpack_ffsl2')
    timer_unpack_ffsl3  = new_timer('tr_gp_ffsl:unpack_ffsl3')
    timer_gp2ffsl       = new_timer('tr_gp_ffsl:gp2ffsl')
    
    timer_pack_ffsl3    = new_timer('tr_gp_ffsl:pack_ffsl3')
    timer_unpack_gp3    = new_timer('tr_gp_ffsl:unpack_gp3')
    timer_unpack_gp4    = new_timer('tr_gp_ffsl:unpack_gp4')
    timer_ffsl2gp       = new_timer('tr_gp_ffsl:ffsl2gp')
    
    IF (PRESENT(default_buffer_size)) THEN
      tracer_max = default_buffer_size
    ELSE
      tracer_max = 5
    ENDIF

    d_nprocs = ldc%d_nprocs
    nproca   = ldc%nproca
    nprocb   = ldc%nprocb
    IF (nproca*nprocb /= d_nprocs) CALL finish(context,'nproca*nprocb /= d_nprocs')
    IF (nprocb /= SIZE(ldc%mapmesh,1)) CALL finish(context,'nprocb /= SIZE(ldc%mapmesh,1)')
    IF (nproca /= SIZE(ldc%mapmesh,2)) CALL finish(context,'nproca /= SIZE(ldc%mapmesh,2)')
    IF (SIZE(gdc) < p_nprocs) CALL finish(context, 'SIZE(gdc) < p_nprocs')

    key   = 0
    color = ldc%spe 
    CALL mpi_comm_split(p_all_comm, color, key, dcomm, ierror)
    IF (ierror /= MPI_SUCCESS) CALL finish(context, 'mpi_comm_split failed')

    CALL MPI_COMM_SIZE(dcomm, dcomm_size, ierror)
    IF (ierror /= MPI_SUCCESS) CALL finish(context, 'MPI_COMM_SIZE failed')
    IF (dcomm_size /= d_nprocs) CALL finish(context, 'dcomm_size /= d_nprocs')

    CALL MPI_COMM_RANK(dcomm, my_pe, ierror)
    IF (ierror /= MPI_SUCCESS) CALL finish(context, 'MPI_COMM_RANK failed')
    IF (lcheck) CALL check(__LINE__, my_pe, 0, d_nprocs-1)
!!$    IF (my_pe == 0) WRITE(0,*) context//': create dcomm of size ', dcomm_size

    IF (ldc%epe - ldc%spe + 1 /= d_nprocs) CALL finish(context, 'bad epe-spe diff')
    pe_shift=ldc%pe-my_pe
    IF (pe_shift /= ldc%spe-1) CALL finish(context, 'bad pe_shift (1)')

    IF (pe_shift < 0 .OR. pe_shift > p_nprocs-1) CALL finish(context, 'bad pe_shift (2)')

    pe_shift_vec = HUGE(pe_shift)
    ierror = 0 
    CALL MPI_GATHER(pe_shift, 1, mpi_integer, &
         & pe_shift_vec, 1, mpi_integer, &
         & 0, dcomm, ierror) 
    IF (ierror /= MPI_SUCCESS) THEN
      error = ierror; ierror = 0;
      CALL mpi_error_string(error, err_str, err_len, ierror)
      WRITE(0,*) 'err_str =', TRIM(ADJUSTL(err_str))
      CALL finish(context, 'MPI_GATHER failed')
    ENDIF
    IF (my_pe == 0) THEN
      IF (ANY(pe_shift_vec /= pe_shift) ) CALL finish(context, 'bad pe_shift_vec')
    ENDIF
    
    ALLOCATE(pe2idx(0:dcomm_size-1))

    DO ip2=1, nproca
      DO ip1=1, nprocb
        idx = ldc%mapmesh(ip1,ip2)
        IF (lcheck) CALL check(__LINE__, idx, 1, p_nprocs)
        pe  = gdc(idx)%pe - pe_shift
        IF (lcheck) CALL check(__LINE__, pe, 0, d_nprocs-1)
        pe2idx(pe) = idx
      ENDDO
    ENDDO
    my_idx = pe2idx(my_pe)

    DO pe = 0, d_nprocs-1
      idx = pe2idx(pe)
      IF (gdc(idx)%pe - pe_shift /= pe) CALL finish(context, 'bad pe2idx')
    ENDDO

    ALLOCATE(kstack(ldc%nlev,0:d_nprocs-1), kstack_n(0:d_nprocs-1), stat = iastat)
    IF (iastat /= 0) CALL finish(context,'alloc failed')
    
  END SUBROUTINE init_tr_gp_ffsl

  !------------------------------------------------------------------------------

  SUBROUTINE split_blocks(n1, nb1, nz1, ncross, n2, nb2, nz2)
    INTEGER, INTENT(out) :: n1, nb1, nz1,  ncross, n2, nb2, nz2

    INTEGER :: na, nglon, nglh(2)
    
    ! todo: results are const., move up to module variables

    na     = ldc%nproma
    nglon  = ldc%nglon
    nglh   = ldc%nglh

    
    ! 1. segment:
    n1     = nglon*nglh(1)
    nb1    = (n1-1)/na+1   ! number of blocks
    nz1    = n1-(nb1-1)*na ! last blocksize, if nz1 < na then we have an equator crossing block 
    IF (lcheck) CALL check(__LINE__, nz1, 1, na)
    IF (ldebug .AND. my_pe == 0) WRITE(0,*) 'split_blocks: na, n1, nb, nz =',na, n1, nb1, nz1


    ! 2. segment:
    ncross  = MIN( na-nz1, nglon*nglh(2) ) ! 2. part of equator crossing block 
    n2     = nglon*nglh(2) - ncross
    IF (n2 > 0) THEN
      nb2    = (n2-1)/na+1 ! this fails for n2 == 0      
      nz2    = n2-(nb2-1)*na
      IF (lcheck) CALL check(__LINE__, nz2, 1, na)
    ELSE
      nb2=0
      nz2=0
    ENDIF
    IF (ldebug .AND. my_pe == 0) WRITE(0,*) 'split_blocks: ncross, n2, nb2, nz2 =',ncross, n2, nb2, nz2

  END SUBROUTINE split_blocks

  !------------------------------------------------------------------------------

  SUBROUTINE pack_gp3(x_gp)
    REAL(dp), INTENT(in) :: x_gp(:,:,:) ! grid point field (nproma, nlev, nblock) 

    INTEGER :: idest, dest_pe, p0, ir
    INTEGER :: kk, k
    INTEGER :: n1, n2, nb1, nb2, na, nz1, nz2, ncross

    IF (require_init) CALL init_tr_gp_ffsl

    IF (ltimer) CALL timer_start(timer_pack_gp3)

    IF (kstack_max == 0) CALL prepare_tr

!!  IF (gp_pack_p0 > 300) &      !!++mgs 05/04/2011
!!       CALL finish('mo_tr_gp_ffsl:pack_gp3','unexpected huge buffer pos')
    IF (gp_pack_p0 >= tracer_max) &
         CALL finish('mo_tr_gp_ffsl:pack_gp3','dynamic sbuf not supported yet')

    IF (tr_dim == 0) THEN
      tr_dim = 3
    ELSE
      IF (tr_dim /= 3) &
           CALL finish('mo_tr_gp_ffsl:pack_gp3','cannot mix 2d and 3d packing yet')
    ENDIF

    na = ldc%nproma
    !CALL split_blocks(n1, nb1, nz1, ncross, n2, nb2, nz2)
    n1     = bsplit_n1    
    nb1    = bsplit_nb1   
    nz1    = bsplit_nz1   
    ncross = bsplit_ncross
    n2     = bsplit_n2    
    nb2    = bsplit_nb2   
    nz2    = bsplit_nz2   

    DO idest = 1, dest_n
      dest_pe = dest(idest)%xpe
      p0      = gp_pack_p0*kstack_n(dest_pe)
      ir = dest(idest)%ir

      IF (ir == 1) THEN
        CALL my_tr1
      ELSEIF (ir == 2) THEN
        CALL my_tr2
      ELSEIF (ir == 3) THEN
        p0 = 2*p0
        CALL my_tr1
        p0 = p0+kstack_n(dest_pe)
        CALL my_tr2
      ELSE
        CALL finish ('mo_tr_gp_ffsl:pack_gp3','unexpected case 1',1)
      ENDIF
    ENDDO

    gp_pack_p0 = gp_pack_p0+1

    IF (ltimer) CALL timer_stop(timer_pack_gp3)

  CONTAINS

    SUBROUTINE my_tr1

      INTEGER :: i0, ia, ib

      DO kk = 1, kstack_n(dest_pe)
        k = kstack(kk,dest_pe)

        i0 = 0
        DO ib = 1, nb1-1
          DO ia = 1, na
            gp_msg(idest)%buf(i0+ia,p0+kk) = x_gp(ia,k,ib)
          ENDDO
          i0 = i0+na
        ENDDO
        ib = nb1
        DO ia = 1, nz1
          gp_msg(idest)%buf(i0+ia,p0+kk) = x_gp(ia,k,ib)
        ENDDO

      ENDDO

    END SUBROUTINE my_tr1

    SUBROUTINE my_tr2

      INTEGER :: i0, ia, ib

      DO kk = 1, kstack_n(dest_pe)
        k = kstack(kk,dest_pe)

        i0 = 0
        DO ia = 1, ncross
          gp_msg(idest)%buf(i0+ia,p0+kk) = x_gp(nz1+ia,k,nb1)
        ENDDO
        i0 = i0 + ncross
        IF (nb2>0) THEN
          DO ib = nb1+1, nb1+nb2-1
            DO ia = 1, na
              gp_msg(idest)%buf(i0+ia,p0+kk) = x_gp(ia,k,ib)
            ENDDO
            i0 = i0+na
          ENDDO
          ib = nb1+nb2
          DO ia = 1, nz2
            gp_msg(idest)%buf(i0+ia,p0+kk) = x_gp(ia,k,ib)
          ENDDO
        ENDIF

      ENDDO

    END SUBROUTINE my_tr2

  END SUBROUTINE pack_gp3

  !------------------------------------------------------------------------------

  SUBROUTINE pack_gp2(x_gp, with_exp)
    REAL(dp), INTENT(in)           :: x_gp(:,:) ! grid point field (nproma, nblock) 
    LOGICAL,  INTENT(in), OPTIONAL :: with_exp  ! apply exp-function in gp-space

    ! TODO: check if vector-exp is available

    INTEGER :: idest, dest_pe, p0, ir
    INTEGER :: n1, n2, nb1, nb2, na, nz1, nz2, ncross

    LOGICAL :: use_exp

    IF (require_init) CALL init_tr_gp_ffsl

    IF (ltimer) CALL timer_start(timer_pack_gp2)

    IF (kstack_max == 0) CALL prepare_tr

!!  IF (gp_pack_p0 > 100) &      !!++mgs 05/04/2011
!!       CALL finish('mo_tr_gp_ffsl:pack_gp3','unexpected huge buffer pos')
    IF (gp_pack_p0 >= tracer_max) &
         CALL finish('mo_tr_gp_ffsl:pack_gp3','dynamic sbuf not supported yet')

    IF (tr_dim == 0) THEN
      tr_dim = 2
    ELSE
      IF (tr_dim /= 2) &
           CALL finish('mo_tr_gp_ffsl:cannot mix 2d and 3d packing yet')
    ENDIF

    IF (PRESENT(with_exp)) THEN
      use_exp = with_exp
    ELSE
      use_exp = .FALSE.
    ENDIF

    na = ldc%nproma
    !CALL split_blocks(n1, nb1, nz1, ncross, n2, nb2, nz2)
    n1     = bsplit_n1    
    nb1    = bsplit_nb1   
    nz1    = bsplit_nz1   
    ncross = bsplit_ncross
    n2     = bsplit_n2    
    nb2    = bsplit_nb2   
    nz2    = bsplit_nz2   


    DO idest = 1, dest_n
      dest_pe = dest(idest)%xpe
      p0      = gp_pack_p0
      ir = dest(idest)%ir

      IF (ir == 1) THEN
        CALL my_tr1
      ELSEIF (ir == 2) THEN
        CALL my_tr2
      ELSEIF (ir == 3) THEN
        p0 = 2*p0
        CALL my_tr1
        p0 = p0+1
        CALL my_tr2
      ELSE
        CALL finish ('mo_tr_gp_ffsl:pack_gp3>','unexpected case 1',1)
      ENDIF
    ENDDO

    gp_pack_p0 = gp_pack_p0+1

    IF (ltimer) CALL timer_stop(timer_pack_gp2)

  CONTAINS

    SUBROUTINE my_tr1

      INTEGER :: i0, ia, ib

      i0 = 0
      IF (use_exp) THEN
        DO ib = 1, nb1-1
          DO ia = 1, na
            gp_msg(idest)%buf(i0+ia,p0+1) = EXP(x_gp(ia,ib))
          ENDDO
          i0 = i0+na
        ENDDO
        ib = nb1
        DO ia = 1, nz1
          gp_msg(idest)%buf(i0+ia,p0+1) = EXP(x_gp(ia,ib))
        ENDDO
      ELSE
        DO ib = 1, nb1-1
          DO ia = 1, na
            gp_msg(idest)%buf(i0+ia,p0+1) = x_gp(ia,ib)
          ENDDO
          i0 = i0+na
        ENDDO
        ib = nb1
        DO ia = 1,nz1
          gp_msg(idest)%buf(i0+ia,p0+1) = x_gp(ia,ib)
        ENDDO
      ENDIF

    END SUBROUTINE my_tr1

    SUBROUTINE my_tr2

      INTEGER :: i0, ia, ib

      i0 = 0
      IF (use_exp) THEN
        DO ia = 1, ncross
          gp_msg(idest)%buf(i0+ia,p0+1) = EXP(x_gp(nz1+ia,nb1))
        ENDDO
        i0 = i0 + ncross
        IF (nb2>0) THEN
          DO ib = nb1+1 ,nb1+nb2-1
            DO ia = 1, na
              gp_msg(idest)%buf(i0+ia,p0+1) = EXP(x_gp(ia,ib))
            ENDDO
            i0 = i0+na
          ENDDO
          ib = nb1+nb2
          DO ia = 1,nz2
            gp_msg(idest)%buf(i0+ia,p0+1) = EXP(x_gp(ia,ib))
          ENDDO
        ENDIF
      ELSE
        DO ia = 1, ncross
          gp_msg(idest)%buf(i0+ia,p0+1) = x_gp(nz1+ia,nb1)
        ENDDO
        i0 = i0 + ncross
        IF (nb2>0) THEN
          DO ib = nb1+1, nb1+nb2-1
            DO ia = 1, na
              gp_msg(idest)%buf(i0+ia,p0+1) = x_gp(ia,ib)
            ENDDO
            i0 = i0+na
          ENDDO
          ib = nb1+nb2
          DO ia = 1, nz2
            gp_msg(idest)%buf(i0+ia,p0+1)=x_gp(ia,ib)
          ENDDO
        ENDIF
      ENDIF

    END SUBROUTINE my_tr2

  END SUBROUTINE pack_gp2

  !------------------------------------------------------------------------------

  SUBROUTINE unpack_ffsl3(x_ffsl)
    REAL(dp), INTENT(out) :: x_ffsl (:,:,:) ! ffsl decomposition

    INTEGER :: p0, i0, j0, i, j, k, kk
    INTEGER :: nglon, ffsl_nlat, ffsl_nlatx, ffsl_rlat(2), jr, idx, src_pe, isrc

    IF (ltimer) CALL timer_start(timer_unpack_ffsl3)

    IF (kstack_max == 0) &
         CALL finish('mo_tr_gp_ffsl:unpack_ffsl3','bad case (1)')

    IF (ffsl_unpack_p0 > gp_pack_p0) &
         CALL finish('mo_tr_gp_ffsl:unpack_ffsl3','buffer_offset too big')


    !nglat2 = ldc%nglat/2
    ffsl_nlat  = ldc%ffsl%nlat
    ffsl_nlatx = ldc%ffsl%nlatx
    ffsl_rlat(1) = ffsl_nlat-ffsl_nlatx
    ffsl_rlat(2) = ffsl_nlatx


    !IF (ffsl_nlatx /= nglat2) &
    !     CALL finish('mo_tr_gp_ffsl:unpack_ffsl3','bad case nlatx')

    DO isrc = 1, src_n
      jr = src(isrc)%jr
      IF (jr == 0) &
           CALL finish('mo_tr_gp_ffsl:unpack_ffsl3','bad case')

      src_pe = src(isrc)%xpe
      p0     = ffsl_unpack_p0*kstack_n(my_pe)        
      idx    = src(isrc)%idx
      nglon  = gdc(idx)%nglon

      IF (jr == 1) THEN
        CALL my_unpack1(src(isrc)%ir)
      ELSEIF (jr == 2) THEN
        CALL my_unpack2(src(isrc)%ir)
      ELSEIF (jr == 3) THEN
        p0 = 2*p0;
        CALL my_unpack1(1)
        p0 = p0+kstack_n(my_pe)
        CALL my_unpack2(2)
      ELSE
        CALL finish('mo_tr_gp_ffsl:unpack_ffsl3','unexpected case')
      ENDIF

    ENDDO

    ffsl_unpack_p0 = ffsl_unpack_p0+1

    IF (ffsl_unpack_p0 == gp_pack_p0) THEN
      gp_pack_p0     = 0
      ffsl_unpack_p0 = 0
      tr_dim         = 0
    ENDIF

    IF (ltimer) CALL timer_stop(timer_unpack_ffsl3)

  CONTAINS

    SUBROUTINE my_unpack1(ir)
      INTEGER, INTENT(in) :: ir ! remote gp-region

      INTEGER :: nglat

      IF (ir > 2) &
           & CALL finish('mo_tr_gp_ffsl:unpack_ffsl3:my_unpack1','bad ir')

      nglat=gdc(idx)%nglh(ir)

      i0 = gdc(idx)%glons(ir)-1
      j0 = ffsl_nlat+1

      DO kk = 1, kstack_n(my_pe)
        k = kstack(kk,my_pe)
        DO j = 1, nglat
          DO i = 1, nglon
            x_ffsl(i0+i,j0-j,kk) = ffsl_msg(isrc)%buf(i+(j-1)*nglon,p0+kk)
          ENDDO
        ENDDO

      ENDDO

    END SUBROUTINE my_unpack1

    SUBROUTINE my_unpack2(ir)
      INTEGER, INTENT(in) :: ir

      INTEGER :: nglat

      IF (ir > 2) &
           & CALL finish('mo_tr_gp_ffsl:unpack_ffsl3:my_unpack2','bad ir')

      nglat=gdc(idx)%nglh(ir)

      i0 = gdc(idx)%glons(ir)-1
      j0 = nglat+1

      DO kk = 1, kstack_n(my_pe)
        k = kstack(kk,my_pe)
        DO j = 1, nglat
          DO i = 1, nglon
            x_ffsl(i0+i,j0-j,kk) = ffsl_msg(isrc)%buf(i+(j-1)*nglon,p0+kk)
          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE my_unpack2

  END SUBROUTINE unpack_ffsl3

  !------------------------------------------------------------------------------

  SUBROUTINE unpack_ffsl2(x_ffsl)
    REAL(dp), INTENT(out) :: x_ffsl (:,:) ! (nlon,ffsl_nlat)

    INTEGER :: p0, i0, j0, i, j, nglon, ffsl_nlat, ffsl_nlatx, jr, idx, src_pe, isrc

    IF (ltimer) CALL timer_start(timer_unpack_ffsl2)

    IF (kstack_max == 0) &
         CALL finish('mo_tr_gp_ffsl:unpack_ffsl2','bad case (1)')

    IF (ffsl_unpack_p0 > gp_pack_p0) &
         CALL finish('mo_tr_gp_ffsl:unpack_ffsl2','buffer_offset too big')

    IF (tr_dim /= 2) &
         CALL finish('mo_tr_gp_ffsl:cannot mix 2d and 3d packing yet')

    ffsl_nlat  = ldc%ffsl%nlat
    ffsl_nlatx = ldc%ffsl%nlatx

    DO isrc = 1, src_n
      jr = src(isrc)%jr
      IF (jr == 0) &
           CALL finish('mo_tr_gp_ffsl:unpack_ffsl2','bad case')

      src_pe = src(isrc)%xpe
      p0     = ffsl_unpack_p0        
      idx    = src(isrc)%idx
      nglon  = gdc(idx)%nglon

      IF (jr == 1) THEN
        CALL my_unpack1(src(isrc)%ir)
      ELSEIF (jr == 2) THEN
        CALL my_unpack2(src(isrc)%ir)
      ELSEIF (jr == 3) THEN
        p0 = 2*p0;
        CALL my_unpack1(1)
        p0 = p0+1
        CALL my_unpack2(2)
      ELSE
        CALL finish('mo_tr_gp_ffsl:unpack_ffsl2','unexpected case')
      ENDIF

    ENDDO

    ffsl_unpack_p0 = ffsl_unpack_p0+1

    IF (ffsl_unpack_p0 == gp_pack_p0) THEN
      gp_pack_p0     = 0
      ffsl_unpack_p0 = 0
      tr_dim         = 0
    ENDIF

    IF (ltimer) CALL timer_stop(timer_unpack_ffsl2)

  CONTAINS

    SUBROUTINE my_unpack1(ir)
      INTEGER, intent(in) :: ir

      INTEGER :: nglat

      IF (ir > 2) &
           & CALL finish('mo_tr_gp_ffsl:unpack_ffsl2:my_unpack1','bad ir')

      nglat=gdc(idx)%nglh(ir)

      i0 = gdc(idx)%glons(ir)-1
      j0 = ffsl_nlat+1
      DO j = 1, nglat
        DO i = 1, nglon
          x_ffsl(i0+i,j0-j) = ffsl_msg(isrc)%buf(i+(j-1)*nglon,p0+1)
        ENDDO
      ENDDO

    END SUBROUTINE my_unpack1

    SUBROUTINE my_unpack2(ir)
      INTEGER, intent(in) :: ir

      INTEGER :: nglat

      IF (ir > 2) &
           & CALL finish('mo_tr_gp_ffsl:unpack_ffsl2:my_unpack2','bad ir')

      nglat=gdc(idx)%nglh(ir)

      i0 = gdc(idx)%glons(ir)-1
      j0 = nglat+1

      DO j = 1,nglat
        DO i = 1,nglon
          x_ffsl(i0+i,j0-j) = ffsl_msg(isrc)%buf(i+(j-1)*nglon,p0+1)
        ENDDO
      ENDDO

    END SUBROUTINE my_unpack2

  END SUBROUTINE unpack_ffsl2

  !------------------------------------------------------------------------------

  SUBROUTINE pack_ffsl3(x_ffsl)
    REAL(dp), INTENT(in) :: x_ffsl (:,:,:) ! ffsl decomposition

    INTEGER :: p0, i0, i, j, k, kk, ffsl_nlat, ffsl_nlatx, jr, idx, src_pe, isrc
    INTEGER :: nglon

    ! reverse meaning of src, dest

    IF (require_init) CALL init_tr_gp_ffsl

    IF (kstack_max==0) CALL prepare_tr

    IF (ltimer) CALL timer_start(timer_pack_ffsl3)

!!  IF (ffsl_pack_p0 > 100) &    !!++mgs 05/04/2011
!!       CALL finish('mo_tr_gp_ffsl:pack_ffsl3','unexpected huge buffer pos')
!!  write(*,*) 'pack_ffsl3: ffsl_pack_p0, tracer_max = ',ffsl_pack_p0,tracer_max
    
    IF (ffsl_pack_p0 >= tracer_max) &
         CALL finish('mo_tr_gp_ffsl:pack_ffsl3','dynamic sbuf not supported yet')

    IF (tr_dim == 0) THEN
      tr_dim = 3
    ELSE
      IF (tr_dim /= 3) &
           CALL finish('mo_tr_gp_ffsl:pack_ffsl3:cannot mix 2d and 3d packing yet')
    ENDIF

    ffsl_nlat  = ldc%ffsl%nlat
    ffsl_nlatx = ldc%ffsl%nlatx

    DO isrc = 1, src_n
      jr = src(isrc)%jr
      IF (jr == 0) &
           CALL finish('mo_tr_gp_ffsl:pack_ffsl3','bad case')

      src_pe = src(isrc)%xpe
      p0     = ffsl_pack_p0*kstack_n(my_pe)        
      idx    = src(isrc)%idx
      nglon  = gdc(idx)%nglon

      IF (jr == 1) THEN
        CALL my_pack1(src(isrc)%ir)
      ELSEIF (jr == 2) THEN
        CALL my_pack2(src(isrc)%ir)
      ELSEIF (jr == 3) THEN
        p0 = 2*p0
        CALL my_pack1(1)
        p0 = p0+kstack_n(my_pe)
        CALL my_pack2(2)
      ELSE
        CALL finish('mo_tr_gp_ffsl:pack_ffsl3','unexpected case')
      ENDIF

    ENDDO

    ffsl_pack_p0 = ffsl_pack_p0+1

    IF (ltimer) CALL timer_stop(timer_pack_ffsl3)

  CONTAINS

    SUBROUTINE my_pack1(ir)
      INTEGER, INTENT(in) :: ir

      INTEGER :: nglat

      i0 = gdc(idx)%glons(ir)-1
      nglat=gdc(idx)%nglh(ir)

      DO kk = 1, kstack_n(my_pe)
        k = kstack(kk,my_pe)
        DO j = 1, nglat
          DO i = 1, nglon
            ffsl_msg(isrc)%buf(i+(j-1)*nglon,p0+kk) = x_ffsl(i0+i,ffsl_nlat+1-j,kk)
          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE my_pack1

    SUBROUTINE my_pack2(ir)
      INTEGER, INTENT(in) :: ir

      INTEGER :: nglat

      i0 = gdc(idx)%glons(ir)-1
      nglat=gdc(idx)%nglh(ir)

      DO kk = 1, kstack_n(my_pe)
        k = kstack(kk,my_pe)
        DO j = 1, nglat
          DO i = 1, nglon
            ffsl_msg(isrc)%buf(i+(j-1)*nglon,p0+kk) = x_ffsl(i0+i,nglat+1-j,kk)
          ENDDO
        ENDDO
      ENDDO
    END SUBROUTINE my_pack2

  END SUBROUTINE pack_ffsl3

  !------------------------------------------------------------------------------
  SUBROUTINE pack_ffsl4(x_ffsl)
    REAL(dp), INTENT(in) :: x_ffsl (:,:,:,:) ! ffsl decomposition

    INTEGER :: nt,itracer

    nt = SIZE(x_ffsl,4)
    DO itracer = 1, nt
      CALL pack_ffsl3(x_ffsl(:,:,:,itracer))
    ENDDO

  END SUBROUTINE pack_ffsl4

  !------------------------------------------------------------------------------

  SUBROUTINE unpack_gp3(x_gp)
    REAL(dp),INTENT(out) :: x_gp(:,:,:) ! grid point field (nproma, nlev, nblock) 

    INTEGER :: idest, dest_pe, p0, ir
    INTEGER :: kk, k
    INTEGER :: na, ncross, n1, n2, nb1, nz1, nb2, nz2

    ! reverse meaning of src, dest

    IF (ltimer) CALL timer_start(timer_unpack_gp3)

    IF (kstack_max == 0) CALL prepare_tr

    IF (gp_unpack_p0 > ffsl_pack_p0) &
         CALL finish('mo_tr_gp_ffsl:unpack_gp3','buffer_offset too big')

    na      = ldc%nproma
    !CALL split_blocks(n1, nb1, nz1, ncross, n2, nb2, nz2)
    n1     = bsplit_n1    
    nb1    = bsplit_nb1   
    nz1    = bsplit_nz1   
    ncross = bsplit_ncross
    n2     = bsplit_n2    
    nb2    = bsplit_nb2   
    nz2    = bsplit_nz2   

    DO idest = 1, dest_n
      dest_pe = dest(idest)%xpe
      p0      = gp_unpack_p0*kstack_n(dest_pe)
      ir = dest(idest)%ir

      IF (ir == 1) THEN
        CALL my_tr1
      ELSEIF (ir == 2) THEN
        CALL my_tr2
      ELSEIF (ir == 3) THEN
        p0=2*p0
        CALL my_tr1
        p0 = p0+kstack_n(dest_pe)
        CALL my_tr2
      ELSE
        CALL finish ('mo_tr_gp_ffsl:unpack_gp3','unexpected case 1',1)
      ENDIF

    ENDDO


    gp_unpack_p0 = gp_unpack_p0+1
    IF (gp_unpack_p0 == ffsl_pack_p0) THEN
      gp_unpack_p0 = 0
      ffsl_pack_p0 = 0
      tr_dim       = 0
    ENDIF

    IF (ltimer) CALL timer_stop(timer_unpack_gp3)

  CONTAINS

    SUBROUTINE my_tr1
      INTEGER :: i0, ia, ib
      DO kk = 1, kstack_n(dest_pe)
        k = kstack(kk,dest_pe)
        i0 = 0
        DO ib = 1, nb1-1
          IF (lcheck) CALL check(__LINE__, i0+na, 1, SIZE(gp_msg(idest)%buf,1))
          DO ia = 1, na
            x_gp(ia,k,ib) = gp_msg(idest)%buf(i0+ia,p0+kk)
          ENDDO
          i0 = i0+na
        ENDDO
        ib = nb1
        IF (lcheck) CALL check(__LINE__, i0+nz1, 1, SIZE(gp_msg(idest)%buf,1))
        DO ia = 1, nz1
          x_gp(ia,k,ib) = gp_msg(idest)%buf(i0+ia,p0+kk)
        ENDDO
      ENDDO
    END SUBROUTINE my_tr1

    SUBROUTINE my_tr2
      INTEGER :: i0, ia, ib
      DO kk = 1,kstack_n(dest_pe)
        k = kstack(kk,dest_pe)
        i0 = 0
        IF (lcheck) CALL check(__LINE__, i0+ncross, 0, SIZE(gp_msg(idest)%buf,1))
        DO ia = 1, ncross
          x_gp(nz1+ia,k,nb1) = gp_msg(idest)%buf(i0+ia,p0+kk)
        ENDDO
        i0 = i0 + ncross 
        IF (nb2>0) THEN
          DO ib = nb1+1,nb1+nb2-1
            IF (lassert) CALL check(__LINE__, i0+na, 1, SIZE(gp_msg(idest)%buf,1))
            DO ia = 1,na
              x_gp(ia,k,ib) = gp_msg(idest)%buf(i0+ia,p0+kk)
            ENDDO
            i0 = i0+na
          ENDDO
          ib = nb1+nb2
          IF (lassert) CALL check(__LINE__, i0+nz2, 1, SIZE(gp_msg(idest)%buf,1))
          DO ia = 1,nz2
            x_gp(ia,k,ib) = gp_msg(idest)%buf(i0+ia,p0+kk)
          ENDDO
          DO ia = nz2+1,na
            x_gp(ia,k,ib) = 0.0_dp
          ENDDO
        ELSE
          DO ia = nz1+ncross+1,na
            x_gp(ia,k,nb1) = 0.0_dp
          ENDDO
        ENDIF
      ENDDO
    END SUBROUTINE my_tr2

  END SUBROUTINE unpack_gp3

  !------------------------------------------------------------------------------

  SUBROUTINE unpack_gp4(x_gp)
    REAL(dp), INTENT(out) :: x_gp(:,:,:,:) ! grid point field (nproma, nlev, ntracer, nblock) 

    INTEGER :: idest, dest_pe, p0, ir
    INTEGER :: kk, k
    INTEGER :: na,ncross,n1,n2,nb1,nb2,nz1,nz2,itracer,ntracer

    ! reverse meaning of src, dest

    IF (ltimer) CALL timer_start(timer_unpack_gp4)

    IF (kstack_max == 0) CALL prepare_tr

    IF (gp_unpack_p0 > ffsl_pack_p0) &
         CALL finish('unpack_gp4','buffer_offset too big')

    ntracer = SIZE(x_gp,3)
    na = ldc%nproma
    !CALL split_blocks(n1, nb1, nz1, ncross, n2, nb2, nz2)
    n1     = bsplit_n1    
    nb1    = bsplit_nb1   
    nz1    = bsplit_nz1   
    ncross = bsplit_ncross
    n2     = bsplit_n2    
    nb2    = bsplit_nb2   
    nz2    = bsplit_nz2   

    DO itracer = 1, ntracer

      DO idest = 1, dest_n
        dest_pe = dest(idest)%xpe
        p0      = gp_unpack_p0*kstack_n(dest_pe)
        ir = dest(idest)%ir

        IF (ir  == 1) THEN
          CALL my_tr1
        ELSEIF (ir == 2) THEN
          CALL my_tr2
        ELSEIF (ir == 3) THEN
          p0 = 2*p0
          CALL my_tr1
          p0 = p0+kstack_n(dest_pe)
          CALL my_tr2
        ELSE
          CALL finish ('unpack_gp4','unexpected case 1',1)
        ENDIF
      ENDDO

      gp_unpack_p0 = gp_unpack_p0+1

    ENDDO

    IF (gp_unpack_p0  ==  ffsl_pack_p0) THEN
      gp_unpack_p0 = 0
      ffsl_pack_p0 = 0
      tr_dim       = 0
    ENDIF

    IF (ltimer) CALL timer_stop(timer_unpack_gp4)

  CONTAINS

    SUBROUTINE my_tr1
      INTEGER :: i0, ia, ib
      DO kk = 1, kstack_n(dest_pe)
        k = kstack(kk,dest_pe)
        i0 = 0
        DO ib = 1,nb1-1
          IF (lcheck) CALL check(__LINE__, i0+na, 1, SIZE(gp_msg(idest)%buf,1))
          DO ia = 1,na
            x_gp(ia,k,itracer,ib) = gp_msg(idest)%buf(i0+ia,p0+kk)
          ENDDO
          i0 = i0+na
        ENDDO
        ib = nb1
        IF (lcheck) CALL check(__LINE__, i0+nz1, 1, SIZE(gp_msg(idest)%buf,1))
        DO ia = 1,nz1
          x_gp(ia,k,itracer,ib) = gp_msg(idest)%buf(i0+ia,p0+kk)
        ENDDO
      ENDDO
    END SUBROUTINE my_tr1


    SUBROUTINE my_tr2
      INTEGER :: i0, ia, ib
      DO kk = 1, kstack_n(dest_pe)
        k = kstack(kk,dest_pe)
        i0 = 0
        IF (lcheck) CALL check(__LINE__, i0+ncross, 0, SIZE(gp_msg(idest)%buf,1))
        DO ia = 1, ncross
          x_gp(nz1+ia,k,itracer,nb1) = gp_msg(idest)%buf(i0+ia,p0+kk)
        ENDDO
        i0 = i0 + ncross
        IF (nb2>0) THEN
          DO ib = nb1+1,nb1+nb2-1
            IF (lcheck) CALL check(__LINE__, i0+na, 1, SIZE(gp_msg(idest)%buf,1))
            DO ia = 1,na
              x_gp(ia,k,itracer,ib) = gp_msg(idest)%buf(i0+ia,p0+kk)
            ENDDO
            i0 = i0+na
          ENDDO
          ib = nb1+nb2
          IF (lcheck) CALL check(__LINE__, i0+nz2, 1, SIZE(gp_msg(idest)%buf,1))
          DO ia = 1,nz2
            x_gp(ia,k,itracer,ib) = gp_msg(idest)%buf(i0+ia,p0+kk)
          ENDDO
          DO ia = nz2+1,na
            x_gp(ia,k,itracer,ib) = 0.0_dp
          ENDDO
        ELSE
          DO ia = nz1+ncross+1,na
            x_gp(ia,k,itracer,nb1) = 0.0_dp
          ENDDO
        ENDIF
      ENDDO
    END SUBROUTINE my_tr2

  END SUBROUTINE unpack_gp4

  !------------------------------------------------------------------------------

  SUBROUTINE unpack_gp4_slow(x_gp)
    REAL(dp), INTENT(out) :: x_gp(:,:,:,:) ! grid point field

    INTEGER :: nt, itracer

    nt = SIZE(x_gp,3)
    DO itracer = 1, nt
      CALL unpack_gp3(x_gp(:,:,itracer,:)) ! slow
    ENDDO

  END SUBROUTINE unpack_gp4_slow

  !------------------------------------------------------------------------------

  SUBROUTINE gp2ffsl
    INTEGER :: nlat2, nt, nglon

    IF (gp_pack_p0  ==  0) RETURN ! nothing to do

    IF (ffsl_unpack_p0 /=  0) &
         CALL finish('tr_gp_ffsl_exchange: rbuf not empty')

    IF (ltimer) CALL timer_start(timer_gp2ffsl)

    nlat2  = ldc%nglat/2
    nglon  = ldc%nglon
    nt     = gp_pack_p0

#ifndef NOMPI
    IF (use_mpi_collectives) THEN
      CALL my_alltoallv
    ELSE
      CALL my_isend_recv
    ENDIF
#else
    CALL my_serial
#endif

    IF (ltimer) CALL timer_stop(timer_gp2ffsl)

  CONTAINS

    SUBROUTINE my_serial
      INTEGER :: idest, dest_pe, isrc, src_pe, nk, nr,  ilon, ilat, ikt, i0

      idest   = 1
      dest_pe = dest(idest)%xpe
      nr      = dest(idest)%nr

      isrc    = 1
      src_pe  = src(isrc)%xpe
      IF (tr_dim == 3) THEN
        nk = kstack_n(my_pe)
      ELSE
        nk = 1
      ENDIF
      
      DO ikt = 1, nr*nk*nt
        DO ilat = 1, nlat2
          i0 = (ilat-1)*nglon
          DO ilon = 1, nglon
            ffsl_msg(isrc)%buf(ilon+(ilat-1)*nglon,ikt) = gp_msg(idest)%buf(i0+ilon,ikt)
          ENDDO
        END DO
      ENDDO
      
    END SUBROUTINE my_serial

#ifndef NOMPI
    SUBROUTINE my_alltoallv
      INTEGER :: idest, dest_pe, isrc, src_pe, nk, nr

      INTEGER :: sc(0:ldc%d_nprocs-1), sd(0:ldc%d_nprocs-1), st
      INTEGER :: rc(0:ldc%d_nprocs-1), rd(0:ldc%d_nprocs-1), rt
      INTEGER :: ierror
      INTEGER :: idx

      sc(:)  = 0
      sd(:)  = 0
      st     = p_real_dp

      DO idest = 1, dest_n
        dest_pe = dest(idest)%xpe
        IF (tr_dim == 3) THEN
          nk = kstack_n(dest_pe)
        ELSE
          nk = 1
        ENDIF
        nr = dest(idest)%nr
        sc(dest_pe) = nlat2*nglon*nk*nt*nr
        sd(dest_pe) = gp_msg(idest)%offset
      ENDDO

      rc(:) = 0
      rd(:) = 0
      rt = p_real_dp
      IF (tr_dim == 3) THEN
        nk = kstack_n(my_pe)
      ELSE
        nk = 1
      ENDIF

      DO isrc = 1, src_n
        src_pe = src(isrc)%xpe
        nr     = src(isrc)%nr
        idx    = src(isrc)%idx
        rc(src_pe) = (gdc(idx)%nglat/2)*gdc(idx)%nglon*nk*nt*nr
        rd(src_pe) = ffsl_msg(isrc)%offset
      ENDDO


      CALL MPI_ALLTOALLV(gp_msg(1)%buf(1,1),     sc, sd, st, &
                         ffsl_msg(1)%buf(1,1),   rc, rd, rt, &
                         dcomm, ierror)

    END SUBROUTINE my_alltoallv
#endif

    SUBROUTINE my_isend_recv
      INTEGER :: idest, dest_pe, isrc, src_pe, nk, nr
      INTEGER :: idx, ir, hsize, nglh(3)

      nglh(1:2) = ldc%nglh
      nglh(3)   = ldc%nglat

      DO idest = 1, dest_n
        dest_pe = dest(idest)%xpe
        IF (tr_dim == 3) THEN
          nk = kstack_n(dest_pe)
        ELSE
          nk = 1
        ENDIF
        ir = dest(idest)%ir
        hsize = nglon * nglh(ir)
        CALL p_isend(gp_msg(idest)%buf(1,1), dest_pe, gp2ffsl_tag, hsize*nk*nt, dcomm)
      ENDDO

      IF (tr_dim == 3) THEN
        nk = kstack_n(my_pe)
      ELSE
        nk = 1
      ENDIF
      DO isrc = 1, src_n
        src_pe = src(isrc)%xpe
        nr     = src(isrc)%nr
        ir     = src(isrc)%ir
        idx    = src(isrc)%idx
        nglh(1:2) = gdc(idx)%nglh
        nglh(3)   = gdc(idx)%nglat
        hsize  = gdc(idx)%nglon * nglh(ir)
        ffsl_msg(isrc)%buf(:,:)=-4.0_dp
        CALL  p_recv(ffsl_msg(isrc)%buf(1,1), src_pe, gp2ffsl_tag, hsize*nk*nt, dcomm)
      ENDDO

      CALL p_wait
    END SUBROUTINE my_isend_recv

  END SUBROUTINE gp2ffsl

  !------------------------------------------------------------------------------

  SUBROUTINE ffsl2gp
    INTEGER :: nlat2, nglon, nt

    ! reverse meaning of src, dest

    IF (ffsl_pack_p0  ==  0) RETURN ! nothing to do

    IF (gp_unpack_p0 /=  0) &
         CALL finish('tr_gp_ffsl_exchange: rbuf not empty')

    IF (tr_dim /= 3) &
         CALL finish('ffsl2gp:unexpected tr_dim state')

    IF (ltimer) CALL timer_start(timer_ffsl2gp)

    nlat2  =  ldc%nglat/2
    nglon  =  ldc%nglon
    nt     = ffsl_pack_p0


#ifndef NOMPI
    IF (use_mpi_collectives) THEN
      CALL my_alltoallv
    ELSE
      CALL my_isend_recv
    ENDIF
#else
    CALL my_serial
#endif



    IF (ltimer) CALL timer_stop(timer_ffsl2gp)

  CONTAINS

    SUBROUTINE my_serial
      INTEGER :: idest, dest_pe, isrc, src_pe,  nk, nr
      INTEGER :: ilon, ilat, ikt, i0

      idest   = 1
      dest_pe = dest(idest)%xpe
      nr      = dest(idest)%nr
      isrc    = 1
      src_pe  = src(isrc)%xpe
      nk      = kstack_n(my_pe)

      DO ikt = 1, nr*nk*nt
        DO ilat = 1, nlat2
          i0 = (ilat-1)*ldc%nglon
          DO ilon = 1, ldc%nglon
            gp_msg(isrc)%buf(i0+ilon,ikt) = ffsl_msg(idest)%buf(ilon+(ilat-1)*ldc%nglon,ikt)
          ENDDO
        END DO
      ENDDO
      
    END SUBROUTINE my_serial

#ifndef NOMPI
    SUBROUTINE my_alltoallv
      INTEGER :: idest, dest_pe, isrc, src_pe,  nk

      INTEGER :: sc(0:ldc%d_nprocs-1), sd(0:ldc%d_nprocs-1), st
      INTEGER :: rc(0:ldc%d_nprocs-1), rd(0:ldc%d_nprocs-1), rt
      INTEGER :: ierror
      INTEGER :: ir, idx, hsize, nglh(3)

      sc(:)  = 0
      sd(:)  = 0
      st     = p_real_dp

      nglh(1:2) = ldc%nglh
      nglh(3)   = ldc%nglat

      DO idest = 1, dest_n
        dest_pe     = dest(idest)%xpe
        nk          = kstack_n(dest_pe)
        ir = dest(idest)%ir
        hsize = nglon * nglh(ir)
        sc(dest_pe) = hsize*nk*nt
        sd(dest_pe) = gp_msg(idest)%offset
      ENDDO

      rc(:)  = 0
      rd(:)  = 0
      rt     = p_real_dp

      DO isrc = 1, src_n
        src_pe     = src(isrc)%xpe
        nk         = kstack_n(my_pe)
        ir     = src(isrc)%ir
        idx    = src(isrc)%idx
        nglh(1:2) = gdc(idx)%nglh
        nglh(3)   = gdc(idx)%nglat
        hsize  = gdc(idx)%nglon * nglh(ir)
        rc(src_pe) = hsize*nk*nt
        rd(src_pe) = ffsl_msg(isrc)%offset
      ENDDO


      CALL MPI_ALLTOALLV(ffsl_msg(1)%buf(1,1), rc, rd, rt, &
           &             gp_msg(1)%buf(1,1),   sc, sd, st, &
           &             dcomm, ierror)

    END SUBROUTINE my_alltoallv
#endif

    SUBROUTINE my_isend_recv
      INTEGER :: idest, dest_pe, isrc, src_pe,  nk
      INTEGER :: idx, ir, hsize, nglh(3)
      ! todo: use new communicator
      nglh(1:2) = ldc%nglh
      nglh(3)   = ldc%nglat

      DO idest = 1, dest_n
        dest_pe = dest(idest)%xpe
        nk      = kstack_n(dest_pe)
        ir = dest(idest)%ir
        hsize = nglon * nglh(ir)
        CALL p_irecv(gp_msg(idest)%buf(1,1), dest_pe, ffsl2gp_tag, hsize*nk*nt, dcomm)
      ENDDO

      DO isrc = 1, src_n
        src_pe = src(isrc)%xpe
        nk     = kstack_n(my_pe)
        ir     = src(isrc)%ir
        idx    = src(isrc)%idx
        nglh(1:2) = gdc(idx)%nglh
        nglh(3)   = gdc(idx)%nglat
        hsize  = gdc(idx)%nglon * nglh(ir)
        CALL p_send(ffsl_msg(isrc)%buf(1,1), src_pe, ffsl2gp_tag, hsize*nk*nt, dcomm)
      ENDDO

      CALL p_wait

    END SUBROUTINE my_isend_recv

  END SUBROUTINE ffsl2gp


  !------------------------------------------------------------------------------

  SUBROUTINE check(line, i, imin, imax)
    INTEGER, INTENT(in) :: line, i, imin, imax

    CHARACTER(len=*), PARAMETER :: context = 'mo_tr_gp_ffsl:check'    
    INTEGER, SAVE :: state = 0
    CHARACTER(len=16) :: line_str, state_str, i_str, imin_str, imax_str
    CHARACTER(len=80) :: m

    state = state + 1
    
    IF (i<imin .OR. i>imax) THEN

      WRITE(line_str,'(I8)') line
      WRITE(state_str,'(I8)') state
      WRITE(i_str,'(I8)') i
      WRITE(imin_str,'(I8)') imin
      WRITE(imax_str,'(I8)') imax
      
      m =   'line='//TRIM(ADJUSTL(line_str))//&
           &', state='//TRIM(ADJUSTL(state_str))//&
           &', i='//TRIM(ADJUSTL(i_str))//&
           &', imin='//TRIM(ADJUSTL(imin_str))//&
           &', imax='//TRIM(ADJUSTL(imax_str))
           
      CALL finish(context, TRIM(ADJUSTL(m)))

    ENDIF
    
  END SUBROUTINE check

  !------------------------------------------------------------------------------

  SUBROUTINE assert(line, condition)
    INTEGER, INTENT(in) :: line
    LOGICAL, INTENT(in) :: condition

    CHARACTER(len=*), PARAMETER :: context = 'mo_tr_gp_ffsl:assert'    
    INTEGER, SAVE :: state = 0
    CHARACTER(len=16) :: line_str, state_str
    CHARACTER(len=80) :: m

    state = state + 1
    
    IF (.NOT. condition) THEN

      WRITE(line_str,'(I8)') line
      WRITE(state_str,'(I8)') state

      m =   'line '//TRIM(ADJUSTL(line_str))//&
           &', state='//TRIM(ADJUSTL(state_str))//&
           &', assert failed'

      CALL finish(context, TRIM(ADJUSTL(m)))

    ENDIF
    
  END SUBROUTINE assert

  !------------------------------------------------------------------------------

  SUBROUTINE prepare_tr
    CHARACTER(len=*), PARAMETER :: context = 'mo_tr_gp_ffsl:prepare_tr'

    INTEGER :: nglon, nglat, nglev, n, na, nb
    INTEGER :: ia, k, p
    INTEGER :: pe, idx, dest_set_b, src_pe, dest_pe
    INTEGER :: gp_lat1, gp_lat2, gp_latn, i, j
    INTEGER :: recv_ir(0:d_nprocs-1), recv_nr(0:d_nprocs-1)
    INTEGER :: recv_jr(0:d_nprocs-1)
    INTEGER :: send_ir(0:d_nprocs-1), send_nr(0:d_nprocs-1)
    INTEGER :: send_jr(0:d_nprocs-1)
    INTEGER :: send_hsize(0:d_nprocs-1), recv_hsize(0:d_nprocs-1)
    INTEGER :: k2dest_idx(2,ldc%nlev)
    INTEGER :: kflag(ldc%nlev,0:d_nprocs-1)
    INTEGER :: nglh(2), hsize(2), gp2ffsl_rmap(2), ffsl2gp_rmap(2)

    IF (require_init) CALL finish(context,'expected initialization done')

    ! prepare block splitting:
    CALL split_blocks(bsplit_n1, bsplit_nb1, bsplit_nz1, bsplit_ncross, bsplit_n2, bsplit_nb2, bsplit_nz2)
        

    ! In the following we have the transposition from gp to ffsl in mind.
    ! For the inverse transposition (ffsl->gp) we will later just exchange
    ! the meaning of src and dest, so we don't need to handle that direction
    ! here at all.

    nglon           = ldc%nglon
    nglat           = ldc%nglat
    nglh            = ldc%nglh
    nglev           = ldc%nlev
    na              = ldc%nproma
    nb              = ldc%ngpblks

    ! vertical analysis:

    ! from here on we put ourself in the role of a sender:

    k2dest_idx(:,:) = 0
    kstack(:,:)     = 0
    kstack_n(:)     = 0
    kstack_max      = 0

    DO k = 1, ldc%nlev
      !
      ! for PEs with odd set_a:
      !   Northern Hemisphere sent to same set_a
      !   Southern Hemisphere sent to set_a+1 (unless set_a  ==  nproca)
      ! for PEs with even set_a:
      !   Northern Hemisphere sent to set_a-1
      !   Southern Hemisphere sent to same set_a
      !
      dest_set_b  =  MOD(k-1,nprocb)+1 ! set_b of receiving PE for this k
      IF (lcheck) CALL check(__LINE__, dest_set_b, 1, nprocb)
      IF( MOD(ldc%set_a,2)  ==  1 ) THEN
        IF (lcheck) CALL check(__LINE__, ldc%set_a, 1, nproca)
        k2dest_idx(1,k) = ldc%mapmesh(dest_set_b,ldc%set_a) ! N-region goes to my nprocb proc-family
        ia = MIN(ldc%set_a+1,ldc%nproca)
        IF (lcheck) CALL check(__LINE__, ia, 1, nproca)
        k2dest_idx(2,k) = ldc%mapmesh(dest_set_b,ia) ! S-region goes preferably to northern nprocb proc-family
      ELSE
        IF (lcheck) CALL check(__LINE__, ldc%set_a-1, 1, nproca)  
        k2dest_idx(1,k) = ldc%mapmesh(dest_set_b,ldc%set_a-1) ! N-region goes to southern nprocb proc-family
        IF (lcheck) CALL check(__LINE__, ldc%set_a, 1, nproca)
        k2dest_idx(2,k) = ldc%mapmesh(dest_set_b,ldc%set_a) ! S-region goes to my nprocb proc-family
      ENDIF
    ENDDO

    kflag(:,:) = 0
    DO k = 1, ldc%nlev
      DO p = 1,2
        idx              = k2dest_idx(p,k)
        pe               = gdc(idx)%pe - pe_shift
        IF (lcheck) CALL check(__LINE__, pe, 0, d_nprocs-1)
        kflag(k,pe)     = kflag(k,pe)+1
      ENDDO
    ENDDO

    kstack_n(:) = 0       ! how many k-levels do we send to a pe
    kstack(:,:) = 0       ! which levels ..
    DO pe = 0, d_nprocs-1 ! loop over dest pes
      n = 0
      DO k = 1, ldc%nlev
        IF (kflag(k,pe) > 0 ) THEN
          n = n+1
          IF (lcheck) CALL check(__LINE__, n, 1, ldc%nlev)
          kstack(n,pe) = k
        ENDIF
      ENDDO
      kstack_n(pe) = n
      IF (n > kstack_max) kstack_max = n
    ENDDO

    ! check if my ffsl region overlaps with my own gp-data
    IF (lassert) CALL assert(__LINE__, kstack_n(my_pe) > 0)

    ! horizontal analysis:
    
    ! from here on we put ourself in the role of a receiver in the gp->ffsl transfer:

    ! local ffsl-region:
    gp_lat1 = ldc%nlat+1-ldc%ffsl%latn  !translate ffsl region boundaries to gp latitudes
    gp_lat2 = ldc%nlat+1-ldc%ffsl%lats  ! 
    IF (lassert) CALL assert(__LINE__, gp_lat1 <= gp_lat2 )
    gp_latn = gp_lat2 - gp_lat1 +1
    src_n = 0

    DO src_pe = 0, d_nprocs-1
      idx=pe2idx(src_pe)

      recv_ir(src_pe)     = 0 ! which regions do we get? 1=first, 2=second, 3=both
      recv_jr(src_pe)     = 0 ! where does it go? to 1=first, 2=second, 3=both ffsl-substribes
      recv_hsize(src_pe)  = 0 ! horizontal size of one single region (must be valid for both regions) 
      recv_nr(src_pe)     = 0 ! number of regions we get from src_pe

      hsize(:)            = 0 ! size of first and second gp region
      gp2ffsl_rmap(:)     = 0 ! region map: gp_region to ffsl_substripe
      ffsl2gp_rmap(:)     = 0 ! inverse region map

      ! i-loop over remote gp-regions:
      DO i = 1, 2
        IF (gp_lat1 <= gdc(idx)%glats(i) .AND. gdc(idx)%glate(i) <= gp_lat2) THEN
          recv_nr(src_pe)         = recv_nr(src_pe) + 1
          hsize(i)                = gdc(idx)%nglon*gdc(idx)%nglh(i) ! remote region size
          IF ( gp_lat1 == gdc(idx)%glats(i) ) THEN
            gp2ffsl_rmap(i)       = 1
            ffsl2gp_rmap(1)       = i
          ELSEIF ( gp_lat2 == gdc(idx)%glate(i)) THEN
            gp2ffsl_rmap(i)       = 2
            ffsl2gp_rmap(2)       = i
          ELSE
            CALL finish('mo_tr_gp_ffsl:prepare_tr','unexpected region map')                        
          ENDIF
        ENDIF
      ENDDO

      
      IF ( recv_nr(src_pe) > 0 ) THEN
        IF (lassert) CALL assert(__LINE__, hsize(1) == 0 .OR. hsize(2) == 0 .OR. hsize(1) == hsize(2))
        recv_hsize(src_pe)   = MAXVAL(hsize)
        recv_ir(src_pe)      = MERGE(1, 0, gp2ffsl_rmap(1) > 0) + MERGE(2, 0, gp2ffsl_rmap(2) > 0)
        recv_jr(src_pe)      = MERGE(1, 0, ffsl2gp_rmap(1) > 0) + MERGE(2, 0, ffsl2gp_rmap(2) > 0)
        IF (recv_nr(src_pe) == 2) CALL assert(__LINE__, gp2ffsl_rmap(1) == 1 .AND. gp2ffsl_rmap(2) == 2)
        src_n = src_n + 1
      ENDIF

    ENDDO

    ALLOCATE(src(src_n))

    j = 0
    DO src_pe = 0, d_nprocs-1
      IF (recv_nr(src_pe) > 0) THEN
        j = j+1
        IF (lcheck) CALL check(__LINE__, j, 1, src_n)
        src(j)%xpe   = src_pe
        src(j)%idx   = pe2idx(src_pe)
        src(j)%ir    = recv_ir(src_pe) 
        src(j)%nr    = recv_nr(src_pe)
        src(j)%jr    = recv_jr(src_pe)
        src(j)%hsize = recv_hsize(src_pe)
      ENDIF
    ENDDO

    IF (ldebug) THEN
      ! check if we receive all the volume data we need (k dimsion is meaningless)
      n=0
      DO j=1, src_n
        n = n + src(j)%hsize*src(j)%nr * kstack_n(my_pe)
      ENDDO
      CALL assert(__LINE__, n == ldc%nlon * ldc%ffsl%nlat * kstack_n(my_pe) )
    ENDIF

    ! from here on we put ourself again in the role of a sender in the gp->ffsl transfer:

    dest_n = 0
    DO dest_pe = 0, d_nprocs-1
      idx = pe2idx(dest_pe)
      ! remote ffsl-region:
      gp_lat1              = gdc(idx)%nlat+1-gdc(idx)%ffsl%latn ! translate fssl-lat to gp-lat
      gp_lat2              = gdc(idx)%nlat+1-gdc(idx)%ffsl%lats ! ''
      IF (lassert) CALL assert(__LINE__, gp_lat1 <= gp_lat2)
      send_ir(dest_pe)     = 0 ! which gp region to send
      send_jr(dest_pe)     = 0 ! on which substripe to put the data at the receiver
      send_hsize(dest_pe)  = 0 ! horizontal size of one my regions (they must have equal size)
      send_nr(dest_pe)     = 0 ! number of regions to transmit

      hsize(:)             = 0 ! see above
      gp2ffsl_rmap(:)      = 0 ! 
      ffsl2gp_rmap(:)      = 0 ! 

      ! i-loop over local gp-regions:
      DO i = 1, 2
        IF (gp_lat1 <= ldc%glats(i) .AND. ldc%glate(i) <= gp_lat2) THEN
          send_nr(dest_pe)      = send_nr(dest_pe) + 1
          hsize(i)              = nglon*nglh(i) ! local region size
          IF ( gp_lat1 == ldc%glats(i) ) THEN
            gp2ffsl_rmap(i)       = 1
            ffsl2gp_rmap(1)       = i
          ELSEIF ( gp_lat2 == ldc%glate(i)) THEN
            gp2ffsl_rmap(i)       = 2
            ffsl2gp_rmap(2)       = i
          ELSE
            CALL finish('mo_tr_gp_ffsl:prepare_tr','unexpected region map')
          ENDIF
        ENDIF          
      ENDDO
      IF (send_nr(dest_pe) > 0) THEN
        IF (lassert) CALL assert(__LINE__, hsize(1) == 0 .OR. hsize(2) == 0 .OR. hsize(1) == hsize(2))
        send_hsize(dest_pe)  = MAXVAL(hsize)
        send_ir(dest_pe)     = MERGE(1, 0, gp2ffsl_rmap(1) > 0) + MERGE(2, 0, gp2ffsl_rmap(2) > 0)
        send_jr(dest_pe)     = MERGE(1, 0, ffsl2gp_rmap(1) > 0) + MERGE(2, 0, ffsl2gp_rmap(2) > 0)
        IF (lassert .AND. send_nr(dest_pe) == 2) CALL assert(__LINE__, gp2ffsl_rmap(1) == 1 .AND. gp2ffsl_rmap(2) == 2)
        dest_n = dest_n+1
      ENDIF
      
    ENDDO

    ALLOCATE(dest(dest_n))    

    j = 0
    DO dest_pe = 0, d_nprocs-1
      IF (send_nr(dest_pe) > 0) THEN
        j = j+1
        IF (lcheck) CALL check(__LINE__, j, 1, dest_n)
        dest(j)%xpe    = dest_pe
        dest(j)%idx    = pe2idx(dest_pe)
        dest(j)%ir     = send_ir(dest_pe)
        dest(j)%nr     = send_nr(dest_pe)
        dest(j)%jr     = send_jr(dest_pe)
        dest(j)%hsize  = send_hsize(dest_pe)
      ENDIF
    ENDDO


    IF (ldebug) THEN
      ! check if send all the volume data we have
      n=0
      DO j=1, dest_n
        dest_pe=dest(j)%xpe
        n = n + dest(j)%hsize * dest(j)%nr * kstack_n(dest_pe)
      ENDDO
      CALL assert(__LINE__, n == nglon*nglat*ldc%nlev )
    ENDIF


    CALL my_alloc

  CONTAINS

    SUBROUTINE my_alloc
#ifndef NOMPI
      INTEGER(kind=MPI_ADDRESS_KIND) :: abase, addr, adiff
      INTEGER :: j, ierror, idiff
#endif

      IF (.NOT. ALLOCATED(gp_msg)) THEN

        ALLOCATE(gp_msg(dest_n))

        DO j=1, dest_n
          ALLOCATE( gp_msg(j)%buf( dest(j)%hsize, dest(j)%nr*kstack_max*tracer_max) ) ! with local region size
          gp_msg(j)%buf(:,:) = 0.0_dp
          gp_msg(j)%offset   = 0
        ENDDO

        ALLOCATE(ffsl_msg(src_n))

        DO j=1, src_n
          ALLOCATE( ffsl_msg(j)%buf( src(j)%hsize, src(j)%nr*kstack_max*tracer_max) ) ! with remote region size
          ffsl_msg(j)%buf(:,:) = 0.0_dp
          ffsl_msg(j)%offset   = 0
        ENDDO

#ifndef NOMPI
        IF (use_mpi_collectives) THEN
          ! calculate buffer offsets
          CALL MPI_GET_ADDRESS(gp_msg(1)%buf(1,1), abase, ierror)          
          IF (ierror>0) CALL finish('mo_tr_gp_ffsl:prepare_tr:my_alloc','bad case (0a)')

          DO j=1, dest_n
            CALL MPI_GET_ADDRESS(gp_msg(j)%buf(1,1), addr, ierror)
            IF (ierror>0) CALL finish('mo_tr_gp_ffsl:prepare_tr:my_alloc','bad case (1a)')
            adiff=addr - abase
            IF ( ABS(adiff/8) > HUGE(idiff) ) CALL finish('mo_tr_gp_ffsl:prepare_tr:my_alloc','bad offset')
            idiff=INT(adiff/8)
            IF ( idiff*8 /= adiff  ) CALL finish('mo_tr_gp_ffsl:prepare_tr:my_alloc','bad alignment')
            gp_msg(j)%offset=idiff
          ENDDO
          
          CALL MPI_GET_ADDRESS(ffsl_msg(1)%buf(1,1), abase, ierror)
          IF (ierror>0) CALL finish('mo_tr_gp_ffsl:prepare_tr:my_alloc','bad case (0b)')
          
          DO j=1, src_n
            CALL MPI_GET_ADDRESS(ffsl_msg(j)%buf(1,1), addr, ierror)
            IF (ierror>0) CALL finish('mo_tr_gp_ffsl:prepare_tr:my_alloc','bad case (1b)')
            adiff=addr - abase
            IF ( ABS(adiff/8) > HUGE(idiff) ) CALL finish('mo_tr_gp_ffsl:prepare_tr:my_alloc','bad offset')
            idiff=INT(adiff/8)
            IF ( idiff*8 /= adiff  ) CALL finish('mo_tr_gp_ffsl:prepare_tr:my_alloc','bad alignment')
            ffsl_msg(j)%offset=idiff
          ENDDO
        ENDIF
#endif

      END IF


    END SUBROUTINE my_alloc


  END SUBROUTINE prepare_tr

END MODULE mo_tr_gp_ffsl
