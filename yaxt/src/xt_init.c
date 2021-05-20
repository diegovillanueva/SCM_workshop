/**
 * @file xt_init.c
 *
 * @copyright Copyright  (C)  2012 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://doc.redmine.dkrz.de/yaxt/html/
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include "core/core.h"
#include "xt_idxlist_internal.h"
#include "xt_idxstripes_internal.h"
#include "xt_idxsection_internal.h"
#include "xt_idxempty_internal.h"
#include "xt_exchanger.h"
#include "xt_mpi_internal.h"
#include "instr.h"

INSTR_DEF(instr,"YAXT_lifetime")

static enum {
  xt_lib_pre_init,
  xt_lib_initialized,
  xt_lib_finalized,
} xt_lib_state = xt_lib_pre_init;

void
xt_initialize(MPI_Comm default_comm)
{
  Xt_default_comm = default_comm;
  xt_mpi_init();
  xt_idxempty_init();
  xt_idxstripes_initialize();
  xt_idxsection_initialize();
  xt_idxlist_intersection_init();
#if ! defined TLS && defined HAVE_PTHREAD
  xt_exchanger_init();
#endif
#ifdef INSTR_WITH_SCT
  sct_init(32, "YAXT", Xt_default_comm);
  setenv("SCT_PROC_CHOICE", "SCT_REDUCE_ALL", 0);
  setenv("SCT_CALLSTATS", "0", 0);
  INSTR_START(instr);
#endif
  xt_lib_state = xt_lib_initialized;
}

void
xt_finalize(void)
{
  if (xt_lib_state == xt_lib_initialized)
  {
    xt_idxsection_finalize();
    xt_idxstripes_finalize();
    xt_idxempty_finalize();
#if ! defined TLS && defined HAVE_PTHREAD
    xt_exchanger_finalize();
#endif
    xt_mpi_finalize();
#ifdef INSTR_WITH_SCT
    INSTR_STOP(instr);
    INSTR_REPORT(SCT_GETENV, SCT_GETENV, SCT_GETENV);
    INSTR_FINALIZE();
#endif
    xt_lib_state = xt_lib_finalized;
  }
}

int
xt_initialized(void)
{
  return xt_lib_state > xt_lib_pre_init;
}

int
xt_finalized(void)
{
  return xt_lib_state == xt_lib_finalized;
}
