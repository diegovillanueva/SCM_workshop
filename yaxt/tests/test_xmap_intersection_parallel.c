/**
 * @file test_xmap_intersection_parallel.c
 *
 * @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
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

#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>

#include <yaxt.h>

#include "tests.h"
#include "xt/xt_xmap_intersection.h"

struct test_message {

  int rank;       // rank of communication partner
  const int *pos;   // positions to be sent/received
  int num_pos; // number of positions
};

static void test_xmap(
  Xt_xmap xmap,
  int num_sends, const struct test_message send_messages[num_sends],
  int num_recvs, const struct test_message recv_messages[num_recvs]);

static Xt_xmap (*xmi_new)(
  int nsrc_com, const struct Xt_com_list src_com[nsrc_com],
  int ndst_com, const struct Xt_com_list dst_com[ndst_com],
  Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist, MPI_Comm comm);

static void
parse_options(int *argc, char ***argv);

int main(int argc, char **argv)
{

  // init mpi
  xt_mpi_call(MPI_Init(&argc, &argv), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  int my_rank, comm_size;

  xt_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &my_rank), MPI_COMM_WORLD);
  xt_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &comm_size), MPI_COMM_WORLD);

  if (comm_size != 3) {

    xt_finalize();
    MPI_Finalize();
    return 77;
  }

  parse_options(&argc, &argv);

  { // simple test (round robin)

    // setup

    Xt_int src_index = (Xt_int)my_rank;
    Xt_int dst_index = (Xt_int)((my_rank + 1)%comm_size);
    Xt_idxlist src_idxlist = xt_idxvec_new(&src_index, 1);
    Xt_idxlist dst_idxlist = xt_idxvec_new(&dst_index, 1);
    int num_src_intersections = 1;
    struct Xt_com_list src_com = {.list = src_idxlist,
                                  .rank = (my_rank+1)%comm_size};
    int num_dst_intersections = 1;
    struct Xt_com_list dst_com = {.list = dst_idxlist,
                                  .rank = (my_rank+comm_size-1)%comm_size};

    Xt_xmap xmap = xmi_new(num_src_intersections, &src_com,
                           num_dst_intersections, &dst_com,
                           src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // test

    int send_pos = 0;
    int num_sends = 1;
    struct test_message send_messages[1] = {{.rank = (my_rank+1)%comm_size,
                                             .pos = &send_pos, .num_pos = 1}};
    int recv_pos = 0;
    int num_recvs = 1;
    struct test_message recv_messages[1] = {{.rank = (my_rank+comm_size-1)%comm_size,
                                             .pos = &recv_pos, .num_pos = 1}};

    test_xmap(xmap, num_sends, send_messages, num_recvs, recv_messages);

    // cleanup

    xt_xmap_delete(xmap);
    xt_idxlist_delete(dst_idxlist);
    xt_idxlist_delete(src_idxlist);
  }

  { // rank 0 receives the same point from rank 1 and 2

    // setup

    Xt_int src_index = 0;
    Xt_int dst_index = 0;
    Xt_idxlist src_idxlist = (my_rank == 0)?xt_idxempty_new():xt_idxvec_new(&src_index, 1);
    Xt_idxlist dst_idxlist = (my_rank == 0)?xt_idxvec_new(&dst_index, 1):xt_idxempty_new();
    int num_src_intersections = (my_rank == 0)?0:1;
    struct Xt_com_list src_com = {.list = src_idxlist,
                                  .rank = 0};
    int num_dst_intersections = (my_rank == 0)?2:0;
    struct Xt_com_list dst_com[2] = {{.list = dst_idxlist, .rank = 1},
                                     {.list = dst_idxlist, .rank = 2}};

    Xt_xmap xmap = xmi_new(num_src_intersections, &src_com,
                           num_dst_intersections, dst_com,
                           src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // test

    int send_pos = 0;
    int num_sends = (my_rank == 1)?1:0;
    struct test_message send_messages[1] = {{.rank = 0, .pos = &send_pos,
                                             .num_pos = 1}};
    int recv_pos = 0;
    int num_recvs = (my_rank == 0)?1:0;
    struct test_message recv_messages[1] = {{.rank = 1, .pos = &recv_pos,
                                             .num_pos = 1}};

    test_xmap(xmap, num_sends, send_messages, num_recvs, recv_messages);

    // cleanup

    xt_xmap_delete(xmap);
    xt_idxlist_delete(dst_idxlist);
    xt_idxlist_delete(src_idxlist);
  }

  { // all ranks can receive data from the others
    // rank               |  0  |  1  |  2  |
    // source indices     | 1,2 | 2,0 | 0,1 |
    // destination indice |  0  |  1  |  2  |

    // setup

    Xt_int src_indices[2] = { (Xt_int)((my_rank+1)%comm_size),
                              (Xt_int)((my_rank+2)%comm_size) };
    Xt_int dst_index = (Xt_int)my_rank;
    Xt_idxlist src_idxlist = xt_idxvec_new(src_indices, 2);
    Xt_idxlist dst_idxlist = xt_idxvec_new(&dst_index, 1);
    int num_src_intersections[3] = {2, 1, 0};
    Xt_idxlist src_intersection_idxlist[2] = {xt_idxvec_new(src_indices, 1),
                                              xt_idxvec_new(src_indices+1, 1)};
    struct Xt_com_list src_com[2] = {{.list = src_intersection_idxlist[0],
                                      .rank = 1},
                                     {.list = src_intersection_idxlist[1],
                                      .rank = (my_rank == 0)?2:0}};
    int num_dst_intersections = 1;
    struct Xt_com_list dst_com = {.list = dst_idxlist, .rank = (my_rank == 0)?1:0};

    Xt_xmap xmap = xmi_new(num_src_intersections[my_rank], src_com + my_rank,
                           num_dst_intersections, &dst_com,
                           src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // test

    if (my_rank == 0) {

      int send_pos[2] = {0,1};
      int num_sends = 2;
      struct test_message send_messages[2] = {{.rank = 1, .pos = send_pos+0,
                                               .num_pos = 1},
                                              {.rank = 2, .pos = send_pos+1,
                                               .num_pos = 1}};
      int recv_pos = 0;
      int num_recvs = 1;
      struct test_message recv_messages[1] = {{.rank = 1, .pos = &recv_pos,
                                               .num_pos = 1}};

      test_xmap(xmap, num_sends, send_messages, num_recvs, recv_messages);

    } else if (my_rank == 1) {

      int send_pos = 1;
      int num_sends = 1;
      struct test_message send_messages[1] = {{.rank = 0, .pos = &send_pos,
                                               .num_pos = 1}};
      int recv_pos = 0;
      int num_recvs = 1;
      struct test_message recv_messages[1] = {{.rank = 0, .pos = &recv_pos,
                                               .num_pos = 1}};

      test_xmap(xmap, num_sends, send_messages, num_recvs, recv_messages);

    } else {

      int num_sends = 0;

      int recv_pos = 0;
      int num_recvs = 1;
      struct test_message recv_messages[1] = {{.rank = 0, .pos = &recv_pos,
                                               .num_pos = 1}};

      test_xmap(xmap, num_sends, NULL, num_recvs, recv_messages);
    }

    // cleanup

    xt_xmap_delete(xmap);
    xt_idxlist_delete(src_intersection_idxlist[1]);
    xt_idxlist_delete(src_intersection_idxlist[0]);
    xt_idxlist_delete(dst_idxlist);
    xt_idxlist_delete(src_idxlist);
  }

  { // all ranks can receive data from the others
    // rank               |         0         |         1         |         2         |
    // source indices     |     0,1,2,3,4     |     3,4,5,6,7     |     6,7,8,0,1     |
    // destination indice | 0,1,2,3,4,5,6,7,8 | 0,1,2,3,4,5,6,7,8 | 0,1,2,3,4,5,6,7,8 |

    // setup

    static const Xt_int src_indices[3][5]
      = {{0,1,2,3,4}, {3,4,5,6,7}, {6,7,8,0,1}};
    static const Xt_int dst_indices[9] = {0,1,2,3,4,5,6,7,8};
    Xt_idxlist src_idxlist = xt_idxvec_new(src_indices[my_rank], 5);
    Xt_idxlist dst_idxlist = xt_idxvec_new(dst_indices, 9);

    struct Xt_com_list src_com[3] = {{.list = src_idxlist, .rank = 0},
                                     {.list = src_idxlist, .rank = 1},
                                     {.list = src_idxlist, .rank = 2}};
    struct Xt_com_list dst_com[3] =
      {{.list = xt_idxvec_new(src_indices[0], 5), .rank = 0},
       {.list = xt_idxvec_new(src_indices[1], 5), .rank = 1},
       {.list = xt_idxvec_new(src_indices[2], 5), .rank = 2}};

    Xt_xmap xmap = xmi_new(3, src_com, 3, dst_com,
                           src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // test

    static const int send_pos[3][5] = {{0,1,2,3,4}, {2,3,4}, {2}};
    static const int num_send_pos[3] = {5, 3, 1};
    struct test_message send_messages[3] = {{.rank = 0,
                                             .pos = send_pos[my_rank],
                                             .num_pos = num_send_pos[my_rank]},
                                            {.rank = 1,
                                             .pos = send_pos[my_rank],
                                             .num_pos = num_send_pos[my_rank]},
                                            {.rank = 2,
                                             .pos = send_pos[my_rank],
                                             .num_pos = num_send_pos[my_rank]}};
    static const int recv_pos[3][5] = {{0,1,2,3,4}, {5,6,7}, {8}};
    static const int num_recv_pos[3] = {5, 3, 1};
    struct test_message recv_messages[3] = {{.rank = 0, .pos = recv_pos[0],
                                             .num_pos = num_recv_pos[0]},
                                            {.rank = 1, .pos = recv_pos[1],
                                             .num_pos = num_recv_pos[1]},
                                            {.rank = 2, .pos = recv_pos[2],
                                             .num_pos = num_recv_pos[2]}};

    test_xmap(xmap, 3, send_messages, 3, recv_messages);

    // cleanup

    xt_xmap_delete(xmap);
    xt_idxlist_delete(dst_com[2].list);
    xt_idxlist_delete(dst_com[1].list);
    xt_idxlist_delete(dst_com[0].list);
    xt_idxlist_delete(dst_idxlist);
    xt_idxlist_delete(src_idxlist);
  }

  { // one rank receives data from the other two, that have duplicated indices
    // (this provokes a bug found by Joerg Behrens)
    // rank               |  0  |  1  |   2   |
    // source indices     | 0,2 | 1,2 |       |
    // destination indice |     |     | 0,1,2 |

    // setup

    Xt_int src_indices[2][2] = {{0,2}, {1,2}};

    struct Xt_com_list * src_com, * dst_com;
    int num_src_intersections, num_dst_intersections;
    Xt_idxlist src_idxlist, dst_idxlist;

    if (my_rank == 2) {

      src_com = NULL;
      num_src_intersections = 0;

      dst_com = malloc(2 * sizeof(*dst_com));
      num_dst_intersections = 2;
      dst_com[0].list = xt_idxvec_new(src_indices[0], 2);
      dst_com[0].rank = 0;
      dst_com[1].list = xt_idxvec_new(src_indices[1], 2);
      dst_com[1].rank = 1;

      Xt_int dst_indices[3] = {0,1,2};

      src_idxlist = xt_idxempty_new();
      dst_idxlist = xt_idxvec_new(dst_indices, 3);

    } else {

      src_com = malloc(1 * sizeof(*src_com));
      src_com->list = xt_idxvec_new(src_indices[my_rank], 2);
      src_com->rank = 2;
      num_src_intersections = 1;

      dst_com = NULL;
      num_dst_intersections = 0;

      src_idxlist = xt_idxvec_new(src_indices[my_rank], 2);
      dst_idxlist = xt_idxempty_new();
    }

    Xt_xmap xmap = xmi_new(num_src_intersections, src_com,
                           num_dst_intersections, dst_com,
                           src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // test

      if (my_rank == 2) {

        static const int recv_pos[2][2] = {{0,2}, {1}};
        static const struct test_message recv_messages[2]
          = {{.rank = 0, .pos = recv_pos[0], .num_pos = 2},
             {.rank = 1, .pos = recv_pos[1], .num_pos = 1}};

        test_xmap(xmap, 0, NULL, 2, recv_messages);

      } else {

        static const int send_pos[2][2] = {{0,1}, {0}};
        static const int num_send_pos[2] = {2, 1};
        struct test_message send_messages = {.rank = 2,
                                             .pos = send_pos[my_rank],
                                             .num_pos = num_send_pos[my_rank]};

        test_xmap(xmap, 1, &send_messages, 0, NULL);
      }

    // cleanup

    xt_xmap_delete(xmap);
    xt_idxlist_delete(dst_idxlist);
    xt_idxlist_delete(src_idxlist);
    for (int i = 0; i < num_dst_intersections; ++i)
      xt_idxlist_delete(dst_com[i].list);
    free(dst_com);
    for (int i = 0; i < num_src_intersections; ++i)
      xt_idxlist_delete(src_com[i].list);
    free(src_com);
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void test_xmap_iter(Xt_xmap_iter iter, int num_msgs,
                           const struct test_message msgs[num_msgs]) {

  if (num_msgs == 0) {

    if (iter != NULL)
      PUT_ERR("ERROR: xt_xmap_get_*_iterator (iter should be NULL)\n");

  } else if (iter == NULL) {

    PUT_ERR("ERROR: xt_xmap_get_*_iterator (iter should not be NULL)\n");

  } else {

    int i = 0;

    do {

      if (xt_xmap_iterator_get_rank(iter) != msgs[i].rank)
        PUT_ERR("ERROR: xt_xmap_iterator_get_rank\n");

      if (xt_xmap_iterator_get_num_transfer_pos(iter) != msgs[i].num_pos)
        PUT_ERR("ERROR: xt_xmap_iterator_get_num_transfer_pos\n");

      const int *restrict pos = xt_xmap_iterator_get_transfer_pos(iter);

      for (int j = 0; j < msgs[i].num_pos; ++j)
        if (pos[j] != msgs[i].pos[j])
          PUT_ERR("ERROR: xt_xmap_iterator_get_transfer_pos\n");

      ++i;

    } while (xt_xmap_iterator_next(iter));

    if (i != num_msgs)
      PUT_ERR("ERROR: xt_xmap_iterator_next (wrong number of message)\n");
  }
}

static void test_xmap(
  Xt_xmap xmap,
  int num_sends, const struct test_message send_messages[num_sends],
  int num_recvs, const struct test_message recv_messages[num_recvs]) {

  if (xt_xmap_get_num_destinations(xmap) != num_sends)
    PUT_ERR("ERROR: xt_xmap_get_num_destinations\n");
  if (xt_xmap_get_num_sources(xmap) != num_recvs)
    PUT_ERR("ERROR: xt_xmap_get_num_sources\n");

  Xt_xmap_iter send_iter = xt_xmap_get_out_iterator(xmap);
  Xt_xmap_iter recv_iter = xt_xmap_get_in_iterator(xmap);

  test_xmap_iter(send_iter, num_sends, send_messages);
  test_xmap_iter(recv_iter, num_recvs, recv_messages);

  if (recv_iter != NULL) xt_xmap_iterator_delete(recv_iter);
  if (send_iter != NULL) xt_xmap_iterator_delete(send_iter);
}

static void
parse_options(int *argc, char ***argv)
{
  xmi_new = xt_xmap_intersection_new;
  int opt;
  while ((opt = getopt(*argc, *argv, "m:")) != -1) {
    switch (opt) {
    case 'm':
      if (!strcmp(optarg, "xt_xmap_intersection_new"))
        xmi_new = xt_xmap_intersection_new;
      else if (!strcmp(optarg, "xt_xmap_intersection_ext_new"))
        xmi_new = xt_xmap_intersection_ext_new;
      else
      {
        fprintf(stderr, "Unknown xmap intersection constructor requested %s\n",
                optarg);
        exit(EXIT_FAILURE);
      }
    }
  }
}
