#ifndef PGA_H
#define PGA_H

#include <stdio.h>
#include "addr.h"
// the type of migrant chrom
#define MIG_ELITE 1
#define MIG_RANDOM 2
//void create_mpi_chrom_type(Chrom chrom, MPI_Datatype *t);
int send_emi();
int recv_imi();
int fill_imi();
int fetch_imi(int *buffer, int *origin);
void psearch(void);
int get_best_dim(int np, int *rowSize, int *colSize);
int imi_find_origin(long long fitv);

extern int pop_size;
extern int ga_done; // indicate if ga search() is done
extern IDTYPE globalId;
extern FILE * myout;
// performance study: immigrant history buffer
#define IMI_HIST_BUFFER_SIZE 20
#endif
