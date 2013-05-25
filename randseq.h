#ifndef RANDSEQ_H
#define RANDSEQ_H

/* generate a random sequence of size n, the problem size */
int randseq_init(int **randseq_holder, int size);
int randseq_shuffle(int *randseq, int size);
int randseq_finalize(int *randseq, int size);
int randseq_verify(int *randseq, int size);
#endif
