#ifndef RANDSEQ_C
#define RANDSEQ_C
/* generate a random sequence of size n, the problem size */
#include <stdlib.h>
#include <string.h>
#include "randseq.h"
#include "ga.h"

// init the random sequence memory
int randseq_init(int **randseq_holder, int size) {
	int i, j;
	char ihash[size];
	int *randseq; 
	if (*randseq_holder == NULL) {
		randseq = (int *)malloc(sizeof(int) * size);
		*randseq_holder = randseq;
	}
	memset(ihash, 0, sizeof(char) * size);
	j = 0;
	while (j<size) {
		i = MYRANDI(size);
		if (ihash[i] == 0) { // this number not used before
			randseq[j] = i;
			ihash[i] = 1;
			j++;
		}
	}
	return 1;
}
// shuffle the random sequence memory
int randseq_shuffle(int *randseq, int size)
{
	if (randseq == NULL || size < 2) return 0;
	int num, p, q, i, tmp;
	// decide number of exchanges
	num = MYRANDI(size/2) + 1;
	for (i=0; i<num; i++) {
		// select two indices randomly
		p = MYRANDI(size);
		do {
			q = MYRANDI(size);
		} while (p == q);
		tmp = randseq[q];
		randseq[q] = randseq[p];
		randseq[p] = tmp;
#ifdef DEBUG_RANDSEQ
		printf("%d<->%d, ", p, q);
#endif
	}
	return 1;
}
int randseq_finalize(int *randseq, int size)
{
	if (randseq != NULL) free(randseq);
	return 1;
}
int randseq_verify(int *randseq, int size)
{
	int i; char ihash[size];
	memset(ihash, 0, sizeof(int) * size);
	for (i=0; i<size; i++) {
		if (ihash[randseq[i]] > 0) return 0;
		else ihash[randseq[i]] ++;
	}
	return 1;
}

#endif
