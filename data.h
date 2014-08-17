#ifndef DATA_H
#define DATA_H

/* read data into memory
 * data format conforms OR-LIB format
 */

typedef int VTYPE; // type of value
typedef int WTYPE; // type of weight

// matrix access: size of a row is implicitly decided as n
#define M(matrix, i, j) (matrix[n * i + j])
// transposed matrix access
#define MT(matrix, i, j) (matrix[m * i + j])
#define Mr(matrix, i) (matrix[n * i])

void readData(char *fileName);

#endif
