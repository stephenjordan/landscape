#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct {
  int K;               //number of terms
  uint64_t *vminus;    //positive terms
  uint64_t *vplus;     //negative terms
}instance;

void fprintinstance(FILE *fp, instance sss) {
  int k;
  fprintf(fp, "%d\n", sss.K);
  for(k = 0; k < sss.K; k++)
    fprintf(fp, "%llu\t%llu\n", (unsigned long long)sss.vplus[k], (unsigned long long)sss.vminus[k]);
}

uint64_t rand64() {
  return ((uint64_t)rand())^((uint64_t)rand()<<16)^((uint64_t)rand()<<32)^((uint64_t)rand()<<48);
}

int main(int argc, char *argv[]) {
  instance sss;         //the instance of subset sum
  unsigned int seed;    //the seed for the rng
  int k;                //counter variable for subset sum numbers
  int maxpow;           //the biggest allowed number is 2^maxpow
  uint64_t maxval;      //the biggest allowed number is 2^maxpow
  FILE *fp;
  if(argc != 4) {
    printf("Usage: generate K maxpow filename.sss\n");
    return 0;
  }
  sss.K = atoi(argv[1]);
  maxpow = atoi(argv[2]);
  seed = time(NULL);
  srand(seed);
  fp = fopen(argv[3], "w");
  if(fp == NULL) {
    printf("Unable to create file %s.\n", argv[3]);
    return 0;
  }
  sss.vminus = malloc(sss.K*sizeof(uint64_t));
  sss.vplus = malloc(sss.K*sizeof(uint64_t));
  if(sss.vminus == NULL || sss.vplus == NULL) {
    printf("Unable to allocate instance.\n");
    return 0;
  }
  //generate a random instance
  maxval = 1;
  maxval <<= maxpow;
  for(k = 0; k < sss.K; k++) {
    sss.vminus[k] = rand64()%maxval;
    sss.vplus[k] = rand64()%maxval;
  }
  fprintinstance(fp, sss);
  free(sss.vminus);
  free(sss.vplus);
  fclose(fp);
  return 0;
}
