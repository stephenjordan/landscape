#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

//K is the number of numbers in the subset sum instance. Due to my hardcoding,
//we must restrict it to 64*BINS. Put differently, we should normally set BINS
//to the smallest integer i such that 64*i >= K. W is the number of walkers.
//W = 20 was good for MAX3SAT.
#define W 20
int K;
int BINS;

//I'll use the convention that w is an index for walkers and k is an index for bits.

//RAND_MAX is guaranteed to be always at least 2^16. Sometimes it is more but I will
//not assume this. On my machine it is actually 2^31.

typedef struct {
  int64_t residual;
  uint64_t *m; //maximum value of K is BINS*64
}walker;

typedef struct {
  uint64_t *vminus;
  uint64_t *vplus; 
}instance;

uint64_t rand64() {
  return ((uint64_t)rand())^((uint64_t)rand()<<16)^((uint64_t)rand()<<32)^((uint64_t)rand()<<48);
}

//return the kth bit of m
int mval(uint64_t *m, int k) {
  int bin, shift;
  bin = k >> 6;
  if(bin > BINS-1) {
    printf("Bin %i out of range", bin);
    return 0;
  }
  shift = k - (bin<<6);
  return (m[bin]>>shift)&1LLU;
}

//flip the ith bit of m
void flip(uint64_t *m, int k) {
  uint64_t newval;
  int bin;
  int shift;
  bin = k>>6;
  if(bin > BINS-1) {
    printf("Bin %i out of range", bin);
    return;
  }
  shift = k - (bin<<6);
  newval = (m[bin])^(1LLU<<shift);
  m[bin] = newval;
}

void printstring(uint64_t *m) {
  int k;
  for(k = K-1; k >= 0; k--) printf("%i", mval(m,k));
  printf("\n");
}

void evaluate(walker *w, instance *sss) {
  int k;
  w->residual = 0;
  for(k = 0; k < K; k++) {
    if(mval(w->m, k)) w->residual += (int64_t)sss->vplus[k];
    else w->residual -= (int64_t)sss->vminus[k];
  }
}

void initpop(walker *pop, instance *sss) {
  int w, k;
  for(w = 0; w < W; w++) {
    pop[w].m = malloc(BINS*sizeof(uint64_t));
    if(pop[w].m == NULL) printf("Unable to allocate bins!\n");
  }
  for(w = 0; w < W; w++) {
    for(k = 0; k < BINS; k++) pop[w].m[k] = rand64();
    evaluate(&pop[w], sss);
  }
}

void freepop(walker *pop) {
  int w;
  for(w = 0; w < W; w++) free(pop[w].m);
}

void printinstance(instance sss) {
  int k;
  for(k = 0; k < K; k++) {
printf("vplus[%i] = %llu \n", k, (unsigned long long)sss.vplus[k]);
printf("vminus[%i] = -%llu \n", k, (unsigned long long)sss.vminus[k]);
  }
}

void printwalker(walker w) {
  printf("%lld\n", (long long)w.residual);
  printstring(w.m);
}

//Hop to a random neighbor by flipping one bit.
void hop(walker *cur, walker *pro, instance *sss) {
  int bflip; //the index of the bit that gets flipped
  int bin;
  bflip = rand()%K;
  //printf("flipped %i\n", bflip);
  for(bin = 0; bin < BINS; bin++) pro->m[bin] = cur->m[bin];
  flip(pro->m, bflip);
  if(mval(cur->m, bflip)) pro->residual = cur->residual - sss->vplus[bflip] - sss->vminus[bflip];
  else pro->residual = cur->residual + sss->vplus[bflip] + sss->vminus[bflip];
}

void sit(walker *cur, walker *pro) {
  int bin;
  for(bin = 0; bin < BINS; bin++) pro->m[bin] = cur->m[bin];
  pro->residual = cur->residual;
}

//ternary Bernoulli random variable
//return 0 with probability p0, 1 with probability p1, 2 otherwise 
int tern(double p0, double p1) {
  int threshold1, threshold2;
  int r;
  if(p0 < 0 || p1 < 0 || p0+p1 > 1) {
    printf("Invalid probabilities %f %f in tern.\n", p0, p1);
    return 0;
  }
  threshold1 = (int)(p0*(double)RAND_MAX);
  threshold2 = (int)((p0+p1)*(double)RAND_MAX);
  r = rand();
  if(r < threshold1) return 0;
  if(r < threshold2) return 1;
  return 2;
}

int64_t pot(int64_t residual) {
  if(residual >= 0) return residual;
  else return -residual;
}

void walk(double duration, double vscale, instance *sss) {
  int w;                //counter variable for walkers
  walker pop1[W];       //a population of W walkers
  walker pop2[W];       //a population of W walkers
  walker *cur;          //the current locations of walkers
  walker *pro;          //the locations in progress
  walker *tmp;          //temporary holder for pointer swapping
  double s;             //current value of s
  double phop;          //probability of hopping to a neighboring vertex
  double ptel;          //probability of teleporting to another walker's location
  int action;           //0 = hop, 1 = teleport, 2 = sit
  int64_t umin, umax;   //the min&max residual among the occupied locations
  int winners;          //number of times a walker hits zero residual
  double dt;            //the adjustable timestep
  double ttot;          //the total time evolution elapsed
  int dest;             //destination walker
  int stepcount;        //total number of timesteps
  printf("K = %i\n", K);
  printf("W = %i\n", W);
  printf("BINS = %i\n", BINS);
  printf("duration = %e\n", duration);
  printf("vscale = %e\n", vscale);
  //initialize the walkers to random locations
  initpop(pop1, sss);
  initpop(pop2, sss);
  cur = pop1;
  pro = pop2;
  winners = 0;
  ttot = 0;
  stepcount = 0;
  do {
    s = ttot/duration;
    //calculate the minimum potential among currently occupied locations
    umin = pot(cur[0].residual);
    umax = umin;
    for(w = 0; w < W; w++) {
      if(pot(cur[w].residual) < umin) umin = pot(cur[w].residual);
      if(pot(cur[w].residual) > umax) umax = pot(cur[w].residual);
    }
    dt = 0.99/(1.0-s+s*vscale*(double)(umax-umin)); //this ensures we have no negative probabilities
    //printf("dt = %e\n", dt);
    w = rand()%W;
    phop = (1.0-s)*dt;
    dest = 0;
    do {
      //subtracting umin yields invariance under uniform potential change
      ptel = dt*s*vscale*(double)(pot(cur[w].residual)-umin); //here we subtract the offset
      //printf("phop = %f\tptel = %f\n", phop, ptel);
      action = tern(phop, ptel);
      if(action == 2) { //sit
        sit(&cur[w], &pro[dest]);
        dest++;
      }
      //if(action == 1) walker dies, do nothing
      if(action == 0) { //hop
        hop(&cur[w], &pro[dest], sss);
        dest++;
      }
      w = (w+1)%W;
    }while(dest < W);
    //swap pro with cur
    tmp = pro;
    pro = cur;
    cur = tmp;
    stepcount++;
    for(w = 0; w < W; w++) if(pot(cur[w].residual) == 0) winners++;
    ttot += dt;
  }while(ttot < duration && winners == 0);
  if(winners > 0) {
    if(winners == 1) printf("Found 1 solution:\n");
    else printf("Found %i solutions:\n", winners);
    for(w = 0; w < W; w++) if(pot(cur[w].residual) == 0) printstring(cur[w].m);
  }
  //if no satisfying assignments were found, print the best ones------------------
  else {
    umin = pot(cur[0].residual);
    umax = umin;
    for(w = 0; w < W; w++) {
      if(pot(cur[w].residual) < umin) umin = pot(cur[w].residual);
      if(pot(cur[w].residual) > umax) umax = pot(cur[w].residual);
    }
    printf("Best solutions found have %llu potential.\n", (long long)umin);
    for(w = 0; w < W; w++) if(pot(cur[w].residual) == umin) printstring(cur[w].m);
  }
  printf("stepcount: %i\n", stepcount);
  freepop(pop1);
  freepop(pop2);
  //for(w = 0; w < W; w++) printwalker(cur[w]);
}

int main() {
  instance sss;         //the instance of subset sum
  unsigned int seed;    //the seed for the rng
  double duration;      //physical time
  int k;                //counter variable for subset sum numbers
  int maxpow;           //the biggest allowed number is 2^maxpow
  uint64_t maxval;      //the biggest allowed number is 2^maxpow
  double vscale;        //scales the residual
  //seed = time(NULL);
  seed = 1481847052;    //fixed seed for testing
  srand(seed);
  printf("seed = %u\n", seed);
  duration = 1000;
  K = 20;
  BINS = 1;
  sss.vminus = malloc(K*sizeof(uint64_t));
  sss.vplus = malloc(K*sizeof(uint64_t));
  if(sss.vminus == NULL || sss.vplus == NULL) {
    printf("Unable to allocate instance.\n");
    return 0;
  }
  //generate a random instance
  maxpow = 11; //max value is 2^maxpow
  maxval = 1;
  maxval <<= maxpow;
  printf("maxval = 2^%d\n", maxpow);
  for(k = 0; k < K; k++) {
    sss.vminus[k] = rand64()%maxval;
    sss.vplus[k] = rand64()%maxval;
  }
  printinstance(sss);
  vscale = (double)100.0/(double)maxval; //a guess, really
  walk(duration, vscale, &sss);
  free(sss.vminus);
  free(sss.vplus);
  return 0;
}
