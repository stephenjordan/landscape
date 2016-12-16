#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//m is an array of size BINS. m stores the ones and zeros specifying the subset.
//They are packed into 64-bit blocks. Residual is the value of the subset sum
//instance.
typedef struct {
  int64_t residual;
  uint64_t *m;
}walker;

//K is the number of numbers in the subset sum instance. Due to my hardcoding,
//we must restrict it to 64*BINS. Put differently, we should normally set BINS
//to the smallest integer such that 64*i >= k. bincount() will compute this
//number.
typedef struct {
  int K; 
  int BINS;
  uint64_t *vminus;
  uint64_t *vplus; 
}instance;

int bincount(int K) {
  if(K%64 == 0) return K>>6;
  return (K>>6)+1;
}

//RAND_MAX is guaranteed to be always at least 2^16. Sometimes it is more but I will
//not assume this. On my machine it is actually 2^31.
uint64_t rand64() {
  return ((uint64_t)rand())^((uint64_t)rand()<<16)^((uint64_t)rand()<<32)^((uint64_t)rand()<<48);
}

//return the kth bit of m
int mval(uint64_t *m, int k) {
  int bin, shift;
  bin = k >> 6;
  shift = k - (bin<<6);
  return (m[bin]>>shift)&1LLU;
}

//flip the ith bit of m
void flip(uint64_t *m, int k) {
  uint64_t newval;
  int bin;
  int shift;
  bin = k>>6;
  shift = k - (bin<<6);
  newval = (m[bin])^(1LLU<<shift);
  m[bin] = newval;
}

void printstring(uint64_t *m, int K) {
  int k;
  for(k = K-1; k >= 0; k--) printf("%i", mval(m,k));
  printf("\n");
}

void evaluate(walker *w, instance *sss) {
  int k;
  w->residual = 0;
  for(k = 0; k < sss->K; k++) {
    if(mval(w->m, k)) w->residual += (int64_t)sss->vplus[k];
    else w->residual -= (int64_t)sss->vminus[k];
  }
}

void initpop(walker *pop, int W, instance *sss) {
  int w, k;
  for(w = 0; w < W; w++) {
    pop[w].m = malloc(sss->BINS*sizeof(uint64_t));
    if(pop[w].m == NULL) printf("Unable to allocate bins!\n");
  }
  for(w = 0; w < W; w++) {
    for(k = 0; k < sss->BINS; k++) pop[w].m[k] = rand64();
    evaluate(&pop[w], sss);
  }
}

void freepop(walker *pop, int W) {
  int w;
  for(w = 0; w < W; w++) free(pop[w].m);
}

void printinstance(instance sss) {
  int k;
  printf("k\tv+\tv-\n");
  for(k = 0; k < sss.K; k++)
    printf("%i\t%llu\t-%llu \n", k, (unsigned long long)sss.vplus[k], (unsigned long long)sss.vminus[k]);
}

//Hop to a random neighbor by flipping one bit.
void hop(walker *cur, walker *pro, instance *sss) {
  int bflip; //the index of the bit that gets flipped
  int bin;
  bflip = rand()%sss->K;
  //printf("flipped %i\n", bflip);
  for(bin = 0; bin < sss->BINS; bin++) pro->m[bin] = cur->m[bin];
  flip(pro->m, bflip);
  if(mval(cur->m, bflip)) pro->residual = cur->residual - sss->vplus[bflip] - sss->vminus[bflip];
  else pro->residual = cur->residual + sss->vplus[bflip] + sss->vminus[bflip];
}

//hop to a neigbor or next-neighbor by flipping one or two bits
void hop2(walker *cur, walker *pro, instance *sss) {
  int type; //0 = 1 hop, 2 = 2 hops
  type = rand()%2;
  hop(cur, pro, sss);
  if(type) hop(cur, pro, sss);
}

void sit(walker *cur, walker *pro, instance *sss) {
  int bin;
  for(bin = 0; bin < sss->BINS; bin++) pro->m[bin] = cur->m[bin];
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

int walk(double duration, double vscale, int W, instance *sss) {
  int w;                //counter variable for walkers
  walker *pop1;         //a population of W walkers
  walker *pop2;         //a population of W walkers
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
  //initialize the walkers to random locations
  pop1 = (walker *)malloc(W*sizeof(walker));
  pop2 = (walker *)malloc(W*sizeof(walker));
  initpop(pop1, W, sss);
  initpop(pop2, W, sss);
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
        sit(&cur[w], &pro[dest], sss);
        dest++;
      }
      //if(action == 1) walker dies, do nothing
      if(action == 0) { //hop
        //hop2(&cur[w], &pro[dest], sss); //doesn't seem to really help
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
    for(w = 0; w < W; w++) if(pot(cur[w].residual) == 0) printstring(cur[w].m, sss->K);
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
    for(w = 0; w < W; w++) if(pot(cur[w].residual) == umin) printstring(cur[w].m, sss->K);
  }
  printf("stepcount: %i\n", stepcount);
  freepop(pop1, W);
  freepop(pop2, W);
  free(pop1);
  free(pop2);
  return winners;
}

int main(int argc, char *argv[]) {
  int W;                //the number of walkers
  instance sss;         //the instance of subset sum
  unsigned int seed;    //the seed for the rng
  double duration;      //physical time
  int k;                //counter variable for subset sum numbers
  double vscale;        //scales the residual
  FILE *fp;             //the file containing the instance
  long long int vp, vm; //for loading vplus and vminus
  uint64_t maxval;      //the biggest value in the instance
  int trials;           //fresh trials of Markov Chain
  int t;                //trial counter
  int success;          //flag set to one if exact solution found
  if(argc != 5) {
    printf("usage: landscape walkers trials duration input.sss\n");
    return 0;
  }
  W = atoi(argv[1]);
  trials = atoi(argv[2]);
  duration = atof(argv[3]);
  fp = fopen(argv[4], "r");
  if(fp == NULL) {
    printf("Unable to read instance file %s\n", argv[4]);
    return 0;
  }
  fscanf(fp, "%i", &sss.K);
  sss.BINS = bincount(sss.K);
  sss.vminus = malloc(sss.K*sizeof(uint64_t));
  sss.vplus = malloc(sss.K*sizeof(uint64_t));
  if(sss.vminus == NULL || sss.vplus == NULL) {
    printf("Unable to allocate instance.\n");
    return 0;
  }  
  for(k = 0; k < sss.K; k++) {
    fscanf(fp, "%lld %lld", &vp, &vm);
    sss.vplus[k]  = (uint64_t)vp;
    sss.vminus[k] = (uint64_t)vm;
  }
  printf("W = %d\n", W);
  printf("trials = %d\n", trials);
  printf("duration = %lf\n", duration);
  printf("infile= %s\n", argv[4]);
  seed = time(NULL);
  //seed = 1481847052;    //fixed seed for testing
  srand(seed);
  //for(k = 0; k < 160; k++) rand(); //for testing against earlier versions
  printf("seed = %u\n", seed);
  printinstance(sss);
  maxval = 0;
  for(k = 0; k < sss.K; k++) {
    if(sss.vminus[k] > maxval) maxval = sss.vminus[k];
    if(sss.vplus[k] > maxval) maxval = sss.vplus[k];
  }
  printf("maxval = %llu\n", (unsigned long long)maxval);
  //vscale = (double)100.0/(double)2048; //for testing against earlier versions
  vscale = (double)100.0/(double)maxval; //a guess, really
  printf("vscale = %e\n", vscale);
  t = 0;
  do {
    printf("trial %i------------------------------------\n", t+1);
    success = walk(duration, vscale, W, &sss);
    t++;
  }while(t < trials && !success);
  free(sss.vminus);
  free(sss.vplus);
  fclose(fp); 
  return 0;
}
