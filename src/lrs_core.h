#ifndef LRSCORE_H_
#define LRSCORE_H_

#include <pthread.h>
#include "lrs_opts.h"
#include "lrs_subs.h"
#include "sw/ssw_cpp.h"

struct Locate {
  int beg, end;
  char typ; // F:5+ f:5- T:3+ t:3-
  bool operator<(const Locate& ob) const {
    return beg < ob.beg;
  }
};

struct Anchor {
  char typ;
  int score;
};

void _segm();

static map<Locate, char> pis;
static map<int, vector<Locate> > loc;
static pthread_mutex_t mutex;
static map<int, Anchor> anc;

extern Option opt;
extern string seq;
extern vector<int> pos;

#endif
