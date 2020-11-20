#ifndef LRSOPTS_H_
#define LRSOPTS_H_

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <algorithm>
#include <sstream>

using namespace std;

struct Option {
  string inp, out, pm5, pm3, fmt;
  int umi, cln, len, thd;
  bool gtail, help;
};

#endif
