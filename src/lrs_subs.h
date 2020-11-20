#ifndef LRSSUBS_H_
#define LRSSUBS_H_

#include <iostream>
#include <getopt.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <algorithm>
#include <sstream>
#include "lrs_opts.h"

using namespace std;

extern Option opt;

void _usage();
void _getopt(int argc, char* argv[]);
string _revc(string &str, int flag);

#endif

