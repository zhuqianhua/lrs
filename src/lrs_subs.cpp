#include "lrs_subs.h"

void _usage()
{
  printf ("\nProgram: lrs (Long Reads Segmentation)\n");
  printf ("Contact: Qianhua ZHU <zhuqianhua@bgi.com>\n\n");
  printf ("Usage  : lrs [options]\n\n");
  printf ("Options: -i, --input  *<s>  input file *.bam or *.fa\n");
  printf ("         -f, --format *<s>  input file format, b: bam and f: fasta\n");
  printf ("         -o, --out     <s>  prefix of output, default ./lrs\n");
  printf ("         -u, --umi     <i>  length of umi, default 8\n");
  printf ("         -a, --primer5 <s>  5' primer sequence, default AAGCAGTGGTATCAACGCAGAGTACGGGGG\n");
  printf ("         -b, --primer3 <s>  3' primer sequence, default AAGCAGTGGTATCAACGCAGAGTAC\n");
  printf ("         -G, --pm5g    <b>  5' primer with GGGGG and 3' primer without GGGGG\n");
  printf ("         -L, --cln     <i>  minimum continuously align length of primer, default 5\n");
  printf ("         -l, --len     <i>  minimum align length for primer, default 16\n");
  printf ("         -p, --process <i>  number of threads, default 6\n");
  printf ("         -h, --help    <b>  print this information\n\n");
  exit(0);
}

void _getopt(int argc, char* argv[])
{
  char const * shortOpt = "i:f:o:u:a:b:GL:l:p:h";
  struct option longOpt[] = {
    {"input", 1, NULL, 'i'},
    {"format", 1, NULL, 'f'},
    {"out", 1, NULL, 'o'},
    {"umi", 1, NULL, 'u'},
    {"primer5", 1, NULL, 'a'},
    {"primer3", 1, NULL, 'b'},
    {"pm5g", 1, NULL, 'G'},
    {"cln", 1, NULL, 'L'},
    {"len", 1, NULL, 'l'},
    {"process", 1, NULL, 'p'},
    {"help", 0, NULL, 'h'},
    {NULL, 0, NULL, 0},
  };
  int nextOpt;
  while ((nextOpt = getopt_long(argc, argv, shortOpt, longOpt, NULL)) != -1) {
    switch (nextOpt) {
      case 'i':
        opt.inp = optarg;
        break;
      case 'f':
        opt.fmt = optarg;
        break;
      case 'o':
        opt.out = optarg;
        break;
      case 'u':
        opt.umi = atoi(optarg);
        break;
      case 'a':
        opt.pm5 = optarg;
        break;
      case 'b':
        opt.pm3 = optarg;
        break;
      case 'G':
        opt.gtail = true;
        break;
      case 'L':
        opt.cln = atoi(optarg);
        break;
      case 'l':
        opt.len = atoi(optarg);
        break;
      case 'p':
        opt.thd = atoi(optarg);
        break;
      case 'h':
        opt.help = true;
        break;
    }
  }
  if (opt.help == true or opt.inp.compare("") == 0 or opt.fmt.compare("") == 0)
    _usage();
  if (opt.umi < 0) 
    opt.umi = 0;
  transform(opt.pm5.begin(), opt.pm5.end(), opt.pm5.begin(), ::toupper);
  transform(opt.pm3.begin(), opt.pm3.end(), opt.pm3.begin(), ::toupper);
}

void _info(string des)
{
  string str = "[%Y-%m-%d %X] " + des;
  time_t t = time(0);
  char tmp[64];
  strftime (tmp, sizeof(tmp), str.c_str(), localtime(&t));
  puts(tmp);
}

string _revc(string &str, int flag)
{
  // reverse and complementary
  map<char, char> base;
  base['A'] = 'T'; base['T'] = 'A'; base['C'] = 'G'; base['G'] = 'C'; base['U'] = 'A';
  string rev(str);
  reverse(rev.begin(), rev.end());
  if (flag == 1) {
    for (int i = 0; i < str.size(); i++)
      if (base.find(rev[i]) != base.end())
        rev[i] = base[rev[i]];
  }
  return rev;
}

