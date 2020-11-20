#include "lrs.h"

#define polya 0.6
#define gtail 0.8
#define glen 10
#define rtTtail 6
#define NELE(x) (sizeof(x)/sizeof(x[0]))

static const char init_header[] = "@HD\tVN:1.5\tSO:unknown\n";
htsFile *ofl, *onfl;
ofstream ofs, opa;
bam_hdr_t *header;
int _obam(string &nam, bam1_t *r);

template <class T>
string _itos(T n)
{
  stringstream ss;
  ss << n;
  return ss.str();
}

void ibam()
{
  samFile *ibam = sam_open(opt.inp.c_str(), "r");
  header = sam_hdr_read(ibam);
  // write header
  if (sam_hdr_write(ofl, header) < 0 || sam_hdr_write(onfl, header) < 0) {
    cerr << "header writing error: " << endl;
    exit(-1);
  }
  // read bam
  bam1_t *r = bam_init1(); int mline = 0;
  while (sam_read1(ibam, header, r) > 0)
  {
    char *r_qname = bam_get_qname(r);
    uint8_t *r_seq = bam_get_seq(r);
    seq = "";
    for (int i = 0; i < r->core.l_qseq; i++)
      seq += seq_nt16_str[bam_seqi(r_seq, i)];
    pos.clear(); //cout << r_qname << endl; 
    _segm();
    string qname = r_qname;
    if (pos.size() > 0)
      _obam(qname, r);
    mline ++;
    //if (mline >= 1)
      //break;
  }
  sam_close(ibam);
}

int _fasta(string &str)
{
  int flag = 0;
  if (str.find_first_of('>') != string::npos) {
    int end = str.length();
    if (str.find_first_of(' ') != string::npos)
      end = str.find_first_of(' ');
    str = str.substr(str.find_first_of('>') + 1, end - str.find_first_of('>') - 1);
    return 0;
  } else {
    transform(str.begin(), str.end(), str.begin(), ::toupper);
    return 1;
  }
}

void _fas(string &nam)
{
  pos.clear(); 
  _segm();
  bam1_t *r;
  if (pos.size() > 0)
    _obam(nam, r);
}

int _obam(string &nam, bam1_t *r)
{
  /* pos: typ, beg, end, rev
          typ = [0: fl, 1: pm5 only, 2: pm3 only, 3: pm5 and pm3, but not fl]
          beg, end: 0-start
          rev = [0: fwd, 1: rev]
  */
  uint8_t *r_qual = bam_get_qual(r);
  /*uint8_t *np = bam_aux_get(r, "np"); //np:i
  uint8_t *rq = bam_aux_get(r, "rq"); //rq:f
  uint8_t *zm = bam_aux_get(r, "zm"); //zm:i
  uint8_t *sn = bam_aux_get(r, "sn"); //sn:B
  uint8_t *RG = bam_aux_get(r, "RG"); //RG:Z
  double asn[bam_auxB_len(sn)];
  for (int tt = 0; tt <= 4; tt ++) 
    asn[tt] = bam_auxB2f(sn, tt);
  char * trg = bam_aux2Z(RG); */
  int stat[4] = {0, 0, 0, 0};
  for (int i = 0; i <= pos.size() - 4; i += 4) { //cout << pos[i] << " " << pos[i+1] << " " << pos[i+2] << " " << pos[i+3] << endl;
    if (pos[i+2] < pos[i+1]) {
      cerr << "something wrong with pos[] in " << nam << " with seq length: " << seq.size()
            << " at: " << pos[i+1] << " - " << pos[i+2] << endl;
      exit(0);
    }
    // write bam
    bam1_t *n = bam_init1();
    if (pos[i+2] - pos[i+1] < 20)
      continue;
    int send = pos[i+2];
    string qname, umi = "_U", sseq = seq.substr(pos[i+1], pos[i+2] - pos[i+1] + 1);
    if (pos[i+3] == 1)
      sseq = _revc(sseq, 1); 
    // umi
    umi += sseq.substr(sseq.size() - opt.umi, opt.umi);
    sseq = sseq.substr(0, sseq.size() - opt.umi);
    send = send - opt.umi;
    // detect gtail for full polyA 12G
    size_t end = sseq.size() - 1;
    int pgg = end;
    while (1) {
      size_t pgs = sseq.find_last_of('G', end);
      size_t npg = sseq.find_last_not_of('G', pgs);
      if (pgs == string::npos || npg == string::npos)
        break;
      if (end - npg > 2) {
        float rate = static_cast<float>(pgs-npg) / static_cast<float>(end-npg);
        if (rate < gtail)
          break;
      }
      pgg = npg + 1;
      end = npg;
    }
    if (sseq.size() - pgg >= glen) {
      send = send - (sseq.size() - pgg);
      sseq = sseq.substr(0, pgg);
    } 
    // detect polyA tail
    end = sseq.size() - 1;
    int ppa = end;
    while (1) {
      size_t pas = sseq.find_last_of('A', end);
      size_t npa = sseq.find_last_not_of('A', pas);
      if (pas == string::npos || npa == string::npos)
        break;
      float rate = static_cast<float>(pas-npa) / static_cast<float>(end-npa);
      if (rate < polya) 
        break;
      ppa = npa + 1;
      end = npa;
    }
    string spolya;
    if (sseq.size() - ppa >= rtTtail) {
      spolya = sseq.substr(ppa, sseq.size() - ppa);
      send = send - spolya.size();
      sseq = sseq.substr(0, ppa);
    }

    // output bam: qname = Sp_Ep_Pp_Us
    int p_nam = nam.find_last_of('/');
    qname = nam.substr(0, p_nam+1) + "L" + _itos(sseq.size()) + umi + "_ccs";
    if (spolya.size() >= rtTtail)
      opa << qname << "\t" << spolya << endl;

	n->core.tid = -1;
    n->core.pos = -1;
    n->core.mtid = -1;
    n->core.mpos = -1;
    n->l_data = qname.size() + 1 + (int)(1.5 * sseq.size() + (sseq.size() % 2 != 0));
    if (n->m_data < n->l_data) {
      n->m_data = n->l_data;
      kroundup32(n->m_data);
      n->data = (uint8_t*)realloc(n->data, n->m_data);
    }
    n->core.flag = BAM_FMUNMAP;
    n->core.l_qname = qname.size() + 1;
    n->core.l_qseq = sseq.size();
    n->core.n_cigar = 0;
    memcpy(n->data, (void *)qname.c_str(), qname.size() + 1);
    uint8_t *n_seq = bam_get_seq(n);
    uint8_t *n_qual = bam_get_qual(n);
    for (int j = pos[i+1]; j <= send; j ++) {
      bam_set_seqi(n_seq, j - pos[i+1], seq_nt16_table[sseq[j-pos[i+1]]]);
      if (opt.fmt.compare("b") == 0)
        n_qual[j - pos[i+1]] = r_qual[j];
      else 
        n_qual[j - pos[i+1]] = 72;
    }
    /*// add aux
    if (rq != 0)
      bam_aux_update_float(n, "rq", bam_aux2f(rq));
    if (zm != 0)
      bam_aux_update_int(n, "zm", bam_aux2i(zm));
    if (np != 0)
      bam_aux_update_int(n, "np", bam_aux2i(np));
    if (RG != 0)
      bam_aux_update_str(n, "RG", strlen(trg) + 1, trg);
    if (sn != 0)
      bam_aux_update_array(n, "sn", 'f', 4, asn);*/

    if (pos[i] == 0) {
      sam_write1(ofl, header, n);
      stat[0] ++;
    } else {
      sam_write1(onfl, header, n);
      if (pos[i] == 1) 
        stat[1] ++;
      else if (pos[i] == 2)
        stat[2] ++;
      else
        stat[3] ++;
    }
    bam_destroy1(n);
  }
  ofs << stat[0] << "\t" << stat[1] << "\t" << stat[2] << "\t" << stat[3] << endl;
}

void ifasta()
{
  //write header
  header = bam_hdr_init();
  header->l_text = strlen(init_header);
  header->text = strdup(init_header);
  header->n_targets = 0;
  if (sam_hdr_write(ofl, header) < 0 || sam_hdr_write(onfl, header) < 0) {
    cerr << "header writing error: " << endl;
    exit(-1);
  } 

  ifstream ifs(opt.inp.c_str());
  string str(""), nam("");
  while (!ifs.eof())
  {
    getline(ifs, str);
    int flag = _fasta(str);
    if (flag == 0) {
      if (nam.compare("") != 0)
        _fas(nam);
      nam = str;
      seq = "";
    } else {
      seq += str;
    }
  }
  ifs.close();
  _fas(nam);
}

int main(int argc, char* argv[])
{
  _getopt(argc, argv);
  string fl = opt.out + "_fl.bam", nfl = opt.out + "_nfl.bam", fls = opt.out + "_stat.xls", pa = opt.out + "_polya.xls";
  ofl = hts_open(fl.c_str(), "wb");
  onfl = hts_open(nfl.c_str(), "wb");
  ofs.open(fls.c_str());
  ofs << "FL\tOnly primer5\tOnly primer3\tChimeric\n";
  opa.open(pa.c_str());
  //opa << "ID\tSeq\n";
  if (opt.fmt.compare("b") == 0)
    ibam();
  else if (opt.fmt.compare("f") == 0)
    ifasta();
  else { 
    cerr << "illegal parameter for opt.fmt, pls check it" << endl; 
    exit(-1);
  }
  sam_close(ofl);
  sam_close(onfl);
  ofs.close();
  opa.close();
  return 0;
}

