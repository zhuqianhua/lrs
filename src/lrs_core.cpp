#include "lrs_core.h"

#define link 6

void _loc(string &prm, char typ)
{
  Locate tag;
  for (int i = 0; i < prm.size(); i ++) { 
    if (prm.size() - i < opt.cln)
      continue;
    string anc = prm.substr(i, opt.cln);
    size_t loc = -1;
    while (1) {
      loc = seq.find(anc, loc + 1);
      if (loc == seq.npos)
        break;
      tag.beg = loc - i;
      tag.end = loc + prm.size() - i - 1;
      if (tag.beg < 0)
        tag.beg = 0;
      if (tag.end > seq.length() - 1)
        tag.end = seq.length() - 1;
      tag.typ = typ;
      if (tag.end - tag.beg + 1 < opt.len)
        continue;
      if (opt.gtail == true) {
        int bg = 0, ba = 0;
        if (typ == 'F') {
          for (int z = tag.end - 4; z <= tag.end; z++)
            if (seq[z] == 'G')
              bg ++;
          for (int z = tag.end + 20; z < tag.end + 25; z++) {
            if (z >= seq.size() - 1) 
              break;
            if (seq[z] == 'T')
              ba ++;
          }
        } else if (typ == 'f') {
          for (int z = tag.beg; z <= tag.beg + 4; z++)
            if (seq[z] == 'C')
              bg ++;
          for (int z = tag.beg - 25; z < tag.beg - 20; z++) {
            if (z < 0)
              break;
            if (seq[z] == 'A')
              ba ++;
          }
        } else if (typ == 'T') {
          for (int z = tag.end + 1; z <= tag.end + 5; z++)
            if (seq[z] == 'G') 
              bg ++;
          for (int z = tag.end + 20; z < tag.end + 25; z++) {
            if (z >= seq.size() - 1);
              break;
            if (seq[z] == 'T')
              ba ++;
          }
        } else if (typ == 't') {
          for (int z = tag.beg - 5; z < tag.beg; z ++)
            if (seq[z] == 'C')
              bg ++;
          for (int z = tag.beg - 25; z < tag.beg - 20; z++) {
            if (z < 0)
              break;
            if (seq[z] == 'A')
              ba ++;
          }
        }
        if (((typ == 'F' || typ == 'f') && bg >= 3 && ba < 4) || ((typ == 'T' || typ == 't') && bg < 3 && ba >= 4))
          pis[tag] = ' ';
      } else {
        pis[tag] = ' ';
      }
    }
  }
}

void _locate()
{
  loc.clear();
  string fr = _revc(opt.pm5, 1), tr = _revc(opt.pm3, 1);
  _loc(opt.pm5, 'F');
  _loc(fr, 'f');
  _loc(opt.pm3, 'T');
  _loc(tr, 't');
  int i = 0;
  for (map<Locate, char>::iterator im = pis.begin(); im != pis.end(); im++) {
    //cout << im->first.beg << " " << im->first.end << " " << im->first.typ << endl;
    loc[i % opt.thd].push_back(im->first);
    i ++;
  }
  pis.clear();
}

void *_segmentation(void *thread) 
{
  int THD = *(int *) thread;
  if (loc.find(THD) == loc.end())
    return 0;
  for (vector<Locate>::iterator iv = loc[THD].begin(); iv != loc[THD].end(); iv++) {
    string _read = seq.substr(iv->beg, iv->end - iv->beg + 1);
    int32_t mask = _read.size() / 2;
    // Declares a default Aligner
    StripedSmithWaterman::Aligner aligner;
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result
    StripedSmithWaterman::Alignment align;
    string _primer;
    if (iv->typ == 'F') {
      _primer = opt.pm5;
    } else if (iv->typ == 'T') {
      _primer = opt.pm3;
    } else if (iv->typ == 'f') {
      _primer = _revc(opt.pm5, 1);
    } else if (iv->typ == 't') {
      _primer = _revc(opt.pm3, 1);
    } else {
      cerr << "illegal primer type" << endl;
      exit(-1);
    }
    aligner.Align(_read.c_str(), _primer.c_str(), _primer.size(), filter, &align, mask);
    string dgt;
    int match = 0, clip = 0;
    for (int i = 0; i < align.cigar_string.size(); i++) {
      if (isdigit(align.cigar_string[i]))
        dgt += align.cigar_string[i];
      else {
        if (align.cigar_string[i] == '=')
          match += atoi(dgt.c_str());
        dgt = "";
      }
    }
    if (match < opt.len)
      continue;
    pthread_mutex_lock(&mutex); 
    if (anc.find(align.ref_begin + iv->beg) == anc.end()) {
      Anchor tanc;
      tanc.typ = iv->typ;
      tanc.score = align.sw_score;
      anc[iv->beg] = tanc; //anc[align.ref_begin + iv->beg] = tanc;
    }
    pthread_mutex_unlock(&mutex);
  }
}

int _adjust()
{
  /* anc pos_type: 0-start, F f T t  
  pos: fl_nfl,begin,end,reverse_complementary */
  if (anc.size() == 0)
    return 0;
  map <char, int> mlen;
  mlen['F'] = mlen['f'] = opt.pm5.size();
  mlen['T'] = mlen['t'] = opt.pm3.size();
  map <int, char> mseg;
  struct Linker {
    int beg, score;
    char typ;
  };
  Linker M, F, T, f, t;
  M.score = F.score = T.score = f.score = t.score = 0;
  f.beg = t.beg = F.beg = T.beg = -1;
  char last = NULL;
  for (map<int, Anchor>::iterator it = anc.begin(); it != anc.end(); it ++) {
    //cout << "anc: " << it->first << " " << it->second.typ << " " << it->second.score << endl;
    map<int, Anchor>::iterator in = it;
    ++ in;
    if (it->second.score > M.score) {
      M.beg = it->first;
      M.score = it->second.score;
      M.typ = it->second.typ;
    }
    if (it->second.typ == 'F' && mseg.size() == 0 && it->second.score > F.score) {
      F.beg = it->first; 
      F.score = it->second.score; 
      F.typ = it->second.typ;
    } else if (it->second.typ == 'T' && mseg.size() == 0 && it->second.score > T.score) {
      T.beg = it->first; 
      T.score = it->second.score; 
      T.typ = it->second.typ;
    } else if (it->second.typ == 'f' && it->second.score > f.score) {
      f.beg = it->first; 
      f.score = it->second.score; 
      f.typ = it->second.typ;
    } else if (it->second.typ == 't' && it->second.score > t.score) {
      t.beg = it->first; 
      t.score = it->second.score; 
      t.typ = it->second.typ;
    }
    if (in == anc.end() || in->first - it->first - mlen[it->second.typ] > link || in->first - it->first < mlen[it->second.typ])
      continue;
    if ((it->second.typ=='t' && (in->second.typ=='F'||in->second.typ=='T')) || (it->second.typ=='f' && in->second.typ=='T')) {
      mseg[it->first] = it->second.typ;
      mseg[in->first] = in->second.typ;
      last = in->second.typ;
      f.score = t.score = 0;
      it ++;
    }
  }
  anc.clear();  
  /* pos: typ, beg, end, rev
          typ = [0: fl, 1: pm5 only, 2: pm3 only, 3: pm5 and pm3, but not fl]
          beg, end: 0-start
          rev = [0: fwd, 1: rev] */
  int _typ = 0, _beg = 0, _end = seq.size() - 1, _rev = 0;
  if (mseg.size() == 0) {
    if (M.typ == 'F') {
      _typ = 1;
      _beg = M.beg + mlen[M.typ];
      if (t.beg > M.beg) {
        _end = t.beg - 1;
        _typ = 0;
      }
    } else if (M.typ == 'f') {
      _typ = 1;
      _end = M.beg - 1;
      _rev = 1;
      if (T.beg + mlen['T'] < M.beg && T.beg != -1) {
        _beg = T.beg + mlen['T'];
        _typ = 0;
      }
    } else if (M.typ == 'T') {
      _typ = 2;
      _beg = M.beg + mlen[M.typ];
      _rev = 1;
      if (f.beg > M.beg) {
        _end = f.beg - 1;
        _typ = 0;
      }
    } else {
      _typ = 2;
      _end = M.beg - 1;
      if (F.beg != -1 && F.beg + mlen['F'] < M.beg) {
        _beg = F.beg + mlen['F'];
        _typ = 0;
      }
    }
    if (_end > _beg) {
      pos.push_back(_typ); 
      pos.push_back(_beg); 
      pos.push_back(_end); 
      pos.push_back(_rev);
    }
  } else {
    if (mseg.size() % 2 != 0) {
      cerr << "something wrong with the mseg, pls debug..." << endl;
      exit(0);
    }
    for (map<int, char>::iterator im = mseg.begin(); im != mseg.end(); im ++) {
      map<int, char>::iterator ib = im, ia = im;
      ib --; ia ++;
      if (ib == mseg.end()) {
        if (im->second == 't') {
          _rev = 0; 
          _end = im->first - 1;
          if (F.beg != -1 && F.beg + mlen['F'] < im->first) {
            _beg = F.beg + mlen['F'];
            _typ = 0;
          } else {
            _beg = 0; 
            _typ = 2;
          }
        } else if (im->second == 'f') {
          _rev = 1; 
          _end = im->first - 1;
          if (T.beg != -1 && T.beg + mlen['T'] < im->first) {
            _beg = T.beg + mlen['T']; 
            _typ = 0;
          } else {
            _beg = 0; 
            _typ = 1;
          }
        }
        if (_end <= _beg)
          continue;
        pos.push_back(_typ); 
        pos.push_back(_beg); 
        pos.push_back(_end); 
        pos.push_back(_rev);
      } else if (ia == mseg.end()) {
        if (im->second == 'F') {
          _rev = 0; 
          _beg = im->first + mlen['F'];
          if (t.beg > im->first) {
            _end = t.beg - 1; 
            _typ = 0;
          } else {
            _end = seq.size() - 1; 
            _typ = 1;
          }
        } else if (im->second == 'T') {
          _rev = 1; _beg = im->first + mlen['T'];
          if (f.beg > im->first) {
            _end = f.beg - 1; 
            _typ = 0;
          } else {
            _end = seq.size() - 1;
            _typ = 2;
          }
        }
        if (_end <= _beg)
          continue;
        pos.push_back(_typ); 
        pos.push_back(_beg); 
        pos.push_back(_end); 
        pos.push_back(_rev);
      } else {
        _beg = im->first + mlen[im->second];
        _end = ia->first - 1;
        if (im->second == 'F') {
          _rev = 0;
          if (ia->second == 't')
            _typ = 0;
          else 
            _typ = 3;
        } else if (im->second == 'T') {
          _rev = 1;
          if (ia->second = 'f') 
            _typ = 0;
          else 
            _typ = 3;
        }
        im++;
        if (_end <= _beg)
          continue;
        pos.push_back(_typ); 
        pos.push_back(_beg); 
        pos.push_back(_end); 
        pos.push_back(_rev);
      }
    }
  }
}

void _segm()
{
  _locate();
  int THREAD[opt.thd];
  pthread_t thread[opt.thd];
  for (int i = 0; i < opt.thd; i++) {
    THREAD[i] = i;
    int res = pthread_create(&thread[i], NULL, _segmentation, (void *) &THREAD[i]);
    if (res != 0) {
      printf("Thread create %d failed\n", i);
      exit(-1);
    }
  }
  while (1)
  {
    int mak = 0;
    void *status;
    for (int i = 0; i < opt.thd; i++)
      if (pthread_join(thread[i], &status) != 0)
        mak ++;
    if (mak == 0)
      break;
  }
  loc.clear();
  _adjust();
}
