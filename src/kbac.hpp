//: src:kbac.hpp
// kbac test
// Copyright 2011 Gao Wang
#ifndef _KBAC_HPP   
#define _KBAC_HPP
#include <iostream>
#include <list>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <vector>
#include <iomanip>

class KbacTest 
{

public:
  KbacTest(int* nn, int* qq, double* aa, double* mafUpper,
	   double* xdatIn, double* ydatIn, double* mafIn, int* xcol, int* ylen);
  ~KbacTest();

  void calcKbacP(double* pvalue, int* sided);

private:
  std::vector< std::vector<double> > m_xdat;
  std::vector<double> m_ydat;
  std::vector<double> m_observedMafs;
  double m_mafLower;
  double m_mafUpper;
  bool m_quiet;
  double m_alpha;
  double m_nPermutations;
  unsigned int m_adaptive;

  double m_checkAdaptivePvalue(unsigned int permcount1, unsigned int permcount2,
			       unsigned int currentIdx, unsigned int checkPoint,
			       unsigned int alternative) const;   
  void m_trimXdat() {
    std::vector< std::vector<double> > xdat = m_xdat;
    m_xdat.clear();
    m_xdat.resize(xdat.size());

    for (unsigned int j = 0; j != m_observedMafs.size(); ++j) {
      if (m_observedMafs[j] <= m_mafLower || m_observedMafs[j] > m_mafUpper) 
	continue;

      else {
	for (unsigned int i = 0; i != xdat.size(); ++i)
	  m_xdat[i].push_back(xdat[i][j]);
      }
    }  
    return;
  }
};

namespace std {

  //!- Dump a vector to screen
  template<class T> ostream & operator<<(ostream & out, const vector<T> & vec)
  {
    if (!vec.empty()) {
      typename vector<T>::const_iterator it = vec.begin();
      out << *it;
      for (++it; it != vec.end(); ++it)
        out << " " << *it ;
    }
    return out;
  }
}
#endif
