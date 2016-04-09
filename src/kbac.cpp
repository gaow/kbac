//: src:kbac.cpp
// kbac test
// Copyright 2011 Gao Wang
#include "kbac.hpp"
#include "ranker.hpp"
#include "gsl_cdf.h"
#include "gsl_randist.h"
namespace {
  const double AFFECTED = 1.0, UNAFFECTED = 0.0, HOMO_ALLELE = 2.0, 
        MINOR_ALLELE = 1.0, MAJOR_ALLELE = 0.0, MISSING_ALLELE = -9.0;
}

KbacTest::KbacTest (int* nn, int* qq, double* aa, double* mafUpper,
		    double* xdatIn, double* ydatIn, double* mafIn, int* xcol, int* ylen) 
{

  this->m_nPermutations = *nn;
  this->m_quiet = *qq;
  this->m_alpha = *aa;
  this->m_mafUpper = *mafUpper;
  this->m_mafLower = 0.0;
  
  if (m_alpha >= 1.0) {
    m_adaptive = 0;
  }
  else {
    m_adaptive = 5000;
  }

  for (unsigned int i = 0; i != *ylen; ++i) { 
    m_ydat.push_back(ydatIn[i]);
  }
  for (unsigned int i = 0; i != *xcol; ++i) { 
    m_observedMafs.push_back(mafIn[i]);
  }

  int start = 0;
  int nInvalidGenotypes = 0;
  for (unsigned int i = 0; i != *ylen; ++i) {
    std::vector<double> tmpvct(0);
    for (unsigned int j = 0; j != *xcol; ++j) {
      if (xdatIn[start+j] != MAJOR_ALLELE && xdatIn[start+j] != MINOR_ALLELE && xdatIn[start+j] != HOMO_ALLELE) {
        ++nInvalidGenotypes;
        tmpvct.push_back(MAJOR_ALLELE);
      }
      else {
        tmpvct.push_back(xdatIn[start+j]);    
      }
    }
    m_xdat.push_back(tmpvct);
    start += *xcol;
  }
     
  if (nInvalidGenotypes > 0) {
    std::cout << "**Warning: Invalid genotype codings detected on " << nInvalidGenotypes << " variant sites (codings must be 0 = wild-type, 1 = heterozygous, or 2 = homozygous. KBAC will treat them as 0 = wild-type" << std::endl;
  }

  if (m_quiet == false) {
    std::cout << "Read data Y length " << m_ydat.size() << " and X dimension " << m_xdat.size() << " by " << m_xdat[0].size() << std::endl;
  }

  bool isInputOk = (start == *xcol * *ylen && m_xdat.size() != 0 && m_ydat.size() == m_xdat.size() 
      && m_observedMafs.size() == m_xdat[0].size() && m_mafLower >= 0.0 
      && m_mafUpper <= 1.0 && m_mafUpper > m_mafLower && m_alpha > 0.0);

  if (isInputOk);
  else {
    std::cerr << "**Error** Input data problem in KbacTest::KbacTest(). Now Quit." << std::endl;
    exit(-1);
  }
}


KbacTest::~KbacTest() {};


void KbacTest::calcKbacP(double* pvalue, int* sided)
{
  /*! * the KBAC Statistic: sum of genotype pattern frequencies differences, weighted by hypergeometric kernel. <br>
   * * It is a permutation based two-sided test.  <br>
   * * See <em> Liu DJ 2010 PLoS Genet. </em> <br><br>
   *Implementation:
   */
  if (*sided == 0 || *sided > 2)
    *sided = 1;

  //!- trim data by mafs upper-lower bounds
  m_trimXdat();


  unsigned int sampleSize = m_ydat.size();
  // sample size
  unsigned int regionLen = m_xdat[0].size();
  // candidate region length
  unsigned int nCases = 0; 
  // case size

  for(unsigned int i = 0; i !=  m_ydat.size(); ++i) {
    if(m_ydat[i] == AFFECTED) 
      ++nCases;
  }
  unsigned int nCtrls = sampleSize - nCases;

  if (m_quiet == false) {
    std::cout << "Sample size: " << sampleSize << std::endl;
    std::cout << "Number of cases: " << nCases << std::endl;
    std::cout << "Candidate region length: " << regionLen << std::endl;
  }

  std::vector<double> genotypeId(sampleSize);   
  //!-Compute unique genotype patterns (string) as ID scores (double) 
  bool hasWt = false;
  for (unsigned int i = 0; i != sampleSize; ++i) {

    double vntIdL = 0.0; 
    double vntIdR = 0.0;
    const double ixiix= pow(9.0, 10.0);
    unsigned int lastCnt = 0;
    unsigned int tmpCnt = 0;

    for (unsigned int j = 0; j != regionLen; ++j) { 

      if (m_xdat[i][j] != MISSING_ALLELE && m_xdat[i][j] != MAJOR_ALLELE) 
        vntIdR += pow(3.0, 1.0 * (j - lastCnt)) * m_xdat[i][j];
      else 
        continue;
      if (vntIdR >= ixiix) {
        vntIdL = vntIdL + 1.0;
        vntIdR = vntIdR - ixiix;
        lastCnt = lastCnt + tmpCnt + 1;
        tmpCnt = 0;
        continue;
      }
      else { 
        ++tmpCnt; 
        continue; 
      }
    }

    // one-to-one "ID number" for a genotype pattern
    genotypeId[i] = vntIdL + vntIdR * 1e-10;
    if (!(genotypeId[i] > MAJOR_ALLELE) && !hasWt)
      hasWt = true;
  }
  rank(genotypeId, genotypeId, "default");

  // remove wildtype and get unique genotype patterns 
  //!- Unique genotype patterns that occur in the sample
  std::vector<double> uniquePattern = genotypeId;
  std::sort(uniquePattern.begin(), uniquePattern.end());
  std::vector<double>::iterator it = std::unique(uniquePattern.begin(), uniquePattern.end());
  uniquePattern.resize(it - uniquePattern.begin()); 
  if (hasWt) uniquePattern.erase(uniquePattern.begin());
  if (uniquePattern.size() == 0) {
    std::cout << "**Warning** non-wildtype genotype data is empty. KBAC has nothing to work on. Return p-value 1.0" << std::endl;
    *pvalue = 1.0;
    return;
  }

  // count number of sample individuals for each genotype pattern
  unsigned int uniquePatternCounts[uniquePattern.size()];
  for (unsigned int u = 0; u != uniquePattern.size(); ++u) 
    uniquePatternCounts[u] = 0;

  for (unsigned int i = 0; i != sampleSize; ++i) {
    // for each sample, identify/count its genotype pattern

    for (unsigned int u = 0; u != uniquePattern.size(); ++u) {

      if (genotypeId[i] == uniquePattern[u]) {
        // genotype pattern identified
        ++uniquePatternCounts[u];
        // count this genotype pattern
        break;
      }
      else;
      // genotype pattern not found -- move on to next pattern
    }
  }

  if (m_quiet == false) {
    std::cout << "\nGenotype pattern loading ranks: " << std::endl;
    std::cout << genotypeId << std::endl;
    std::cout << "\nNumber of each unique individual genotype patterns (totaling " << uniquePattern.size() << " patterns excluding wildtype): " << std::endl;
    for (unsigned int u = 0; u != uniquePattern.size(); ++u) std::cout << uniquePatternCounts[u] << ", "; 
    std::cout << std::endl;
  }

  unsigned int iPermutation = 0;
  unsigned int permcount1 = 0, permcount2 = 0;
  double observedStatistic = 0.0;
  while (iPermutation <= m_nPermutations) {

    // the KBAC statistic. Will be of length 1 or 2
    std::vector<double> kbacStatistics(0);
    // two models
    for (unsigned int s = 0; s != *sided; ++s) {

      //!- count number of sample cases (for the 1st model, or ctrls for the 2nd model) for each genotype pattern
      unsigned int uniquePatternCountsSub[uniquePattern.size()];
      for (unsigned int u = 0; u != uniquePattern.size(); ++u) 
        uniquePatternCountsSub[u] = 0;
      // genotype pattern counts in cases (for the 1st model, or ctrls for the 2nd model) 

      for (unsigned int i = 0; i != sampleSize; ++i) {
        if ( m_ydat[i] == (AFFECTED - 1.0 * s) ) {
          // for each "case (for the 1st model, or ctrls for 2nd model)", identify/count its genotype pattern
          for (unsigned int u = 0; u != uniquePattern.size(); ++u) {
            if (genotypeId[i] == uniquePattern[u]) {
              // genotype pattern identified in cases (for the 1st model, or ctrls for 2nd model)
              ++uniquePatternCountsSub[u];
              // count this genotype pattern
              break;
            }
            else;
            // genotype pattern not found -- move on to next pattern
          }
        }
        else;
      }

      //!- KBAC weights
      std::vector<double> uniquePatternWeights(uniquePattern.size());
      // genotype pattern weights, the hypergeometric distribution cmf
      for (unsigned int u = 0; u != uniquePattern.size(); ++u) 
        uniquePatternWeights[u] = 0.0;

      for (unsigned int u = 0; u != uniquePattern.size(); ++u) {
        if (s == 0) 
          uniquePatternWeights[u] = gsl_cdf_hypergeometric_P(uniquePatternCountsSub[u], uniquePatternCounts[u], sampleSize - uniquePatternCounts[u], nCases);
          //uniquePatternWeights[u] = gw_hypergeometric_cmf(uniquePatternCountsSub[u], uniquePatternCounts[u], sampleSize - uniquePatternCounts[u], nCases);
        else
          uniquePatternWeights[u] = gsl_cdf_hypergeometric_P(uniquePatternCountsSub[u], uniquePatternCounts[u], sampleSize - uniquePatternCounts[u], nCtrls);
          //uniquePatternWeights[u] = gw_hypergeometric_cmf(uniquePatternCountsSub[u], uniquePatternCounts[u], sampleSize - uniquePatternCounts[u], nCtrls);
      }

      if (m_quiet == false && iPermutation == 0) {
	std::cout << "\nUnique genotype patterns weights (model " << s+1 << "):" << std::endl;
	std::cout << uniquePatternWeights << std::endl;
      }

      //!- KBAC statistic: sum of genotype pattern frequencies differences in cases vs. controls, weighted by the hypergeometric distribution kernel
      double kbac = 0.0;
      for (unsigned int u = 0; u != uniquePattern.size(); ++u) { 
        if (s == 0)
          kbac = kbac + ( (1.0 * uniquePatternCountsSub[u]) / (1.0 * nCases) - (1.0 * (uniquePatternCounts[u] - uniquePatternCountsSub[u])) / (1.0 * nCtrls) ) *  uniquePatternWeights[u];
        else
          kbac = kbac + ( (1.0 * uniquePatternCountsSub[u]) / (1.0 * nCtrls) - (1.0 * (uniquePatternCounts[u] - uniquePatternCountsSub[u])) / (1.0 * nCases) ) *  uniquePatternWeights[u];
      }

      //std::cout << kbac << std::endl;

      //FIXME
      //gw_round(kbac, 0.0001);
      kbacStatistics.push_back(kbac);
    }

    double statistic = 0.0;
    //!- one model statistic
    if (kbacStatistics.size() == 1) {
      statistic = kbacStatistics[0];
    }
    //!- two model statistic
    else if (kbacStatistics.size() == 2) {
      statistic = fmax(kbacStatistics[0], kbacStatistics[1]);
    }
    else {
      std::cerr << "**Error KBAC statistic (Error code -5)" << std::endl;
      exit(-1);
    }

    if (iPermutation == 0) 
      observedStatistic = statistic;
    else {
      if (statistic >= observedStatistic) 
        ++permcount1;
      if (statistic <= observedStatistic)
        ++permcount2;
      if (m_adaptive != 0)
        *pvalue = m_checkAdaptivePvalue(permcount1, permcount2, iPermutation, m_adaptive, 0);
    }
    if (*pvalue <= 1.0)
      break;
    //!- Permutation
    random_shuffle(m_ydat.begin(), m_ydat.end());
    ++iPermutation;
  }

  if (*pvalue <= 1.0);
  else {
    *pvalue = (1.0 * permcount1 + 1.0) / (1.0 * m_nPermutations + 1.0);
  }
  return;
}


double KbacTest::m_checkAdaptivePvalue(unsigned int permcount1, unsigned int permcount2, unsigned int currentIdx, unsigned int checkPoint, unsigned int alternative) const
{
  if (currentIdx % checkPoint == 0 && checkPoint > 5) {
    //!- adaptive p-value calculation, at an interval of #checkPoint permutations 
    // apply the "six-sigma" rule

    double adaptivePvalue = 1.0;
    if (alternative == 1 || alternative == 0) { 
      adaptivePvalue = (1.0 * permcount1 + 1.0) / (1.0 * currentIdx + 1.0);
    }
    else {
      double permcount = (permcount1 < permcount2) ? permcount1 : permcount2;
      adaptivePvalue = (2.0 * permcount + 1.0) / (1.0 * currentIdx + 1.0);
    }

    double sd = sqrt(adaptivePvalue * (1.0 - adaptivePvalue) / (1.0 * currentIdx));
    double sixsigma = adaptivePvalue - 6.0 * sd;

    if (m_quiet == false) {
      std::cout << "nPerm" << "\t" << "Current.P" << "\t" << "St.Error" << "\t" << "6-sigma.lower.bound" << std::endl;
      std::cout << currentIdx << "\t" << adaptivePvalue << "\t" << sd << "\t" << sixsigma << std::endl;
    }

    if (sixsigma > m_alpha) 
      return adaptivePvalue;
    else
      return 9.0;
  }
  else 
    return 9.0;
}
