// build: g++ -O2 -std=c++17 -I ../../src bench_compiled.cxx -o bench_compiled
#include "analyzer_pvnew.h"
#include <chrono>
#include <cstdio>
#include <fstream>
using namespace FCCAnalyses::AlephPVNew;
int main(int argc, char** argv) {
  // input: flat binary written by the python driver: nev, then per event
  // [n, n*5 par doubles, n*21 cov doubles]
  std::ifstream f(argv[1], std::ios::binary);
  int nev; f.read((char*)&nev, 4);
  std::vector<detail::TrackSet> evs;
  for (int e = 0; e < nev; ++e) {
    int n; f.read((char*)&n, 4);
    std::vector<double> par(5*n), cov(21*n);
    f.read((char*)par.data(), 8*5*n);
    f.read((char*)cov.data(), 8*21*n);
    detail::TrackSet ts;
    for (int i = 0; i < n; ++i)
      detail::add_track(ts, par[5*i], par[5*i+1], par[5*i+2], par[5*i+3], par[5*i+4], cov.data()+21*i);
    evs.push_back(ts);
  }
  BeamSpot bs; FitConfig cfg;
  double s = 0;
  auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < 3; ++r)
    for (auto& ts : evs) { auto res = detail::fit_core(ts, &bs, nullptr, cfg); s += res.chi2; }
  auto t1 = std::chrono::steady_clock::now();
  printf("compiled fit_core: %.3f ms/event, checksum %.6g\n",
         std::chrono::duration<double, std::milli>(t1-t0).count()/(evs.size()*3.0), s);
  double s2 = 0;
  t0 = std::chrono::steady_clock::now();
  for (auto& ts : evs) { auto r2 = detail::select_core(ts, &bs, 5.0, cfg, 2, false); s2 += r2.kept.size(); }
  t1 = std::chrono::steady_clock::now();
  printf("compiled select_core: %.3f ms/event, kept-sum %.0f\n",
         std::chrono::duration<double, std::milli>(t1-t0).count()/evs.size(), s2);
  return 0;
}
