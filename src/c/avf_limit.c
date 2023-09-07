#include "strdata.h"

void avf_limit(double *avf, int n, double max_limit) {
  
  double avfsum = 0.0;
  
  for (int v = 0; v < n; ++v) {
    if (avf[v] <= AMIN / max_limit) {
      avf[v] = AMIN / max_limit;
    } else if (avf[v] >= (1.0 - (AMIN / max_limit))) {
      avf[v] = 1.0 - (AMIN / max_limit);
    }
    avfsum += avf[v];
  }
  if (avfsum > (1.0 - (AMIN / max_limit))) {
    double avf_dif = avfsum - (1.0 - (AMIN / max_limit));
    for (int j = 0; j < n; ++j) {
      avf[j] -= (avf[j] / avfsum) * avf_dif;
    }
  }
  avf[n] = 1.0 - avfsum;
}