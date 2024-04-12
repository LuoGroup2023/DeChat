#ifndef _TOOL_utils_HPP_
#define _TOOL_utils_HPP_
#include <ctime>

class ProgressManager {
  public:
    ProgressManager(long long totalValue, const char* unitName);
    ProgressManager(long long totalValue);
    void updateProgress(long long nbDone);

  private:
    long long nbDone;
    long long totalValue;
    long long onePc;
    float pcDone;
    long long timestampBegin;
    long long lastDisplayElapsedSeconds;
    long long lastElapsed[10];
    float lastPercents[10];
    float lastDone[10];
    char durationString[20];

    struct timespec timeTmp;
    const char* unitName;

    char* formatDuration(int nbSec);
    int getTermWidth();
    void printProgress(float speed=0);
};

#endif
