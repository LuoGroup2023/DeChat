#include <sys/ioctl.h>
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <ctime>
using namespace std;

#include "utils.hpp"

ProgressManager::ProgressManager(long long totalValue, const char *unitName) : totalValue(totalValue), unitName(unitName) {
  this->lastDisplayElapsedSeconds = -10;
  this->nbDone = 0;
  this->pcDone = 0;
  this->onePc = 0.01 * totalValue;
  struct timespec timeBegin;
  clock_gettime(CLOCK_MONOTONIC, &timeBegin);
  this->timestampBegin = timeBegin.tv_sec;
  for (int i=0; i < 10; i++) {
    lastPercents[i] = 0;
    lastElapsed[i] = 0;
    lastDone[i] = 0;
  }
}

ProgressManager::ProgressManager(long long totalValue) : ProgressManager(totalValue, "u") {}

char * ProgressManager::formatDuration(int nbSec) {
  int time = nbSec;

  int hour = time/3600;
  time = time%3600;
  int min = time/60;
  time = time%60;
  int sec = time;

  sprintf(this->durationString, "%.2d:%.2d:%.2d", hour, min, sec);
}

void ProgressManager::updateProgress(long long nbDone) {
  clock_gettime(CLOCK_MONOTONIC, &timeTmp);
  long long elapsedSeconds = (timeTmp.tv_sec - timestampBegin);
  float pcDone = (long double)nbDone * 100 / (long double)totalValue;

  // we display progress if
  // last progress has not been displayed yet
  // OR it's finished NOW
  // OR if it's been more than 5 seconds and percent has changed
  if (   lastPercents[9] < 100
      && (nbDone >= totalValue
          || (   elapsedSeconds - this->lastDisplayElapsedSeconds >= 1
              && pcDone != this->pcDone)
         )
  ) {
    this->nbDone = nbDone;
    this->pcDone = pcDone;
    this->lastDisplayElapsedSeconds = elapsedSeconds;

    float speed = 0;

    // are we done ?
    if (nbDone < totalValue) {
      // estimate remaining time only if something has already be done
      if (nbDone > 0) {
        long long estimatedRemaining;
        if (pcDone <= 10) {
          // before 10% : estimate with elapsed from beginning
          estimatedRemaining = (float)elapsedSeconds * (float)(totalValue - nbDone) / (float)nbDone;
          if (elapsedSeconds > 0) {
            speed = (float)nbDone / (float)elapsedSeconds;
          }
        }
        else {
          // after 10% : estimate with last 10 measures
          long long timeLast10Measures = elapsedSeconds - lastElapsed[0];
          float pcDiff = pcDone - lastPercents[0];
          estimatedRemaining = (float)timeLast10Measures * (float)(100 - pcDone) / pcDiff;
          long long nbDiff = nbDone - lastDone[0];
          if (timeLast10Measures > 0) {
            speed = (float)nbDiff / (float)timeLast10Measures;
          }
        }
        // estimation precision : +-10
        estimatedRemaining = estimatedRemaining - (estimatedRemaining%10) + 10;

        // generate estimated time string
        formatDuration(estimatedRemaining);
      }
      // nbDone is 0, no estimation
      else {
        sprintf(durationString, "--:--:--");
      }
    }
    else {
      formatDuration(elapsedSeconds);
      if (elapsedSeconds > 0) {
        speed = (long double)nbDone / (long double)elapsedSeconds;
      }
    }

    // manage array of previous elapsed
    for (int i=0; i < 9; i++) {
      lastElapsed[i] = lastElapsed[i+1];
      lastPercents[i] = lastPercents[i+1];
      lastDone[i] = lastDone[i+1];
    }
    lastElapsed[9] = elapsedSeconds;
    lastPercents[9] = pcDone;
    lastDone[9] = nbDone;

    // actually display progress
    this->printProgress(speed);
  }
}

void ProgressManager::printProgress(float speed) {
  int progressBarWidth = getTermWidth() - 60;
  int progressPos = progressBarWidth * this->pcDone / 100;
  cout << "[";
  for (int i = 0; i < progressBarWidth; ++i) {
    if (i < progressPos) cout << "=";
    else if (i == progressPos) cout << ">";
    else cout << " ";
  }
  // are we done ?
  if (this->nbDone < this->totalValue) {
    cout << "] ";
    printf("%.2f", this->pcDone);
    cout << "%";
    cout << " [rmng: " << this->durationString << "]";
    cout << " [spd: ";
    printf("%.2f %s", speed, this->unitName);
    cout << "/s]";
    cout << " " << '\r' << std::flush;
  }
  else {
    cout << "] " << "100% [Done in " << durationString << "] [avg spd: ";
    printf("%.2f %s", speed, this->unitName);
    cout << "/s]\n" << std::flush;
  }
}

int ProgressManager::getTermWidth() {
  struct winsize w;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
  return w.ws_col;
}
