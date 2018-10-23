// ReadFileFindPeaks.cpp

#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "ReadDRS4.h"
#include "FindPeaks.h"

using namespace std;


int main(int argc, const char * argv[]){

    string filename;
    int last_event_n;
    double thr_to_find_peaks;

    bool display = false;

    if (argc > 3){
        filename = argv[1];
        last_event_n = atoi(argv[2]);
        thr_to_find_peaks = atof(argv[3]);
    }
    else {
       printf("Usage: ReadFileFindPeaks <filename> <last_event_n> <thr_to_find_peaks>\n");
       return 0;
    }

    ReadDRS4(filename, last_event_n, thr_to_find_peaks, display);

    return 0;


}
