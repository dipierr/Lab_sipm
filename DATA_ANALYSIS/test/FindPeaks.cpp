// FindPeaks.cpp

#include <iostream>
#include <vector>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "FindPeaks.h"

using namespace std;

std::vector <std::vector<double> > trace_DLED{{},{}};
int gap_between_peaks;


//______________________________________________________________________________
void FindPeaks(std::vector <std::vector<double> > &trace, double thr_to_find_peaks) {

    int dleddt = 5;
    gap_between_peaks = 2*dleddt;

    DLED(trace, dleddt);

    FindPeaksRisingFalling(trace_DLED, thr_to_find_peaks, dleddt*10,2,4);

}

//______________________________________________________________________________
void DLED(std::vector <std::vector<double> > &trace, int dleddt){

    int trace_DLED_length = trace[0].size() - dleddt;

    // trace_DLED
    for(int i=0; i<trace[0].size(); i++){
        trace_DLED[0].push_back(trace[0][i + dleddt]);
        trace_DLED[1].push_back(trace[1][i + dleddt]-trace[1][i]);
    }

}


//______________________________________________________________________________
void FindPeaksRisingFalling(std::vector <std::vector<double> > &t, double thr, int max_peak_width, int rising_points, int falling_points){
    int i=1+rising_points;
    int index_old = -1;
    int index_new = -1;
    int peak_start = 0;
    int peak_end = 0;
    int peak_width = max_peak_width;
    double time_delay = 0;
    bool rising = true;
    bool falling = true;
    int half_rising_points = 0;
    int half_falling_points = 0;
    int diff = 0;



    int length = trace_DLED[0].size();

    bool find_rising_edge = true;

    while(i<length-1){ // loop on trace t

        if(find_rising_edge){ // find rising edge
            if(t[1][i] > thr){ // loop > thr

                // find if rising edge or not
                rising = true;
                half_rising_points = (int)(rising_points/2);
                diff = rising_points%2;
                // cout<<endl<<"rising_points = "<<rising_points<<", i = "<<i<<endl;
                for(int j=i-half_rising_points; j<i+half_rising_points+diff; j++){ // rising?
                    // cout<<j<<", ";
                    if( t[1][j] >= t[1][j+1] ){
                        rising = false;
                    }
                } // end rising?

                // getchar();

                if(rising){ // rising edge
                    peak_start = i-1;
                    i++;
                    find_rising_edge = false;
                }
                else{
                    i++;
                } // end rising edge

            }else{
                i++;
            } // end loop > thr
        } // end find rising edge
        else{ // wait until the peak ends

            // find if falling edge or not
            falling = true;
            half_falling_points = (int)(falling_points/2);
            diff = falling_points%2;
            // cout<<endl<<"falling_points = "<<falling_points<<", i = "<<i<<endl;
            for(int j=i-half_falling_points; j<i+half_falling_points+diff; j++){ // falling?
                // cout<<j<<", ";
                if( t[1][j] <= t[1][j+1] ){
                    falling = false;
                }
            } // end falling?

            if(falling){ // falling edge
                peak_end = i+1;
                peak_width = peak_end-peak_start;

                // I have found a falling edge. The following time I will look for a new peak:
                find_rising_edge = true;


                // I have to check if this can be a peak:
                if(peak_width < max_peak_width){ // peak_width < max_peak_width
                    // I have found a new peak

                    //Now I look for the maximum value in that window
                    index_new = find_peak_fix_time(peak_start, peak_end);

                    // PEAK FOUND
                    cout<<t[0][index_new]<<"\t"<<t[1][index_new]<<endl;

                    // now index_new will be index_old, for the following loop cycle
                    index_old = index_new;

                    i+=gap_between_peaks;

                } // end peak_width < max_peak_width
                else{
                    // this is not a peak, too wide in time
                    // cout<<"*******   PEAK TOO WIDE   *******"<<endl;
                    i++;
                }
            } // end falling edge
            else{
                // I have not found a falling edge
                i++;
            }
        } //end wait until the peak ends

    } // end loop on trace t

}

//______________________________________________________________________________
int find_peak_fix_time(int mintp, int maxtp){
    double max_func = -10000;
    int index_func=-1;
    for(int i=mintp; i<maxtp; i++){
        //cout<<mintp<<"\t"<<maxtp<<endl;
        if(trace_DLED[1][i]>max_func){
            max_func=trace_DLED[1][i];
            index_func=i;
        }
    }

    return index_func;
}
