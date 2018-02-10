'''
Goal: to extract information from SiPM traces digitized by the AGILENT MS06054A.
Steps: 
        to find the peaks (at a fixed time with LED on and external trig and with Self trigger. At a random time with external trigger and LED off, for DCR)
                to evaluate a baseline value relative to each peak (still to be added)
        to measure amplitude of each peak
                to measure the charge relative to of each peak
        optionally display traces and the found peaks (baselines and amplitudes)
                to produce amplitude spectrum
                to produce charge spectrum
                to analyze spectra: gaussian fit on peaks, cross-talk evaluation, gain

TEMPORARY VERSION

'''

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.mlab as mlab
import numpy as np
from numpy import array
from scipy import stats
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import pandas as pd
import argparse
import statistics
import peakutils #Peak detection utilities for 1D data


#--inputfile: file.txt containing the traces

__description__ = 'Producing the azimuth and zenith angle for sim_telarray'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__, formatter_class=formatter)
PARSER.add_argument('-f', '--input_file', type=str, required=False, default='wave0.txt', help='ascii file produced by oscilloscope.')
PARSER.add_argument('-flist', '--input_filelist', type=str, required=False, default='', help='list of file names to be analyzed')
PARSER.add_argument('-d', '--display', action='store_true', required=False, default=False, help='display the selected event trace')
PARSER.add_argument('-1plot', '--superimpose', action='store_true', required=False, default=False, help='display all events on the same figure')
PARSER.add_argument('-all', '--all_events', action='store_true', required=False, default=False, help='Analyze all events in the file (max=100000)')
PARSER.add_argument('-first', '--first_event', type=int, required=False, default=0, help='First event number to analyze')
PARSER.add_argument('-last', '--last_event',type=int, required=False, default=100000, help='Last event number to analyze')
PARSER.add_argument('-miny', '--miny', type=float, required=False, default=-0.01, help='Amplitude signal scale (y), minimum')
PARSER.add_argument('-maxy', '--maxy', type=float, required=False, default= 0.2, help='Amplitude signal scale (y), maximum')
PARSER.add_argument('-mintp', '--min_time_peak', type=float, required=False, default=100, help='Starting time interval [index] for peak search')
PARSER.add_argument('-maxtp', '--max_time_peak', type=float, required=False, default=160, help='Ending time interval [index] for peak search')
PARSER.add_argument('-minoff', '--min_ind_offset', type=float, required=False, default=0, help='Starting index for offset search')
PARSER.add_argument('-maxoff', '--max_ind_offset', type=float, required=False, default=79, help='Ending index for offset search')
PARSER.add_argument('-estoff', '--estimated_offset', type=float, required=False, default=0, help='Estimated offset (for plots)')
PARSER.add_argument('-noise', '--noise_level', type=float, required=False, default=0, help='Noise level')
PARSER.add_argument('-fit', '--fit_hist', type=bool, required=False, default=False, help='Fit histogram?')


max_a = 3
bins_Volt = 180
bins_Time = 79

def gaussian(x, mean, amplitude, standard_deviation):
    return amplitude * np.exp( - ((x - mean) / standard_deviation) ** 2)

def find_offset(trace_selected):
        lift = 100
        base_vect = peakutils.baseline(trace_selected+lift, 1) #if i have a negative baseline the code crashes (maybe)
        base = statistics.mean(base_vect) - lift
        return base

def read_next_event(infile, offset_found, offset, min_index_find_offset,max_index_find_offset):
        trace = []
        trace.append([])
        trace.append([])
        trace_length = 0
        for i,line in enumerate(infile):
                if line == '':
                        break
                        
                for j,dum in enumerate(line.split()):
                        if (i == 1) and (j == 0):
                                trace_length=int(dum)+3
                        if (i >= 3) and (j == 0):
                                trace[0].append(float(dum))
                        if (i >= 3) and (j == 1):       
                                trace[1].append(-float(dum))                             
                if (trace_length != 0 and i >= trace_length-1):
                        break
        
        trace_np = np.array(trace)
        
        if ((trace_np[1]!=[])):# and (offset_found==False)):
                offset_local = find_offset(trace_np[1][min_index_find_offset:max_index_find_offset])
                
        else:
                offset_local = offset
        return trace_np, offset_local

# with self trigger or with external trigger (with light) the peak is right after the trigger at fixed time. Only 1 peak/trace
def find_peak_fixtime(trace,min_ind,max_ind, noise_level):
        first_ind = int(min_ind*0.5)
        #print(first_ind)
        x=trace[0][first_ind:]
        y=trace[1][first_ind:]
        indexes = peakutils.indexes(y, thres=0., min_dist=50)
        #indexes = peakutils.interpolate(x, y, ind=indexes) #Tries to enhance the resolution of the peak detection by using Gaussian fitting, centroid computation or an arbitrary function on the neighborhood of each previously detected peak index
        sel_indexes = []
        found = False
        
        #now I look for the triggered peak, between min_ind and max_ind
        for i in range (0, np.size(indexes)):
            if(y[indexes[i]]>=noise_level):
                sel_indexes.append(indexes[i])
                if((indexes[i]>=min_ind-first_ind) and (indexes[i]<=max_ind-first_ind) and (found==False)):
                    peak = np.array([x[indexes[i]], y[indexes[i]]]) #in the [min_ind, max_ind] iterval I choose only the first peak (which is the triggered one)
                    found = True
                    
        peaks = np.array([x[sel_indexes], y[sel_indexes]]) #all the peaks, not only the triggered one
        if(found==False): #if i don't find a triggered peak)
            #print('ERROR: peak not found; check mintp and maxtp')
            peak = np.array([-100, -100])
                        
        return peak, peaks, found

def show_trace(trace,peak,miny,maxy, points_layout):            
        plt.figure(num=None, figsize=(24, 6), dpi=80, facecolor='w', edgecolor='k')
        #plt.ylim([miny,maxy])
        plt.plot(trace[0], trace[1], 'black')
        plt.plot(peak[0],peak[1], str(points_layout))
        plt.ylabel("V")
        plt.xlabel("s")
        plt.grid(True)
        plt.axvline(x=trace[0][100])
        plt.axvline(x=trace[0][150])
        plt.show()

def fit(bin_centers,bin_heights,bin_borders,fit_start,fit_end):
        found_fit_start=False
        found_fit_end=False
        for i in range (0,bins_Volt):
            if((found_fit_start==False) and (bin_centers[i]>fit_start)):
                fit_start_ind=i
                found_fit_start=True
                #print fit_start_ind
            elif((found_fit_end==False) and (bin_centers[i]>fit_end)):
                fit_end_ind=i
                found_fit_end=True
                #print fit_end_ind
            
        popt, pcov = curve_fit(gaussian, bin_centers[fit_start_ind:fit_end_ind], bin_heights[fit_start_ind:fit_end_ind], p0=[1., 0., 1.])
        
        x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 10000)
        plt.plot(x_interval_for_fit, gaussian(x_interval_for_fit, *popt), 'r',label='fit')
        
        mean = popt[0]
        err_mean = pcov[0][0]
        standard_deviation = popt[2]
        return mean, err_mean, standard_deviation

def main(**kwargs):
        
        min_index_find_offset = kwargs['min_ind_offset']
        max_index_find_offset = kwargs['max_ind_offset']
        estimated_offset = kwargs['estimated_offset']
        noise_level = kwargs['noise_level']
        offset_found=False
        offset=0.
        first_event_n = kwargs['first_event']
        last_event_n = kwargs['last_event']
        if kwargs['input_filelist'] != '':
                kwargs['all_events'] = True
        if kwargs['all_events']:
                if kwargs['display'] == True:
                        kwargs['display'] = False
                        print("Warning: Analyzing all events, --display option is not possible.")
                first_event_n = 0
                last_event_n = 100000   
        if not(kwargs['all_events']):
                print("Reading from event "+str(first_event_n)+" to event "+str(last_event_n))
        if kwargs['superimpose']:
                trace_all = []
                trace_all.append([])
                trace_all.append([])
        
        #if only the peak amplitude is needed, to add one dimension in case also peak time is relavant
        peaks_all = []
        signal_all = []
        offset_all = []
        peak_not_found=0
        if kwargs['input_filelist'] != '':
                cnt=0
                with open(kwargs['input_filelist'], "r") as infilelist:
                        offset_found=False
                        for f in infilelist:
                                f = f.rstrip('\n')
                                cnt=cnt+1
                                if(cnt%128 == 0):
                                        print("Reading file ",f)
                                with open(f, "r") as infile:
                                        for i in range(first_event_n,last_event_n): 
                                                if (i%500 == 0):
                                                        print("Read event "+str(i))
                                                if(i==first_event_n+1):
                                                        offset_found=True
                                                trace, offset = read_next_event(infile, offset_found, offset, min_index_find_offset,max_index_find_offset)
                                                if (len(trace[0]) < 10) :
                                                        print("Reached EOF...exiting")
                                                        break   
                                                peak, peaks, found = find_peak_fixtime(trace,kwargs['min_time_peak'],kwargs['max_time_peak'], noise_level)
                                                if(found==True): #I consider the peaks only if I have found the triggered peak
                                                        peaks_all.append(peak[1])
                                                        signal_all.append(peaks[1][:])
                                                        offset_all.append(offset)
                                                else:
                                                        peak_not_found = peak_not_found+1
                                infile.close()
                infilelist.close()
        else:
                with open(kwargs['input_file'], "r") as infile:
                        for i in range(first_event_n,last_event_n):
                                if (i%500 == 0):
                                        print("Read event "+str(i))
                                if(i==first_event_n+1):
                                        offset_found=True
                                trace, offset = read_next_event(infile, offset_found, offset, min_index_find_offset,max_index_find_offset)
                                if len(trace[0]) < 10 :
                                        print("Reached EOF...exiting")
                                        break   
                                if kwargs['superimpose']:
                                        for jj in range(0,len(trace[0])-1):
                                                trace_all[0].append(trace[0][jj])
                                                trace_all[1].append(trace[1][jj])
                                peak, peaks, found = find_peak_fixtime(trace,kwargs['min_time_peak'],kwargs['max_time_peak'],noise_level)
                                if(found==True):  #I consider the peaks only if I have found the triggered peak
                                        peaks_all.append(peak[1])
                                        signal_all.append(peaks[1][:])
                                        offset_all.append(offset)
                                else:
                                        peak_not_found = peak_not_found+1
                                if(kwargs['display'] and not(kwargs['superimpose'])):
                                        show_trace(trace,peak,kwargs['miny'] + estimated_offset,kwargs['maxy'] + estimated_offset, 'ro') #in order to consider the offset
                                        #show_trace(trace,peaks,kwargs['miny'] + estimated_offset,kwargs['maxy'] + estimated_offset, 'bo')
                                elif(kwargs['display'] and (kwargs['superimpose'])):
                                        if i == first_event_n :
                                                plt.figure(num=None, figsize=(24, 6), dpi=80, facecolor='w', edgecolor='k')
                                                plt.ylim([kwargs['miny']+estimated_offset,kwargs['maxy']+estimated_offset])  #in order to consider the offset
                                                plt.ylabel("V")
                                                plt.xlabel("s")
                                                plt.grid(True)
                                        if i == last_event_n -1 :
                                                plt.hist2d(trace_all[0], trace_all[1], bins=[bins_Time,bins_Volt], norm=LogNorm())
                                                plt.colorbar()          
                                                plt.show()
                                                
                infile.close()  
        
        print('Number of peaks not found = ' + str(peak_not_found))
        
        
        
        plt.figure(num=None, figsize=(24, 6), dpi=80, facecolor='w', edgecolor='k')
        offset_mean = statistics.mean(offset_all)
        print('Offset = ' + str(offset_mean))
        peaks_all_np = np.array(peaks_all)
        peaks_all_np = peaks_all_np - offset_mean
        plt.grid(True)
        #plt.yscale('log')
        plt.xlabel('V')
        #plt.ylim(1,1000)
        
        bin_heights, bin_borders, _ = plt.hist(peaks_all_np, bins = bins_Volt, range = (0, kwargs['maxy']), facecolor='g', edgecolor='g',histtype='step', stacked=True, fill=False)
        bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
        
        if(kwargs['fit_hist']==True):
            if(kwargs['input_file']=='20180209_DARK_34_04.txt' or kwargs['input_file']=='20180209_DARK_34_03.txt'):
                fit_range = [0, 0.02,  0.035,   0.05,   0.065,  0.08,   0.095,   0.11]
            elif(kwargs['input_file']=='20180209_DARK_35_01.txt' or kwargs['input_file']=='20180209_DARK_35_02.txt'):
                fit_range = [0, 0.025,  0.04,   0.06,   0.075,  0.09,   0.11,   0.125] 
            elif(kwargs['input_file']=='20180209_DARK_36_02.txt' or kwargs['input_file']=='20180209_DARK_36_01.txt'):
                fit_range = [0, 0.025,  0.045,  0.065,  0.085,  0.105,  0.125,  0.145]
            else:
                print('ERROR: SPECIFY fit_range')
                quit()
            peaks_num = len(fit_range)-1
            mean = np.linspace(0,0,peaks_num)
            err_mean = np.linspace(0,0,peaks_num)
            standard_deviation = np.linspace(0,0,peaks_num)
            diff_mean = np.linspace(0,0,peaks_num-1)
            sum_err_mean=0
            for i in range(0,peaks_num):
                fit_start = fit_range[i]
                fit_end = fit_range[i+1]
                mean[i], err_mean[i], standard_deviation[i] = fit(bin_centers,bin_heights,bin_borders,fit_start,fit_end)
                sum_err_mean = sum_err_mean + err_mean[i]
                #print(mean[i])
            for i in range(0,peaks_num-1):
                diff_mean[i] = mean[i+1]-mean[i]
            print('Diff mean\t'+str(diff_mean))
            mean_diff_mean = statistics.mean(diff_mean)
            print('Mean diff mean\t'+str(mean_diff_mean))
            print('Err Mean\t'+str(sum_err_mean))
            
            
        plt.show()
        
        
if __name__ == '__main__':
        args = PARSER.parse_args()
        main(**args.__dict__)
