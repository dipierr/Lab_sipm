#!/usr/bin/python

'''
Goal: to extract information from SiPM traces digitized by the CAEN Digitizer.
Steps: 
	to find the peaks
	to evaluate a baseline value relative to each peak
	to measure amplitude of each peak
	to measure the charge relative to of each peak
	optionally display traces and the found peaks (baselines and amplitudes)
	to produce amplitude spectrum
	to produce charge spectrum
'''

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import pandas as pd
import argparse


#ottengo le opzioni: 
#--inputfile: file.txt containing the traces

__description__ = 'Producing the azimuth and zenith angle for sim_telarray'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('-f', '--input_file', type=str, required=True, default='wave0.txt', help='ascii file produced by wavedump sw reading the CAEN Digitizer')
PARSER.add_argument('-d', '--display', action='store_true', required=False, default=False, help='display the selected event trace')
PARSER.add_argument('-all', '--all_events', action='store_true', required=False, default=False, help='Analyze all events in the file (max=100000)')
PARSER.add_argument('-first', '--first_event', type=int, required=False, default=0, help='First event number to analyze')
PARSER.add_argument('-last', '--last_event',type=int, required=False, default=10, help='Last event number to analyze')
#PARSER.add_argument('-eve', '--events', required=False, default=[0,2], help='First and last event numbers')

def gaus(x,a,x0,sigma):
	return a*exp(-(x-x0)**2/(2*sigma**2))

def read_trace_length(infile):
	trace_length = 0
	for i,line in enumerate(infile):
		for j,dum in enumerate(line.split()):
			if (i == 0) and (j == 2):
				trace_length=int(dum)
		if i == (6 + trace_length):
			break
	return trace_length
	
def read_event(infile,trace_length,event_n):
	trace = []	
	for i,line in enumerate(infile):
		if (i >= event_n*(trace_length+7)+7) and (i < (event_n+1)*(trace_length+7)-1):
			trace.append(float(line))
		elif (i >= (event_n+1)*(trace_length+7)-1):
			break						
	trace_np = np.array(trace)	
	return trace_np

def read_next_event(infile,trace_length):
	trace = []	
	for i,line in enumerate(infile):
		if line == '':
			break
		if (i >= 7 and i < trace_length+6 ):
			trace.append(float(line))
		elif (i >= trace_length+6 ):
			break						
	trace_np = np.array(trace)	
	return trace_np	

#######
# To be used and called in case the baseline on a single trace is unstable
def find_run_baseline(inputfile,nevents_baseline):
	for n in range(1,nevents_baseline):
		trace = read_event(inputfile,n,False)
		trace_ave = np.sum(trace) / trace.size
		print(trace_ave,np.mean(trace))	
	baseline = np.mean(trace)	
	return baseline
#######
	
def find_trace_baseline(trace):
	bl = stats.mode(trace) # ADC count distribution's mode (is an array [0]=value of the mode,[1] = counts of tyhe mode) as Baseline	
	#I want to fit the histogram of ADC counts with a gaussian	
	#first simple way using the scipy.stats.norm.fit method, but not flexible parameter bounds (I want to fix the mean value of the gaussian)
	#bl_m, bl_s = stats.norm.fit(trace) # Sigma of gaussian fit of the data	
	
	#second way: histogram and curve_fit
	# first, produce the histogram with plt.hist
	#y, binx, patches = plt.hist(trace, bins = int(np.amax(trace)-np.amin(trace))+1, range = (np.amin(trace)-0.5, np.amax(trace)+0.5), normed=0, facecolor='green', alpha=0.75)
	y, binx = np.histogram(trace, bins = int(np.amax(trace)-np.amin(trace))+1, range = (np.amin(trace)-0.5, np.amax(trace)+0.5), normed=0)
	x = []
	for i in range(0,binx.size-1):
		x.append( binx[i]+(binx[i+1] - binx[i])/2)#bin center			
	# second, fit the histogram with the optimize.curve_fit method (coincident bounds = fixed parameter not allowed!)
	popt,pcov = curve_fit(gaus,x,y,p0=[bl[1],bl[0],1],bounds=([-np.inf,bl[0],-np.inf],[np.inf,bl[0]+0.000001,np.inf]))		
	return bl[0], popt[2]

def find_peaks(trace,display):	
	bl,bl_s = find_trace_baseline(trace)
	
	#peak begins when trace > 3 sigma and ends when trace<peak_max/2, minimum 3 time slices.
	peaks = [[2]]#first index is peak_bin, second peak_amplitude	
	#trace = trace - bl
	peak_num=0
	peaks[peak_num]=[0,0]
	peak_ready = True
	for i in range(0,trace.size):
		if (trace[i] > bl[0] + 2*bl_s):
			if (peak_ready):
				istart = i
				peaks[peak_num] = [i,trace[i]]
				peak_ready=False
			if (trace[i] > peaks[peak_num][1]):
				peaks[peak_num] = [i,trace[i]]
		if (trace[i] < (peaks[peak_num][1]-bl[0])/2+bl[0]) and (i-istart > 3) and (not(peak_ready)) and (trace[i] < trace[i-1]) and (trace[i] < trace[i-2]):
			peak_ready = True
			peak_num=peak_num+1
			peaks.append([0,0])
	
	peaks_np = np.asarray(peaks)
	peaks_np = np.delete(peaks_np,len(peaks_np)-1,0) 	
	#print(peaks_np)
	
	if display:
		xpeaks = peaks_np[:,0]
		ypeaks = peaks_np[:,1]
		plt.figure(num=None, figsize=(24, 6), dpi=80, facecolor='w', edgecolor='k')
		plt.plot(trace, 'black')
		plt.plot(xpeaks,ypeaks,'ro')
		plt.ylabel("ADC counts")
		plt.grid(True)			
		plt.show()	
			
	return peaks_np,bl[0]
	
def main(**kwargs):
	
	#first_event_n = 0
	#last_event_n = 1500
	#print(kwargs['events'])
	#first_event_n = kwargs['events'][0]
	#last_event_n = kwargs['events'][1]
	first_event_n = kwargs['first_event']
	last_event_n = kwargs['last_event']
	if kwargs['all_events']:
		kwargs['display'] = False
		print("Warning: Analyzing all events, --display option is not possible.")
		first_event_n = 0
		last_event_n = 100000
	
	print("Reading from event ",first_event_n," to event ",last_event_n)
		
	with open(kwargs['input_file'], "r") as infile:
		trace_length = read_trace_length(infile)
		for i in range(first_event_n,last_event_n):
			if (i%500 == 0):
				print("Read event ",i)
			trace = read_next_event(infile,trace_length)
			if len(trace) < trace_length-1 :
				print("Event ",i-1)
				print("Reached EOF...exiting")
				break				
			peaks,bl = find_peaks(trace,kwargs['display'])
			peaks[:,1] = peaks[:,1] - bl
			if (i == first_event_n):
				hist0,bin_edges0 = np.histogram(peaks[:,1],bins=30,range=(-0.5,29.5))
			else:
				hist,bin_edges = np.histogram(peaks[:,1],bins=30,range=(-0.5,29.5))
				hist0 = hist0 + hist
				#print(hist,hist0)
		bin_cen = []
		bin_half_width = (bin_edges[2] - bin_edges[1])/2.
		for i in range(0,len(bin_edges)-1):
			bin_cen.append(bin_edges[i] + bin_half_width)
			#print(i,bin_cen)
		#print(bin_cen,hist0)
		plt.hist(bin_cen, bins = 30, range = (-0.5,29.5), normed=0, weights = hist0)
		#plt.yscale('log')	
		plt.show()
			
	infile.close()

if __name__ == '__main__':
	args = PARSER.parse_args()
	main(**args.__dict__)

