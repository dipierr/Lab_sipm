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

def gaus(x,a,x0,sigma):
	return a*exp(-(x-x0)**2/(2*sigma**2))

def read_trace_length(inputfile):
	with open(inputfile, "r") as infile:
		for i,line in enumerate(infile):
			for j,dum in enumerate(line.split()):
				if (i == 0) and (j == 2):
					trace_length=dum
			if i == 7:
				break
	infile.close()	
	return int(trace_length)
	
def read_event(inputfile,event_n):
	trace_length = read_trace_length(inputfile)
	trace = []
	with open(inputfile, "r") as infile:
		for i,line in enumerate(infile):
			if (i >= event_n*(trace_length+7)+7) and (i < (event_n+1)*(trace_length+7)):
				#trace.append(float(line.split()[0]))
				trace.append(float(line))
				#print(line,float(line.split()[0]))
			elif (i < event_n*(trace_length+7)):
				continue
			elif (i >= (event_n+1)*(trace_length+7)):
				break						
	infile.close()
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
	y, binx, patches = plt.hist(trace, bins = int(np.amax(trace)-np.amin(trace))+1, range = (np.amin(trace)-0.5, np.amax(trace)+0.5), normed=0, facecolor='green', alpha=0.75)
	
	x = []
	for i in range(0,binx.size-1):
		x.append( binx[i]+(binx[i+1] - binx[i])/2)#bin center			
	# second, fit the histogram with the optimize.curve_fit method (coincident bounds = fixed parameter not allowed!)
	popt,pcov = curve_fit(gaus,x,y,p0=[bl[1],bl[0],1],bounds=([-np.inf,bl[0],-np.inf],[np.inf,bl[0]+0.000001,np.inf]))
	
	#plt.clf()
	plt.plot(x,gaus(x,*popt),'r:',label='fit')
	plt.legend()
	plt.show()
		
	return bl[0], popt[2]

def find_peaks(trace,display):
	
	bl,bl_s = find_trace_baseline(trace)
	print(bl[0],bl_s)
	
	#peak begins when trace > 3 sigma and ends when trace<peak_max/2, minimum 3 time slices.
	peaks = [[2]]#first index is peak_bin, second peak_amplitude	
	#trace = trace - bl
	peak_num=0
	peaks[peak_num]=[0,0]
	peak_ready = True
	for i in range(0,trace.size):
		print(i,trace[i],peak_ready,peak_num)
		if (trace[i] > bl[0] + 2*bl_s):
			if (peak_ready):
				istart = i
				peaks[peak_num] = [i,trace[i]]
				peak_ready=False
				print("Trace above threshold",i,trace[i])
			if (trace[i] > peaks[peak_num][1]):
				peaks[peak_num] = [i,trace[i]]
		if (trace[i] < (peaks[peak_num][1]-bl[0])/2+bl[0]) and (i-istart > 3) and (not(peak_ready)) and (trace[i] < trace[i-1]) and (trace[i] < trace[i-1]):
			print("Peak found ",peaks[peak_num])
			peak_ready = True
			peak_num=peak_num+1
			peaks.append([0,0])
	
	peaks_np = np.asarray(peaks)
	peaks_np = np.delete(peaks_np,len(peaks_np)-1,0) 	
	print(peaks_np)
	
	if display:
		xpeaks = peaks_np[:,0]
		ypeaks = peaks_np[:,1]
		plt.figure(num=None, figsize=(24, 6), dpi=80, facecolor='w', edgecolor='k')
		plt.plot(trace, 'black')
		plt.plot(xpeaks,ypeaks,'ro')
		plt.ylabel("ADC counts")		
		plt.show()	
			
	return peaks_np

	
def read_all_events(inputfile):
	trace_length = read_trace_length(inputfile)

	
def main(**kwargs):
	for i in range(1,5):
		trace = read_event(kwargs['input_file'],i)
		peaks = find_peaks(trace,kwargs['display'])

if __name__ == '__main__':
	args = PARSER.parse_args()
	main(**args.__dict__)

