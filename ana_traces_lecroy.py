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
from matplotlib.colors import LogNorm
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
PARSER.add_argument('-f', '--input_file', type=str, required=False, default='wave0.txt', help='ascii file produced by Lecroy osccilloscope.')
PARSER.add_argument('-flist', '--input_filelist', type=str, required=False, default='', help='list of file names to be analyzed')
PARSER.add_argument('-d', '--display', action='store_true', required=False, default=False, help='display the selected event trace')
PARSER.add_argument('-1plot', '--superimpose', action='store_true', required=False, default=False, help='display all events on the same figure')
PARSER.add_argument('-all', '--all_events', action='store_true', required=False, default=False, help='Analyze all events in the file (max=100000)')
PARSER.add_argument('-first', '--first_event', type=int, required=False, default=0, help='First event number to analyze')
PARSER.add_argument('-last', '--last_event',type=int, required=False, default=10, help='Last event number to analyze')
PARSER.add_argument('-miny', '--miny', type=float, required=False, default=-0.002, help='Amplitude signal scale (y), minimum')
PARSER.add_argument('-maxy', '--maxy', type=float, required=False, default= 0.010, help='Amplitude signal scale (y), maximum')
#PARSER.add_argument('-eve', '--events', required=False, default=[0,2], help='First and last event numbers')
PARSER.add_argument('-mintp', '--min_time_peak', type=float, required=False, default=0., help='Starting time interval [s] for peak search')
PARSER.add_argument('-maxtp', '--max_time_peak', type=float, required=False, default=3.0e-9, help='Starting time interval [s] for peak search')

def read_next_event_lecroy(infile):
	trace = []
	trace.append([])
	trace.append([])
	trace_length = 0
	for i,line in enumerate(infile):
		if line == '':
			break

		for j,dum in enumerate(line.split()):		
			if (i == 0) and (j == 2):
				trace_length=int(dum)			
			if (i < 6) and (j == 4):
				trace[0].append(float(dum))
			if (i < 6) and (j == 5):				
				trace[1].append(float(dum))			
			if (i >= 6) and (j == 0):
				trace[0].append(float(dum))
			if (i >= 6) and (j == 1):	
				trace[1].append(float(dum))
				
		if (trace_length != 0 and i >= trace_length-1):
			break		
	trace_np = np.array(trace)	
	return trace_np

# with self trigger the peak is right after the trigger at fixed time. Only 1 peak/trace
def find_peak(trace,mintp,maxtp):	
	peak = np.array([0.,0.]) #first index is peak_time, second peak_amplitude	
	for i in range(0,trace[0].size-1):
		if (trace[0][i] > mintp and trace[0][i] < maxtp):
			if(trace[1][i] > peak[1]):
				peak[0] = trace[0][i]
				peak[1] = trace[1][i]			
	return peak
###Visto che trace e' un numpy array, va pythonizzato! prendo il max in un certo intervallo!!!			

def show_trace(trace,peak,miny,maxy):		
	plt.figure(num=None, figsize=(24, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.ylim([miny,maxy])
	plt.plot(trace[0], trace[1], 'black')
	plt.plot(peak[0],peak[1],'ro')
	plt.ylabel("V")
	plt.xlabel("s")
	plt.grid(True)			
	plt.show()			
	
def main(**kwargs):	
	
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
		print("Reading from event ",first_event_n," to event ",last_event_n)
	if kwargs['superimpose']:
		trace_all = []
		trace_all.append([])
		trace_all.append([])
	
	peaks = []
	if kwargs['input_filelist'] != '':
		with open(kwargs['input_filelist'], "r") as infilelist:
			for f in infilelist:
				f = f.rstrip('\n')
				print("Reading file ",f)
				with open(f, "r") as infile:
					for i in range(first_event_n,last_event_n): 
						if (i%500 == 0):
							print("Read event ",i)
						trace = read_next_event_lecroy(infile)
						if len(trace[0]) < 10 :
							print("Reached EOF...exiting")
							break	
						peak = find_peak(trace,kwargs['min_time_peak'],kwargs['max_time_peak'])
						peaks.append(peak[1])	
				infile.close()									
		infilelist.close()
	else:
		with open(kwargs['input_file'], "r") as infile:
			for i in range(first_event_n,last_event_n):
				if (i%500 == 0):
					print("Read event ",i)
				trace = read_next_event_lecroy(infile)
				if len(trace[0]) < 10 :
					print("Reached EOF...exiting")
					break	
				if kwargs['superimpose']:
					for jj in range(0,len(trace[0])-1):
						trace_all[0].append(trace[0][jj])
						trace_all[1].append(trace[1][jj])
				peak = find_peak(trace,kwargs['min_time_peak'],kwargs['max_time_peak'])
				peaks.append(peak[1])
				if(kwargs['display'] and not(kwargs['superimpose'])):
					show_trace(trace,peak,kwargs['miny'],kwargs['maxy'])
				elif(kwargs['display'] and (kwargs['superimpose'])):
					if i == first_event_n :
						plt.figure(num=None, figsize=(24, 6), dpi=80, facecolor='w', edgecolor='k')
						plt.ylim([kwargs['miny'],kwargs['maxy']])
						plt.ylabel("V")
						plt.xlabel("s")
						plt.grid(True)
						#plt.plot(trace[0], trace[1], 'black')
					#else:
						#plt.plot(trace[0], trace[1], 'black')	
					if i == last_event_n -1 :
						plt.hist2d(trace_all[0], trace_all[1], bins=100, norm=LogNorm())
						plt.colorbar()		
						plt.show()									
		infile.close()	
	
	plt.grid(True)
	plt.hist(peaks, bins = 40, range = (0,20e-3), facecolor='g', edgecolor='g',histtype='step', stacked=True, fill=False)
	plt.yscale('log')
	plt.xlabel('V')	
	plt.show()
	
if __name__ == '__main__':
	args = PARSER.parse_args()
	main(**args.__dict__)

