from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

def Jurkevich(times, mags, minperiod=10, num_points=200, nbins=10):
	'''Performs a Jurkevich periodicity analysis as given in Jurkevich 1971.
	INPUTS:
	times: 1d array containing times of datapoints in lightcurve. Times do not need to be in ascending order.

	mags: magnitudes associated with respective times.

	minperiod: minimum period to calculate Jurkevich statistic for. 

	num_points: number of datapoints used in sampling the period and Jurkevich statistics.

	nbins: number of bins to split the phased light currve into. 

	OUTPUTS:
	trial_period_range: gives an array of length num_points, between minperiod and maxperiod, of the times used to fold the light curve.

	array_Jurkevich_statistic: the corresponding Jurkevich statistic calculated for each period in trial_period_range.
	'''
	timeinds = times.argsort()
	times = times[timeinds]
	mags = mags[timeinds]


	times = times - times[0] # aking time array start at t=0 to make things easier later on.
	baseline_length = times[-1]
	maxperiod = baseline_length #cannot calculate Jurkevich statistic for a time greater than the length of the light curve.

	array_Jurkevich_statistics = np.zeros(num_points)

	trial_period_range = np.linspace(minperiod, maxperiod, num_points)
	for (p, trialperiod) in enumerate(trial_period_range):
		trial_freq = 1/trialperiod
		phases = times*trial_freq
		phases = phases - np.floor(phases) 

		inds = phases.argsort()
		phases = phases[inds]
		new_mags = mags[inds]

		binned_phases = np.array_split(phases, nbins)
		binned_mags = np.array_split(new_mags, nbins)

		sum_square_dev = np.zeros(nbins)
		
		for (k,maggroup) in enumerate(binned_mags):
			num_datapoints = len(maggroup)
			mean_mag = np.mean(maggroup)
			square_vals = maggroup**2

			if len(maggroup) > 1:
				V_k = np.sum(square_vals) - num_datapoints*mean_mag**2
				#V_k = np.var(maggroup)*num_datapoints #variance of points within a bin times by the number of points in the bin
				sum_square_dev[k] = V_k
			else:
				raise Exception('BIN IS UNDERSAMPLED - Decrease the number of bins!')
				#V_k = num_datapoints 
				#sum_square_dev[k] = V_k #
		print sum_square_dev
		Jurkevich_statistic = np.sum(sum_square_dev)
		array_Jurkevich_statistics[p] = Jurkevich_statistic

	overall_sum_square_dev = np.var(mags)*len(mags)
	array_Jurkevich_statistics = array_Jurkevich_statistics/overall_sum_square_dev
	print array_Jurkevich_statistics

	return (trial_period_range, array_Jurkevich_statistics)
