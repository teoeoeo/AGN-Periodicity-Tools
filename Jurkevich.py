from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

def weighted_mean(values, errors):
	weightings = 1/errors**2
	sum_weightings = sum(weightings)
	weighted_mean = sum(weightings*values)/sum_weightings

	weighted_err  = np.sqrt(1/sum_weightings)

	return weighted_mean, weighted_err

def Jurkevich(times, mags, magerrs, minperiod=100, maxperiod=3000, num_points=200, nbins=5, type='normal'):
	'''Performs a Jurkevich periodicity analysis as given in Jurkevich 1971.
	INPUTS:
	times: 1d array containing times of datapoints in lightcurve. Times do not need to be in ascending order.

	mags: magnitudes associated with respective times.

	magerrs; the errors on the respevtive data points.

	minperiod: minimum period to calculate Jurkevich statistic for. 

	maxperiod: maximum period to calculate Jurkevich statistic for. 
	
	num_points: number of datapoints used in sampling the period and Jurkevich statistics.

	nbins: number of bins to split the phased light currve into. 

	OUTPUTS:
	trial_period_range: gives an array of length num_points, between minperiod and maxperiod, of the times used to fold the light curve.

	array_Jurkevich_statistic: the corresponding Jurkevich statistic calculated for each period in trial_period_range.
	'''
	timeinds = times.argsort()
	times = times[timeinds] 
	mags = mags[timeinds]


	times = times - times[0] # taking time array start at t=0 to make things easier later on.
	baseline_length = times[-1]
	#cannot calculate Jurkevich statistic for a time greater than the length of the light curve.

	array_Jurkevich_statistics = np.zeros(num_points)


	if type == 'normal':

		trial_period_range = np.linspace(minperiod, maxperiod, num_points)
		for (p, trialperiod) in enumerate(trial_period_range):
			#print trialperiod
			trial_freq = 1/trialperiod
			phases = times*trial_freq
			phases = phases - np.floor(phases) 

			inds = phases.argsort()
			phases = phases[inds]
			new_mags = mags[inds]

			#binned_phases = np.array_split(phases, nbins)
			#binned_mags = np.array_split(new_mags, nbins)
			#THIS SHIT IT WRONG


			binned_phases = []
			binned_mags = []
			binned_magerrs = []

			for i in range(nbins):
				phaseinds = np.where((phases > i/nbins) & (phases < (i+1)/nbins))
				thisbinphases = phases[phaseinds]
				thisbinmags = new_mags[phaseinds]

				binned_phases.append(np.array(thisbinphases))
				binned_mags.append(np.array(thisbinmags))

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
					#raise Exception('BIN IS UNDERSAMPLED - Decrease the number of bins!')
					#print str(trialperiod) + 'is a bit undersampled'
					pass
					#V_k = num_datapoints 
					#sum_square_dev[k] = V_k #
			#print sum_square_dev
			Jurkevich_statistic = np.sum(sum_square_dev)
			array_Jurkevich_statistics[p] = Jurkevich_statistic

		overall_sum_square_dev = np.var(mags)*len(mags)
		array_Jurkevich_statistics = array_Jurkevich_statistics/overall_sum_square_dev


	elif type == 'weighted':

		trial_period_range = np.linspace(minperiod, maxperiod, num_points)
		for (p, trialperiod) in enumerate(trial_period_range):
			trial_freq = 1/trialperiod
			phases = times*trial_freq
			phases = phases - np.floor(phases) 

			inds = phases.argsort()
			phases = phases[inds]
			new_mags = mags[inds]
			new_magerrs = magerrs[inds]

			#binned_phases = np.array_split(phases, nbins)
			#binned_mags = np.array_split(new_mags, nbins)
			#binned_magerrs = np.array_split(new_magerrs, nbins) THIS SHIT IT WRONG

			binned_phases = []
			binned_mags = []
			binned_magerrs = []

			for i in range(nbins):
				phaseinds = np.where((phases > i/nbins) & (phases < (i+1)/nbins))
				thisbinphases = phases[phaseinds]
				thisbinmags = new_mags[phaseinds]
				thisbinmagerrs = new_magerrs[phaseinds]

				binned_phases.append(np.array(thisbinphases))
				binned_mags.append(np.array(thisbinmags))
				binned_magerrs.append(np.array(thisbinmagerrs))

			sum_square_dev = np.zeros(nbins)
			
			for (k,maggroup) in enumerate(binned_mags):
				num_datapoints = len(maggroup)
				maggrouperrs = binned_magerrs[k]
				mean_mag, mean_err = weighted_mean(maggroup, maggrouperrs)
				#square_vals = maggroup**2
				resids = maggroup - mean_mag


				if len(maggroup) > 1:
					V_k = num_datapoints * np.sum(resids**2/maggrouperrs**2)/np.sum(1/maggrouperrs**2) 
					#V_k = np.var(maggroup)*num_datapoints #variance of points within a bin times by the number of points in the bin
					sum_square_dev[k] = V_k
				else:
					#raise Exception('BIN IS UNDERSAMPLED - Decrease the number of bins!')
					print str(trialperiod) + 'is a bit undersampled'
					pass
					#V_k = num_datapoints 
					#sum_square_dev[k] = V_k #
			#print sum_square_dev
			Jurkevich_statistic = np.sum(sum_square_dev)
			array_Jurkevich_statistics[p] = Jurkevich_statistic

		overall_mean_mag, overall_mean_err = weighted_mean(mags, magerrs)
		overall_resids = mags - overall_mean_mag
		overall_sum_square_dev = len(mags) * np.sum(overall_resids**2/magerrs**2)/np.sum(1/magerrs**2)
		array_Jurkevich_statistics = array_Jurkevich_statistics/overall_sum_square_dev

	else:
		print 'YOU HAVE NOT PICKED A VALID VERSION OF JURKEVICH'

	return (trial_period_range, array_Jurkevich_statistics)

def Jurk_best_curve(times, mags, magerrs, best_per, minperiod=100, maxperiod=3000, nbins=5):
	'''dont actually need min and max period as variables in the function.'''
	timeinds = times.argsort()
	times = times[timeinds] 
	mags = mags[timeinds]
	magerrs = magerrs[timeinds]

	times = times - times[0] # taking time array start at t=0 to make things easier later on.
	baseline_length = times[-1]

	best_freq = 1/best_per
	phases = times*best_freq
	phases = phases - np.floor(phases)

	inds = phases.argsort()
	phases = phases[inds]
	new_mags = mags[inds]
	new_magerrs = magerrs[inds]

	#binned_phases = np.array_split(phases, nbins)
	#binned_mags = np.array_split(new_mags, nbins)
	#binned_magerrs = np.array_split(new_magerrs, nbins)

	binned_phases = []
	binned_mags = []
	binned_magerrs = []

	for i in range(nbins):
		phaseinds = np.where((phases > i/nbins) & (phases < (i+1)/nbins))
		thisbinphases = phases[phaseinds]
		thisbinmags = new_mags[phaseinds]
		thisbinmagerrs = new_magerrs[phaseinds]

		binned_phases.append(np.array(thisbinphases))
		binned_mags.append(np.array(thisbinmags))
		binned_magerrs.append(np.array(thisbinmagerrs))

	best_mags = np.zeros(nbins)
	best_magerrs = np.zeros(nbins)
	for (i,themags) in enumerate(binned_mags):
		mean_mag, mean_magerr = weighted_mean(themags, binned_magerrs[i]) 
		best_mags[i] = mean_mag
		best_magerrs[i] = mean_magerr

	phased_bin_centres = np.linspace(1/(2*nbins),1-1/(2*nbins),nbins) #Centres of the bins in terms of phase

	time_bin_centres = phased_bin_centres*best_per #converting phases to times
	temp_time_bin_centres = time_bin_centres
	num_periods = int(np.ceil(baseline_length/best_per)) #getting how many periods are in the full length of the light curve (rounded up to nearest int)
	for i in range(num_periods-1):
		time_bin_centres = np.append(time_bin_centres, temp_time_bin_centres+(i+1)*best_per) #repeating the phased light curve through the entire time span of the data (and a bit more due to rounding up of num_periods)
	best_mags = np.tile(best_mags, num_periods)
	best_magerrs = np.tile(best_magerrs, num_periods)

	phased_bin_edgesplus = phased_bin_centres + 0.499 * 1/nbins 
	phased_bin_edgesminus = phased_bin_centres - 0.499* 1/nbins
	phased_bin_edges = np.append(phased_bin_edgesplus, phased_bin_edgesminus)
	phased_bin_edges = np.sort(phased_bin_edges)

	time_bin_edges = phased_bin_edges*best_per
	temp_time_bin_edges = time_bin_edges
	for i in range(num_periods-1):
		time_bin_edges = np.append(time_bin_edges, temp_time_bin_edges+(i+1)*best_per)
	#time_bin_edges = np.tile(time_bin_edges, num_periods)

	best_mags_edges = np.repeat(best_mags,2)
	best_magerrs_edges = np.repeat(best_magerrs,2) #need two copies of each entry in order as need to plot it both to the left and right of the real point.

	return (phased_bin_centres, phased_bin_edges, time_bin_centres, time_bin_edges, best_mags, best_magerrs, best_mags_edges, best_magerrs_edges)


if __name__ == '__main__':

	name = 'PKS0157+011'
	data = np.loadtxt('../Optical_data/Catalina/'+ name + '.txt', delimiter=',', comments='#',skiprows=1)

	MJD = data[:,5]
	mag = data[:,1]
	magerr = data[:,2]

	inds = MJD.argsort()
	sortedMJD = MJD[inds]
	sorted_mag = mag[inds]
	sorted_magerr = magerr[inds]
	
	min_MJD = min(MJD)
	shifted_MJD = sortedMJD - min_MJD

	#shifted_MJD = np.linspace(0,1000)
	#sorted_mag = np.sin(2*np.pi*shifted_MJD/450)

	Jurk_period, Jurk_stat = Jurkevich(shifted_MJD, sorted_mag)

	plt.figure()
	plt.subplot(211)
	plt.plot(shifted_MJD, sorted_mag, 'ko', markersize=2)
	#plt.errorbar(shifted_MJD, sorted_mag, sorted_magerr,capsize=2,linestyle='', color='black')
	plt.xlabel('MJD - ' + str(min_MJD))
	plt.ylabel('Magnitude')

	plt.subplot(212)
	plt.plot(Jurk_period, Jurk_stat)
	plt.xlabel('Period')
	plt.ylabel('$V^2_m$')

	plt.show()