plt.hist(scipy.stats.zscore(abs(this_pearson.flatten())), alpha=.5, label='zscore(abs(Pearson))')
plt.hist(scipy.stats.zscore(this_mi.flatten()), alpha=.5, label='zscore(MI)', bins=25)
plt.xlabel('Pearson r = {:.02f}\np-value = {}'.format(scipy.stats.pearsonr(scipy.stats.zscore(abs(this_pearson.flatten())),\
                                     scipy.stats.zscore(this_mi.flatten()))[0], scipy.stats.pearsonr(scipy.stats.zscore(abs(this_pearson.flatten())),\
                                     scipy.stats.zscore(this_mi.flatten()))[1]))
plt.legend()
plt.show()