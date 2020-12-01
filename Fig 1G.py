import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

## define functions for calculation
def Normallization(trace, preframe):
    ts = np.array(trace)
    return ts/np.min(ts[:preframe])

def Mean_after(trace, preframe):
    nts = Normallization(trace, preframe)
    return (np.mean(nts[:preframe]), np.mean(nts[preframe:]))

def Mean_after_p(alltrace, preframe):
    from scipy import stats
    N = len(alltrace.columns) # number of samples in the data frame
    mbefore = []
    mafter = []
    for i in range(N):
        ts = alltrace.iloc[:,i]
        mean_before, mean_after = Mean_after(ts, preframe)
        mbefore.append(mean_before)
        mafter.append(mean_after)
        del ts, mean_before, mean_after
    # Mann-Whitney test to see whether there is a difference between groups
    #_, p = stats.mannwhitneyu(mbefore, mafter, use_continuity=True, alternative=None)
	# paired t test to see whether the avg decreases
    p = stats.ttest_rel(mbefore, mafter)
    outcome = {}
    tmp1 = pd.DataFrame( data = mbefore, columns = ['normalized mean (A.U.)'] )
    tmp1['type'] = 'before'
    tmp2 = pd.DataFrame( data = mafter, columns = ['normalized mean (A.U.)'] )
    tmp2['type'] = 'after'
    outcome['raw'] = pd.concat([tmp1, tmp2])
    outcome['p'] = p
    return outcome
        

medium = pd.read_excel('Medium, BMS536924.xlsx', sheet_name = '+medium', header = 1)
bms_plusglu = pd.read_excel('Medium, BMS536924.xlsx', sheet_name = '+glucose, +bms', header = 1)
bms_minusglu = pd.read_excel('Medium, BMS536924.xlsx', sheet_name = '-glucose, +bms', header = 1)


time_glu = np.array(range(241)) # time in seconds

print(Mean_after_p(medium, 120)['p'])
print(Mean_after_p(bms_plusglu, 120)['p'])
print(Mean_after_p(bms_minusglu, 120)['p'])



pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)


plt.figure(1)

plt.subplot(222)
sns.boxplot(x = 'type', y = 'normalized mean (A.U.)', data = Mean_after_p(medium, 120)['raw'],  palette="Greys")
sns.swarmplot(x = 'type', y = 'normalized mean (A.U.)', data = Mean_after_p(medium, 120)['raw'],  color="k")
plt.ylim([0.9,2.2])
plt.xlabel('N = 22')
plt.title('+imaging solution')
plt.ylabel('normalized mean \n intensity (A.U.)', fontsize = 12)
plt.plot([0, 0, 1, 1], [1.7, 1.9, 1.9, 1.7], lw=1.5, c='k')
plt.text(.5, 1.9, "p = 0.76", ha='center', va='bottom', color='k')
plt.yticks([1,2], ['1', '2'])

plt.figure(2)

plt.subplot(132)
sns.boxplot(x = 'type', y = 'normalized mean (A.U.)', data = Mean_after_p(bms_plusglu, 120)['raw'],  palette="Greys")
sns.swarmplot(x = 'type', y = 'normalized mean (A.U.)', data = Mean_after_p(bms_plusglu, 120)['raw'],  color="k")
plt.ylim([0.9,4.7])
plt.title('+ IR inhibitor, \n high glucose')
plt.ylabel('normalized mean intensity (A.U.)', fontsize = 14)
plt.plot([0, 0, 1, 1], [3, 4, 4, 2.5], lw=1.5, c='k')
plt.text(.5, 4, "p = 0.02", ha='center', va='bottom', color='k')
plt.xlabel('N = 12')
plt.yticks([1,2,3,4], ['1', '2', '3', '4'])


plt.subplot(133)
sns.boxplot(x = 'type', y = 'normalized mean (A.U.)', data = Mean_after_p(bms_minusglu, 120)['raw'],  palette="Greys")
sns.swarmplot(x = 'type', y = 'normalized mean (A.U.)', data = Mean_after_p(bms_minusglu, 120)['raw'],  color="k")
plt.ylim([0.9,4.7])
plt.title('+ IR inhibitor, \n low glucose')
plt.plot([0, 0, 1, 1], [2, 2.5, 2.5, 2], lw=1.5, c='k')
plt.text(.5, 2.5, "p = 0.11", ha='center', va='bottom', color='k')
plt.ylabel('')
plt.xlabel('N = 12')
plt.yticks([1,2,3,4], ['1', '2', '3', '4'])


plt.show()



