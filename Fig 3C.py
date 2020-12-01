#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


## define functions for calculation
def Normallization(trace, preframe):
    ts = np.array(trace)
    return ts/np.min(ts[:preframe])

## read in data 
# glu_bms12 = pd.read_excel('12min_response.xlsx', sheet_name = 'BMS', header = 1)

# a trial to lump data together, modified 7/8/20
glu_bms12_1 = pd.read_excel('12min_response.xlsx', sheet_name = 'BMS', header = 1)
glu_bms12_2 = pd.read_excel('12min_response.xlsx', sheet_name = '+cilium, +glucose & bms', header = 1)
glu_bms12_3 = pd.read_excel('12min_response.xlsx', sheet_name = '-cilium, +glucose & bms', header = 1)
glu_bms12 = pd.concat([glu_bms12_1, glu_bms12_2, glu_bms12_3], axis = 1)

glu_dmso12 = pd.read_excel('12min_response.xlsx', sheet_name = 'DMSO', header = 1)


# In[3]:


# normalized calcium concentration over smaller time points

def Mean_over_time(trace, preframe, frameinterval):
    # the frame interval is the number of frames that is considered in analyzing mean over a period of time
    num_interval = int(np.floor(len(trace)/ frameinterval)) # number of the time intervals the trace can be devided into
    nts = np.array(Normallization(trace, preframe))
    mean_overtime = []
    for j in range(num_interval):
        startt = j * frameinterval
        endt = startt + frameinterval
        meanpertime = np.mean(nts[startt:endt])
        mean_overtime.append(meanpertime)
        del startt, endt, meanpertime
    # print(num_interval)
    return np.array(mean_overtime)

def Mean_overtime(alltrace, preframe, frameinterval):
    from scipy import stats
    N = len(alltrace.columns) # number of samples in the data frame
    M = int( np.floor(len(alltrace.iloc[:,0])/ frameinterval) )
    mtime = np.zeros((N, M))
    for i in range(N):
        ts = alltrace.iloc[:,i]
        mtime[i, :] = Mean_over_time(ts, preframe, frameinterval)
        del ts
    return mtime


def Plot_avg_mean_overtime(alltrace, preframe, frameinterval, plotFigure = True):
    mtime = Mean_overtime(alltrace, preframe, frameinterval)
    tmp_mean = np.mean(mtime, axis = 0)
    tmp_sd = np.std(mtime, axis = 0)/ np.sqrt(len(mtime))
    
    time = np.array(range(mtime.shape[-1])) * frameinterval + frameinterval/ 2
    
    if plotFigure:
    
        fig = plt.figure(figsize = (10, 5))
        ax = fig.add_subplot(111)
        ax.errorbar(time, tmp_mean, tmp_sd, linestyle = '')
        ax.bar(time, tmp_mean, width = 8)
        ax.axvline(x=120, color = 'k', linestyle = ':')
        plt.xlim([0, 240])
        plt.xlabel('time (s)')
        plt.ylabel('normalized R-GECO intensity (A.U.)')
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["font.size"] = "16"
        
    return time, tmp_mean, tmp_sd


# In[13]:


time, bms_avg, bms_sd = Plot_avg_mean_overtime(glu_bms12, 120, 60, False)

_, dmso_avg, dmso_sd = Plot_avg_mean_overtime(glu_dmso12, 120, 60, False)


plt.figure(figsize = (5, 4))

plt.bar(time-10, bms_avg, width = 20, color = 'r', alpha = 0.5, label = "+BMS536924 (n = 72)")
plt.errorbar(time-10, bms_avg, bms_sd, linestyle = '', color = 'r')


plt.bar(np.array(range(119, 721, 1)), height = 0.1, width = 1, bottom = 2.95, color = 'k')
plt.text(130, 2.8, '+Glucose')
plt.ylim([0, 3])
# plt.text(200, 2.3, '*', fontsize=16) # add star to label statistical significance


plt.bar(time+10, dmso_avg, width = 20, color = 'g', alpha = 0.5, label = "+DMSO (n = 52)")
plt.errorbar(time+10, dmso_avg, dmso_sd, linestyle = '', color = 'g')
plt.xlim([0, 721])
plt.xlabel('time (min)')
plt.ylabel('normalized mean intensity (A.U.)')
#plt.rcParams["font.family"] = "Times New Roman"
#plt.rcParams["font.size"] = "16"
plt.xticks([0, 120, 240, 480, 720], [-2, 0, 2, 6, 10]) # convert time into minute



plt.legend(loc = 'upper right')

plt.show()


# In[17]:


# run statistics side-by-side for each minute

from scipy import stats

i = 1

for i in range(1, 12):
    print(i)
    print(stats.ttest_ind_from_stats(dmso_avg[i], dmso_sd[i]*np.sqrt(57), 52, bms_avg[i], bms_sd[i]*np.sqrt(72), 72, equal_var=False))
    
