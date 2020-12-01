import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statistics

## define functions for calculation
def Normallization(trace, preframe):
    ts = np.array(trace)
    return ts/statistics.median(ts[:preframe])
    #return ts/np.min(ts[:preframe])
   


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
        plt.xlabel('time (s)')
        plt.ylabel('normalized R-GECO intensity (A.U.)')
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["font.size"] = "16"
        
    return time, tmp_mean, tmp_sd




## read in +glucose data 
glu_ciliated = pd.read_excel('12min_response.xlsx', sheet_name = '+cilium, +glucose', header = 1)

glu_nocilia = pd.read_excel('12min_response.xlsx', sheet_name = '-cilium, +glucose', header = 1)

glu_pooled = pd.concat([glu_ciliated, glu_nocilia], axis = 1)

print(len(glu_pooled))

## read in + glucose & insuln data
gluins_ciliated = pd.read_excel('12min_response.xlsx', sheet_name = '+cilium, +glucose & insulin', header = 1)

gluins_nocilia = pd.read_excel('12min_response.xlsx', sheet_name = '-cilium, +glucose & insulin', header = 1)

gluins_pooled = pd.concat([gluins_ciliated, gluins_nocilia], axis = 1)



# calculate avg normalized calcium change by minute

time, glu_avg, glu_sd = Plot_avg_mean_overtime(glu_pooled, 120, 60, False)

_, gluins_avg, gluins_sd = Plot_avg_mean_overtime(gluins_pooled, 120, 60, False)


# see the statistical differences between the first 2 min and the rest

from scipy import stats

for i in range(0, 12):
    print(i)
    print(stats.ttest_ind_from_stats(glu_avg[i], glu_sd[i]*np.sqrt(75), 75, gluins_avg[i], gluins_sd[i]*np.sqrt(43), 43, equal_var=False))


plt.figure(figsize = (6, 5))

plt.bar(time-10, glu_avg, width = 20, color = 'dimgray', alpha = 0.5, label = "+glucose (n = 75)")
plt.errorbar(time-10, glu_avg, glu_sd, linestyle = '', color = 'dimgray')

plt.bar(time+10, gluins_avg, width = 20, color = 'lightgray', alpha = 0.5, label = "+glucose & insulin (n = 43)", hatch = '\\')
plt.errorbar(time+10, gluins_avg, gluins_sd, linestyle = '', color = 'lightgray')

# * for p < 0.001
plt.text( 4*60-40, 2.2, '*', fontsize = 20 )
plt.text( 5*60-40, 2.3, '*', fontsize = 20 )
plt.text( 6*60-40, 2.3, '*', fontsize = 20 )
plt.text( 7*60-40, 2.2, '*', fontsize = 20 )
plt.text( 8*60-40, 2.1, '*', fontsize = 20 )
plt.text( 9*60-40, 2.1, '*', fontsize = 20 )
plt.text( 10*60-40, 2, '*', fontsize = 20 )
plt.text( 11*60-40, 2, '*', fontsize = 20 )
plt.text( 12*60-40, 2, '*', fontsize = 20 )


plt.bar(np.array(range(119, 721, 1)), height = 0.1, width = 1, bottom = 2.96, color = 'k')
plt.text(130, 2.8, '+stimulus')
plt.ylim([0, 3])
plt.xlim([0, 721])
plt.xlabel('time (min)')
plt.ylabel('normalized mean intensity (A.U.)')
plt.xticks([0, 120, 240, 480, 720], [-2, 0, 2, 6, 10]) # convert time into minute
plt.legend(loc = 'upper right')
plt.show()






