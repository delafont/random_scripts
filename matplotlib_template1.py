from bbcflib.track import track
from bbcflib.gfminer.common import unroll
import numpy
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from os.path import basename


def smooth(x,window_len=100,window='hanning'):
    # From Matplotlib's cookbook: http://www.scipy.org/Cookbook/SignalSmooth
    if x.ndim != 1:         raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len: raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window must be one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s = numpy.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': # moving average
        w = numpy.ones(window_len,'d')
    else:
        w = eval('numpy.'+window+'(window_len)')
    y = numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1]


def plot(chromosomes,region,nbins,format='png'):

    k = 0 # subplot index
    ncols = 2 # number of columns
    nrows = 4 # number of rows
    sub = ncols*nrows # number of subplots per page

    for chr in chromosomes:
        if k%sub == 0: # new figure
            fig = plt.figure(figsize=[10,15])
            fig.suptitle('Suptitle', fontsize=14)

        x = range(region)
        y = [z**2 for z in x]

        ax = fig.add_subplot(nrows,ncols,k%sub+1)
        ax.set_title('Title', ha="left")
        ax.set_xticks(range(region))
        ax.set_xticklabels([str(z) for z in x])
        ax.set_yticks(range(0,100,10))
        ax.set_ylim([0,100])
        ax.set_yscale('log')

        ax.plot(x,y,color='red')
        ax.fill_between(range(region),y,1,where=y>0,color='blue')
        ax.fill_between(range(region),y,1,where=y<0,color='blue')

        w = smooth(y,window_len=region/nbins)
        ax.plot(range(region),w,color='blue',linewidth=2.)

        ## TEST
        if 0:
            fig.savefig('test.png',format='png')
            plt.close(fig)
            break

        if (k+1)%sub==0 or (k+1)==len(chromosomes):
            fig.subplots_adjust(left=0.08, right=0.95, bottom=0.08, top=0.95, hspace=0.3, wspace=0.09)
            name = 'figname'+'.'+format
            fig.savefig(name,format=format)
            plt.close(fig)
        k += 1
