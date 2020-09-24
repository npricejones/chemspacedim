import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
fs = 20

def histogram2d(x,y,bins=100):
    """
    Create a 2D histogram of data represented by the two dimensions x and y

    x:        Array of data values in 'x'
    y:        Array of data values in 'y'
    bins:     Number of bins in which to divide each axis

    """
    # Create histogram
    mask = (np.invert(np.isnan(x))) & (np.invert(np.isnan(y)))
    H,xedges,yedges = np.histogram2d(x[mask],y[mask],bins=bins)
    # Reorient appropriately
    H = np.rot90(H)
    H = np.flipud(H)
    # Mask where bins are empty
    Hmasked = np.ma.masked_where(H==0,H)
    return xedges,yedges,Hmasked

def plothist2d(fig,ax,histdata,clabel=False,vmin=0,vmax=110,norm='lin',
               default_cmap='viridis'):
    """
    Create a 2D histogram of data represented by the two dimensions x and y

    fig:      Figure to plot in
    ax:       Subplot object to plot in
    histdata: containes information on how to plot the histogram
    clabel:   Label for the colourbar - no colourbar is plotted if not given
    vmin:     Minimum value of the histogram
    vmax:     Maximum value of the histogram

    """
    xedges,yedges,Hmasked = histdata
    # Plot histogram
    if norm == 'lin':
        im = ax.pcolormesh(xedges,yedges,Hmasked,
                           cmap = plt.get_cmap(default_cmap),
                           vmin=vmin,vmax=vmax,rasterized=True)
    elif norm == 'log':
        im = ax.pcolormesh(xedges,yedges,Hmasked,
                           cmap = plt.get_cmap(default_cmap),
                           norm=LogNorm(vmin=vmin,vmax=vmax),rasterized=True)
    # If colourbar is desired, plot and label it
    if clabel:
        cbar=fig.colorbar(im,pad = 0)
        cbar.set_label(label=clabel,fontsize=fs)
        cbar.ax.tick_params(labelsize=fs)
    return im
