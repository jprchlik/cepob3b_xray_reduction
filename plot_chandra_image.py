import pyfits
import matplotlib.pyplot as plt
import numpy as np

class format_image:

#initialize stuff to format image
    def __init__(self,fname,vmin=0.8,vmax=1000.,scale='log'):

        self.fname = fname
        self.scale = scale
        self.vmin = vmin
        self.vmax = vmax

        self.fits = pyfits.open(fname)[0]


#Define new extent
    def plot_boundary(self):
        self.origin = 'lower'
        xs, ys = np.shape(self.fits.data)
        xrang = np.array([0,xs])
        yrang = np.array([0,ys])
        

        xvals = (xrang-self.fits.header['CRPIX1'])*self.fits.header['CDELT1']+self.fits.header['CRVAL1']
        yvals = (yrang-self.fits.header['CRPIX2'])*self.fits.header['CDELT2']+self.fits.header['CRVAL2']
    
        self.extent = [xvals[0],xvals[1],yvals[0],yvals[1]]

    
#Plot new figure
    def make_plot(self,figsize=(8,8),dpi=200,cmap=plt.cm.binary,starreg=False,backreg=False):
       self.dpi = dpi
       self.figsize = figsize

       self.plot_boundary()

       self.fig, self.ax = plt.subplots(figsize=figsize,dpi=dpi)
       self.ax.imshow(self.fits[1].data,cmap=cmap,extent=self.extent,
                      vmin=self.vmin,vmax=self.vmax,origin=self.origin)
       if starreg:
           self.plot_starreg()
       if backreg:
           self.plot_backreg()

       self.save_fig()

#Save figure
    def save_fig(self):
        self.fig.savefig(fname.replace('fits','png'),bbox_pad=.1,bbox_inches='tight',dpi=self.dpi)
        self.fig.savefig(fname.replace('fits','eps'),bbox_pad=.1,bbox_inches='tight',dpi=self.dpi)
