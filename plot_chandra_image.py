import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
from shapely.wkt  import dumps, loads
from matplotlib.patches import Ellipse
import pyregion
from astropy.wcs import WCS
from wcsaxes import WCSAxes


class format_image:

#initialize stuff to format image
    def __init__(self,epoch,vmin=0.8,vmax=100.,scale='log'):

        self.sdir = '{0:4d}/'.format(epoch)
      
        self.fname = self.sdir+'obs{0:4d}_img.fits'.format(epoch)
        self.scale = scale
        self.vmin = vmin
        self.vmax = vmax

        self.fits = pyfits.open(self.fname)[0]

    def convert_pxy(self,x,y):
        xvals = (x-self.fits.header['CRPIX1'])*self.fits.header['CDELT1']+self.fits.header['CRVAL1']
        yvals = (y-self.fits.header['CRPIX2'])*self.fits.header['CDELT2']+self.fits.header['CRVAL2']
        return xvals,yvals

    def convert_exy(self,x,y):
#from http://cxc.harvard.edu/ciao/threads/regions/
#DOES NOT WORK
        xvals = (187.2768)+np.tan((-0.000136667)*(x-(self.fits.header['CRVAL1P'])))
        yvals = (+2.0542 )+np.tan((+0.000136667)*(y-(self.fits.header['CRVAL2P'])))

        return xvals,yvals



#Define new extent
    def plot_boundary(self):
        self.origin = 'lower'
        xs, ys = np.shape(self.fits.data)
        xrang = np.array([0,xs])
        yrang = np.array([0,ys])
        

#covert pixel values to RA, DEC coordinates
        xvals,yvals = self.convert_pxy(xrang,yrang)
    
        self.extent = [xvals[0],xvals[1],yvals[0],yvals[1]]

    
#Plot new figure
    def make_plot(self,figsize=(8,8),dpi=200,cmap=plt.cm.gray,starreg=False,backreg=False):
       self.dpi = dpi
       self.figsize = figsize

       self.plot_boundary()

#       self.fig, self.ax = plt.subplots(figsize=figsize,dpi=dpi)
#Set up WCS coordinate axes
       self.fig = plt.figure()
       wcs = WCS(self.fits.header)
       self.ax = WCSAxes(self.fig,[0.1,0.1,0.8,0.8],wcs=wcs)
       self.fig.add_axes(self.ax)
       self.ax.imshow(np.log10(self.fits.data),cmap=cmap,
 
                      vmin=np.log10(self.vmin),vmax=np.log10(self.vmax),origin=self.origin)
       if starreg:
           self.plot_starreg()
       if backreg:
           self.plot_backreg()

       self.ax.set_xlabel('RA')
       self.ax.set_ylabel('Dec')

       self.save_fig()
#parse ellipse
    def parse_ellipse(self,i):
        fstr = i[0].coord_list
        x,y = self.convert_pxy(float(fstr[0]),float(fstr[1]))
        w,h = float(fstr[2])*self.fits.header['CDELT1'],float(fstr[3])*self.fits.header['CDELT2']
        a   = float(fstr[4])
        return x,y,np.abs(w),np.abs(h),a
#add patches to plot
    def add_patches(self,r):
        patch_list, artist_list = r.get_mpl_patches_texts()
        # ax is a mpl Axes object
        for p in patch_list:
            self.ax.add_patch(p)
        for t in artist_list:
           self.ax.add_artist(t)

#plot source regions
    def plot_starreg(self):
        regs = np.loadtxt(self.sdir+'sources/wav.src.reg',dtype='str')
        for i in regs:
            try:
                r = pyregion.parse(i).as_imagecoord(self.fits.header)
                r[0].attr[1]['color'] = 'green'
 
                self.add_patches(r)

            except ValueError:
                print 'Failed for now'
                i = i.split('&')[0]
                r = pyregion.parse(i).as_imagecoord(self.fits.header)
                r[0].attr[1]['color'] = 'teal'
                self.add_patches(r)
 
#create an exclusion region
            m = i.replace('Ellipse(','').replace(')','').split(',')
            m[2] = float(m[2])*2.
            m[3] = float(m[3])*2.
            n = 'Ellipse({0},{1},{2:6.5f},{3:6.5f},{4})'.format(m[0],m[1],m[2],m[3],m[4])
            r = pyregion.parse(n).as_imagecoord(self.fits.header)
            
            r[0].attr[1]['color'] = 'red'
            self.add_patches(r)
  

    def plot_backreg(self):
        reg = np.loadtxt('ascii_region.reg')
        print 'Not implemented'

#Save figure
    def save_fig(self):
        self.fig.savefig(self.fname.replace('fits','png'),bbox_pad=.1,bbox_inches='tight',dpi=self.dpi)
        self.fig.savefig(self.fname.replace('fits','eps'),bbox_pad=.1,bbox_inches='tight',dpi=self.dpi)

#The image is projected Gnomomically so we need to correct for that to get the circles in the right place
    def fix_project(self,ra,dec):
        a = np.cos(np.radians(dec))*np.cos(np.radians(ra-self.fits.header['CRVAL1']))
        f = 1./self.fits.header['CDELT2']*(180./np.pi)/(np.sin(np.radians(self.fits.header['CRVAL2']))*np.sin(np.radians(dec))+a*np.cos(np.radians(self.fits.header['CRVAL2'])))
        decn = -f*(np.cos(np.radians(self.fits.header['CRVAL2']))*np.sin(np.radians(dec))-a*np.sin(np.radians(self.fits.header['CRVAL2'])))
        ran = -f*np.cos(np.radians(dec))*np.sin(np.radians(ra-self.fits.header['CRVAL1']))

        ran = self.fits.header['CRVAL1']+(ran)*self.fits.header['CDELT1']
        decn = self.fits.header['CRVAL2']-(decn)*self.fits.header['CDELT2']    
        return ran,decn
