import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
from shapely.wkt  import dumps, loads
from matplotlib.patches import Ellipse
import pyregion


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
#        xvals = (x-self.fits.header['CRPIX1']*3.)*self.fits.header['CDELT1']+self.fits.header['CRVAL1']
#        yvals = (y-self.fits.header['CRPIX2']*3.)*self.fits.header['CDELT2']+self.fits.header['CRVAL2']
#        xvals = (x-self.fits.header['CRPIX1P']-self.fits.header['CRVAL1P'])*self.fits.header['CDELT1P']
#        yvals = (y-self.fits.header['CRPIX2P']-self.fits.header['CRVAL2P'])*self.fits.header['CDELT2P']
#        xvals, yvals = self.convert_pxy(xvals,yvals)
#        xvals, yvals = self.fix_project(xvals,yvals)
#from http://cxc.harvard.edu/ciao/threads/regions/
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

       self.fig, self.ax = plt.subplots(figsize=figsize,dpi=dpi)
       self.ax.imshow(np.log10(self.fits.data),cmap=cmap,extent=self.extent,
                      vmin=np.log10(self.vmin),vmax=np.log10(self.vmax),origin=self.origin)
       if starreg:
           self.plot_starreg()
       if backreg:
           self.plot_backreg()

       self.save_fig()
#parse ellipse
    def parse_ellipse(self,i):
        fstr = i.replace('Ellipse','').replace('(','').replace(')','').split(',')
        x,y = self.convert_exy(float(fstr[0]),float(fstr[1]))
        w,h = float(fstr[2])*self.fits.header['CDELT1'],float(fstr[3])*self.fits.header['CDELT2']
        a   = float(fstr[4])
        return x,y,w,h,a

#plot source regions
    def plot_starreg(self):
        regs = np.loadtxt(self.sdir+'sources/wav.src.reg',dtype='str')
        regs = regs[2:3]
        for i in regs:
#            try:
#                x,y,w,h,a = self.parse_ellipse(i)
                print 'HERE'
                print i
                r = pyregion.parse('fk5;'+i).as_imagecoord(self.fits.header)
                
#                print x,y,w,h,a
                ell = Ellipse(xy=(x,y),width=w,height=h,angle=a,color='green',linewidth=1.0,fill=None,zorder=50)
                self.ax.add_patch(ell)
#            except ValueError:
#                print 'Failed for now'
       

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
