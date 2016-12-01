import matplotlib as mpl
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
        self.epoch = epoch
      
        self.fname = self.sdir+'obs{0:4d}_img.fits'.format(epoch)
        self.scale = scale
        self.vmin = vmin
        self.vmax = vmax
# format for add exclusion regions in background
        self.outreg = '\n-Ellipse({0:6.2f},{1:6.2f},{2:6.4f},{3:6.4f},{4:3.3f})'

        self.fits = pyfits.open(self.fname)[0]
#convert to physical coordinates
    def convert_pxy(self,x,y):
        xvals = (x-self.fits.header['CRPIX1P'])*self.fits.header['CDELT1P']+self.fits.header['CRVAL1P']
        yvals = (y-self.fits.header['CRPIX2P'])*self.fits.header['CDELT2P']+self.fits.header['CRVAL2P']
        return xvals,yvals

#    def convert_exy(self,x,y):
##from http://cxc.harvard.edu/ciao/threads/regions/
##DOES NOT WORK
#        xvals = (187.2768)+np.tan((-0.000136667)*(x-(self.fits.header['CRVAL1P'])))
#        yvals = (+2.0542 )+np.tan((+0.000136667)*(y-(self.fits.header['CRVAL2P'])))
#
#        return xvals,yvals



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
    def make_plot(self,figsize=(8,8),dpi=200,cmap=plt.cm.gray,starreg=False,backreg=False,writreg=False):
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
       if writreg:
           self.write_reg()

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

#write in registract background regions with sources excluded
    def write_reg(self):
        import matplotlib.path as mplPath
#counter for regions
        f = 1
        for j,i in enumerate(self.backreg):
            try:
                if i.name == 'polygon':
                    outfile = open(self.sdir+'sources/e{0:5d}_background_I{1:1d}.reg'.format(self.epoch,f).replace(' ','0'),'w')
                    c = i.coord_list 
#in physical coordinates
                    xvals,yvals = self.convert_pxy(np.array(c[::2]),np.array(c[1::2]))
                    bbPath = mplPath.Path(np.array([[c[0],c[1]],
                                          [c[2],c[3]],
                                          [c[4],c[5]],
                                          [c[6],c[7]]]))
                    inbox, = np.where(bbPath.contains_points(self.starreg_coor[:,:2]))
#start with the Polygon
                    polyfmt = 'polygon({0:8.4f},{1:8.4f},{2:8.4f},{3:8.4f},{4:8.4f},{5:8.4f},{6:8.4f},{7:8.4f})'
                    outfile.write(polyfmt.format(xvals[0],yvals[0],xvals[1],yvals[1],xvals[2],yvals[2],xvals[3],yvals[3]).replace(' ',''))
#write all exclusion regions to file
                    for k in inbox:
#write out in physical cooridnates
                        ta = self.starreg_phys[k,:]
                        outfile.write(self.outreg.format(ta[0],ta[1],ta[2],ta[3],ta[4]).replace(' ',''))
                    outfile.close()
#increase counter for regions              
                    f += 1
                                          
#find ASCI-I regions
#                print i.attr[1]['text']
#                if i.attr[1]['text'][0] == 'I':
#                    print i.attr[1]['text'],i.coord_list,i.name
            except KeyError:
                continue
    

    def rot_data(self,p,rotation):
        t_start = self.ax.transData
        t = mpl.transforms.Affine2D().rotate_deg(rotation)
        t_end = t_start + t
        p.set_transformation(t_end)

        return p 

#add patches to plot
    def add_patches(self,r,rotation=0):
        patch_list, artist_list = r.get_mpl_patches_texts()

            

        # ax is a mpl Axes object
        for p in patch_list:
    #add rotation for asci region
            if rotation != 0:
                p = self.rot_data(p,rotation)
            self.ax.add_patch(p)
        for t in artist_list:
    #add rotation for asci region text
            if rotation != 0:
                t = self.rot_data(t,rotation)
            self.ax.add_artist(t)

#plot source regions
    def plot_starreg(self):
        regs = np.loadtxt(self.sdir+'sources/wav.src.reg',dtype='str')
        self.starreg = regs
#add to check whether coords are region
#add reject region in image    coordinates
        self.starreg_coor = np.zeros((regs.size,5))-9999.9
#add reject region in physical coordinates
        self.starreg_phys = np.zeros((regs.size,5))-9999.9
        for j,i in enumerate(regs):
            try:
                r = pyregion.parse(i).as_imagecoord(self.fits.header)
                r[0].attr[1]['color'] = 'green'
#output in physical cooridinates
#                i = i.replace('Ellipse(','').replace(')','')
#                p = np.array(i.split(',')).astype('float')
                
                self.add_patches(r)
            except ValueError:
                print 'Failed for now'
                i = i.split('&')[0]
                r = pyregion.parse(i).as_imagecoord(self.fits.header)
                r[0].attr[1]['color'] = 'teal'
#output in physical cooridnates
#                i = i.replace('Ellipse(','').replace(')','')
#                p = np.array(i.split(',')).astype('float')
                self.add_patches(r)
 
#create an exclusion region
            m = i.replace('Ellipse(','').replace(')','').split(',')
            m[2] = float(m[2])*2.
            m[3] = float(m[3])*2.
#create region with double width and height 
            n = 'Ellipse({0},{1},{2:6.5f},{3:6.5f},{4})'.format(m[0],m[1],m[2],m[3],m[4])
#in image coordinates
            r = pyregion.parse(n).as_imagecoord(self.fits.header)
#in image physical coordinates
            p = pyregion.parse(n)
            
            r[0].attr[1]['color'] = 'red'
#add reject regions to plots
            self.add_patches(r)
#add reject region in image coordinates
            self.starreg_coor[j,:] = np.array(r[0].coord_list)
#add reject region in physical coordinates
            self.starreg_phys[j,:] = np.array(p[0].coord_list)
  

    def plot_backreg(self):
#Doesn't work for now
#        asciir = 'ascii_default.reg'
#        r = pyregion.open(asciir).as_imagecoord(self.fits.header)
#        self.add_patches(r,rotation=self.fits.header['ROLL_PNT'])

        asciir = self.sdir+'boxes_e{0:5d}_asci.reg'.format(self.epoch).replace(' ','0')
#in image coordinates
        r = pyregion.open(asciir).as_imagecoord(self.fits.header)
#in image physical coordinates
        p = pyregion.open(asciir)

#in image coordinates
        self.backreg = r
#in image physical coordinates
        self.backimg = p
        self.add_patches(r)

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
