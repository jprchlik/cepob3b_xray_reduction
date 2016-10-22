from astropy.io import ascii

from fancy_plot import *
from sherpa.astro.ui import *
from astropy.table import Table,hstack,join
import pyfits
import matplotlib.pyplot as plt
import numpy as np
import glob
from astropy import units as u
from astropy.coordinates import SkyCoord
import os 
import pyfits

dirs = ['9919','9920','10809','10810','10811','10812']

odir = '/home/jakub/wandajune_home/Ancillary_Projects/tol_Halpha/new_xray/'
mine = 0.5
maxe = 8.0

for i in dirs:
    print i

    try:
        os.mkdir(i+'/cstat_fits')
    except OSError:
        print 'Directory Already Exists'

    bigtab= ascii.read('reduced_xrays_obsid{0:5d}.ipac'.format(int(i)).replace(' ','0'),format='ipac')
    keep, = np.where(bigtab['akT'] > 0.)
    bigtab = bigtab[keep]

    for j in np.arange(1,bigtab['RA'].size+1):
        try:
            fig,ax = plt.subplots()
            load_pha(i+'/srcflux_output/sers_{0:4d}.pi'.format(j).replace(' ','0'))
            load_grouping(i+'/srcflux_output/sers_{0:4d}_group_10.pi'.format(j).replace(' ','0'))
    
            load_arf(i+'/srcflux_output/sers_{0:4d}.arf'.format(j).replace(' ','0'))
            load_rmf(i+'/srcflux_output/sers_{0:4d}.rmf'.format(j).replace(' ','0'))
    #        group_counts(10)
            
            ufit = pyfits.open(i+'/srcflux_output/sers_{0:4d}_cstat_umod.pi'.format(j).replace(' ','0'))
            afit = pyfits.open(i+'/srcflux_output/sers_{0:4d}_cstat_amod.pi'.format(j).replace(' ','0'))
            data = get_data()
    #        data.notice(mine,maxe)
            data.group_counts(10)
     
    #        data = pyfits.open(i+'/srcflux_output/sers_{0:4d}_group_10.pi'.format(j).replace(' ','0'))
    #        dat = []
    #        for m in np.arange(data[1].data['ENERG_HI'].size):
    #            dat.append(sum(data[1].data['MATRIX'][m]))
            plotd = data.to_plot()
            ax.scatter(plotd[0],plotd[1],color='black',label='source')
            ax.plot((ufit[1].data['XHI']+ufit[1].data['XLO'])/2.,ufit[1].data['MODEL'],'--',color='red',label='unabs')
            ax.plot((afit[1].data['XHI']+afit[1].data['XLO'])/2.,afit[1].data['MODEL'],'--',color='blue',label='abs')
            ax.legend(loc='best')
            ax.set_title('unabs kT = {0:3.2f}eV, abs kT = {1:3.2f}'.format(float(bigtab['ukT'][j]),float(bigtab['akT'][j])))
         
    
            fancy_plot(ax)
            ax.set_xlim([0.,8.])
            ax.set_yscale('log')
#            ax.set_xscale('log')
            ax.set_ylim([0.000001,1.0])
            
            ax.set_xlabel('Energy [eV]')
            ax.set_ylabel('Flux [dn/s/cm$^2$]')
            
            
            fig.savefig('spectrum_plots/cstat/e_{1:5d}_src_{0:4d}_fit.png'.format(j,int(i)).replace(' ','0'),bbox_pad=.1,bbox_inches='tight')
            plt.clf()
            plt.close()
        except:
            print '' 
