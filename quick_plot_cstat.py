from astropy.io import ascii
from fancy_plot import *
from astropy.table import Table,hstack,join
import pyfits
import matplotlib.pyplot as plt
import numpy as np
import glob
from astropy import units as u
from astropy.coordinates import SkyCoord


def match_stars(ss,xs,cut=2.0):
    ras = ss['RA']
    des = ss['Dec']
    xarray = []
    sarray = []
    for i in np.arange(xs['RA'].size):
        rax = xs['RA'][i]
        dex = xs['Dec'][i]
        ddr = (rax-ras)*np.cos(np.radians(dex))
        ddd = dex-des
        dis = np.sqrt((ddd**2.+ddr**2))
#convert to arcseconds
        dis = dis*3600.
        found, = np.where(dis < cut)
        if found.size == 1:
            sarray.append(found[0])
            xarray.append(i)
        if found.size > 1:
            mindis, = np.where(dis == dis.min())
            sarray.append(mindis[0])
            xarray.append(i)
    sarray = np.array(sarray)
    xarray = np.array(xarray)
    return sarray,xarray



dirs = ['9919','9920','10809','10810','10811','10812']

odir = '/home/jakub/wandajune_home/Ancillary_Projects/tol_Halpha/new_xray/'

fig1,ax1 = plt.subplots()
fig2,ax2 = plt.subplots()
fig3,ax3 = plt.subplots()
fig4,ax4 = plt.subplots()
fig5,ax5 = plt.subplots()
for i in dirs:
    print i

    ocounts = ascii.read(odir+'{0:5d}_cstat_all.dat'.format(int(i)).replace(' ','0'))
    counts1 = pyfits.open(i+'/obsid_{0:5d}_exp_src.fits'.format(int(i)).replace(' ','0'))
    counts = Table(counts1[1].data)
    counts = counts['RA','DEC','NET_COUNTS','NET_COUNTS_ERR']
    counts.rename_column('DEC','Dec')
    fluxes = ascii.read(i+'/prelim_out_cstat_{0:5d}.dat'.format(int(i)).replace(' ','0'))
    bigtab = hstack([counts,fluxes])
    bigtab.write('reduced_xrays_obsid{0:5d}.ipac'.format(int(i)).replace(' ','0'),format='ascii.ipac')

    nxarray,oxarray = match_stars(counts,ocounts,cut=2.0)
    good, = np.where(fluxes['uflux'] > 0)
    
    ax1.scatter(counts['NET_COUNTS'][good],100.*(fluxes['uflux_err']/fluxes['uflux'])[good],color='black')
    ax2.scatter(counts['NET_COUNTS'][good],100*(fluxes['aflux_err']/fluxes['aflux'])[good],color='black')
    ax3.scatter(counts['NET_COUNTS'][nxarray],ocounts['net_cnts'][oxarray],color='black')
    ax4.scatter(np.log10(bigtab['aflux'][nxarray]),np.log10(ocounts['abs.flux'][oxarray]),color='black')

# difference in counts
    diffc = counts['NET_COUNTS'][nxarray]-ocounts['net_cnts'][oxarray]
    bad, = np.where(np.abs(diffc) > 3.*diffc.std())
    print '{0:^5}{1:^5}{2:^7}{3:^7}{4:^15}{5:^15}'.format('epoch','id','ncnt','ocnt','ra','dec')  
    for j in bad:
        ra = float(bigtab['RA'][nxarray][j])
        dec = float(bigtab['Dec'][nxarray][j])
        c = SkyCoord(ra=ra*u.degree,dec=dec*u.degree)
        ra = c.ra.hms
        dec = c.dec.dms
        inc = c.to_string('hmsdms')
        print '{0:^5d}{1:^5d}{2:^7.1f}{3:^7.1f}{4:^30}'.format(int(i),int(bigtab['src'][nxarray][j]),float(bigtab['NET_COUNTS'][nxarray][j]),float(ocounts['net_cnts'][oxarray][j]),inc)

ax4.set_xlim([-16,-11])
ax4.set_ylim([-16,-11])
ax5.set_xlim([-16,-11])
ax5.set_ylim([-16,-11])

xs4 = np.array(ax4.get_xlim())
ys4 = np.array(ax4.get_ylim())
xs5 = np.array(ax5.get_xlim())
ys5 = np.array(ax5.get_ylim())

ax3.plot([-100,1000],[-100,1000],'--',color='red',linewidth=3)
ax4.plot([-100,1000],[-100,1000],'--',color='red',linewidth=3)
ax5.plot([-100,1000],[-100,1000],'--',color='red',linewidth=3)

fancy_plot(ax1)
fancy_plot(ax2)
fancy_plot(ax3)
fancy_plot(ax4)
fancy_plot(ax5)

ax1.set_ylim([-10.0,100.0])
ax1.set_xlim([-10.0,300.0])
ax1.set_xlabel('Net Counts')
ax1.set_ylabel('unabs Flux Error [%]')

ax2.set_ylim([-10.0,100.0])
ax2.set_xlim([-10.0,300.0])
ax2.set_xlabel('Net Counts')
ax2.set_ylabel('abs Flux Error [%]')

ax3.set_ylim([-10.0,300.0])
ax3.set_xlim([-10.0,300.0])
ax3.set_xlabel('Net Counts (CIAO)')
ax3.set_ylabel('Net Counts (ACIS Extract)')

ax4.set_ylim(xs4)
ax4.set_xlim(ys4)
ax4.set_xlabel('Jakub Sherpa abs flux')
ax4.set_ylabel('Brad Sherpa abs flux')
ax4.set_title('Cstat')

ax5.set_ylim(xs5)
ax5.set_xlim(ys5)
ax5.set_xlabel('Jakub Sherpa abs flux')
ax5.set_ylabel('Brad Sherpa abs flux')
ax5.set_title('Cstat')

fig1.savefig('test_cstat_errors_unabs.eps',bbox_pad=.1,bbox_inches='tight')
fig1.savefig('test_cstat_errors_unabs.png',bbox_pad=.1,bbox_inches='tight')

fig2.savefig('test_cstat_errors_abs.eps',bbox_pad=.1,bbox_inches='tight')
fig2.savefig('test_cstat_errors_abs.png',bbox_pad=.1,bbox_inches='tight')

fig3.savefig('test_cstat_errors_new_cnts.eps',bbox_pad=.1,bbox_inches='tight')
fig3.savefig('test_cstat_errors_new_cnts.png',bbox_pad=.1,bbox_inches='tight')

fig4.savefig('test_cstat_abs_flux.eps',bbox_pad=.1,bbox_inches='tight')
fig4.savefig('test_cstat_abs_flux.png',bbox_pad=.1,bbox_inches='tight')

fig4.savefig('test_cstat_unabs_flux.eps',bbox_pad=.1,bbox_inches='tight')
fig4.savefig('test_cstat_unabs_flux.png',bbox_pad=.1,bbox_inches='tight')
