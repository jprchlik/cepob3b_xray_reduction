from astropy.io import ascii
from fancy_plot import *
from astropy.table import Table,hstack,join
import pyfits
import matplotlib.pyplot as plt
import numpy as np
import glob
import pandas


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

physprop = ascii.read('../tol_Halpha/all_spectra/Combined_hectospec.ipac',format='ipac')
physprop.rename_column('DEC','Dec_1')
physprop.rename_column('RA','RA_1')


#new RA and DEC
astr = pandas.read_pickle('../tol_Halpha/Ha_photometry/cep_astr_struc_astr.pic')

#Add ra and dec
physprop['RA'] = astr['J2MRA'][0][physprop['index']]
physprop['Dec'] = astr['J2MDEC'][0][physprop['index']]



dirs = ['9919','9920','10809','10810','10811','10812']

odir = '/home/jakub/wandajune_home/Ancillary_Projects/tol_Halpha/new_xray/'

fig1,ax1 = plt.subplots()
fig2,ax2 = plt.subplots()
fig3,ax3 = plt.subplots()
fig4,ax4 = plt.subplots()
fig5,ax5 = plt.subplots()
fig6,ax6 = plt.subplots()
fig7,ax7 = plt.subplots()
fig8,ax8 = plt.subplots(nrows=3,figsize=(8.5,11.))

for i in dirs:
    print i

    ocounts = ascii.read(odir+'{0:5d}_cstat_all.dat'.format(int(i)).replace(' ','0'))
    counts1 = pyfits.open(i+'/obsid_{0:5d}_exp_src.fits'.format(int(i)).replace(' ','0'))
    counts = Table(counts1[1].data)
    counts = counts['RA','DEC','NET_COUNTS','NET_COUNTS_ERR','COMPONENT']
    counts.rename_column('COMPONENT','src')
    counts.rename_column('DEC','Dec')
    fluxes = ascii.read(i+'/prelim_out_wstat_{0:5d}_test.dat'.format(int(i)).replace(' ','0'))
    fluxes.rename_column('RA','RA_f')
    fluxes.rename_column('Dec','Dec_f')
#    bigtab = hstack([counts,fluxes])
    bigtab = join(counts,fluxes,join_type='inner',keys=['src'])
    bigtab.write('reduced_xrays_obsid{0:5d}.ipac'.format(int(i)).replace(' ','0'),format='ascii.ipac')
    #Color array
    bigtab['color'] = 'black'
    bigtab['color'][((bigtab['NET_COUNTS'] < 30.) & (bigtab['NET_COUNTS'] >= 10.))] = 'blue'
    bigtab['color'][(bigtab['NET_COUNTS'] < 10.)] = 'red'

    nxarray,oxarray = match_stars(bigtab,ocounts,cut=2.0)
    good, = np.where(bigtab['uflux'] > 0)
######Match Phys prob
    pnxarray,poxarray = match_stars(bigtab,physprop,cut=2.0)
    
    ax1.scatter(bigtab['NET_COUNTS'][good],100.*(bigtab['uflux_err']/bigtab['uflux'])[good],color='black')
    ax2.scatter(bigtab['NET_COUNTS'][good],100*(bigtab['aflux_err']/bigtab['aflux'])[good],color='black')
    ax3.scatter(bigtab['NET_COUNTS'][nxarray],ocounts['net_cnts'][oxarray],color='black')
    ax3.scatter(bigtab['NET_COUNTS'][nxarray],ocounts['net_cnts'][oxarray],color='black')

    ax4.scatter(np.log10(bigtab['b_uflux'][nxarray]),np.log10(ocounts['total.unabs.flux'][oxarray]),color='black')
    bad, = np.where(bigtab['NET_COUNTS'][nxarray] <= 10.)
    ax4.scatter(np.log10(bigtab['b_uflux'][nxarray][bad]),np.log10(ocounts['total.unabs.flux'][oxarray][bad]),color='red',zorder=200)
    print np.median(bigtab['aflux'][nxarray]),np.median(ocounts['abs.flux'][oxarray])
    ax6.scatter(ocounts['net_cnts'][oxarray],np.log10(bigtab['b_uflux'])[nxarray],color=bigtab['color'][nxarray])
    ax7.scatter(bigtab['NET_COUNTS'],np.log10(bigtab['b_uflux']),color=bigtab['color'][nxarray])


    
    use, = np.where((np.isfinite(np.log10(physprop['lbol'][poxarray]))) & (np.isfinite(np.log10(bigtab['b_uflux'][pnxarray]))))
    ax8[0].scatter(np.log10(physprop['lbol'][poxarray][use]),np.log10(bigtab['b_uflux'][pnxarray][use]),color=bigtab['color'][pnxarray][use])
    ax8[1].scatter(np.log10(physprop['lbol'][poxarray][use]),bigtab['ukT'][pnxarray][use],color=bigtab['color'][pnxarray][use])
    ax8[2].scatter(np.log10(physprop['lbol'][poxarray][use]),bigtab['unH'][pnxarray][use],color=bigtab['color'][pnxarray][use])

ax4.set_xlim([-17,-11])
ax4.set_ylim([-17,-11])
ax5.set_xlim([-17,-11])
ax5.set_ylim([-17,-11])

ax6.set_xlim([0,300])
ax7.set_xlim([0,300])
ax6.set_ylim([-17,-11])
ax7.set_ylim([-17,-11])
ax8[0].set_ylim([-17,-11])
ax8[1].set_ylim([-.5,21])
ax8[2].set_ylim([-.5,2.3])

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
fancy_plot(ax6)
fancy_plot(ax7)
fancy_plot(ax7)

ax1.set_ylim([-10.0,300.0])
ax1.set_xlim([-10.0,300.0])
ax1.set_xlabel('Net Counts')
ax1.set_ylabel('unabs Flux Error [%]')

ax2.set_ylim([-10.0,300.0])
ax2.set_xlim([-10.0,300.0])
ax2.set_xlabel('Net Counts')
ax2.set_ylabel('abs Flux Error [%]')

ax3.set_ylim([-10.0,300.0])
ax3.set_xlim([-10.0,300.0])
ax3.set_xlabel('Net Counts (CIAO)')
ax3.set_ylabel('Net Counts (ACIS Extract)')

ax4.set_ylim(xs4)
ax4.set_xlim(ys4)
ax4.set_xlabel('Jakub Sherpa unabs flux')
ax4.set_ylabel('Brad Sherpa unabs flux')
ax4.set_title('Wstat')

ax5.set_ylim(xs5)
ax5.set_xlim(ys5)
ax5.set_xlabel('Jakub Sherpa abs flux')
ax5.set_ylabel('Brad Sherpa abs flux')
ax5.set_title('Wstat')

ax6.set_ylabel('unabs flux [log(ergs/s/cm$^2$)]')
ax7.set_ylabel('unabs flux [log(ergs/s/cm$^2$)]')
ax8[0].set_ylabel('unabs flux [log(ergs/s/cm$^2$)]')
ax8[1].set_ylabel('kT [keV]')
ax8[2].set_ylabel('nH [cm$^2$]')

ax6.set_xlabel('Net Counts (CIAO)')
ax7.set_xlabel('Net Counts (ACIS Extract)')
ax8[0].set_xlabel('L/L$_{\odot}$')
ax8[1].set_xlabel('L/L$_{\odot}$')
ax8[2].set_xlabel('L/L$_{\odot}$')
fancy_plot(ax8[0])
fancy_plot(ax8[1])
fancy_plot(ax8[2])

fig1.savefig('test_wstat_errors_unabs.eps',bbox_pad=.1,bbox_inches='tight')
fig1.savefig('test_wstat_errors_unabs.png',bbox_pad=.1,bbox_inches='tight')

fig2.savefig('test_wstat_errors_abs.eps',bbox_pad=.1,bbox_inches='tight')
fig2.savefig('test_wstat_errors_abs.png',bbox_pad=.1,bbox_inches='tight')

fig3.savefig('test_wstat_errors_new_cnts.eps',bbox_pad=.1,bbox_inches='tight')
fig3.savefig('test_wstat_errors_new_cnts.png',bbox_pad=.1,bbox_inches='tight')

fig4.savefig('test_wstat_abs_flux.eps',bbox_pad=.1,bbox_inches='tight')
fig4.savefig('test_wstat_abs_flux.png',bbox_pad=.1,bbox_inches='tight')

fig4.savefig('test_wstat_unabs_flux.eps',bbox_pad=.1,bbox_inches='tight')
fig4.savefig('test_wstat_unabs_flux.png',bbox_pad=.1,bbox_inches='tight')

fig6.savefig('test_wstat_unabs_flux_ciao_counts.eps',bbox_pad=.1,bbox_inches='tight')
fig6.savefig('test_wstat_unabs_flux_ciao_counts.png',bbox_pad=.1,bbox_inches='tight')

fig7.savefig('test_wstat_unabs_flux_acis_counts.eps',bbox_pad=.1,bbox_inches='tight')
fig7.savefig('test_wstat_unabs_flux_acis_counts.png',bbox_pad=.1,bbox_inches='tight')

fig8.savefig('test_wstat_unabs_flux_sa_counts.eps',bbox_pad=.1,bbox_inches='tight')
fig8.savefig('test_wstat_unabs_flux_sa_counts.png',bbox_pad=.1,bbox_inches='tight')
