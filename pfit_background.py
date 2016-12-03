import numpy as np
from sherpa.astro.ui import *
#from sherpa.astro.chips import *
from sherpa.plot.chips_backend import *
from sherpa_contrib.profiles import *
from sherpa.astro.utils import _charge_e as q # ~1.6e-9 ergs/1 keV photon (nist.gov)
import glob
import sys
import getopt


diri = '/srcflux_output/'

mine = 0.5
maxe = 8.0
AvnH = 1.69*10**21
Av = 2.5


#for p,i in enumerate(dirs):
#def main(argv):
##read in command line
#    try:
#        opts, args = getopt.getopt(argv,"ho:",["obsid="])
#    except getopt.GetoptError:
#        print 'pfit_cstat.py -o <Chandra OBSID>'
#        sys.exit(2)
#
#    for opt,arg in opts:
#        if opt == '-h':
#            print 'pfit_cstat.py -o <Chandra OBSID>'
#            print 'Chandra OBSID is the Chandra OBSID assuming a directory structure contains its name'
#            sys.exit()
##        elif opt in ('-o','--obsid'):
#            i = arg

set_analysis("energy")
#    files = glob.glob(i+diri+'*bkg.pi')
#    files = files[:1]
#    outf = open(i+'/prelim_out_cstat_{0:5d}.dat'.format(int(i)).replace(' ','0'),'w')
#    outf.write('{8:^10}{0:^10}{1:^10}{2:^10}{3:^10}{4:^10}{5:^10}{6:^10}{7:^10}{9:^10}{10:^10}{11:^10}{12:^10}{13:^10}{14:^10}{15:^10}{16:^10}\n'.format('uflux','uflux_err','aflux','aflux_err','cnts',' cnts_err',' ncnts',' ncnts_err','src','ukT','unH','uchi','urstat','akT','anH','achi','arstat'))
#    for j in np.arange(1,len(files)+1):
#    for j in np.arange(1,10):
#       try:
gc = 8
dirs = ['9919','9920','10809','10810','10811','10812']
#set stat to chidvar
set_stat('chi2datavar')
#set_stat("cstat")
for i,j in enumerate(dirs):

 #   try:
            diri = dirs[i]
            exam = '/summed_background_{0}_src'.format(j).replace(' ','0')
            sfile = diri+exam
            #un#subtracted source spectrum
            load_pha(sfile+".pi")
            #instrument response on the source
#            rmf = get_rmf()
#            umod = xsconstant.c*xspowerlaw.p1*xsgaussian.g1
            umod = xsgaussian.g1+xsgaussian.g2+xswabs.a1*xsvapec.b1+xswabs.a1*xspowerlaw.p1*xsconstant.c#xsgaussian.g3+xsgaussian.g4+xsgaussian.g5

            a1.nH = 0.02
            a1.nH.min = 0.0001
            a1.nH.max = 0.30

            b1.kT = 0.18
            b1.kT.min = 0.10
            b1.kT.max = 0.40

            c.factor = 0.00000000001

            p1.PhoIndex=1.40
            p1.PhoIndex.max = 1.60
            p1.PhoIndex.min = 1.00
            p1.norm=10.0
            p1.norm.max = 13.00
            p1.norm.min = 5.00

            g1.LineE=7.400
            g1.LineE.max = 7.600
            g1.LineE.min = 7.200
            g1.Sigma=0.005
            g1.Sigma.max = 0.017
            g1.Sigma.min = 0.003

            g2.LineE=2.200
            g2.LineE.max = 2.600
            g2.LineE.min = 2.000
            g2.Sigma=0.005
            g2.Sigma.max = 0.017
#
#            g3.LineE=0.500
#            g3.LineE.max = 0.800
#            g3.LineE.min = 0.200
#            g3.Sigma=0.005
#            g3.Sigma.max = 0.417
#            g3.Sigma.min = 0.003
#
#            g4.LineE=1.300
#            g4.LineE.max = 1.400
#            g4.LineE.min = 1.200
#            g4.Sigma=0.005
#            g4.Sigma.max = 0.117
#            g4.Sigma.min = 0.003
#
#
#            g5.LineE=1.500
#            g5.LineE.max = 1.600
#            g5.LineE.min = 1.400
#            g5.Sigma=0.005
#            g5.Sigma.max = 0.117
#            g5.Sigma.min = 0.003
#            g1.LineE=7.4
#            g1.LineE.max = 7.6
#            g1.LineE.min = 7.2
#            g1.Sigma=0.05
#            g1.Sigma.max = 0.07
#            g1.Sigma.min = 0.03
#            g1.norm=3.

            #assume an initial extinction
            set_source(umod)
#model for background
            #guess initial model and fit
#            guess(a1)
#            guess(b1)
#set the maximum and minimum range to fit
            notice(mine,maxe)
#do the fit
#            group_counts(10)
#            fit()
            ####calc_stat_info()
#            proj()
#            get_proj()
#            print get_proj_results()

#no counting errors
#            #change the abs fit color to blue
#            plt.set_curve("crv4",["*.color","blue"])
#            #increase label size
#            plt.set_axis(["label.size",20])
#            #put parameters in title
#          plot the unabs fit
#            plt['xlog'] = True
#            plt['ylog'] = True
#            plt['size'] = 20
#            plt["title"] = "unabs kT = {0:4.2f} keV, nH = {1:4.2f}, chi = {2:4.3f}".format(ukt,unh,uchi)
#            plt['color'] = 'blue'
#            plot_fit()
            get_data_plot_prefs()["xlog"] = True
            get_data_plot_prefs()["ylog"] = True
            get_data_plot_prefs()['symbolcolor'] = 'white'
            get_data_plot_prefs()['errcolor'] = 'white'
            plot_data(1,overplot=False)
            #best fit values absorbed values
#            k = raw_input('Press Enter to Continue')
            colors = ['red','blue','green','orange']
            for p in range(4):
                load_pha(p+4,"{0:4d}/background/e{0:5d}_I{1:1d}.fits.pi".format(int(j),p+1).replace(' ','0'))
                notice(mine,maxe)
                get_data_plot_prefs()['symbolcolor'] = colors[p]
                get_data_plot_prefs()['errcolor'] = colors[p]
                plot_data(id=p+4,overplot=True)

#    except:
#        print 'ERROR'
#            g1.LineE=7.4
#            g1.LineE.max = 7.6
#            g1.LineE.min = 7.2
#            g1.Sigma=0.05
#            g1.Sigma.max = 0.07
#            g1.Sigma.min = 0.03
#            g1.norm=3.

            #assume an initial extinction
#            set_source(umod)
#model for background
            #guess initial model and fit
#            guess(a1)
#            guess(b1)
#set the maximum and minimum range to fit
#            notice(mine,maxe)
#do the fit
#            group_counts(10)
#            fit()
            ####calc_stat_info()
#            proj()
#            get_proj()
#            print get_proj_results()

#no counting errors
#            #change the abs fit color to blue
#            plt.set_curve("crv4",["*.color","blue"])
#            #increase label size
#            plt.set_axis(["label.size",20])
#            #put parameters in title
#          plot the unabs fit
#            plt['xlog'] = True
#            plt['ylog'] = True
#            plt['size'] = 20
#            plt["title"] = "unabs kT = {0:4.2f} keV, nH = {1:4.2f}, chi = {2:4.3f}".format(ukt,unh,uchi)
#            plt['color'] = 'blue'
#            get_data_plot_prefs()["xlog"] = True
#            get_data_plot_prefs()["ylog"] = True
#            plot_fit()
            #best fit values absorbed values
            k = raw_input('Press Enter to Continue')
#    except:
#        print 'ERROR'
