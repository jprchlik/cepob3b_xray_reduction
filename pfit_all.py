import numpy as np
from sherpa.astro.ui import *
#CHIPS plotting
#from pychips import *
from sherpa.astro.utils import _charge_e as q # ~1.6e-9 ergs/1 keV photon (nist.gov)
import glob
import sys
import getopt


#load_pha("3c273.pi")
diri = '/srcflux_output/'
#dirs = ['9919']
#dirs = ['10809']

mine = 0.5
maxe = 3.0
AvnH = 1.69*10**21
Av = 2.5
get_data_plot_prefs()["xlog"] = True
get_data_plot_prefs()["ylog"] = True

#for p,i in enumerate(dirs):
def main(argv):
#read in command line
    try:
        opts, args = getopt.getopt(argv,"ho:c:w",["obsid=","chidvar=","wstat"])
    except getopt.GetoptError as err:
        print err
        print 'pfit_chi2datavar.py -o <Chandra OBSID>'
        sys.exit(2)

    chid = False
    wsta = False
    stat = ''
    for opt,arg in opts:
        if opt == '-h':
            print 'pfit_chi2datavar.py -o <Chandra OBSID> '
            print 'Chandra OBSID is the Chandra OBSID assuming a directory structure contains its name'
            sys.exit()
        elif opt in ('-o','--obsid'):
            i = arg
        elif opt in ('-c','--chidvar'):
            chid = True
            stat = 'chi2datavar'
        elif opt in ('-w','--wstat'):
            wsta = True
#changed to cstat to check background issue
            stat = 'wstat'

    set_analysis("energy")
    files = glob.glob(i+diri+'*bkg.pi')
#create an average background




#    files = files[:1]
    outf = open(i+'/prelim_out_{1}_{0:5d}_3keV.dat'.format(int(i),stat).replace(' ','0'),'w')
    outf.write('{8:^10}{0:^10}{1:^10}{2:^10}{3:^10}{4:^10}{5:^10}{6:^10}{7:^10}{9:^10}{10:^10}{11:^10}{12:^10}{13:^10}{14:^10}{15:^10}{16:^10}{17:^10}{18:^10}\n'.format('uflux','uflux_err','aflux','aflux_err','cnts',' cnts_err',' ncnts',' ncnts_err','src','ukT','unH','ukterr','unHerr','uchi','urstat','akT','anH','achi','arstat'))
    for j in np.arange(1,len(files)+1):
#    for j in np.arange(1,10):
       try:
            exam = 'sers_{0:4d}'.format(j).replace(' ','0')
            
            
    
            sfile = i+diri+exam
            #un#subtracted source spectrum
            load_pha(sfile+".pi")
            #instrument response on the source
            load_arf(sfile+".arf")
            #instrument response on the background 
            load_rmf(sfile+".rmf")
            #background spectrum
#            load_bkg(i+"/summed_background_{0}_src.pi".format(i))
#          Update with full background source
            load_bkg(i+"/background/summed_background_{0}_src.pi".format(i))
            rmf = get_rmf()
            arf = get_arf()
            bkg = get_bkg()
            bkg_scal = get_bkg_scale()
            set_backscal(bkg_scal)
#            notice_id(1,mine,maxe)
            if chid:
#Skipping background subtraction for now
                subtract()
#group points in bins of at least 10
                group_counts(8)
#            save_grouping(sfile+"_group_10.pi",ascii=False,clobber=True)
#applying group to backgroud
#            set_grouping(1,get_grouping(1),bkg_id=1)
            
            
            data_sum = calc_data_sum(mine,maxe)  
             
            data_cnt_rate = calc_data_sum()/get_exposure() 
            
            bkg_sum = calc_data_sum(mine,maxe,bkg_id=1)  
            
            
            bkg_cnt_rate = calc_data_sum(bkg_id=1)/get_exposure(bkg_id=1) 

            
            
            
            #Winston 2010 used xsraymond
#           Scott Wolk says use xsvapec because it has more lines (9/1/16)
#            umod = xswabs.a1*xsraymond.b1
            set_stat(stat) 
            umod = xswabs.abs1*xsvapec.b1
            b1.cache=0
            b1.kT = 1.5 #kT
            b1.kT.max = 50.
            b1.kT.min = 0.01
#            b1.norm = 2.1599743e+21 #cm
            
            #use a background model for the background
            #assume an initial extinction
            abs1.cache = 0
            abs1.nH = 0.02
            abs1.nH.min = 1e-4
            abs1.nH.max = 100
#            set_full_model(umod)
            thaw(abs1)
            set_source(umod)
            #guess initial model and fit
#            guess(abs1)
#            guess(b1)
#set the maximum and minimum range to fit
            notice(mine,maxe)
#Set up back ground and background model
#            set_bkg(bkg)
##setup model components from Brad Spitzerbar (constant,power law, and Gaussian)
#            set_bkg_model(abs1*xsconstant.c*xspowerlaw.p1*xsgaussian.g1)
##set up backgroud model parameters
#            c.factor = 0.0
#            p1.PhoIndex=2
#            p1.norm=0.
#            g1.LineE=0.
#            g1.Sigma=.0
#            g1.norm=0.0
#            
#           
##freeze in backgroud parameters
#            freeze(c.factor)
#            freeze(p1.PhoIndex)
#            freeze(p1.norm)
#            freeze(g1.LineE)
#            freeze(g1.Sigma)
#            freeze(g1.norm)

#fit the model
            fit()

#         plot the unabs fit
#            plot_fit()
            #best fit values absorbed values
            unh,ukt,unorm = get_fit_results().parvals
            ####calc_stat_info()
            uchi = get_fit_results().statval
            urst = get_fit_results().rstat
#            ugoodness=get_stat_info()
#            uchi = ugoodness[0].statval
#            urst = ugoodness[0].rstat
            #Calculate Error for model using 100 samples
#            unabs = sample_flux(umod,mine,maxe,num=100)
#Switch from total model to just xsraymond model because of Doug Burke (2016/12/07)
            unabs = sample_flux(b1,mine,maxe,num=1000)
	    #      save best fit model as fits file
            save_model(sfile+"_{0}_umod.pi".format(stat),ascii=False,clobber=True)
            save_data(sfile+"_{0}_data.pi".format(stat),ascii=False,clobber=True)
#get 1sigma errors of parameters
            conf()
            get_conf()
            perr_a = get_conf_results()
            perr = (np.array(perr_a.parmaxes)-np.array(perr_a.parmins))/2.
            unHerr = perr[0]
            ukTerr = perr[1]
            

            
###            #Winston 2010 used xsraymond
###            amod = xsvapec.b1
###            set_stat(stat) 
###            set_source(amod)
###            b1.cache=0
###            b1.kT = 1.5 #kT
###            b1.kT.max = 50.
###            b1.kT.min = 0.01
###            #guess initial model and fit
###            guess(b1)
#setup model components from Brad Spitzerbar (constant,power law, and Gaussian)
#            set_bkg_model(xsconstant.c*xspowerlaw.p1*xsgaussian*g1)
##freeze in previous backgroud parameters
#            freeze(c.factor)
#            freeze(p1.PhoIndex)
#            freeze(p1.norm)
#            freeze(g1.LineE)
#            freeze(g1.Sigma)
#            freeze(g1.norm)

####            fit()
            #plot the abs fit
#            plot_fit(overplot=True)
            #change the abs fit color to blue
#            set_curve("crv4",["*.color","blue"])
#            #best fit values absorbed values
#Commented out seperate unattenuated fit because of Doug Burke (2016/12/07)
###            akt,anorm = get_fit_results().parvals
###            anh = 0.0
###            ####calc_stat_info()
###            achi = get_fit_results().statval
###            arst = get_fit_results().rstat
            akt, anorm = -9999.9,-9999.9
            anh = 0.0
            achi,arst = -9999.9,-9999.9
###            
###            
###            #Calculate Error for model using 100 samples
###            anabs = sample_flux(amod,mine,maxe,num=100)
####      save best fit model as fits file
###            save_model(sfile+"_{0}_amod.pi".format(stat),ascii=False,clobber=True)
####get the total fluxes and errors in a 3D vector 
            unabs = unabs[1]
            anabs = anabs[0]
####total array
            tot_a = anabs
            tot_u = unabs 
### 
####get the average err abs[2] is a negative number
            anabs_err = (anabs[1]-anabs[2])/2.
            unabs_err = (unabs[1]-unabs[2])/2. 
####get the median fluxes
            anabs = anabs[0]
            unabs = unabs[0]

#no counting errors
            data_sum_err = 0.0

            tot_c_err = 0.0
#write output to file
        
            outf.write('{8:^10d}{0:^10.3e}{1:^10.3e}{2:^10.3e}{3:^10.3e}{4:^10.1f}{5:^10.1f}{6:^10.1f}{7:^10.1f}{9:^10.4f}{10:^10.4f}{17:^10.4f}{18:^10.4f}{11:^10.4f}{12:^10.4f}{13:^10.4f}{14:^10.4f}{15:^10.4f}{16:^10.4f}\n'.format(unabs,unabs_err,anabs,anabs_err,data_sum,data_sum_err,data_sum-bkg_sum,tot_c_err,j,ukt,unh,uchi,urst,akt,anh,achi,arst,ukTerr,unHerr))
        
       except:
#########Write out -9999 if cannot find solution
            print 'FAILED'
            outf.write('{8:^10d}{0:^10d}{1:^10d}{2:^10d}{3:^10d}{4:^10d}{5:^10d}{6:^10d}{7:^10d}{9:^10d}{10:^10d}{11:^10d}{12:^10d}{13:^10d}{14:^10d}{15:^10d}{16:10d}{17:10d}{18:10d}\n'.format(-999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,j,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999))
    outf.close()

if __name__ == "__main__":
    main(sys.argv[1:])
