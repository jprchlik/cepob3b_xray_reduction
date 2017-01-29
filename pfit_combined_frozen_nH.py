import numpy as np
from sherpa.astro.ui import *
#CHIPS plotting
#from pychips import *
from sherpa.astro.utils import _charge_e as q # ~1.6e-9 ergs/1 keV photon (nist.gov)
import glob
import sys
import getopt
from crates_contrib.utils import SimpleCoordTransform
import sherpa_contrib.utils as utils

sys.path.append("/home/jakub/anaconda")
sys.path.append("/home/jakub/wandajune_home/MyModules")


#700 pc in cm
dis = 2.16e21
nH=1.##/cm^3
ne=1.##/cm^3
#local numeber of atoms along the line of sight
den = (nH+ne)*2.*np.pi*dis



def set_sources(j,i,num):
    global chid
    print j,i,num+1
    exam = 'sers_{0:4d}'.format(int(j)).replace(' ','0')
    sfile = i+diri+exam
    #un#subtracted source spectrum
    load_pha(num+1,sfile+".pi")
    #background spectrum
    load_bkg(num+1,i+"/background/summed_background_{0}_src.pi".format(i))
    bkg = get_bkg(num+1)
    bkg_scal = get_bkg_scale(num+1)
    set_backscal(num+1,bkg_scal)
    if chid:
#Skipping background subtraction for now
        subtract(num+1)
#group points in bins of at least 10
        group_counts(num+1,8)
    tr = SimpleCoordTransform(i+'/obs'+i+'_img.fits')
    return sfile,tr
            


#load_pha("3c273.pi")
diri = '/srcflux_output/'
#dirs = ['9919']
#dirs = ['10809']

mine = 0.5
maxe = 3.0
AvnH = 1.69*10**21 #classic value
AvnH = 1.5*10**21 #My cep OB3b Value
Av = 2.5
get_data_plot_prefs()["xlog"] = True
get_data_plot_prefs()["ylog"] = True

#for p,i in enumerate(dirs):
def main(argv):
    global chid
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
#all for fitting of the same source across epochs
            if i == 'all': multifit = True
        elif opt in ('-c','--chidvar'):
            chid = True
            stat = 'chi2datavar'
        elif opt in ('-w','--wstat'):
            wsta = True
#changed to cstat to check background issue
            stat = 'wstat'

#    data = ascii.read('compare_plots/combined_source_cepob3b_list.ipac',format='ipac')
    data = np.genfromtxt('compare_plots/combined_source_cepob3b_list.csv',delimiter=',',names=True)

#array of epochs
    epochs = np.array(['9919','9920','10809','10810','10811','10812'])
#only get ones that have lbol and possibly fit
    fitlist, = np.where(((data['e09919_src_num'] > 0.) | (data['e09920_src_num'] > 0.) | (data['e10809_src_num'] > 0.) | (data['e10810_src_num'] > 0.) | (data['e10811_src_num'] > 0. ) | (data['e10812_src_num'] > 0.)) & (data['teff'] > -100.))

    data = data[fitlist]
    print '{0:4d} sources found that have temperatures'.format(fitlist.size)
#    files = files[:1]
    outf = open('combined/prelim_out_{1}_frozen_nH.dat'.format(i,stat).replace(' ','0'),'w')
    outf.write('{8:^10}{0:^10}{1:^10}{2:^10}{3:^10}{4:^10}{5:^10}{6:^10}{7:^10}{9:^10}{10:^10}{11:^10}{12:^10}{13:^10}{14:^10}{15:^10}{16:^10}{17:^10}{18:^10}{19:^15}{20:^15}{21:^10}{22:^10}{23:^10}\n'.format('uflux','uflux_err','aflux','aflux_err','cnts',' cnts_err',' ncnts',' ncnts_err','src','ukT','unH','ukterr','unHerr','uchi','urstat','akT','anH','achi','arstat','RA','Dec','b_uflux','b_aflux','num_fit'))
#    for j in data['ID']:
    for j in data['ID'][-5:]:
#    for j in np.arange(1,3):
       try:
        #Set stat to use in fitting
            set_analysis("energy")
            set_stat(stat) 

#get list of source ids
            j = int(j)
            print '###############################'
            print 'ID NUMBER {0:4d}'.format(j)
            idnum, = np.where(data['ID'] == j)
            print "Number in list"
            print idnum

            idlist = data[idnum]
            idlist = np.array([int(idlist['e09919_src_num']),int(idlist['e09920_src_num']),int(idlist['e10809_src_num']),int(idlist['e10810_src_num']),int(idlist['e10811_src_num']),int(idlist['e10812_src_num'])])
            #get indexes of matches
            print idlist.T
            print np.where(idlist.T > 0)
            fids, = np.where(idlist > 0)
            num = np.arange(fids.size)

#source string
            scrstr = ''
#source list
            srclis = []
#total counts in observation
            data_sum = 0.
#set up loop id variable
            kk = 1
#number in fit 
            numfit = fids.size
#set up sources and ids
            if numfit > 1:
                for jj in num:
                   sfile,tr = set_sources(idlist[fids][jj],epochs[fids][jj],num[jj])
                   scrstr = scrstr+','+str(num[jj])
                   srclis.append(jj)
                   data_sum = data_sum+calc_data_sum(mine,maxe,id=kk)
                   kk +=1  
            else:
                sfile,tr = set_sources(idlist[fids][0],epochs[fids][0],num[0])
                scrstr = ','+str(num[0])
                srclis = num[0]
                data_sum = data_sum+calc_data_sum(mine,maxe,id=kk)
#remove leading ,
            scrstr = scrstr[1:]
#Source ra and dec
            tcr = utils.pycrates.read_file(sfile+"_srcreg.fits[cols x,y]")
            sx1 = utils.pycrates.copy_colvals(tcr,"x")
            sy1 = utils.pycrates.copy_colvals(tcr,"y")
            (ra, dec) = tr.convert("sky","world",sx1[0],sy1[0])
            ra,dec = float(ra),float(dec)


           #old worthless stuff that remains for convience 
            data_cnt_rate = calc_data_sum()/get_exposure() 
            bkg_sum = calc_data_sum(mine,maxe,bkg_id=1)  
            bkg_cnt_rate = calc_data_sum(bkg_id=1)/get_exposure(bkg_id=1) 

#paramter max values
#           guess nH values
            av = data[idnum]['phot_av']
            try:
               av = float(av)
            except ValueError:
               av = 2.5 #cloud average
            if av < -100: av = 2.5
            avg = AvnH*av*1.e-22 #convert av to nH
            print 'Guess nH  = {0:3.2e}, Av = {1:3.2f}'.format(avg*1.e22,av)
            ne = 0.03 #https://arxiv.org/abs/1308.4010 electron density in diffuse molecular clouds
            bnorm = 1.e-14/(4.*np.pi*dis**2.)*dis*2.*np.pi #http://cxc.harvard.edu/sherpa/ahelp/xsapec.html
            psetmax = np.array([5.0, 3.5 ,1./(den*10.**(-14)/(4.*np.pi))])
            psetmin = np.array([ 1e-4, 0.10,1./(den*10.**(-14)/(4.*np.pi))])
            psetgus = np.array([avg, 1.50,1./(den*10.**(-14)/(4.*np.pi))])
        #Winston 2010 used xsraymond
        #      Scott Wolk says use xsvapec because it has more lines (9/1/16)
#            umod = xswabs.a1*xsraymond.b1
            umod = xswabs.abs1*xsvapec.b1 #http://cxc.harvard.edu/sherpa/ahelp/xsvapec.html
            b1.cache=0
            b1.kT = psetgus[1] #kT
            b1.kT.min = psetmin[1]
            b1.kT.max = psetmax[1]
#            b1.norm   = 1.42317e-05

#            b1.norm = psetgus[2]
#            b1.norm.min = psetmin[2]
#            b1.norm.min = psetmax[2]
#              
            #use a background model for the background
            #assume an initial extinction
            abs1.cache = 0
            abs1.nH = psetgus[0]
            abs1.nH.min = psetmin[0]
            abs1.nH.max = psetmax[0] 
            freeze(abs1)


            for p in num:
                set_source(p+1,umod)
#set the maximum and minimum range to fit
            notice(mine,maxe)

#fit the model using a string seperated by commas

            if fids.size == 3:
                fit(1,2,3,outfile='combined/ID_{0:5d}.dat'.format(int(j)).replace(' ','0'),clobber=True)
            elif fids.size == 2:
                fit(1,2,outfile='combined/ID_{0:5d}.dat'.format(int(j)).replace(' ','0'),clobber=True)
            elif fids.size == 1:
                fit(1,outfile='combined/ID_{0:5d}.dat'.format(int(j)).replace(' ','0'),clobber=True)
            elif fids.size == 4:
                fit(1,2,3,4,outfile='combined/ID_{0:5d}.dat'.format(int(j)).replace(' ','0'),clobber=True)
            elif fids.size == 5:
                fit(1,2,3,4,5,outfile='combined/ID_{0:5d}.dat'.format(int(j)).replace(' ','0'),clobber=True)
            elif fids.size == 6:
                fit(1,2,3,4,5,6,outfile='combined/ID_{0:5d}.dat'.format(int(j)).replace(' ','0'),clobber=True)

            combined = calc_energy_flux()
#         plot the unabs fit
#            plot_fit()
            #best fit values absorbed values
            #remove froze nh value from results and replace with calculated
            ukt,unorm = get_fit_results().parvals
            unh = avg
#            ####calc_stat_info()
#            print '########################################'
#            print calc_stat_info()
#            print get_fit_results()
#            print help(get_fit_results())
#            print help(get_fit_results)
#            print '########################################'
            uchi = get_fit_results().statval
            urst = get_fit_results().rstat
#Switch from total model to just xsraymond model because of Doug Burke (2016/12/07)
#Asked the help desk and the response from Aneta is the flux is calculated based on the model parameters. Therefore, sampling should take into account uncertainties correctly.
	    #      save best fit model as fits file
            exam = 'sers_{0:4d}'.format(int(j)).replace(' ','0')
            sfile = 'combined'+diri+exam
            save_model(sfile+"_{0}_umod.pi".format(stat),ascii=False,clobber=True)
#get 1sigma errors of parameters
            conf()
            get_conf()
            perr_a = get_conf_results()
            pmaxs = np.array(perr_a.parmaxes)
            pmins = np.array(perr_a.parmins)
            gmaxs, = np.where([jj is None for jj in pmaxs])
            gmins, = np.where([jj is None for jj in pmins])


#replace none types with the limits
            if gmaxs.size >= 1:
                pmaxs[gmaxs] = psetmax[gmaxs]
            if gmins.size >= 1:
                pmins[gmins] = psetmin[gmins]

            perr = (pmaxs-pmins)/2.
        
#change for frozen nH
            unHerr = 0.
            ukTerr = perr[0]

            unabs = sample_flux(b1,mine,maxe,num=1000,numcores=3)
#            print help(sample_flux)
#            print '####################'
#            print calc_energy_flux(id=2)
#            print '####################'
#            print calc_energy_flux(id=1)
#            print '####################'
            

            
#No longer fit abs parameters
            akt, anorm = -99.9,-99.9
            anh = 0.0
            achi,arst = -99.9,-99.9
####get the total fluxes and errors in a 3D vector 
            anabs = unabs[0]
            unabs = unabs[1]
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
            anabs = combined
####get best fit flux
            set_source(b1)
            best_uflux = calc_energy_flux(mine,maxe)
            set_source(abs1*b1)
            best_aflux = calc_energy_flux(mine,maxe)



#no counting errors
            data_sum_err = np.sqrt(data_sum) 

            tot_c_err = 0.0
#write output to file
        
            outf.write('{8:^10d}{0:^10.3e}{1:^10.3e}{2:^10.3e}{3:^10.3e}{4:^10.1f}{5:^10.1f}{6:^10.1f}{7:^10.1f}{9:^10.4f}{10:^10.4f}{17:^10.4f}{18:^10.4f}{11:^10.4f}{12:^10.4f}{13:^10.4f}{14:^10.4f}{15:^10.4f}{16:^10.4f}{19:^15.8f}{20:^15.8f}{21:^10.3e}{22:^10.3e}{23:^10d}\n'.format(unabs,unabs_err,anabs,anabs_err,data_sum,data_sum_err,data_sum-bkg_sum,tot_c_err,int(j),ukt,unh,uchi,urst,akt,anh,achi,arst,ukTerr,unHerr,ra,dec,best_uflux,best_aflux,numfit))
        
            clean()
       except:
###########Write out -9999 if cannot find solution
            outf.write('{8:^10d}{0:^10d}{1:^10d}{2:^10d}{3:^10d}{4:^10.1f}{5:^10.1f}{6:^10d}{7:^10d}{9:^10d}{10:^10d}{11:^10d}{12:^10d}{13:^10d}{14:^10d}{15:^10d}{16:10d}{17:10d}{18:10d}{19:^15.8f}{20:^15.8f}{21:^10.5e}{22:^10.5e}{23:10d}\n'.format(-999,-9999,-9999,-9999,data_sum,data_sum_err,-9999,-9999,int(j),-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,ra,dec,-9999,-9999,numfit))
#
            clean()
#remove all trailing information
    outf.close()

if __name__ == "__main__":
    main(sys.argv[1:])
