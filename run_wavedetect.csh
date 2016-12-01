#NEED TO RUN dmcopy and mkpsfmap in xgterm
#       /bin/bash
#	dmcopy “flux/flux_0.5-8.0.img[bin x=2580:6605:1.0,y=1832:5645:1.0]” obs11013_img.fits clob+
#	dmcopy flux/flux_0.5-8.0.img obs11013_img.fits clob+
#	mkpsfmap obs11013_img.fits obs11013_psf.fits 1.7 ecf=(0.15)(0.9)(0.95)

#	dmcopy “10809/flux_0.5-8.0.img[bin x=2580:6605:1.0,y=1832:5645:1.0]” 10809/obs10809_img.fits clob+
#	dmcopy 10809/flux_0.5-8.0.img 10809/obs10809_img.fits clob+
#	mkpsfmap 10809/obs10809_img.fits 10809/obs10809_psf.fits 1.7 ecf=0.9 clob+
#        
####	punlearn wavdetect
####	pset wavdetect infile=10809/obs10809_img.fits
####	pset wavdetect outfile=10809/obsid_10809_exp_src.fits
####	pset wavdetect scellfile=10809/obsid_10809_exp_scell.fits
####	pset wavdetect imagefile=10809/obsid_10809_exp_img.fits
####	pset wavdetect defnbkgfile=10809/obsid_10809_exp_nbgd.fits
####	pset wavdetect regfile=10809/obsid_10809_exp_src.reg
####	pset wavdetect expfile=10809/flux_0.5-8.0.expmap
####	pset wavdetect psffile=10809/obs10809_psf.fits 
####	pset wavdetect scales=1.0,1.414,2.0,2.828,4.0,5.657,8.0,11.314,16.0
####	pset wavdetect maxiter=14
####	pset wavdetect iterstop=0.00006
####	pset wavdetect sigthresh=7e-5
####	pset wavdetect expthresh=0.1
####        pset wavdetect clobber=yes
#####	plist wavdetect | grep -v # > wavdetect_merged.10809_options.txt
####	wavdetect
#####
#####	dmcopy “10810/flux_0.5-8.0.img[bin x=2580:6605:1.0,y=1832:5645:1.0]” 10810/obs10810_img.fits clob+
#####	dmcopy 10810/flux_0.5-8.0.img 10810/obs10810_img.fits clob+
#####	mkpsfmap 10810/obs10810_img.fits 10810/obs10810_psf.fits 1.7 ecf=0.9
#####        
####	punlearn wavdetect
####	pset wavdetect infile=10810/obs10810_img.fits
####	pset wavdetect outfile=10810/obsid_10810_exp_src.fits
####	pset wavdetect scellfile=10810/obsid_10810_exp_scell.fits
####	pset wavdetect imagefile=10810/obsid_10810_exp_img.fits
####	pset wavdetect defnbkgfile=10810/obsid_10810_exp_nbgd.fits
####	pset wavdetect regfile=10810/obsid_10810_exp_src.reg
####	pset wavdetect expfile=10810/flux_0.5-8.0.expmap
####	pset wavdetect psffile=10810/obs10810_psf.fits 
####	pset wavdetect scales=1.0,1.414,2.0,2.828,4.0,5.657,8.0,11.314,16.0
####	pset wavdetect maxiter=14
####	pset wavdetect iterstop=0.00006
####	pset wavdetect sigthresh=7e-5
####	pset wavdetect expthresh=0.1
####        pset wavdetect clobber=yes
#####	plist wavdetect | grep -v # > wavdetect_merged.10810_options.txt
####	wavdetect
#####
####	dmcopy “10811/flux_0.5-8.0.img[bin x=2580:6605:1.0,y=1832:5645:1.0]” 10811/obs10811_img.fits clob+
######	dmcopy 10811/flux_0.5-8.0.img 10811/obs10811_img.fits clob+
####	mkpsfmap 10811/obs10811_img.fits 10811/obs10811_psf.fits 1.7 ecf=0.9
####        
####	punlearn wavdetect
####	pset wavdetect infile=10811/obs10811_img.fits
####	pset wavdetect outfile=10811/obsid_10811_exp_src.fits
####	pset wavdetect scellfile=10811/obsid_10811_exp_scell.fits
####	pset wavdetect imagefile=10811/obsid_10811_exp_img.fits
####	pset wavdetect defnbkgfile=10811/obsid_10811_exp_nbgd.fits
####	pset wavdetect regfile=10811/obsid_10811_exp_src.reg
####	pset wavdetect expfile=10811/flux_0.5-8.0.expmap
####	pset wavdetect psffile=10811/obs10811_psf.fits 
####	pset wavdetect scales=1.0,1.414,2.0,2.828,4.0,5.657,8.0,11.314,16.0
####	pset wavdetect maxiter=14
####	pset wavdetect iterstop=0.00006
####	pset wavdetect sigthresh=7e-5
####	pset wavdetect expthresh=0.1
####        pset wavdetect clobber=yes
#####	plist wavdetect | grep -v # > wavdetect_merged.10811_options.txt
####	wavdetect
####
#####	dmcopy “9920/flux_0.5-8.0.img[bin x=2580:6605:1.0,y=1832:5645:1.0]” 9920/obs9920_img.fits clob+
#####	dmcopy 9920/flux_0.5-8.0.img 9920/obs9920_img.fits clob+
#####	mkpsfmap 9920/obs9920_img.fits 9920/obs9920_psf.fits 1.7 ecf=0.9
#####        
	punlearn wavdetect
	pset wavdetect infile=9920/obs9920_img.fits
	pset wavdetect outfile=9920/obsid_09920_exp_src.fits
	pset wavdetect scellfile=9920/obsid_09920_exp_scell.fits
	pset wavdetect imagefile=9920/obsid_09920_exp_img.fits
	pset wavdetect defnbkgfile=9920/obsid_09920_exp_nbgd.fits
	pset wavdetect regfile=9920/obsid_09920_exp_src.reg
	pset wavdetect expfile=9920/flux_0.5-8.0.expmap
	pset wavdetect psffile=9920/obs9920_psf.fits 
	pset wavdetect scales=1.0,1.414,2.0,2.828,4.0,5.657,8.0,11.314,16.0
	pset wavdetect maxiter=14
	pset wavdetect iterstop=0.00006
	pset wavdetect sigthresh=7e-5
	pset wavdetect expthresh=0.10
        pset wavdetect clobber=yes
#####	plist wavdetect | grep -v # > wavdetect_merged.9920_options.txt
	wavdetect
#####
#####	dmcopy “10812/flux_0.5-8.0.img[bin x=2580:6605:1.0,y=1832:5645:1.0]” 10812/obs10812_img.fits clob+
#####	dmcopy 10812/flux_0.5-8.0.img 10812/obs10812_img.fits clob+
#####	mkpsfmap 10812/obs10812_img.fits 10812/obs10812_psf.fits 1.7 ecf=0.9
####        
####	punlearn wavdetect
####	pset wavdetect infile=10812/obs10812_img.fits
####	pset wavdetect outfile=10812/obsid_10812_exp_src.fits
####	pset wavdetect scellfile=10812/obsid_10812_exp_scell.fits
####	pset wavdetect imagefile=10812/obsid_10812_exp_img.fits
####	pset wavdetect defnbkgfile=10812/obsid_10812_exp_nbgd.fits
####	pset wavdetect regfile=10812/obsid_10812_exp_src.reg
####	pset wavdetect expfile=10812/flux_0.5-8.0.expmap
####	pset wavdetect psffile=10812/obs10812_psf.fits 
####	pset wavdetect scales=1.0,1.414,2.0,2.828,4.0,5.657,8.0,11.314,16.0
####	pset wavdetect maxiter=14
####	pset wavdetect iterstop=0.00006
####	pset wavdetect sigthresh=7e-5
####	pset wavdetect expthresh=0.1
####        pset wavdetect clobber=yes
#####	plist wavdetect | grep -v # > wavdetect_merged.10812_options.txt
####	wavdetect
####
#####	dmcopy “9919/flux_0.5-8.0.img[bin x=2580:6605:1.0,y=1832:5645:1.0]” 9919/obs9919_img.fits clob+
#####	dmcopy 9919/flux_0.5-8.0.img 9919/obs9919_img.fits clob+
#####	mkpsfmap 9919/obs9919_img.fits 9919/obs9919_psf.fits 1.7 ecf=0.9
####        
####	punlearn wavdetect
####	pset wavdetect infile=9919/obs9919_img.fits
####	pset wavdetect outfile=9919/obsid_09919_exp_src.fits
####	pset wavdetect scellfile=9919/obsid_09919_exp_scell.fits
####	pset wavdetect imagefile=9919/obsid_09919_exp_img.fits
####	pset wavdetect defnbkgfile=9919/obsid_09919_exp_nbgd.fits
####	pset wavdetect regfile=9919/obsid_09919_exp_src.reg
####	pset wavdetect expfile=9919/flux_0.5-8.0.expmap
####	pset wavdetect psffile=9919/obs9919_psf.fits 
####	pset wavdetect scales=1.0,1.414,2.0,2.828,4.0,5.657,8.0,11.314,16.0
####	pset wavdetect maxiter=14
####	pset wavdetect iterstop=0.00006
####	pset wavdetect sigthresh=7e-5
####	pset wavdetect expthresh=0.1
####        pset wavdetect clobber=yes
#####	plist wavdetect | grep -v # > wavdetect_merged.9919_options.txt
####	wavdetect
