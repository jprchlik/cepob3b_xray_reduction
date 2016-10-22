# run through CIAO’s srcflux to extract fluxes for the wav detect  results
#
# taken and edited from CIAO website:


# Use roi to create source regions with background regions and no overlaps:
	mkdir 10809/sources
	punlearn roi
 	pset roi group=exclude 
#	pset roi bkgfunction=mul bkgfactor=3
	pset roi radiusmode=mul  bkgradius=5
	pset roi outsrcfile=10809/sources/src%03d.fits
	pset roi targetbkg=target
#        pset roi st
	pset roi evtfile=10809/repro/acisf10809_repro_evt2.fits
        pset roi clobber=yes
	pset roi fovregion=10809/repro/acisf10809_repro_fov1.fits	# would not parse fov1 files from retro. or remade with skyfov
        pset roi infile=10809/obsid_10809_exp_src.fits 
	roi 10809/obsid_10809_exp_src.fits  

#  Create a list of source and background regions from the output of roi:
	splitroi 10809/sources/src\*fits 10809/sources/wav
#
#
#  Run srcflux to extract the flux
	mkdir 10809/srcflux_output
	punlearn srcflux
	pset srcflux bands=0.5:8:1.7
	pset srcflux model=xsraymond.ray1 
	pset srcflux paramvals=ray1.kT=1.0
	pset srcflux absmodel=xsphabs.abs1  
	pset srcflux absparams=abs1.nH=0.0 
	pset srcflux psfmethod=arfcorr
        pset srcflux clobber=yes
        pset srcflux parallel=yes
        pset srcflux fovfile=10809/repro/acisf10809_repro_fov1.fits
        pset srcflux psffile=10809/ob10809_psf.fits
        pset srcflux infile=10809/repro/acisf10809_repro_evt2.fits
        pset srcflux outroot=10809/srcflux_output/sers
        pset srcflux pos=10809/obsid_10809_exp_src.fits
	pset srcflux srcreg=@-10809/sources/wav.src.reg 
	pset srcflux bkgreg=@-10809/sources/wav.bg.reg 
	srcflux



