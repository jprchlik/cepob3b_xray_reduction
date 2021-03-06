#ls 9919/srcflux_output/sers_*bkg.pi > bkg9919.lis
#combine_spectra @bkg9919.lis 9919/summed_background_9919
#
#ls 9920/srcflux_output/sers_*bkg.pi > bkg9920.lis
#combine_spectra @bkg9920.lis 9920/summed_background_9920

#ls 10809/srcflux_output/sers_*bkg.pi > bkg10809.lis
#combine_spectra @bkg10809.lis 10809/summed_background_10809
#
#ls 10810/srcflux_output/sers_*bkg.pi > bkg10810.lis
#combine_spectra @bkg10810.lis 10810/summed_background_10810
#
#ls 10811/srcflux_output/sers_*bkg.pi > bkg10811.lis
#combine_spectra @bkg10811.lis 10811/summed_background_10811
#
#ls 10812/srcflux_output/sers_*bkg.pi > bkg10812.lis
#combine_spectra @bkg10812.lis 10812/summed_background_10812

ls 9919/background/*[0-4].pi > bkg_ext9919.lis
combine_spectra @bkg_ext9919.lis 9919/background/summed_background_9919

ls 9920/background/*[0-4].pi > bkg_ext9920.lis
combine_spectra @bkg_ext9920.lis 9920/background/summed_background_9920

ls 10809/background/*[0-4].pi > bkg_ext10809.lis
combine_spectra @bkg_ext10809.lis 10809/background/summed_background_10809

ls 10810/background/*[0-4].pi > bkg_ext10810.lis
combine_spectra @bkg_ext10810.lis 10810/background/summed_background_10810

ls 10811/background/*[0-4].pi > bkg_ext10811.lis
combine_spectra @bkg_ext10811.lis 10811/background/summed_background_10811

ls 10812/background/*[0-4].pi > bkg_ext10812.lis
combine_spectra @bkg_ext10812.lis 10812/background/summed_background_10812
