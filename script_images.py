import plot_chandra_image as pci
from multiprocessing import Pool

def make_image(epoch):
    chan = pci.format_image(epoch)
    
    chan.make_plot(starreg=True,backreg=True)

    

epochs = [9919,9920,10809,10810,10811,10812]

pool = Pool(processes=4)
out = pool.map(make_image,epochs)
out.close()

