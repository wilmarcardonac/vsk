import numpy as np
import pylab as py
import healpy as hp

vmap = hp.read_map('../zero_mean_vsk_maps/vmap_smica.fits')

smap = hp.read_map('../zero_mean_vsk_maps/smap_smica.fits')

kmap = hp.read_map('../zero_mean_vsk_maps/kmap_smica.fits')

py.figure(1)

hp.mollview(vmap,unit = r'$\mu K^2$', title = 'V-MAP INPAINTED SMICA')

py.savefig('vmap-inpainted-smica.pdf')

py.close(1)

py.figure(2)

hp.mollview(smap, title = 'S-MAP INPAINTED SMICA')

py.savefig('smap-inpainted-smica.pdf')

py.close(2)

py.figure(3)

hp.mollview(kmap, title = 'K-MAP INPAINTED SMICA')

py.savefig('kmap-inpainted-smica.pdf')

py.close(3)

exit()
