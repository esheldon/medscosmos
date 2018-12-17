import numpy as np
import fitsio
import esutil as eu
from . import hst_psfs

FWHM_FAC = 2*np.sqrt(2*np.log(2))

COSMOS_PIXEL_SCALE = 0.03 # arcsec/pixel

class CatMaker(object):
    """
    make the catalog, calculating the box sizes in arcsec
    and matching to the PSFs
    """
    def __init__(self, cat_path, sigma_fac, maxrad):
        self.cat_path = cat_path
        self.sigma_fac = sigma_fac
        self.maxrad = maxrad

        self._set_cat()
        self._set_box_sizes()
        self._set_iso_radius()
        self._match_psf()

    def get_cat(self):
        """
        get the catalog
        """
        return self.matched_cat

    def get_plot(self):
        """
        get the catalog
        """
        return self.plt

    def _set_box_sizes(self):
        """
        get box sizes, assuming pixel scale of 0.03''
        """
        print('setting box sizes in arcsec')

        # element 1 is half light
        fwhm_arcsec = 2*self.cat['flux_radius'][:,1]*COSMOS_PIXEL_SCALE

        sigma_arcsec = fwhm_arcsec/FWHM_FAC
        radius_arcsec = sigma_arcsec*self.sigma_fac

        # box size is twice radius size
        box_size_arcsec = radius_arcsec * 2

        self.cat['box_size_arcsec'] = box_size_arcsec

    def _set_iso_radius(self):
        """
        set iso radius based on isoarea_image
        """
        print('setting iso_radius_arcsec')
        iso_area_pixels = self.cat['isoarea_image'].clip(min=1)
        iso_radius_pixels = np.sqrt(iso_area_pixels/np.pi)
        self.cat['iso_radius_arcsec'] = iso_radius_pixels * COSMOS_PIXEL_SCALE

    def _match_psf(self):
        """
        match to galsim cosmos PSFs
        """

        matcher = hst_psfs.HSTPSFMatcher(
            self.cat,
            self.maxrad,
        )
        matcher.match_catalogs()

        self.plt=matcher.doplot()

        self.matched_cat = matcher.get_catalog_with_matches()

    def _set_cat(self):
        """
        read the catalog, adding some extra fields
        """
        self.cat = read_original_cat(self.cat_path)


def read_original_cat(fname):
    cols=[
        'number',
        'flags',
        'alpha_j2000',
        'delta_j2000',
        'x_image',
        'y_image',
        'mag_iso',
        'isoarea_image',
        'mag_auto',
        'flux_auto',
        'fluxerr_auto',
        'flux_radius',

        # only used for cuts
        'unique',
        'nearstar',
        'masked',
        'mu_class',
        'mask',
    ]

    print('reading:',fname)
    data = fitsio.read(
        fname,
        columns=cols,
        lower=True,
    )

    # clean with blends
    w , = np.where(
        (data['unique'] == 1)
        & (data['nearstar'] == 1)
        & (data['masked'] == 1)
        & (data['mu_class'] < 3)
        & (data['mask'] == 0)
        & (data['mag_iso'] < 26.5)
    )
    print('after cuts %d/%d %g%%' % (w.size,data.size,w.size/data.size*100))
    data=data[w]

    newnames=[]
    for name in data.dtype.names:
        if name=='alpha_j2000':
            name = 'ra'
        elif name=='delta_j2000':
            name = 'dec'
        newnames.append(name)

    data.dtype.names = newnames

    add_dt = [
        ('id','i8'),
        ('box_size_arcsec','f4'),
        ('iso_radius_arcsec','f4'),
    ]
    new_data = eu.numpy_util.add_fields(data, add_dt)
    new_data['id'] = np.arange(new_data.size)
    return new_data

