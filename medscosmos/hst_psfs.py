import os
import numpy as np
import fitsio
import esutil as eu

class HSTPSF(object):
    def __init__(self, cat):
        self.cat=cat
        self._load_galsim_catalog()

    def get_psf(self, index):
        """
        get the psf image for the input index
        """
        gs_index = self.cat['gscosmos_index'][index]
        fname = os.path.join(
            self.psf_dir,
            str(self.gscat[gs_index]['psf_filename'],'utf-8').strip(),
        )
        return fitsio.read(fname)

    def _load_galsim_catalog(self):
        import galsim
        rgcat = galsim.RealGalaxyCatalog()
        self.psf_dir = rgcat.image_dir

        fname = rgcat.file_name
        print('loading galsim cat:',fname)
        self.gscat = fitsio.read(fname,lower=True)


class HSTPSFMatcher(object):
    def __init__(self, cat, maxrad):
        self.cat=cat
        self.maxrad=maxrad
        self._load_galsim_catalog()

    def get_catalog_with_matches(self):
        """
        get a version of the catalog with the match ids included
        """
        add_dt = [
            ('gscosmos_index','i8'),
            ('gscosmos_sep_arcsec','f8'),
        ]
        new_cat = eu.numpy_util.add_fields(self.cat, add_dt)
        new_cat['gscosmos_index'] = -9999
        new_cat['gscosmos_sep_arcsec'] = 9999
        new_cat['gscosmos_index'][self.m_cat] = self.m_gscat
        new_cat['gscosmos_sep_arcsec'][self.m_cat] = self.sep_arcsec

        return new_cat

    def match_catalogs(self):
        import smatch

        print('matching')

        # why does order matter?  We get proper matching when
        # putting the smaller catalog first
        matches = smatch.match(
            self.gscat['ra'],
            self.gscat['dec'],
            self.maxrad,
            self.cat['alpha_j2000'], 
            self.cat['delta_j2000'], 
            #self.cat['ra'], 
            #self.cat['dec'], 
            nside=1024,
        )

        #_, ind = np.unique(matches['i1'], return_index=True)
        _, ind = np.unique(matches['i2'], return_index=True)

        self.m_cat = matches['i2'][ind]
        self.m_gscat = matches['i1'][ind]
        cosdist = matches['cosdist'][ind]
        cosdist.clip(max=1.0, out=cosdist)

        dist_radians = np.arccos(cosdist)
        self.sep_arcsec = np.rad2deg(dist_radians)*3600.0
        
        print('matched %d/%d' % (self.m_cat.size, self.cat.size))


    def doplot(self, show=False):
        import biggles

        plt=biggles.plot(
            self.gscat['ra'],
            self.gscat['dec'],
            type='dot',
            visible=show,
            xlabel='RA',
            ylabel='DEC',
        )

        biggles.plot(
            self.cat['alpha_j2000'], 
            self.cat['delta_j2000'], 
            color='blue',
            type='dot',
            plt=plt,
            visible=False,
        )

        biggles.plot(
            self.cat['alpha_j2000'], 
            self.cat['delta_j2000'], 
            color='red',
            type='dot',
            plt=plt,
            visible=show,
        )

        return plt

    def _match_catalogs_alt(self):

        htm = eu.htm.HTM()
        m1, m2, d12 = htm.match(
            self.cat['ra'],
            self.cat['dec'],
            self.gscat['ra'],
            self.gscat['dec'],
            self.maxrad,
        )

        _, ind = np.unique(m1, return_index=True)

        self.m_cat = m1[ind]
        self.m_gscat = m2[ind]
        
        print('matched %d/%d' % (self.m_cat.size, self.cat.size))


    def _load_galsim_catalog(self):
        import galsim
        rgcat = galsim.RealGalaxyCatalog()
        self.psf_path = rgcat.image_dir

        fname = rgcat.file_name
        print('loading galsim cat:',fname)
        self.gscat = fitsio.read(fname,lower=True)



