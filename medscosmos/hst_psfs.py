import numpy as np
import fitsio
import esutil as eu

class HSTPSF(object):
    def __init__(self, cat, maxrad, nside):
        self.cat=cat
        self.maxrad=maxrad
        self.nside=nside
        self._load_galsim_catalog()
        self._match_catalogs()

    def doplot(self):
        import biggles
        plt=biggles.plot(
            self.gscat['ra'],
            self.gscat['dec'],
            color='yellow',
            type='dot',
            visible=True,
        )
        biggles.plot(
            self.cat['ra'],
            self.cat['dec'],
            color='orange',
            type='dot',
            plt=plt,
            visible=False,
        )

        biggles.plot(
            self.cat['ra'][self.m_cat],
            self.cat['dec'][self.m_cat],
            color='red',
            type='dot',
            plt=plt,
        )



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

    def _match_catalogs(self):
        import smatch

        matches = smatch.match(
            self.gscat['ra'],
            self.gscat['dec'],
            self.maxrad,
            self.cat['ra'], 
            self.cat['dec'], 
            nside=self.nside,
        )

        #_, ind = np.unique(matches['i1'], return_index=True)
        _, ind = np.unique(matches['i2'], return_index=True)

        self.m_cat = matches['i2'][ind]
        self.m_gscat = matches['i1'][ind]

        
        print('matched %d/%d' % (self.m_cat.size, self.cat.size))


    def _load_galsim_catalog(self):
        import galsim
        rgcat = galsim.RealGalaxyCatalog()
        self.psf_path = rgcat.image_dir

        fname = rgcat.file_name
        print('loading galsim cat:',fname)
        self.gscat = fitsio.read(fname,lower=True)



