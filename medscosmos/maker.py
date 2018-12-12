import numpy as np
import meds
import esutil as eu
import fitsio

from . import hst_psfs

FWHM_FAC = 2*np.sqrt(2*np.log(2))

class CosmosMEDSMaker(meds.MEDSMaker):
    def __init__(self, config_path, catname, flistname, **kw):

        config = eu.io.read(config_path)
        self.update(config)

        image_info = self._make_image_info(flistname)
        obj_data = self._make_obj_data(catname, image_info)

        super(CosmosMEDSMaker,self).__init__(
            obj_data,
            image_info,
            config=config,
            **kw
        )

    def _get_full_obj_data(self, obj_data):
        return obj_data

    def _build_meds_layout(self):
        """
        build the object data, filling in the stub we read

        note position offsets appear nowhere in this function
        """


        nim  = self.image_info.size
        nobj = self.obj_data.size

        trim_to_coadd = self.get('trim_to_coadd',False)
        if trim_to_coadd:
            print('    trimming to coadd')
            coadd_wcs, coadd_pos, coadd_bnds, coadd_q = \
                self._get_pos_and_bounds(self.obj_data, 0)
            in_bnds = coadd_bnds.contains_points(coadd_pos['zrow'], coadd_pos['zcol'])
            w_in_bnds, = np.where(in_bnds == True)
            assert w_in_bnds.size > 0,"none found in coadd"

            w_in_bnds = coadd_q[w_in_bnds]
            self.obj_data = self.obj_data[w_in_bnds]

        # we will want something else for DES data
        self.psf_data = hst_psfs.HSTPSF(cat=self.obj_data)

        # box sizes are even
        half_box_size = self.obj_data['box_size']//2

        for file_id in range(nim):

            wcs, pos, bnds, q = self._get_pos_and_bounds(self.obj_data, file_id)

            # do the test
            in_bnds = bnds.contains_points(pos['zrow'], pos['zcol'])
            q_rc, = np.where(in_bnds == True)
            print('    second cut: %6d of %6d objects' % (len(q_rc),len(q)))

            # now make sure everything is there
            if file_id == 0 and len(self.obj_data['ra']) != len(q_rc):
                raise MEDSCreationError('Not all objects were found in first image for '
                                        'MEDS making (which is the coadd/detection '
                                        'image by convention).')
            # compose them
            q = q[q_rc]

            # fill in the object_data structure

            # note q_rc since pos was created using obj_data[q]
            qrow = pos['zrow'][q_rc]
            qcol = pos['zcol'][q_rc]

            icut = self.obj_data['ncutout'][q]
            self.obj_data['file_id'][q,icut] = file_id
            self.obj_data['orig_row'][q,icut] = qrow
            self.obj_data['orig_col'][q,icut] = qcol

            # this results in the object center being close to
            # the natural center (dim-1.)/2.
            ostart_row = qrow.astype('i4') - half_box_size[q] + 1
            ostart_col = qcol.astype('i4') - half_box_size[q] + 1
            crow       = qrow - ostart_row
            ccol       = qcol - ostart_col

            self.obj_data['orig_start_row'][q,icut] = ostart_row
            self.obj_data['orig_start_col'][q,icut] = ostart_col
            self.obj_data['cutout_row'][q,icut]     = crow
            self.obj_data['cutout_col'][q,icut]     = ccol

            # do jacobian, in original, not-offset coords
            # note q_rc since pos was created using self.obj_data[q]
            jacob = wcs.get_jacobian(
                x=pos['wcs_col'][q_rc],
                y=pos['wcs_row'][q_rc])

            # jacob is a tuple of arrays
            self.obj_data['dudcol'][q,icut] = jacob[0]
            self.obj_data['dudrow'][q,icut] = jacob[1]
            self.obj_data['dvdcol'][q,icut] = jacob[2]
            self.obj_data['dvdrow'][q,icut] = jacob[3]

            # increment
            self.obj_data['ncutout'][q] += 1

        self.obj_data = self._make_resized_data(self.obj_data)
        #self.obj_data = self.obj_data[0:10]
        self._set_start_rows_and_pixel_count()

        self._set_psf_layout_hst()

    def _write_psf_cutouts(self):
        self._write_psf_cutouts_hst()

    def _write_psf_cutouts_hst(self):
        """
        write the cutouts for the specified type
        """

        print('writing psf cutouts')

        obj_data=self.obj_data
        psf_data=self.psf_data

        nfile=self.image_info.size
        nobj=obj_data.size

        cutout_hdu = self.fits['psf']

        for file_id in range(nfile):

            for iobj in range(nobj):
                if (iobj+1) % 100 == 0:
                    print('    %d/%d' % (iobj+1,obj_data.size))

                psf_im = self.psf_data.get_psf(iobj)

                ncut=obj_data['ncutout'][iobj]

                for icut in range(ncut):
                    expected_psf_shape = (
                        obj_data['psf_row_size'][iobj,icut],
                        obj_data['psf_col_size'][iobj,icut],
                    )

                    file_id = obj_data['file_id'][iobj, icut]

                    row = obj_data['orig_row'][iobj, icut]
                    col = obj_data['orig_col'][iobj, icut]
                    start_row = obj_data['psf_start_row'][iobj, icut]

                    if psf_im.shape != expected_psf_shape:
                        raise ValueError("psf size mismatch, expected %s "
                                         "got %s" % (expected_psf_shape, psf_im.shape))

                    cutout_hdu.write(psf_im, start=start_row)


    def _set_psf_layout_hst(self):
        """
        set the box sizes and start row for each psf image
        """

        print('setting psf layout for HST')
        obj_data=self.obj_data

        total_psf_pixels = 0
        psf_start_row = 0

        for iobj in range(obj_data.size):
            if (iobj+1) % 100 == 0:
                print('    %d/%d' % (iobj+1,obj_data.size))
            # note assuming same psf for all "epochs"
            # we would only have extra epoch to fake it
            psf_im = self.psf_data.get_psf(iobj)

            psf_shape = psf_im.shape
            psf_npix = psf_im.size

            box_size = max(psf_shape)
            cen = (np.array([box_size]*2)-1.0)/2.0

            # we will expand the psfs

            for icut in range(obj_data['ncutout'][iobj]):

                obj_data['psf_row_size'][iobj,icut] = psf_shape[0]
                obj_data['psf_col_size'][iobj,icut] = psf_shape[1]
                obj_data['psf_cutout_row'][iobj,icut] = cen[0]
                obj_data['psf_cutout_col'][iobj,icut] = cen[1]
                obj_data['psf_start_row'][iobj,icut] = psf_start_row

                psf_start_row += psf_npix
                total_psf_pixels += psf_npix

        self.total_psf_pixels = total_psf_pixels

    def _set_psf_layout_psfex(self):
        """
        set the box sizes and start row for each psf image
        """
        if self.psf_data is None:
            raise ValueError("_set_psf_layout called "
                             "with no psf data set")

        obj_data=self.obj_data
        psf_data=self.psf_data

        total_psf_pixels = 0

        #psf_npix = psf_size*psf_size

        psf_start_row = 0
        psf_shape=None
        for iobj in range(obj_data.size):
            for icut in range(obj_data['ncutout'][iobj]):

                row = obj_data['orig_row'][iobj, icut]
                col = obj_data['orig_col'][iobj, icut]
                file_id = obj_data['file_id'][iobj,icut]

                p = psf_data[file_id]

                pim = p.get_rec(row,col)
                cen = p.get_center(row,col)
                try:
                    sigma = p.get_sigma(row,col)
                except:
                    sigma = p.get_sigma()

                if psf_shape is None:
                    psf_shape = pim.shape
                    psf_npix = psf_shape[0]**2
                    obj_data['psf_box_size'] = psf_shape[0]
                else:
                    tpsf_shape = pim.shape
                    if tpsf_shape != psf_shape:
                        raise ValueError("currently all psfs "
                                         "must be same size")


                obj_data['psf_cutout_row'][iobj,icut] = cen[0]
                obj_data['psf_cutout_col'][iobj,icut] = cen[1]
                obj_data['psf_sigma'][iobj,icut] = sigma
                obj_data['psf_start_row'][iobj,icut] = psf_start_row

                psf_start_row += psf_npix
                total_psf_pixels += psf_npix


        self.total_psf_pixels = total_psf_pixels

    def _get_pos_and_bounds(self, obj_data, file_id):
        nim  = self.image_info.size
        impath=self.image_info['image_path'][file_id].strip()
        position_offset=self.image_info['position_offset'][file_id]

        print("file %4d of %4d: '%s'" % (file_id+1,nim,impath))

        wcs = self._get_wcs(file_id)

        # monkey patching in the position_offset into wcs
        wcs.position_offset=position_offset

        q = self._do_rough_sky_cut(wcs, obj_data['ra'], obj_data['dec'])
        print('    first cut:  %6d of %6d objects' % (q.size,obj_data.size))

        # this is the bottleneck
        pos = self._do_sky2image(wcs,
                                 obj_data['ra'][q],
                                 obj_data['dec'][q])

        # now test if in the actual image space.  Bounds are created
        # in the offset coords
        bnds = self._get_image_bounds(wcs)

        # for coadds add buffer if requested
        if file_id == 0:
            bnds.rowmin -= self['coadd_bounds_buffer_rowcol']
            bnds.rowmax += self['coadd_bounds_buffer_rowcol']
            bnds.colmin -= self['coadd_bounds_buffer_rowcol']
            bnds.colmax += self['coadd_bounds_buffer_rowcol']

        return wcs, pos, bnds, q


    def _read_catalog(self, catname):
        """
        cols=[
            'number',
            'flags',
            'alpha_j2000',
            'delta_j2000',
            'x_image',
            'y_image',
            'isoarea_image',
            'flux_auto',
            'fluxerr_auto',
            'flux_radius',

            # only used for cuts
            'mu_class',
            'mask',
        ]
        with fitsio.FITS(catname,lower=True) as fits:
            #cat = fits[1][cols][100000:110000]
            cat = fits[1].read(columns=cols)
        """
        print('loading catalog:',catname)
        with fitsio.FITS(catname,lower=True) as fits:
            #cat = fits[1][100000:110000]
            cat = fits[1][:]

        w, = np.where(
            (cat['mu_class'] < 3)
            &
            (cat['mask']==0)
            &
            (cat['gscosmos_index'] >= 0)
        )
        print('initial cuts %d/%d %g%%' % (w.size,cat.size,w.size/cat.size*100))

        cat = cat[w]
        return cat

    def _make_obj_data(self, catname, image_info):

        cat = self._read_catalog(catname)

        #ncutout_max=2
        ncutout_max = image_info.size
        if ncutout_max < 2:
            ncutout_max=2
        obj_data = meds.util.get_meds_output_struct(
            cat.size,
            ncutout_max,
            extra_fields=self._get_fields(ncutout_max),
        )

        obj_data['box_size'] = self._get_box_sizes(cat)

        obj_data['id']       = cat[ self['id_name'] ]
        obj_data['number']   = cat['number']
        obj_data['flags']    = cat[ self['flags_name'] ]
        obj_data['ra']       = cat[ self['ra_name'] ]
        obj_data['dec']      = cat[ self['dec_name'] ]
        obj_data['flux']     = cat[ self['flux_name'] ]
        obj_data['flux_err'] = cat[ self['fluxerr_name'] ]

        obj_data['gscosmos_index'] = cat['gscosmos_index']
        obj_data['gscosmos_sep_arcsec'] = cat['gscosmos_sep_arcsec']

        obj_data['iso_radius']  = self._get_iso_radius(cat)

        pos=meds.util.make_wcs_positions(
            cat[ self['row_name'] ],
            cat[ self['col_name'] ],
            self['position_offset'],
        )
        obj_data['input_row'] = pos['zrow']
        obj_data['input_col'] = pos['zcol']

        return obj_data

    def _make_resized_data(self, odata):
        """
        make a new struct with ncutout-sized-arrays based on
        the actual maximum ncutout
        """


        nmax = odata['file_id'].shape[1]
        new_nmax = odata['ncutout'].max()
        if new_nmax < 2:
            new_nmax = 2
        temp_obj_data = odata

        nobj = temp_obj_data.size

        new_data = meds.util.get_meds_output_struct(
            nobj,
            new_nmax,
            extra_fields=self._get_fields(new_nmax),
        )

        tmpst = meds.util.get_meds_output_struct(1, new_nmax)
        required_fields = tmpst.dtype.names

        for name in new_data.dtype.names:
            if name in temp_obj_data.dtype.names:

                shape = new_data[name].shape
                lshape = len(shape)
                #if lshape > 1 and name in required_fields:
                #    new_data[name][:,:] = temp_obj_data[name][:,0:new_nmax]
                #else:
                #    new_data[name][:] = temp_obj_data[name][:]
                if lshape > 1 and shape[1] == new_nmax:
                    new_data[name][:,:] = temp_obj_data[name][:,0:new_nmax]
                else:
                    new_data[name][:] = temp_obj_data[name][:]


        del temp_obj_data

        return new_data


    def _get_iso_radius(self, cat):
        iso_area = cat[self['isoarea_name']].clip(min=1)
        return np.sqrt(iso_area/np.pi)

    def _get_extra_fields(self, obj_data, nmax):
        return []

    def _get_fields(self, ncut):
        return [
            ('number','i8'),
            ('flags','i4'),
            ('iso_radius','f4'),
            ('flux','f4'),
            ('flux_err','f4'),
            ('input_row','f8'),
            ('input_col','f8'),
            ('gscosmos_index','i8'),
            ('gscosmos_sep_arcsec','f8'),
            ('psf_row_size','i4',ncut),
            ('psf_col_size','i4',ncut),
            ('psf_cutout_row','f8',ncut),
            ('psf_cutout_col','f8',ncut),
            ('psf_start_row','i8',ncut),
        ]


    def _get_box_sizes(self, cat):
        """
        get box sizes that are wither 2**N or 3*2**N, within
        the limits set by the user
        """

        box_size = self._get_sigma_size(cat)

        # clip to range
        box_size = box_size.clip(
            min=self['min_box_size'],
            max=self['max_box_size'],
        )

        # now put in fft sizes
        bins = [0]

        bins.extend([sze for sze in self['allowed_box_sizes'] 
                     if sze >= self['min_box_size']
                     and sze <= self['max_box_size']])

        if bins[-1] != self['max_box_size']:
            bins.append(self['max_box_size'])

        bin_inds = np.digitize(box_size,bins,right=True)
        bins = np.array(bins)

        return bins[bin_inds]


    def _get_sigma_size(self, cat):
        """
        "sigma" size, based on flux radius.  There are no ellip
        parameters in catalog
        """

        # note taking element 1 which is half light radius
        sigma = cat['flux_radius'][:,1]*2.0/FWHM_FAC
        drad = sigma*self['sigma_fac']
        drad = np.ceil(drad)
        sigma_size = 2*drad.astype('i4') # sigma size is twice the radius

        return sigma_size


    def _make_image_info(self, flistname):
        """
        won't load any data yet because the files are gzipped and just reading
        the header takes 2.6 G and a long time!

        This means we need to set magzp and scale later when we read
        """

        flist=[]
        with open(flistname) as fobj:
            for line in fobj:
                fname=line.strip()
                flist.append(fname)

        nimage = len(flist)

        path_len = max([len(f) for f in flist])

        try:
            ext_len = len(self['image_ext'])
        except:
            ext_len=None

        image_info = meds.util.get_image_info_struct(
            nimage,
            path_len,
            ext_len=ext_len,
        )
        image_info['position_offset'] = 1
        image_info['image_ext'] = self['image_ext']
        image_info['weight_ext'] = self['weight_ext']

        for i,f in enumerate(flist):
            image_info['image_id'] = i
            image_info['image_path'] = f
            image_info['weight_path'] = f.replace('sci.fits','wht.fits')

        return image_info


