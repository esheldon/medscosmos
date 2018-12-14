"""
TODO:

    - move cuts into catmaker
    - use clean not clean only and apply cuts eric sent
    clean_withblends = (cosmos[‘unique’] == 1) & (cosmos[‘nearstar’] == 1) & (cosmos[‘masked’]==1)
    - put in iso_radius
    - make MEDS box size stuff pre-calculated in arcsec
    - we hope the min box size will make this OK for the DES data too
    - force to common zero point
    - limit to 1024 postage stamp (min 32)
    - just make stamp sizes even, don't pick particular ones

    - implement fake seg map?
        - just a circle outside of iso_radius

"""
import os
import numpy as np
import meds
import esutil as eu
import fitsio

from . import hst_psfs

from . import files
from .files import (
    TempFile,
    StagedOutFile,
)



class CosmosMEDSMaker(meds.MEDSMaker):
    def __init__(self, config_path, catname, flistname, **kw):

        config = eu.io.read(config_path)
        self.update(config)

        image_info = self._make_image_info(flistname)
        
        self.cat_orig = self._read_catalog(catname)
        obj_data = self._make_obj_data(image_info)

        self._setup_fpack()

        if hasattr(self,'psfex_objects'):
            kw['psf_data'] = self.psfex_objects

        super(CosmosMEDSMaker,self).__init__(
            obj_data,
            image_info,
            config=config,
            **kw
        )

    def write(self, filename):
        """
        write compressed meds file

        Parameters
        ----------
        filename: string
            Must end in .fz
        """
        assert filename[-3:]=='.fz','name must end in .fz'

        ucfilename = filename[0:-3]

        with TempFile(ucfilename) as tfile:
            super(CosmosMEDSMaker,self).write(tfile.path)
            self._compress_meds_file(tfile.path, filename)


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

        self._do_psf_setup()

        # box sizes are even
        half_box_size = self.obj_data['box_size']//2

        for file_id in range(nim):

            wcs, pos, bnds, q = self._get_pos_and_bounds(self.obj_data, file_id)

            # do the test
            in_bnds = bnds.contains_points(pos['zrow'], pos['zcol'])
            q_rc, = np.where(in_bnds == True)
            print('    second cut: %6d of %6d objects' % (len(q_rc),len(q)))

            # now make sure everything is there
            if self['check_in_first_image']:
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
        print('setting number field as sequential')
        self.obj_data['number'] = 1+np.arange(self.obj_data.size)

        self._set_start_rows_and_pixel_count()

        if self['survey']=='cosmos':
            self._set_psf_layout_hst()
        else:
            self._set_psf_layout_psfex()


    def _do_psf_setup(self):
        if self['survey']=='cosmos':
            self.psf_data = hst_psfs.HSTPSF(cat=self.obj_data)

    def _set_psfex_objects(self, image_info):
        import psfex
        self.psfex_objects=[]
        for psfex_file in image_info['psfex_path']:
            print('loading:',psfex_file)
            p = psfex.PSFEx(psfex_file)
            self.psfex_objects.append(p)
 
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

                # HST psf is same for every cutout, in fact ncut should always
                # be 1
                try:
                    psf_im = self.psf_data.get_psf(iobj)
                except AttributeError:
                    psf_im = None

                ncut=obj_data['ncutout'][iobj]

                for icut in range(ncut):

                    if psf_im is None:
                        row = obj_data['orig_row'][iobj, icut]
                        col = obj_data['orig_col'][iobj, icut]
                        file_id = obj_data['file_id'][iobj,icut]

                        p = self.psf_data[file_id]

                        psf_im = p.get_rec(row,col)

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

        print('setting psf layout for PSFEx')

        obj_data=self.obj_data
        psf_data=self.psf_data

        total_psf_pixels = 0

        #psf_npix = psf_size*psf_size

        psf_start_row = 0
        for iobj in range(obj_data.size):
            for icut in range(obj_data['ncutout'][iobj]):

                row = obj_data['orig_row'][iobj, icut]
                col = obj_data['orig_col'][iobj, icut]
                file_id = obj_data['file_id'][iobj,icut]

                p = psf_data[file_id]

                pim = p.get_rec(row,col)
                cen = p.get_center(row,col)

                psf_shape = pim.shape
                psf_npix = pim.size

                obj_data['psf_row_size'][iobj,icut] = psf_shape[0]
                obj_data['psf_col_size'][iobj,icut] = psf_shape[1]
                obj_data['psf_cutout_row'][iobj,icut] = cen[0]
                obj_data['psf_cutout_col'][iobj,icut] = cen[1]
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
        read the cosmos catalog
        """
        print('loading catalog:',catname)
        with fitsio.FITS(catname,lower=True) as fits:
            #cat = fits[1][100000:110000]
            if 'object_data' in fits:
                print('reading from MEDS object data')
                ext='object_data'
            else:
                ext=1
            cat = fits[ext][:]

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

    def _add_cat_fields(self, odata, copy=True):
        """
        add fields from the cat

        some will not be in the odata but some will. When copy is True We will
        copy over the ones that are in both, in some cases
        """
        # these are required fileds from get_meds_output_dtype
        # that we have put into the input catalog
        always_copy=[
            'id',
            'ra',
            'dec',
        ]
        cat = self.cat_orig

        add_dt = []
        for d in cat.dtype.descr:
            n = d[0]
            if n not in odata.dtype.names:
                add_dt.append(d)

        obj_data = eu.numpy_util.add_fields(
            odata,
            add_dt,
        )

        if copy:
            for n in always_copy:
                obj_data[n] = cat[n]

            for d in add_dt:
                n = d[0]
                if n in always_copy:
                    continue

                # don't clobber things that should be left at
                # their default values
                if n not in odata.dtype.names:
                    obj_data[n] = cat[n]


        return obj_data

    def _make_obj_data(self, image_info):

        cat = self.cat_orig

        #ncutout_max=2
        ncutout_max = image_info.size
        if ncutout_max < 2:
            ncutout_max=2

        obj_data = meds.util.get_meds_output_struct(
            cat.size,
            ncutout_max,
            extra_fields=self._get_fields(ncutout_max),
        )
        obj_data = self._add_cat_fields(obj_data)

        # convert box size in arcsec to pixels
        obj_data['box_size'] = self._get_box_sizes(image_info, cat)


        pos=meds.util.make_wcs_positions(
            cat[ self['row_name'] ],
            cat[ self['col_name'] ],
            self['position_offset'],
        )
        obj_data['input_row'] = pos['zrow']
        obj_data['input_col'] = pos['zcol']

        return obj_data

    """
        obj_data['id']       = cat[ self['id_name'] ]
        obj_data['number']   = cat['number']
        obj_data['flags']    = cat[ self['flags_name'] ]
        obj_data['ra']       = cat[ self['ra_name'] ]
        obj_data['dec']      = cat[ self['dec_name'] ]
        obj_data['flux']     = cat[ self['flux_name'] ]
        obj_data['flux_err'] = cat[ self['fluxerr_name'] ]

        obj_data['gscosmos_index'] = cat['gscosmos_index']
        obj_data['gscosmos_sep_arcsec'] = cat['gscosmos_sep_arcsec']

        obj_data['box_size_arcsec'] = cat['box_size_arcsec']
    """

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
        new_data = self._add_cat_fields(new_data, copy=False)

        for name in new_data.dtype.names:
            if name in temp_obj_data.dtype.names:

                shape = new_data[name].shape
                lshape = len(shape)

                if lshape > 1 and shape[1] == new_nmax:
                    new_data[name][:,:] = temp_obj_data[name][:,0:new_nmax]
                else:
                    new_data[name][:] = temp_obj_data[name][:]

        del temp_obj_data

        return new_data

    def _get_extra_fields(self, obj_data, nmax):
        return []

    def _get_fields(self, ncut):
        return [
            #('number','i8'),
            #('flags','i4'),
            #('iso_radius','f4'),
            #('flux','f4'),
            #('flux_err','f4'),
            ('input_row','f8'),
            ('input_col','f8'),
            #('gscosmos_index','i8'),
            #('gscosmos_sep_arcsec','f8'),
            ('psf_row_size','i4',ncut),
            ('psf_col_size','i4',ncut),
            ('psf_cutout_row','f8',ncut),
            ('psf_cutout_col','f8',ncut),
            ('psf_start_row','i8',ncut),
            #('box_size_arcsec','f4'),
        ]


    def _get_box_sizes(self, image_info, cat):
        """
        get box sizes that are wither 2**N or 3*2**N, within
        the limits set by the user
        """


        file_id=0
        impath=image_info['image_path'][file_id].strip()
        ext=image_info['image_ext'][file_id]
        wcs_data = fitsio.read_header(impath, ext=ext)
        wcs = eu.wcsutil.WCS(wcs_data)


        jacob = wcs.get_jacobian(100,100)
        dudcol, dudrow, dvdcol, dvdrow = jacob

        det = dvdrow*dudcol - dvdcol*dudrow
        pixel_scale = np.sqrt(abs(det))
        print('found pixel scale:',pixel_scale)
        box_size = cat['box_size_arcsec']/pixel_scale

        # clip to range
        box_size.clip(
            min=self['min_box_size'],
            max=self['max_box_size'],
            out=box_size,
        )
        box_size = box_size.astype('i4')

        w,=np.where( ( (box_size % 2) != 0 ) )
        if w.size > 0:
            box_size[w] += 1

        return box_size

    def _make_image_info(self, flistname):
        survey=self['survey'].lower()

        if survey=='cosmos':
            image_info = self._make_image_info_hst(flistname)
        elif survey=='des':
            image_info = self._make_image_info_des(flistname)
            self._set_psfex_objects(image_info)
        else:
            raise ValueError('bad survey: %s' % self['survey'])

        return image_info

    def _make_image_info_hst(self, flistname):
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

        #image_info = meds.util.get_image_info_struct(
        image_info = get_image_info_struct(
            nimage,
            path_len,
            ext_len=ext_len,
        )
        image_info['position_offset'] = 1
        image_info['image_ext'] = self['image_ext']
        image_info['weight_ext'] = self['weight_ext']

        for i,f in enumerate(flist):
            image_info['image_id'][i] = i
            image_info['image_path'][i] = f
            image_info['weight_path'][i] = f.replace('sci.fits','wht.fits')

        return image_info

    def _make_image_info_des(self, flistname):
        """
        won't load any data yet because the files are gzipped and just reading
        the header takes 2.6 G and a long time!

        This means we need to set magzp and scale later when we read
        """

        flist=[]
        psfex_flist=[]
        with open(flistname) as fobj:
            for line in fobj:
                fname = line.strip()
                flist.append(fname)

                psfex_fname = fname.replace('.fits.fz','_psfcat.psf')
                psfex_flist.append(psfex_fname)

        nimage = len(flist)

        path_len = max([len(f) for f in flist])
        psfex_path_len = max([len(f) for f in psfex_flist])

        try:
            ext_len = len(self['image_ext'])
        except:
            ext_len=None

        extra_dtype = [
            ('psfex_path','U%d' % psfex_path_len),
        ]

        #image_info = meds.util.get_image_info_struct(
        image_info = get_image_info_struct(
            nimage,
            path_len,
            ext_len=ext_len,
            extra_dtype=extra_dtype,
        )
        image_info['position_offset'] = 1
        image_info['image_ext'] = self['image_ext']
        image_info['weight_ext'] = self['weight_ext']

        for i,f in enumerate(flist):
            image_info['image_id'][i] = i
            image_info['image_path'][i] = f
            image_info['weight_path'][i] = f
            image_info['psfex_path'][i] = psfex_flist[i]

        return image_info


    def _setup_fpack(self):
        # -qz 4.0 instead of -q 4.0
        # this means preserve zero pixels
        self['fpack_command'] = \
            'fpack -qz 4.0 -t %d,%d {fname}' % tuple(self['fpack_dims'])

    def _compress_meds_file(self, ucfilename, fzfilename):
        """
        run fpack on the file

        parameters
        ----------
        ucfilename: string
            filename for the uncompressed file
        fzfilename: string
            filename for the compressed file
        """
        from os.path import basename

        tup=(basename(ucfilename),basename(fzfilename))
        print('compressing file: %s -> %s' % tup)
        tpath=files.expandpath(fzfilename)
        if os.path.exists(tpath):
            os.remove(tpath)

        tmpdir = os.path.dirname(ucfilename)
        with StagedOutFile(fzfilename,tmpdir=tmpdir) as sf:
            cmd = self['fpack_command']
            cmd = cmd.format(fname=ucfilename)
            ret=os.system(cmd)

            if ret != 0:
                raise RuntimeError("failed to compress file")

        print('output is in:',fzfilename)


def get_image_info_struct(nimage, path_len,
                          image_id_len=None,
                          wcs_len=None,
                          ext_len=None,
                          extra_dtype=None):
    """
    get the image info structure

    Set default scale to 1.0. The other fields are 0 for
    numbers, or blank for strings

    parameters
    ----------
    nimage: int
        number of images in array
    path_len: int
        length of path strings
    wcs_len: int, optional
        length of wcs strings. If not sent, wcs will not
        be present in the array
    ext_len: int, optional
        If sent, the extension is assumed to be a string
        instead of an integer, and this is the length
    """
    dt = get_image_info_dtype(
        path_len,
        image_id_len=image_id_len,
        wcs_len=wcs_len,
        ext_len=ext_len,
        extra_dtype=extra_dtype,
    )

    data = np.zeros(nimage, dtype=dt)

    data['scale'] = 1.0

    return data

IMAGE_INFO_TYPES = ['image','weight','seg','bmask','bkg']
def get_image_info_dtype(path_len,
                         image_id_len=None,
                         wcs_len=None,
                         ext_len=None,
                         extra_dtype=None):
    """
    get the image_info dtype for the specified path string
    length and wcs string length

    parameters
    ----------
    path_len: int
        length of path strings
    wcs_len: int, optional
        length of wcs strings. If not sent, wcs will not
        be present in data type
    ext_len: int, optional
        If sent, the extension is assumed to be a string
        instead of an integer, and this is the length
    """

    path_fmt = 'U%d' % path_len

    if image_id_len is None:
        image_id_descr = 'i8'
    else:
        image_id_descr = 'U%d' % image_id_len

    if ext_len is not None:
        ext_descr = 'U%d' % ext_len
    else:
        ext_descr = 'i2'
    dt=[]
    for ctype in IMAGE_INFO_TYPES:
        path_name = '%s_path' % ctype
        ext_name  = '%s_ext' % ctype

        dt += [
            (path_name, path_fmt),
            (ext_name,ext_descr),
        ]

    dt += [
        ('image_id', image_id_descr),
        ('image_flags', 'i8'),
        ('magzp', 'f4'),
        ('scale', 'f4'),
        ('position_offset','f8'),
    ]
    if wcs_len is not None:
        wcs_fmt = 'U%d' % wcs_len
        dt += [
            ('wcs',wcs_fmt),
        ]

    if extra_dtype is not None:
        dt += extra_dtype

    return dt


