from __future__ import print_function
import os
import shutil
import tarfile
import yaml
import tempfile

class StagedInFile(object):
    """
    A class to represent a staged file
    If tmpdir=None no staging is performed and the original file path is used

    parameters
    ----------
    fname: string
        original file location
    tmpdir: string, optional
        If not sent, no staging is done.

    examples
    --------
    # using a context for the staged file
    fname="/home/jill/output.dat"
    tmpdir="/tmp"
    with StagedInFile(fname,tmpdir=tmpdir) as sf:
        with open(sf.path) as fobj:
            # read some data

    """
    def __init__(self, fname, tmpdir=None):

        self._set_paths(fname, tmpdir=tmpdir)
        self.stage_in()

    def _set_paths(self, fname, tmpdir=None):
        fname=expandpath(fname)

        self.original_path = fname

        if tmpdir is not None:
            self.tmpdir = expandpath(tmpdir)
        else:
            self.tmpdir = tmpdir

        self.was_staged_in = False
        self._stage_in = False

        if self.tmpdir is not None:
            bdir,bname = os.path.split(self.original_path)
            self.path = os.path.join(self.tmpdir, bname)

            if self.tmpdir == bdir:
                # the user sent tmpdir as the source dir, no
                # staging is performed
                self._stage_in = False
            else:
                self._stage_in = True

    def stage_in(self):
        """
        make a local copy of the file
        """
        import shutil

        if self._stage_in:
            if not os.path.exists(self.original_path):
                raise IOError("file not found:",self.original_path)

            if os.path.exists(self.path):
                print("removing existing file:",self.path)
                os.remove(self.path)
            else:
                makedir_fromfile(self.path)

            print("staging in",self.original_path,"->",self.path)
            shutil.copy(self.original_path,self.path)

            self.was_staged_in = True

    def cleanup(self):
        if os.path.exists(self.path) and self.was_staged_in:
            print("removing temporary file:",self.path)
            os.remove(self.path)
            self.was_staged_in = False

    def __del__(self):
        self.cleanup()

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.cleanup()


class StagedOutFile(object):
    """
    A class to represent a staged file
    If tmpdir=None no staging is performed and the original file
    path is used
    parameters
    ----------
    fname: string
        Final destination path for file
    tmpdir: string, optional
        If not sent, or None, the final path is used and no staging
        is performed
    must_exist: bool, optional
        If True, the file to be staged must exist at the time of staging
        or an IOError is thrown. If False, this is silently ignored.
        Default False.
    examples
    --------

    fname="/home/jill/output.dat"
    tmpdir="/tmp"
    with StagedOutFile(fname,tmpdir=tmpdir) as sf:
        with open(sf.path,'w') as fobj:
            fobj.write("some data")

    """
    def __init__(self, fname, tmpdir=None, must_exist=False):

        self.must_exist = must_exist
        self.was_staged_out = False

        self._set_paths(fname, tmpdir=tmpdir)


    def _set_paths(self, fname, tmpdir=None):
        fname=expandpath(fname)

        self.final_path = fname

        if tmpdir is not None:
            self.tmpdir = expandpath(tmpdir)
        else:
            self.tmpdir = tmpdir

        fdir = os.path.dirname(self.final_path)

        if self.tmpdir is None:
            self.is_temp = False
            self.path = self.final_path
        else:
            if not os.path.exists(self.tmpdir):
                os.makedirs(self.tmpdir)

            bname = os.path.basename(fname)
            self.path = os.path.join(self.tmpdir, bname)

            if self.tmpdir==fdir:
                # the user sent tmpdir as the final output dir, no
                # staging is performed
                self.is_temp = False
            else:
                self.is_temp = True

    def stage_out(self):
        """
        if a tempdir was used, move the file to its final destination
        note you normally would not call this yourself, but rather use a
        context, in which case this method is called for you
        with StagedOutFile(fname,tmpdir=tmpdir) as sf:
            #do something
        """
        import shutil

        if self.is_temp and not self.was_staged_out:
            if not os.path.exists(self.path):
                if self.must_exist:
                    mess = "temporary file not found: %s" % self.path
                    raise IOError(mess)
                else:
                    return

            if os.path.exists(self.final_path):
                print("removing existing file:",self.final_path)
                os.remove(self.final_path)

            makedir_fromfile(self.final_path)

            print("staging out '%s' -> '%s'" % (self.path,self.final_path))
            shutil.move(self.path,self.final_path)

        self.was_staged_out=True

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.stage_out()

class TempFile(object):
    """
    A class to represent a temporary file

    parameters
    ----------
    fname: string
        The full path for file

    examples
    --------

    # using a context for the staged file
    fname="/home/jill/output.dat"
    with TempFile(fname) as sf:
        with open(sf.path,'w') as fobj:
            fobj.write("some data")

            # do something with the file
    """
    def __init__(self, fname):
        self.path = fname

        self.was_cleaned_up = False

    def cleanup(self):
        """
        remove the file if it exists, if not already cleaned up
        """
        import shutil

        if not self.was_cleaned_up:
            if os.path.exists(self.path):
                print("removing:",self.path)
                os.remove(self.path)

            self.was_cleaned_up=True

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.cleanup()


def expandpath(path):
    """
    expand environment variables, user home directories (~), and convert
    to an absolute path
    """
    path=os.path.expandvars(path)
    path=os.path.expanduser(path)
    path=os.path.realpath(path)
    return path


def makedir_fromfile(fname):
    """
    extract the directory and make it if it does not exist
    """
    dname=os.path.dirname(fname)
    try_makedir(dname)

def try_makedir(dir):
    """
    try to make the directory
    """
    if not os.path.exists(dir):
        try:
            print("making directory:",dir)
            os.makedirs(dir)
        except:
            # probably a race condition
            pass

def get_temp_dir():
    """
    get a temporary directory.  Check for batch system specific
    directories in environment variables, falling back to TMPDIR
    """
    tmpdir=os.environ.get('_CONDOR_SCRATCH_DIR',None)
    if tmpdir is None:
        tmpdir=os.environ.get('TMPDIR',None)
        if tmpdir is None:
            tmpdir = tempfile.mkdtemp()
    return tmpdir


def read_yaml(fname):
    with open(fname) as fobj:
        data=yaml.load(fobj)

    return data


#
# specific for the desdm version
#

def get_desdm_file_config(medsconf, tilename, band):
    """
    the desdm version needs a file config

    parameters
    ----------
    medsconf: string
        A name for the meds version or config.  e.g. '013'
        or 'y3a1-v02'
    tilename: string
        e.g. 'DES0417-5914'
    band: string
        e.g. 'i'
    """

    type='fileconf'
    ext='yaml'
    subdir='lists-%s' % band

    return get_meds_datafile_generic(
        medsconf,
        tilename,
        band,
        type,
        ext,
        subdir=subdir,
    )

def get_desdm_finalcut_flist(medsconf, tilename, band):
    """
    the desdm version needs a list

    parameters
    ----------
    medsconf: string
        A name for the meds version or config.  e.g. '013'
        or 'y3a1-v02'
    tilename: string
        e.g. 'DES0417-5914'
    band: string
        e.g. 'i'
    """

    type='finalcut-flist'
    ext='dat'
    subdir='lists-%s' % band

    return get_meds_datafile_generic(
        medsconf,
        tilename,
        band,
        type,
        ext,
        subdir=subdir,
    )


def get_desdm_nullwt_flist(medsconf, tilename, band):
    """
    the desdm version needs a list

    parameters
    ----------
    medsconf: string
        A name for the meds version or config.  e.g. '013'
        or 'y3a1-v02'
    tilename: string
        e.g. 'DES0417-5914'
    band: string
        e.g. 'i'
    """

    type='nullwt-flist'
    ext='dat'
    subdir='lists-%s' % band

    return get_meds_datafile_generic(
        medsconf,
        tilename,
        band,
        type,
        ext,
        subdir=subdir,
    )

def get_coaddinfo_file(medsconf, tilename, band):
    """
    the desdm version needs a list

    parameters
    ----------
    medsconf: string
        A name for the meds version or config.  e.g. '013'
        or 'y3a1-v02'
    tilename: string
        e.g. 'DES0417-5914'
    band: string
        e.g. 'i'
    """

    type='coaddinfo'
    ext='yaml'
    subdir='lists-%s' % band

    return get_meds_datafile_generic(
        medsconf,
        tilename,
        band,
        type,
        ext,
        subdir=subdir,
    )


def get_desdm_seg_flist(medsconf, tilename, band):
    """
    the desdm version needs a list

    parameters
    ----------
    medsconf: string
        A name for the meds version or config.  e.g. '013'
        or 'y3a1-v02'
    tilename: string
        e.g. 'DES0417-5914'
    band: string
        e.g. 'i'
    """

    type='seg-flist'
    ext='dat'
    subdir='lists-%s' % band

    return get_meds_datafile_generic(
        medsconf,
        tilename,
        band,
        type,
        ext,
        subdir=subdir,
    )


def get_desdm_bkg_flist(medsconf, tilename, band):
    """
    the desdm version needs a list

    parameters
    ----------
    medsconf: string
        A name for the meds version or config.  e.g. '013'
        or 'y3a1-v02'
    tilename: string
        e.g. 'DES0417-5914'
    band: string
        e.g. 'i'
    """

    type='bkg-flist'
    ext='dat'
    subdir='lists-%s' % band

    return get_meds_datafile_generic(
        medsconf,
        tilename,
        band,
        type,
        ext,
        subdir=subdir,
    )


def get_desdm_objmap(medsconf, tilename, band):
    """
    the desdm version needs a map

    parameters
    ----------
    medsconf: string
        A name for the meds version or config.  e.g. '013'
        or 'y3a1-v02'
    tilename: string
        e.g. 'DES0417-5914'
    band: string
        e.g. 'i'
    """

    type='objmap'
    ext='fits'
    subdir='lists-%s' % band

    return get_meds_datafile_generic(
        medsconf,
        tilename,
        band,
        type,
        ext,
        subdir=subdir,
    )

def try_remove_timeout(fname, ntry=2, sleep_time=2):
    import time
    
    for i in xrange(ntry):
        try:
            os.remove(fname)
            break
        except:
            if i==(ntry-1):
                raise
            else:
                print("could not remove '%s', trying again "
                      "in %f seconds" % (fname,sleep_time))
                time.sleep(sleep_time)

def try_remove(f):
    try:
        os.remove(f)
        print("removed file:",f)
    except:
        print("could not remove file:",f)


def try_remove_dir(d):
    try:
        shutil.rmtree(d)
        print("removed dir:",d)
    except:
        print("could not remove dir:",d)


def tar_directory(source_dir):
    """
    tar a directory to a tar file called directory.tar.gz
    """
    outfile=source_dir+'.tar.gz'
    print("tarring directory %s -> %s" % (source_dir, outfile))
    with tarfile.open(outfile, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))

