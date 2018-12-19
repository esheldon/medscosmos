from __future__ import print_function
import os
import shutil
import tarfile
import yaml
import tempfile
from glob import glob

def get_cosmos_dir():
    """
    base cosmos directory
    """
    return os.environ['COSMOS_DIR']

def get_meds_dir():
    """
    get the base MEDS_DIR
    """
    return os.path.join(
        get_cosmos_dir(),
        'meds',
    )

def get_meds_run_dir(run):
    """
    get the directory for the given run

    Parameters
    ----------
    run: string
        e.g. cosmos-des01
    """
    md=get_meds_dir()
    return os.path.join(
        md,
        run,
    )

def get_script_dir(run):
    """
    get the directory for the given run scripts

    Parameters
    ----------
    run: string
        e.g. cosmos-des01
    """

    md=get_meds_dir()
    return os.path.join(
        md,
        run,
        'scripts',
    )

def get_script_file(run, tileid):
    """
    get the directory for the given run scripts

    Parameters
    ----------
    run: string
        e.g. cosmos-des01
    """

    d=get_script_dir(run)
    tilestr = get_tilestr(tileid)
    fname='make-%s.sh' % tilestr

    return os.path.join(
        d,
        fname,
    )



def get_tilestr(tileid):
    """
    get zero padded tile string

    Parameters
    ----------
    tileid: int or string
        e.g. 35 or '035'
    """
    tileid=int(tileid)
    return '%03d' % tileid

def get_meds_file_dir(run, tileid):
    """
    get the directory for the given tile

    Parameters
    ----------
    run: string
        e.g. cosmos-des01
    tileid: int
        e.g. 35
    """

    tilestr = get_tilestr(tileid)
    run_dir=get_meds_run_dir(run)
    return os.path.join(
        run_dir,
        tilestr,
    )

def get_meds_file(run, tileid, survey, band):
    """
    get the directory for the given tile

    Parameters
    ----------
    run: string
        e.g. cosmos-des01
    tileid: int
        e.g. 35
    survey: string
        e.g. cosmos or des
    band: string
        band string such as 'g'
    """
    d = get_meds_file_dir(run, tileid)

    tilestr = get_tilestr(tileid)
    fname = 'meds-%(run)s-%(tilestr)s-%(survey)s-%(band)s.fits.fz'
    fname = fname % {
        'run':run,
        'tilestr':tilestr,
        'survey':survey,
        'band':band,
    }

    return os.path.join(
        d,
        fname,
    )


def get_flist_dir():
    """
    get directory holding file lists
    """
    return os.path.join(
        get_cosmos_dir(),
        'flists',
    )

def get_cosmos_flist(tileid):
    """
    get the cosmos coadd file list

    Parameters
    ----------
    tileid: int
        e.g. 35
    """
    tilestr = get_tilestr(tileid)
    fname='cosmos-flist%s.dat' % tilestr
    return os.path.join(
        get_flist_dir(),
        fname,
    )

def get_cosmos_flists():
    """
    get all the cosmos coadd file lists
    """
    d = get_flist_dir()
    fl = glob(d+'/cosmos-flist[0-9][0-9][0-9].dat')
    fl.sort()
    return fl

def get_cosmos_tileids():
    """
    get all the cosmos coadd tile ids
    """
    fl = get_cosmos_flists()
    ids = []

    for f in fl:
        bname = os.path.basename(f)
        tileid = int( bname[-7:].replace('.dat','') ) 
        ids.append(tileid)

    return ids

def get_des_flist(band):
    """
    get the cosmos coadd file list

    Parameters
    ----------
    band: string
        e.g. 'i'
    """
    fname='cosmos-flist-des-%s.dat' % band
    return os.path.join(
        get_flist_dir(),
        fname,
    )


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

