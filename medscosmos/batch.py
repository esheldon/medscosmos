from . import files

class BatchMakerBase(dict):
    """
    base class for batch makers
    """
    def __init__(self, args):
        self.args=args

        self['run'] = args.run
        self['cosmos_config_file'] = args.cosmos_config_file
        self['des_config_file'] = args.des_config_file
        self['cat_file'] = args.cat_file

        self._load_ids()

    def go(self):
        """
        make all the scripts
        """
        for tileid in self.tile_ids:
            self._write_scripts(tileid)

    def _write_scripts(self, tileid):
        raise NotImplementedError("implement in a child class")

    def _load_ids(self):
        self.tile_ids = files.get_cosmos_tileids()

class BatchMakerShell(BatchMakerBase):
    """
    write shell scripts
    """
    def _write_scripts(self, tileid):
        self._write_bash_script(tileid)

    def _write_bash_script(self, tileid):
        self['tileid'] = tileid

        script_file = files.get_script_file(self['run'], tileid)
        print('writing:',script_file)
        #files.makedir_fromfile(script_file)
        
        text = _script_template % self
        #with open(script_file,'w') as fobj:
        #    fobj.write(text)

class BatchMakerWQ(BatchMakerShell):
    """
    write wq submit scripts in addition to shell scripts
    """
    def _write_scripts(self, tileid):
        super(BatchMakerWQ,self)._write_scripts()
        self._write_wq_script(tileid)

    def _write_wq_script(self, tileid):
        pass

_script_template=r"""#!/bin/bash

medscosmos-make-cosmos-des-meds \
    %(run)s \
    %(cosmos_config_file)s  \
    %(des_config_file)s  \
    %(cat_file)s \
    %(tileid)d

"""


