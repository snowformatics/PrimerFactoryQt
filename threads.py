import tempfile
import subprocess
import os
import logging
import platform
import time
from PyQt4 import QtCore, QtGui, QtDeclarative
from types import *

import database_tools

class BowtieThread(QtCore.QThread):
    def __init__(self, sequence_path, db_location, primer_path, missmatches):
        QtCore.QThread.__init__(self)

        self.sequence = sequence_path
        self.db_location = db_location
        self.missmatches = missmatches
        self.primer_path = primer_path

    def run(self):
        """Bowtie search thread."""
        # Start Bowtie analyse
        info_message, bowtie_file_path = self.start_bowtie(self.sequence,
                                                            self.db_location,
                                                            self.primer_path,
                                                            self.missmatches)
        
        self.emit(QtCore.SIGNAL("threadDone(QString, QString)"), info_message, bowtie_file_path)
        
    def start_bowtie(self, sequence, db_location, primer_path, missmatches):
        """Start Bowtie alignment."""

        # Temp file for results
        temp_bowtie_file = tempfile.mkstemp()
        db_name = 'my_query_db'
        database_file_location = sequence

        assert type(missmatches) is StringType
        assert type(database_file_location) is StringType
        assert type(sequence) is StringType
        assert type(primer_path) is StringType
        assert type(temp_bowtie_file[1]) is StringType

        os.chdir(db_location)
        
        try:
            # Create database from query sequence, because length can be > 1000 bp. Input sequence will be all primers.
            process = subprocess.Popen(["bowtie-build", sequence, db_location + db_name])
            process.wait()

            if os.path.exists(db_location + str(db_name) + ".rev.1.ebwt"):

                # Now start Bowtie search
                process = subprocess.Popen(["bowtie", "-a", "-v", missmatches,  "-y", db_location + db_name, "-f",
                                            primer_path, temp_bowtie_file[1]+'.txt'])
                process.wait()
                return 'Search done!', temp_bowtie_file[1]+'.txt'

            else:
                logging.debug(time.strftime("%d.%m.%Y um %H:%M:%S Uhr"))
                logging.debug(str(platform.system()+platform.release()))
                logging.exception('Bowtie DB. Got exception on main handler')
                return "Error, database could not be created!", False

        except (IOError, OSError):
            logging.debug(time.strftime("%d.%m.%Y um %H:%M:%S Uhr"))
            logging.debug(str(platform.system()+platform.release()))
            logging.exception('Bowtie. Got exception on main handler')
            return 'Bowtie was not found! Please install <a href=\"http://bowtie-bio.sourceforge.net/index.shtml\">Bowtie</a>'