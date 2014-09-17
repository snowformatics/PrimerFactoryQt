import time
import os
import logging
import platform
from types import *

    
def delete_databases(db_list, db_location):
    """Deletes all selected databases."""
    bowtie_endings = [".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt"]

    os.chdir(db_location)
    for db in db_list:
        try:
            for extension in bowtie_endings:
                os.remove(db_location + db + extension)
        except (IOError, OSError):
            logging.debug(time.strftime("%d.%m.%Y um %H:%M:%S Uhr"))
            logging.debug(str(platform.system()+platform.release()))
            logging.exception('Got exception on main handler')
            return "Could not delete database!", False