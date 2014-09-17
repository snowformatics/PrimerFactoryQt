import re
import tempfile
from Bio import SeqIO
from Bio.SeqUtils.MeltingTemp import Tm_staluc
from Bio.Emboss.Applications import Primer3Commandline
#from Bio.Application import generic_run
#from Bio.Emboss import Primer3
from math import log10
import os


#=====================================================================
def validate_seq(sequence):
    """Validate DNA sequence."""
    sequence = sequence.strip()
    sequence = sequence.replace(" ", "")
    sequence = sequence.replace("\n", "")
    sequence.upper()
    regex = re.compile('^[ACTGNRYSWKMBDHVEFILPQSXZ]*$', re.I)
    if regex.search(sequence) is not None:
        return True
    else:
        return False


def validate_fasta_seq(sequence):
    """Validate Fasta format."""
    sequence = sequence.replace(" ", "")
    sequence.upper()

    regex = re.compile('>\S*\n[ACTGNRYSWKMBDHVEFILPQSXZ]*', re.MULTILINE)

    if regex.search(sequence) is not None:
    #if regex.search(sequence) is not None and sequence.count('>') == 1:
        return True
    else:
        return False


def add_header(in_file):
    """Adds a header a first line of outfile and return new file."""

    header = "Read name\t" + "Reference strand\t" + "Name of reference sequence\t" \
                 + "Position alignment occurs\t" + "Read sequence\t" + "Read qualities\t" \
                 + "Ceiling\t" + "Mismatch descriptors\n"

    # Temp file for final results including header
    temp_out = tempfile.mkstemp()
    f_in = open(in_file, 'r')
    results = f_in.read()
    f_out = open(temp_out[1] + '.txt', 'w')
    f_out.write(header)
    f_out.write(results)

    f_in.close()
    f_out.close()
    return temp_out[1] + '.txt'


#=====================================================================
# Tm calculations
class Tm(object):
    def __init__(self, natrium_conc, primer_conc):
        self.natrium_conc = float(natrium_conc)
        self.primer_conc = float(primer_conc)

    def import_export_tm_data(self, in_file, method):
        """Import a fasta file, calculate tm and export results."""
        # Open sequence file
        f_fasta = open(in_file, "r")
        temp_tm_file = tempfile.mkstemp()
        f_out = open(temp_tm_file[1] + '.txt', 'w')
        for seq_record in SeqIO.parse(f_fasta, "fasta"):
            id = seq_record.id
            seq = str(seq_record.seq).upper()
            tm = method(seq)
            f_out.write(str(id) + '\t' + seq + '\t' + str(tm) + '\n')
        f_out.close()
        f_fasta.close()
        return temp_tm_file[1] + '.txt'

    def basic_tm(self, sequence):
        """Calculates the TM with basic rules.
        If len(sequence) < 13:
        Tm= (wA+xT) * 2 + (yG+zC) * 4
        If len(sequence) > 13:
        Tm= 64.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC)
        Where w,x,y,z are the number of the bases A,T,G,C in the sequence"""
        # Calculate Tm
        for base in range(4):
            count_a = sequence.count('A')
            count_c = sequence.count('C')
            count_g = sequence.count('G')
            count_t = sequence.count('T')
        if len(sequence) < 14:
            tm_basic = (2 * (count_a + count_t)) + (4 * (count_g + count_c))
        else:
            tm_basic = 64.9 + 41 * (count_g + count_c - 16.4) / (count_a + count_t + count_g + count_c)
        return round(tm_basic, 1)

    def salt_adjusted(self, sequence):
        """Salt adjusted calculation.
        If len(sequence) < 13:
        Tm= (wA+xT)*2 + (yG+zC)*4 - 16.6*log10(0.050) + 16.6*log10([Na+])
        If len(sequence) > 13:
        Tm= 100.5 + (41 * (yG+zC)/(wA+xT+yG+zC)) - (820/(wA+xT+yG+zC)) + 16.6*log10([Na+])
        """
        # Calculate Tm
        for base in range(4):
            count_a = sequence.count('A')
            count_c = sequence.count('C')
            count_g = sequence.count('G')
            count_t = sequence.count('T')

        if len(sequence) < 14:
            tm_salt = (2 * (count_a+count_t)) + (4 * (count_g+count_c)) - 16.6 * (log10(0.050)) + 16.6*log10(self.natrium_conc/1000.0)
        else:
            tm_salt = 100.5 + (41.0 * (count_g + count_c)/(count_a+count_t+count_g+count_c)) - (820.0/(count_a+count_t+count_g+count_c)) + 16.6*log10(self.natrium_conc/1000.0)
        return round(tm_salt, 1)

    def nearest_neighbour(self, sequence):
        """Calculate Tm with rule.
        See http://www.basic.northwestern.edu/biotools/oligocalc.html"""
        tm_nn = Tm_staluc(sequence, dnac=self.primer_conc, saltc=self.natrium_conc)
        return round(tm_nn, 1)

class PrimerDesign(object):
    def __init__(self, parameter_list, db_location, seq_infile, raw_data):

        # Paramters
        self.parameter_list = parameter_list
        self.primer_3_location = "eprimer3.exe"
        self.db_location = db_location
        self.seq_infile = seq_infile
        self.raw_data = raw_data
        # Create result file
        temp_tm_file = tempfile.mkstemp()
        self.out_file = temp_tm_file[1] + '.txt'


    def create_input_file(self, seq_infile):
        """Formats the sequence input file for primer3."""
        seq_in_tmp = tempfile.mkstemp()
        f_in = open(seq_in_tmp[1] + '.txt', 'w')
        for seq_record in SeqIO.parse(seq_infile, "fasta"):
            id = seq_record.id
            seq = str(seq_record.seq).upper()
            f_in.write("PRIMER_SEQUENCE_ID=")
            f_in.write(id + '\n')
            f_in.write("SEQUENCE_TEMPLATE=")
            f_in.write(seq)
            f_in.write("\n=\n")
        f_in.close()
        return seq_in_tmp[1] + '.txt'

    def create_setting_file(self, parameter_list):
        """Formats the setting input file for primer3."""
        setting_tmp = tempfile.mkstemp()
        f_in = open(setting_tmp[1] + '.txt', 'w')
        # Required header
        f_in.write("Primer3 File - http://primer3.sourceforge.net\nP3_FILE_TYPE=settings\n\nP3_FILE_ID=Default settings of primer3 version 2.3.6\n")

        for parameter in parameter_list:
            f_in.write(parameter[0] + '=' + parameter[1] + '\n')
        config_path = self.db_location
        #f_in.write('PRIMER_THERMODYNAMIC_PARAMETERS_PATH=D:\\Users\\lueck\\Google Drive\\Python\\Software projects\\Testing\\Primer3\\release-2.3.6\\primer3_config\\')
        f_in.write('PRIMER_THERMODYNAMIC_PARAMETERS_PATH='+config_path)
        f_in.write('\n=')
        f_in.close()
        return setting_tmp[1] + '.txt'

    def create_parsed_file(self, out_file):
        """Creates a nice to read and csv friendly output file."""
        parsed_tmp = tempfile.mkstemp()
        f_out = open(parsed_tmp[1] + '.txt', 'w')
        f_in = open(out_file, 'r')

        header = 'Primer Name\tPrimer Left\tPrimer Right\tPrimer Left Position\tPrimer Right Position\tPrimer Left Tm\tPrimer Right Tm\tPrimer Left GC\tPrimer Right GC\tProduct size\n'
        f_out.write(header)
        data = f_in.readlines()
        for line in data:

            if line.startswith("PRIMER_SEQUENCE_ID"):
                f_out.write(line.split('=')[1].strip() + '\t')
            if line.split('=')[0].endswith("_SEQUENCE"):
                f_out.write(line.split('=')[1].strip() + '\t')
            if line.split('=')[0].endswith("_TM"):
                f_out.write(line.split('=')[1].strip() + '\t')
            if line.split('=')[0].endswith("_GC_PERCENT"):
                f_out.write(line.split('=')[1].strip() + '\t')
            if line.split('=')[0].endswith("_PRODUCT_SIZE"):
                f_out.write(line.split('=')[1].strip() + '\n')
            if line.startswith("PRIMER_LEFT") and len(line.split('=')[0]) < 15:
                f_out.write(line.split('=')[1].split(',')[0] + '\t')
            if line.startswith("PRIMER_RIGHT") and len(line.split('=')[0]) < 15:
                f_out.write(line.split('=')[1].split(',')[0] + '\t')

        f_in.close()
        f_out.close()
        return parsed_tmp[1] + '.txt'




    def start_primer3(self):
        """Design primers with eprimer3."""
        seq_infile = self.create_input_file(self.seq_infile)
        setting_file = self.create_setting_file(self.parameter_list)

        os.chdir(self.db_location)
        cmd = "primer3_core.exe <" + seq_infile + " >" + self.out_file + "  -p3_settings_file=" + setting_file
        os.system(cmd)
        parsed_file = self.create_parsed_file(self.out_file)
        if self.raw_data:
            return [self.out_file, seq_infile, setting_file]
        else:
            return [parsed_file]




# Testing
# my_tm = Tm(50, 50)
# print my_tm.import_export_tm_data('D:\Users\lueck\Desktop\primer.txt', my_tm.basic_tm)
# print my_tm.import_export_tm_data('D:\Users\lueck\Desktop\primer.txt', my_tm.salt_adjusted)
# print my_tm.import_export_tm_data('D:\Users\lueck\Desktop\primer.txt', my_tm.nearest_neighbour)
