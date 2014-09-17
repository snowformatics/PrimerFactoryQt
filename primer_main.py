import sys
import os
import tempfile
import shutil
import logging
import time
import platform
import webbrowser

import threads

from PyQt4 import QtCore, QtGui, QtDeclarative
from primer_ui import Ui_MainWindow

import database_tools
import sequence_tools
import primer_results

# pyrcc4 primer.qrc > primer_rc.py
# pyuic4 primer.ui > primer_ui.py


class MyMainWindow(QtGui.QMainWindow):
    def __init__(self, *args):
        QtGui.QMainWindow.__init__(self, *args)

        self.ui = Ui_MainWindow()

        # Main paths
        self.data_location = str(QtGui.QDesktopServices.storageLocation(QtGui.QDesktopServices.DataLocation))
        self.home_location = str(QtGui.QDesktopServices.storageLocation(QtGui.QDesktopServices.HomeLocation))
        self.temp_location = str(QtGui.QDesktopServices.storageLocation(QtGui.QDesktopServices.TempLocation))
        self.app_location = self.data_location + '/PrimerFactoryQt/'
        self.db_location = self.app_location + '/databases/'
        self.images_location = self.app_location + '/images/'
        self.primer_path = ''
        self.sequence_length = 0
        
        try:
            self.create_folders()
            self.copying_files()
            # Logging
            self.log_file = self.app_location + '/logging_primerfactory.txt'
            logging.basicConfig(filename=self.log_file, level=logging.DEBUG)

        except (IOError, OSError):
            logging.debug(time.strftime("%d.%m.%Y um %H:%M:%S Uhr"))
            logging.debug(str(platform.system()+platform.release()))
            logging.exception('Got exception on main handler')
            raise

        self.ui.setupUi(self)
        self.create_connects()
        self.ui.plainTextEdit.setFocus()
        #self.ui.commandLinkButton.setFlat(False)

        self.ui.statusbar.showMessage("System Status | Normal | " + self.app_location)

    def create_connects(self):
        """Connect events."""
        
        # Tab Bowtie
        self.ui.commandLinkButton.clicked.connect(self.open_query_file)
        self.ui.pushButton_3.clicked.connect(self.open_primer_file)
        self.ui.pushButton.clicked.connect(self.start_bowtie)
        # Tab Tm
        self.ui.pushButton_6.clicked.connect(self.calculate_primer_tm)
        self.ui.commandLinkButton_5.clicked.connect(self.open_tm_file)
        # Tab Design
        self.ui.pushButton_4.clicked.connect(self.switch_page)
        self.ui.pushButton_5.clicked.connect(self.switch_page)
        self.ui.pushButton_8.clicked.connect(self.switch_page)
        self.ui.pushButton_7.clicked.connect(self.switch_page)
        self.ui.pushButton_2.clicked.connect(self.start_primer_design)
        self.ui.commandLinkButton_2.clicked.connect(self.open_primer_design_file)

        # Menu
        self.ui.actionAbout.triggered.connect(self.show_about_message)
        self.ui.actionDocumentation.triggered.connect(self.show_help)
        self.ui.actionClear_all.triggered.connect(self.clear_all_fields)

    def create_folders(self):
        """Create important folders of PF in the data path."""
        location_folders = [self.app_location, self.db_location, self.images_location]
        for folder in location_folders:
            if not os.path.exists(folder):
                os.mkdir(folder)

    def copying_files(self):
        """Copy important files of PF in the data path."""
        try:
            to_copy_folders = os.listdir(os.getcwd() + '/to_copy/')
            for folder in to_copy_folders:
                to_copy_files = os.listdir(os.getcwd() + '/to_copy/' + folder)
                for files in to_copy_files:
                    if not folder.startswith('.'):
                        if not files.startswith('.'):
                            if not os.path.exists(self.app_location + '/' + folder + '/' + files):
                                shutil.copyfile(os.getcwd() + '/to_copy/' + folder + '/' + files, self.app_location + '/' + folder + '/' + files)
        except (IOError, OSError):
            logging.debug(time.strftime("%d.%m.%Y um %H:%M:%S Uhr"))
            logging.debug(str(platform.system()+platform.release()))
            logging.exception('Got exception on main handler')
            raise

    #================================================================================================================
    ### File open management
    def open_query_file(self):
        """Label the sequence path of Bowtie."""
        self.ui.plainTextEdit.clear()
        self.ui.statusbar.showMessage("System Status | Normal")
        path = str(self.open_sequence_file())
        if path != "None":
            self.ui.label_11.setText(path)
            self.insert_seq(path, self.ui.plainTextEdit)

    def open_primer_file(self):
        """Label the sequence path of Bowtie."""
        self.primer_path = str(self.open_sequence_file())
        if self.primer_path != "None":
            self.ui.label_23.setText(self.primer_path)

    def open_tm_file(self):
        """Label the sequence path of Tm."""
        self.ui.plainTextEdit_5.clear()
        self.ui.statusbar.showMessage("System Status | Normal")
        tm_path = str(self.open_sequence_file())
        if tm_path != "None":
            self.ui.label_3.setText(tm_path)
            self.insert_seq(tm_path, self.ui.plainTextEdit_5)

    def open_primer_design_file(self):
        """Label the sequence path of PrimerDesign."""
        self.ui.plainTextEdit_2.clear()
        self.ui.statusbar.showMessage("System Status | Normal")
        prim_design_path = str(self.open_sequence_file())
        if prim_design_path != "None":
            self.ui.label_14.setText(prim_design_path)
            self.insert_seq(prim_design_path, self.ui.plainTextEdit_2)

    def insert_seq(self, file_location, plain_textedit):
        """Insert a sequence file content into a text edit widget."""
        self.ui.statusbar.showMessage("System Status | Normal")
        f_seq = open(file_location, 'r')
        sequence = f_seq.read()
        plain_textedit.clear()
        plain_textedit.insertPlainText(sequence)
        f_seq.close()
    
    def open_sequence_file(self):
        """Open a sequence file and return the path."""
        sequence_file_location = QtGui.QFileDialog.getOpenFileName(
                                    self,
                                    u"Open a sequence file",
                                    self.home_location,
                                    u"Sequence (*.fasta *.fas *.txt)")
        if not sequence_file_location.isNull():
            if os.path.exists(sequence_file_location):
                return str(sequence_file_location)

    ### File save management
    def save_validate_query_sequence(self, lable, plaintext):
        """Save and validate sequence and fasta format."""

        if os.path.exists(lable.text()):
            sequence_file_location = lable.text()
            f_query = open(sequence_file_location, 'r')
            sequence = f_query.read()
            f_query.close()
        elif str(plaintext.toPlainText()) != '':
            sequence = str(plaintext.toPlainText())
        else:
            sequence = None
            self.show_info_message('Please enter a DNA sequence or upload a file!')
        if sequence is not None:
            validate_seq = sequence_tools.validate_seq(sequence)
            validate_fasta = sequence_tools.validate_fasta_seq(sequence)
            if validate_seq is True or validate_fasta is True:
                temp_seq_file = tempfile.mkstemp()
                f_seq = open(temp_seq_file[1], 'w')
                if validate_seq:
                    f_seq.write('>My_sequence' + '\n')
                f_seq.write(sequence)
                f_seq.close()
                sequence_file_location = temp_seq_file[1]
                return sequence_file_location, len(sequence)
            else:
                self.show_info_message('Please enter a valid DNA sequence, one sequence per run.')

    #===================================================================================================================
    ### Primer Search Tab
    def start_bowtie(self):
        """Starts Bowtie alignment."""

        if self.primer_path != "":
            self.sequence_tmp, self.sequence_length = self.save_validate_query_sequence(self.ui.label_11, self.ui.plainTextEdit)
            if self.sequence_tmp is not None:
                self.setCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
                self.ui.statusbar.showMessage("Starting Search | Please wait...")
                self.searchThread = threads.BowtieThread(str(self.sequence_tmp),
                                                                 str(self.db_location),
                                                                 str(self.primer_path),
                                                                 str(self.ui.comboBox.currentText()))
                self.connect(self.searchThread, QtCore.SIGNAL("threadDone(QString, QString)"), self.thread_bowtie_done)
                self.searchThread.start()
        else:
            self.show_info_message('Please choose a primer file.')
        self.ui.label_11.setText('')

    def thread_bowtie_done(self, info_message, bowtie_file_path):
        """Show message after thread has finished."""
        if os.path.exists(bowtie_file_path):
                statinfo = os.stat(bowtie_file_path)
                if statinfo.st_size > 0:
                    # Add header
                    final_bowtie_file_path = sequence_tools.add_header(bowtie_file_path)
                    self.show_search_resuts(final_bowtie_file_path, self.sequence_length)
                    self.ui.statusbar.showMessage("Search done")
                else:
                    self.ui.statusbar.showMessage("Search failed")
                    self.show_info_message('Sorry, no matches found.')
        # Delete query DB
        database_tools.delete_databases(['my_query_db'], self.db_location)
        self.setCursor(QtGui.QCursor(QtCore.Qt.ArrowCursor))

    def show_search_resuts(self, path, sequence_length):
        """Create a table and a graphic with the primer search result and open a new window."""
        self.result_window = primer_results.PrimerSearchResults(path, sequence_length)
        self.result_window.show()

    #===================================================================================================================
    ### Primer Tm Tab
    def calculate_primer_tm(self):
        """Calculate primer Tm."""
        sequence_tmp, self.sequence_length = self.save_validate_query_sequence(self.ui.label_3, self.ui.plainTextEdit_5)
        if sequence_tmp is not None:
            natrium_conc = self.ui.lineEdit.text()
            primer_conc = self.ui.lineEdit_3.text()
            my_tm = sequence_tools.Tm(natrium_conc, primer_conc)
            if self.ui.radioButton.isChecked():
                f_out = my_tm.import_export_tm_data(sequence_tmp, my_tm.nearest_neighbour)
            elif self.ui.radioButton_3.isChecked():
                f_out = my_tm.import_export_tm_data(sequence_tmp, my_tm.basic_tm)
            elif self.ui.radioButton_4.isChecked():
                f_out = my_tm.import_export_tm_data(sequence_tmp, my_tm.salt_adjusted)

            webbrowser.open(f_out)
    #===================================================================================================================
    ### Primer Design Tab
    def switch_page(self):
        """Switch the setting pages."""
        if self.ui.stackedWidget.currentIndex() == 1:
            self.ui.stackedWidget.setCurrentIndex(0)
        else:
            self.ui.stackedWidget.setCurrentIndex(1)

    def get_all_parameters(self):
        """Returns all parameter in a list."""

        is_not_zero = True
        for number in range(10, 20):
            if eval("self.ui.lineEdit_" + str(number) + ".text()") == '':
                is_not_zero = False

        if is_not_zero:
            all_parm = [('SEQUENCE_EXCLUDED_REGION', self.ui.lineEdit_4.text()), ('SEQUENCE_INCLUDED_REGION', self.ui.lineEdit_5.text()),
                        ('SEQUENCE_TARGET', self.ui.lineEdit_6.text()), ('PRIMER_PRODUCT_SIZE_RANGE', self.ui.lineEdit_7.text()),
                        ('PRIMER_NUM_RETURN', self.ui.lineEdit_20.text()), ('PRIMER_MIN_SIZE', self.ui.lineEdit_10.text()),
                        ('PRIMER_OPT_SIZE', self.ui.lineEdit_8.text()), ('PRIMER_MAX_SIZE', self.ui.lineEdit_12.text()),
                        ('PRIMER_MIN_TM', self.ui.lineEdit_13.text()), ('PRIMER_OPT_TM', self.ui.lineEdit_9.text()),
                        ('PRIMER_MAX_TM', self.ui.lineEdit_14.text()), ('PRIMER_MIN_GC', self.ui.lineEdit_15.text()),
                        ('PRIMER_OPT_GC_PERCENT', self.ui.lineEdit_11.text()), ('PRIMER_MAX_GC', self.ui.lineEdit_16.text()),
                        ('PRIMER_PAIR_MAX_DIFF_TM', self.ui.lineEdit_17.text()), ('PRIMER_SALT_MONOVALENT', self.ui.lineEdit_18.text()),
                        ('PRIMER_DNA_CONC', self.ui.lineEdit_19.text())]
            return all_parm
        else:
            return None

    def start_primer_design(self):
        """Start primer design."""
        parameter_list = self.get_all_parameters()
        if parameter_list:
            sequence_tmp, self.sequence_length = self.save_validate_query_sequence(self.ui.label_14, self.ui.plainTextEdit_2)
            my_design = sequence_tools.PrimerDesign(parameter_list, self.db_location, sequence_tmp, self.ui.checkBox.isChecked())
            f_out = my_design.start_primer3()
            for resulst_file in f_out:
                webbrowser.open(resulst_file)
        else:
            self.show_info_message('Missing advanced settings parameter!\nAll fields are required and must be valid.')


    #===================================================================================================================
    def clear_all_fields(self):
        """Clears all sequence and ID fields."""
        fields = [self.ui.plainTextEdit_2, self.ui.plainTextEdit, self.ui.plainTextEdit_5]

        for field in fields:
            field.clear()

    def show_info_message(self, message):
        QtGui.QMessageBox.information(self,
                    u"Information",
                    message)

    def show_about_message(self):
        """Shows about box."""

        message = """PrimerFactoryQt v. 1.0.2 <br/>
                     (C) 2014 Stefanie Lueck <br/>
                     CC-BY-NC License <br/>
                     If you have any comments or problems please contact the author via labtools[at]ipk-gatersleben.de
                    <br/>
                     <a href=\"http://labtools.ipk-gatersleben.de\">Homepage</a>"""




        about_box = QtGui.QMessageBox.about(self,
                                            "About PrimerFactoryQt",
                                            """<p style="font-family: 'verdana'; font-size:10pt; color:#2E9AFE">"""
                                                + message + """</p>""")\



    def show_help(self):
        """Show help."""

        webbrowser.open("labtools.ipk-gatersleben.de/help/PrimerFactoryQt/PrimerFactoryQt.html")


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)

    QtGui.QApplication.setStyle(QtGui.QStyleFactory.create("gtk"))
    QtGui.QApplication.setPalette(QtGui.QApplication.style().standardPalette())
    myapp = MyMainWindow()
    myapp.show()
    sys.exit(app.exec_())
   

