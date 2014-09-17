from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4 import QtGui
from PyQt4 import QtCore
import pyqtgraph as pg
import operator


class GraphDisplay(QtGui.QMainWindow):
    def __init__(self, found_primer_data, query_name, sequence_length, parent=None):
        QtGui.QMainWindow.__init__(self, parent)

        #self.setGeometry(QtCore.QRect(200, 100, 1500, 300))
        ### Constants for drawing
        # Axe range
        x_range = int(sequence_length) + 50
        y_range = len(found_primer_data)*6
        #y_range = (len(found_primer_data)*0.7) * 12
        #print y_range
        zero_range = -10
        headLen = 10
        tailLen = 10
        primer_pos_offset = headLen + tailLen

        primer_pos_y = 5
        primer_y_offset = 5
        #primer_y_offset = len(found_primer_data)*0.45
        fwd_angle = -180
        rev_angle = -360

        # Layout stuff
        pg.setConfigOption('background', 'w')
        pg.setConfigOption('foreground', 'k')

        self.view = pg.GraphicsView()
        self.lay = pg.GraphicsLayout()

        # Plot
        p = self.lay.addPlot(row=0, col=0)
        p.setMouseEnabled(False, False)
        p.showAxis('left', False)
        p.setLabel('bottom', (str(query_name) +' Sequence length in bp'))
        p.setXRange(zero_range, x_range)
        p.setYRange(zero_range, y_range)

        for hit in found_primer_data:
            # Positions
            primer_pos = int(hit[2]) + primer_pos_offset
            primer_name = str(hit[0])
            primer_length = len(hit[3])

            # Angle
            strand = hit[1]
            if strand == '+':
                angle = fwd_angle
                text_pos = primer_pos-(len(primer_name)*4)
            else:
                angle = rev_angle
                text_pos = primer_pos-(len(primer_name)*3)
                primer_pos -= primer_length

            # Draw arrow
            a = pg.ArrowItem(angle=angle, tipAngle=60, headLen=headLen, tailLen=tailLen, tailWidth=5)
                             #pen={'color': 'w', 'width': 3})
            a.setPos(primer_pos, primer_pos_y)

            # Draw text
            text = pg.TextItem(html=primer_name, anchor=(-0.3, 1.3), color=(0, 0, 0))
            text.setPos(text_pos, primer_pos_y-primer_y_offset)

            p.addItem(a)
            p.addItem(text)

            primer_pos_y += primer_y_offset

        self.view.setCentralItem(self.lay)
        self.setCentralWidget(self.view)


class PrimerSearchResults(QWidget):
    def __init__(self, path, sequence_length):
        QWidget.__init__(self)

        bowtie_result_file = path
        found_primer_data, header, query_name = self.create_data(bowtie_result_file)

        # Set title and size
        self.setWindowTitle('Primers found for query ' + query_name)
        self.resize(1000, 793)

        # Table
        table_window = pg.TableWidget()
        table_window.setData(found_primer_data)
        table_window.setHorizontalHeaderLabels(header)
        table_window.setColumnWidth(0, 200)
        table_window.setColumnWidth(1, 150)
        table_window.setColumnWidth(2, 150)
        table_window.setColumnWidth(3, 250)
        table_window.setColumnWidth(4, 150)
        graphic_window = GraphDisplay(found_primer_data, query_name, sequence_length)
        layout = QVBoxLayout(self)
        layout.addWidget(table_window)
        layout.addWidget(graphic_window)
        self.setLayout(layout)

    def create_data(self, path):
        """Reads a Bowtie file and append results into a list."""
        results = []

        f_bowtie = open(path, 'r')
        bowtie_data = f_bowtie.readlines()

        for line in bowtie_data:
            line = line.replace('\n', '')
            line = line.split('\t')
            # Convert pos to int fro sorting
            if line[3] != 'Position alignment occurs':
                line[3] = int(line[3])

            # Extract query name and remove some stuff we don't need
            query_name = line.pop(2)
            line.pop(4)
            line.pop(4)

            #Append results
            results.append(line)

        header = results[0]
        results = results[1:]
        results.sort(key=operator.itemgetter(2))

        return results, header, query_name



# if __name__ == "__main__":
#     app = QApplication(sys.argv)
#     w = PrimerSearchResults()
#     w.show()
#     sys.exit(app.exec_())