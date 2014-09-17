from distutils.core import setup
import py2exe
# cd D:\Users\lueck\Google Drive\Python\Software projects\Primer Factory Qt
# python setup.py py2exe

# It is important to include the graph widget and items directly. Not nice but work around.

setup(windows=[{"script": "primer_main.py"}],

      options={"py2exe":
                {"includes": ['sip', 'decimal', 'PyQt4.QtSql', "PyQt4.QtCore", "PyQt4.QtGui", "PyQt4.QtNetwork",
                              r'scipy.sparse.csgraph._validation',  r'scipy.special._ufuncs_cxx',
                              'pyqtgraph.graphicsItems.ArrowItem', 'pyqtgraph.graphicsItems.TextItem', 'pyqtgraph.widgets.TableWidget'],
                  #"bundle_files":1,
                  #"optimize": 2,
                  "dll_excludes": ["mswsock.dll", "powrprof.dll"]}
                })

#


