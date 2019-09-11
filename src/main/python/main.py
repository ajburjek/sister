#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Two-Body Orbit"""
import sys
import astropy.units as u
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.plotting import OrbitPlotter2D
from PyQt5.QtWidgets import *
from fbs_runtime.application_context.PyQt5 import ApplicationContext
#matplotlib connections
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
"""Python 3.7
   Simpson Aerospace (c) 2019
   Christopher R. Simpson
   christopher.r.simpson@simpsonaerospace.com"""
#------------------------------------------------------------------------------
class PlotWindow(QDialog):
    def __init__(self, parent=None):
        super(PlotWindow, self).__init__(parent)
        
        def createLeftGroupBox(self):
            self.leftGroupBox = QGroupBox("Classical Orbital Elements")
            
            #connect button to plot method
            def on_button_clicked():
                #Two body calculation
                aPoli = apogee.value() * u.km
                eccPoli = ecc.value() * u.one
                incPoli = inc.value() * u.deg
                raanPoli = raan.value() * u.deg
                aopPoli = aop.value() * u.deg
                nuPoli = nu.value() * u.deg
                ss = Orbit.from_classical(Earth, aPoli, eccPoli, incPoli, raanPoli, aopPoli, nuPoli)
                op = OrbitPlotter2D()
                op.plot(ss, label="Computed Orbit")
                self.canvas.draw()
                
            #Text Entry for Classical Orbital Elements
            apogee = QDoubleSpinBox()
            apogee.setDecimals(6)
            apogee.setMaximum(8.078e+8)
            apogee.setSuffix(' km')
            
            ecc = QDoubleSpinBox()
            ecc.setDecimals(6)
            ecc.setMaximum(1)
            
            inc = QDoubleSpinBox()
            inc.setDecimals(6)
            inc.setMaximum(180.0)
            inc.setMinimum(0)
            inc.setSuffix(' deg')
            
            raan = QDoubleSpinBox()
            raan.setDecimals(6)
            raan.setMaximum(360.0)
            raan.setSuffix(' deg')
            
            aop = QDoubleSpinBox()
            aop.setDecimals(6)
            aop.setMaximum(360.0)
            aop.setSuffix(' deg')
            
            nu = QDoubleSpinBox()
            nu.setDecimals(6)
            nu.setMaximum(360.0)
            nu.setSuffix(' deg')
            
            button = QPushButton('Apply')
            button.clicked.connect(on_button_clicked)
            
            layout = QVBoxLayout()
            layout.addWidget(apogee)
            layout.addWidget(ecc)
            layout.addWidget(inc)
            layout.addWidget(raan)
            layout.addWidget(aop)
            layout.addWidget(nu)
            layout.addWidget(button)
            self.leftGroupBox.setLayout(layout)
            
#        def createRightGroupBox(self):
#            self.rightGroupBox = QGroupBox("Group 2")
#            #figure instance to plot on
#            self.figure = Figure()
#            #canvas widget to display figure
#            self.canvas = FigureCanvas(self.figure)
#            #navigation widget
#            self.toolbar = NavigationToolbar(self.canvas, self)
#            
#            layout = QVBoxLayout()
##            layout.addWidget(self.figure)
#            self.rightGroupBox.setLayout(layout)
            
        #layout
        createLeftGroupBox(self)
#        createRightGroupBox(self)
        
        mainLayout = QGridLayout()
        mainLayout.addWidget(self.leftGroupBox)#, 1, 0)
#        mainLayout.addWidget(self.rightGroupBox, 1, 1)
        self.setLayout(mainLayout)

if __name__ == '__main__':
    appctxt = ApplicationContext()
    gallery = PlotWindow()
    gallery.show()
    sys.exit(appctxt.app.exec_())