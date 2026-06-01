# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2020 replay file
# Internal Version: 2019_09_14-01.49.31 163176
# Run by huawei on Mon Jun  1 13:55:17 2026
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(1.02778, 1.02604), width=151.289, 
    height=101.783)
session.viewports['Viewport: 1'].makeCurrent()
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
execfile('compare_plate4_s4r.py', __main__.__dict__)
#* NameError: name '__file__' is not defined
#* File "compare_plate4_s4r.py", line 30, in <module>
#*     ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), 
#* os.pardir))
