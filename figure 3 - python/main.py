# MEDITATIONS ON GEOMETRY ##############################################################################################
import numpy as np
from utils.spinor_class import *
# MAIN #################################################################################################################

# 1) surface - unitary - quaternion - rotor - path - moving frames - parallel transport
spinor = Spinor()

# 2) Calculate phase grids
sys.stdout.write("Calculating S1 phase grids \n")
spinor.drawProgressBar(0)
ξi, γi, DARBOUX_γi = spinor.phase_grids(principle_axis='i')
spinor.drawProgressBar(33.3)
ξj, γj, DARBOUX_γj = spinor.phase_grids(principle_axis='j')
spinor.drawProgressBar(66.6)
ξk, γk, DARBOUX_γk = spinor.phase_grids(principle_axis='k')
spinor.drawProgressBar(100)

# 3) Plot phase grids
if not path.exists('figures'):
    makedirs('figures')
sys.stdout.write("\nGenerating phase plots \n")
spinor.drawProgressBar(0)
spinor.plot_phase_grids([ξi, γi, DARBOUX_γi], [ξj, γj, DARBOUX_γj], [ξk, γk, DARBOUX_γk],
                        clip_geometricphase=True, progress=0, figsize=(17.8, 20.2))
spinor.plot_phase_grids([ξi, γi, DARBOUX_γi], [ξj, γj, DARBOUX_γj], [ξk, γk, DARBOUX_γk],
                        clip_geometricphase=False, progress=50, figsize=(17.8, 20.2))

# FIN ##################################################################################################################
