# This is a python for the validation test case Beam_Unconfined_Cylinder_Free_Vibration.data

import numpy as np



t, Fpx = np.loadtxt('Beam_Unconfined_Cylinder_Free_Vibration_pb_Force_pression.out', unpack=True, usecols=[0,4])
Fvx = np.loadtxt('Beam_Unconfined_Cylinder_Free_Vibration_pb_Contrainte_visqueuse.out', unpack=True, usecols=[4])
    
Fx = Fvx + Fpx #force per unit length of cylinder

Fpy = np.loadtxt('Beam_Unconfined_Cylinder_Free_Vibration_pb_Force_pression.out', unpack=True, usecols=[5])
Fvy = np.loadtxt('Beam_Unconfined_Cylinder_Free_Vibration_pb_Contrainte_visqueuse.out', unpack=True, usecols=[5])

Fy = Fvy + Fpy #force per unit length of cylinder

Fpz = np.loadtxt('Beam_Unconfined_Cylinder_Free_Vibration_pb_Force_pression.out', unpack=True, usecols=[6])
Fvz = np.loadtxt('Beam_Unconfined_Cylinder_Free_Vibration_pb_Contrainte_visqueuse.out', unpack=True, usecols=[6])

Fz = Fvz + Fpz #force per unit length of cylinder

DataOut = np.column_stack((t,Fx, Fy, Fz))
np.savetxt('Numerical_force.txt', DataOut)


