
import medcoupling as mc
import MEDLoader as ml
import numpy as np
from math import *
import sys, os

def removeDuplicates(interf):
   interf.zipCoords() 
   delta=interf.getCaracteristicDimension()
   eps=0.0001*delta
   interf.mergeNodes(eps) # Supprime les noeuds en doublons pour la joonction de proc!!
   cdg_interf = interf.computeCellCenterOfMass()
   interf.zipCoords() 
   coords = interf.getCoords() 
   x=coords[:,0]
   y=coords[:,1]
   y_TCL, i_TCL  = y.getMinValue()
   elem_TCL = interf.getCellIdsLyingOnNodes([i_TCL],False)
   s0, s1 = interf.getNodeIdsOfCell((int)(elem_TCL))
   if (i_TCL == s1):
      print(f"TCL is the second of vertex of the elem {(int)(elem_TCL)}, so we change the elements orientation") 
      ret = interf.changeOrientationOfCells()
      if ret:
         print("some nodes have been merged")# and hence this field lies on another mesh.")
      pass
   list_bord = interf.findBoundaryNodes()
   if (not interf.isContiguous1D()):
      commons, commons_id = interf.findCommonCells(0)
      cells_to_keep = list(range(interf.getNumberOfCells()))
      if len(commons)>0:
         for i, elem in enumerate(commons_id[:-1]):
            for j, common in enumerate(commons[commons_id[i]:commons_id[i+1]]):
               if (j>0):
                  # p common
                  # <MEDLoader.DataArrayInt32Tuple; proxy of <Swig Object of type 'MEDCoupling::DataArrayInt32Tuple *' at 0x75ea5390ff00> >
                  cells_to_keep.remove((int)(common))
               pass
            pass
         interf = interf.buildPartOfMySelf(cells_to_keep)
   return True, interf
"""
         if (not interf.isContiguous1D()):
            print('Je ne comprends pas pourquoi... ')
         pass
      pass
   
   some_0D_mesh, _d, _dI, _rD, _rDI =  interf.buildDescendingConnectivity()
   _dsi = _rDI.deltaShiftIndex()
   dsii = _dsi.findIdsNotInRange(0,3)
   if (dsii.getNumberOfTuples()):
      raise Exception("GB MEDCouplingUMesh::orderConsecutiveCells1D only work with a mesh being a (piecewise) connected line!")

   ordered_list = interf.orderConsecutiveCells1D().toNumPyArray()
   if (ordered_list[0] != (int)(elem_TCL)): 
      # The list does not start with TCL element/
      print("Reordering list to make it start by TCL element")
      for i, idx in enumerate(ordered_list):
         c0, c1 = interf.getNodeIdsOfCell((int)(idx))
         if (i_TCL == c0) or (i_TCL == c1):
            first_cell = (int)(idx)
            break;
      # The new list should start by TCL FT element: 
      if (i!=0):
         ordered_list = np.hstack([ordered_list[i:],ordered_list[:i]])
   
   interf = dump(interf, ordered_list, 0)
   if (not interf.isContiguous1D()):
      print('Je ne comprends pas pourquoi... ')
   return True, interf
"""

def dump(interf, ordered_list, call=0, fic="mesh"):
   c1pre, _ = interf.getNodeIdsOfCell((int)(ordered_list[0])) # Initialize
   xy = np.zeros((0,2))
   xpre, ypre = interf.getCoordinatesOfNode(c1pre)
   ok = True
   for i, idx in enumerate(ordered_list):
      idx = idx.item()
      c0, c1 = interf.getNodeIdsOfCell(idx)
      x0 , y0 = interf.getCoordinatesOfNode(c0)
      x1 , y1 = interf.getCoordinatesOfNode(c1)
      # Dumps : 
      xy = np.vstack([xy, [x0, y0], [x1, y1]])
      if c1pre == c0:
         # Fine, its contiguous
         c1pre = c1 # Increment it for next passage
      else:
         ok = False
         print("Element: ", idx)
         print("Previous element Edge #1", c1pre)
         print("New element Edge #0", c0)
         fic +=str(call)+".txt"
         np.savetxt(fic, xy)
         print(f"See file {fic}")
         li = [ x.item() for x in ordered_list[:i]]
         # On continue le maillage dans un nouveau fichier : 
         try: 
            new_list = ordered_list[i+1:]
            new_mesh = dump(interf, new_list, call+1)
            del new_mesh # On ne fait rien du nouveau maillage
            if (call>3): 
               raise Exception("Too many calls")
            print(f"Exiting dump() call #{call}")
            
            break  
         except:
            raise Exception("Interfacial mesh is not contiguous") 
   if ok:
      fic +=str(call)+".txt"
      np.savetxt(fic, xy)
   else:
      # interf is modified.
      new_interf = interf.buildPartOfMySelf(li)
      new_interf.zipCoords()
      return new_interf

   # interf is not modified.
   return interf

def reorderingForContinuity(interf, check=True, printer=False):
   c = interf.getNodalConnectivity()
   ctype = type(c)
   cI =interf.getNodalConnectivityIndex()
   ordered_list = interf.orderConsecutiveCells1D()

   # Check if the first som will be the lowest:
   if check:
      coords = interf.getCoords() 
      x=coords[:,0]
      y=coords[:,1]
      y_TCL, i_TCL  = y.getMinValue()
      fa7_TCL = interf.getCellIdsLyingOnNodes([i_TCL],False)[0]
      if interf.getNodeIdsOfCell(fa7_TCL)[0] != i_TCL:
         raise Exception("The first som of the futur First elem is not the TCL! consider something like new orientation")
         interf.changeOrientationOfCells()
         ordered_list = interf.orderConsecutiveCells1D()

   """
   for i, idx in enumerate(ordered_list):
      idx = (int)(idx)
      c0, c1 = interf.getNodeIdsOfCell(idx)
      x0 , y0 = interf.getCoordinatesOfNode(c0)
      x1 , y1 = interf.getCoordinatesOfNode(c1)
      print(idx, c0, c1)
   """

   """ c'est ko ca...
   som0 = c[1::3]
   som1 = c[2::3]
   som0.renumberInPlace(ordered_list)
   som1.renumberInPlace(ordered_list)
   c[1::3] = som0
   c[2::3] = som1
   """
   # A new object (maybe int32 or 64)
   carr = ctype() # ml.DataArrayInt32()
   li = []
   for i, idx in enumerate(ordered_list):
      idx = (int)(idx)
      c0, c1 = interf.getNodeIdsOfCell(idx)
      li.extend([1,c0,c1]) # 1 is for the elem type

   carr.setValues(li)
   interf.setConnectivity(carr, cI)

   ordered_list = [(int)(x) for x in ordered_list ]
   ml.WriteUMesh("interf-gb.med", interf, True)
   if printer:
      for idx in range(interf.getNumberOfCells()):
         c0, c1 = interf.getNodeIdsOfCell(idx)
         print(idx, c0, c1)
   
   return True, interf

def buildDummyMesh(Lx=0.002, Ly=0.004):
    """

    """
    xmin, xmax = 0., Lx
    ymin, ymax = 0., Ly
    Nx = 200 # 200
    Ny = 370 # 400
    coordsX = mc.DataArrayDouble.New()
    v = np.linspace(xmin, xmax, num=Nx + 1, endpoint=True)
    coordsX.setValues(v.tolist())
    coordsY = mc.DataArrayDouble.New()
    v = np.linspace(ymin, ymax, num=Ny + 1, endpoint=True)
    coordsY.setValues(v.tolist())
    mC = mc.MEDCouplingCMesh.New()
    mC.setCoords(coordsX, coordsY)
    mC.setName("Cartesian_Box")
    mU = mC.buildUnstructured()
    mU.setName("Unstructured_Box")
    return mC, mU

def myRemap(srcF,trg):
   # préparer le remapper pour une projection d'un champ au celulles vers un autre champ de celulles
   rem=mc.MEDCouplingRemapper()
   src = srcF.getMesh()
   rem.prepare(src,trg,"P0P0")
   srcF.setNature(mc.IntensiveMaximum)
   # effectuer la projection trgF sera le nouveau champ basé sur le nouveau maillage
   trgF=rem.transferField(srcF,-1)
   return trgF

def readFileInfo(fic="info.txt"):
    data = {}
    with open(fic) as f:
        for line in f.readlines():
            key, value = line.rstrip("\n").split("=")
            data[key] = float(value)
    return data

########################################
# Main 
########################################
if os.path.isfile("info.txt"):
   data = readFileInfo("info.txt")
elif os.path.isfile("../info.txt"):
   data = readFileInfo("../info.txt")
elif os.path.isfile("../../info.txt"):
   data = readFileInfo("../../info.txt")
elif os.path.isfile("../../../info.txt"):
   data = readFileInfo("../../../info.txt")
else:
   raise Exception("File info.txt is missing")
lda = data["lda"]
Lvap = data["Lvap"]

time=sys.argv[1]
pyfile=os.path.basename(sys.argv[0])
fName = f"lata/post_t{time}_0000.med"
euler = mc.ReadUMeshFromFile(fName, "dom")
interf = mc.ReadUMeshFromFile(fName, "INTERFACES") # TODO : only one mesh is read! what about time iterations? 

if (not interf.isContiguous1D()):
   ok, interf = removeDuplicates(interf)
   if (not interf.isContiguous1D()):
      # Still not contiguous, reordering : 
      ok2, interf = reorderingForContinuity(interf)
   ml.WriteUMesh("interf.med", interf, True)

if (not interf.isContiguous1D()):
   raise Exception("Interfacial mesh is not contiguous") 

# Creation of interf2 which is the cut version of interf, 
#    cut by elements from dom (euler)
ml.WriteUMesh("interf-nvo.med", interf, True)

modifMesh = False
delta=interf.getCaracteristicDimension()
eps=0.0001*delta
if modifMesh:
   eps = 1.e-9
   mC, mU = buildDummyMesh()
   euler2, interf2, mp1, mp2 = mU.Intersect2DMeshWith1DLine(mU,interf, eps)
else:
   eps = 10e-12
   euler.mergeNodes(eps)
   euler.zipCoords()
   euler2, interf2, mp1, mp2 = euler.Intersect2DMeshWith1DLine(euler,interf, eps)


# To clean interf2 of unnecessary points:
interf2.zipCoords()
cdg_interf2 = interf2.computeCellCenterOfMass()

mp = 'MPOINT_THERMIQUE_ELEM_dom'
gradT='TEMPERATURE_GRAD_THERMIQUE_ELEM_dom'
ai = 'INTERFACIAL_AREA_ELEM_dom'
aia_part = None
mpaCFD_part = None
iteration = -1 # last iteration
if mp in mc.GetAllFieldNames(fName):
   lst_dt = mc.GetAllFieldIterations(fName,mp) # take last time-step
   fC = mc.ReadFieldCell(fName,"dom",0,mp,lst_dt[iteration][0],lst_dt[iteration][1])
   print(f"[{pyfile}] Selected Time: {lst_dt[iteration][2]}")
   if modifMesh:
      fU = myRemap(fC,mU)
      mpa = fU.getArray()
      cells,cellsIndex = mU.getCellsContainingPoints(cdg_interf2, len(cdg_interf2),eps )
   else:
      mpa = fC.getArray()
      cells,cellsIndex = euler.getCellsContainingPoints(cdg_interf2, len(cdg_interf2),eps )
   #print(cells)
   mpa_part = mpa[cells]
   if ai in mc.GetAllFieldNames(fName):
      fai = mc.ReadFieldCell(fName,"dom",0,ai,lst_dt[iteration][0],lst_dt[iteration][1])
      aia = fai.getArray()
      aia_part = aia[cells]
      pass
   if gradT in mc.GetAllFieldNames(fName):
      fgradT = mc.ReadFieldCell(fName,"dom",0,gradT,lst_dt[iteration][0],lst_dt[iteration][1])
      gradTa = fgradT.getArray()
      mpaCFD_part = lda/(Lvap)*gradTa[cells]
      pass
   # Je ne comprends pas pourquoi c'est seulemnt parfois, mais ils n'ont pas tous la meme taille donc impossible : 
   #mat = np.c_[cdg_interf2.toNumPyArray(),mpa_part.toNumPyArray(),aia_part.toNumPyArray()]
   # J'imaginais en effet que les points de cdg_interf2 qui ne sont pas dans interf sont sur les faces. 
   # Dans ce cas, avec n'importe quelle précision eps, on recupere 2 elem a chaque fois. 
   #np.savetxt("mpoint.txt", mat)
   pass

coord2 = interf2.getCoords()
xminmax, yminmax = interf2.getBoundingBox()
xmin, _ = xminmax
ymin, _ = yminmax

xysm = np.zeros((0,6))
xysm_cells = np.zeros((0,6))
s=0
# Pour ordonner les elements : 
ordered_list = interf2.orderConsecutiveCells1D()
c1pre, _ = interf2.getNodeIdsOfCell(ordered_list[0]) # Initialize
xpre, ypre = interf2.getCoordinatesOfNode(c1pre)
#print("Do we start at lower left corner of the bounding box?")
#print(xmin, xpre)
#print(ymin, ypre)

# for idx in range(interf2.getNumberOfCells()): 
for idx in ordered_list: 
   idx = idx[0]
   c0, c1 = interf2.getNodeIdsOfCell(idx)
   m = mpa_part[idx] # The mpoint of the element
   if c1pre == c0:
      # Fine, its contiguous
      c1pre = c1 # Increment it for next passage
   else:
      print("Element: ", idx)
      print("Previous element Edge #1", c1pre)
      print("New element Edge #0", c0)
      #xysm, xysm_cells = getXYSM(interf2)
      #break
      raise Exception("Interfacial mesh is not contiguous") 
   x0 , y0 = interf2.getCoordinatesOfNode(c0)
   x1 , y1 = interf2.getCoordinatesOfNode(c1)
   xi , yi  = cdg_interf2[idx].getValues() # of the bary in the element idx
   s0 = s
   edge_length = sqrt((x0-x1)**2+(y0-y1)**2)
   ai = 0.
   if aia_part is not None:
      ai = aia_part[idx] # Corresponding interface area of interf2 in the element
   else:
      ai = 2.*pi*xi*edge_length
   m2 = 0.
   if mpaCFD_part is not None:
      m2 = mpaCFD_part[idx] # Corresponding from gradT
   s += edge_length
   s1 = s
   # First segment point:
   tmp1 = np.array((x0,y0,s0,m, ai, m2))
   # second segment point:
   tmp2 = np.array((x1,y1,s1,m, ai, m2))
   # At the center:
   tmpi = np.array((xi,yi,s0+0.5*edge_length,m, ai, m2))
   # Dumps : 
   xysm = np.vstack([xysm, tmp1, tmp2])
   xysm_cells = np.vstack([xysm_cells,tmpi])
   pass

dirTime = f"T{time}"
if not os.path.isdir(dirTime):
   os.makedirs(dirTime, exist_ok=True)
np.savetxt(os.path.join(dirTime, "xysm.dat"), xysm, header="x y s m ai m_from_gradT")
np.savetxt(os.path.join(dirTime, "xysm_cells.dat"), xysm_cells, header="xi yi si m ai m_from_gradT")
# Reloading : 
xysm = np.loadtxt(os.path.join(dirTime, "xysm.dat"))
xysm_cells = np.loadtxt(os.path.join(dirTime, "xysm_cells.dat"))

vx,vy,vs,vm,va,vm2 = xysm.T
vxi,vyi,vsi,vmi,vai,vmi2 = xysm_cells.T

dT = data["DT"]
theta_app = data["theta"]*pi/180
print(f"[{pyfile}] CASE: dT={dT} -- theta_app=", data["theta"])
# Analytical solution from Vadim:
smin = 2e-7
vs_ana = np.linspace(smin,vs.max(),201)
# 1. without interfacial resistance : 
qi_ana = lda*dT/(vs_ana*theta_app)
mp_ana = qi_ana/Lvap
# 2. with interfacial resistance : 
Ri=6.38209e-08
qi_anar = lda*dT/(vs_ana*theta_app+Ri*lda)
mp_anar = qi_anar/Lvap
np.savetxt("sm_ana.dat", np.vstack([vs_ana,mp_ana, mp_anar]).T, header="s m_ana, mp_ana(Ri)")

import matplotlib.pyplot as plt
dirname = os.path.basename(os.getcwd())
if (len(dirname)>1) and (dirname[0] == "R") and (dirname[1:].isdigit()):
   mesh_lvl, tcl_ext = os.getcwd().split('/')[-3:-1]
else:
   mesh_lvl, tcl_ext = os.getcwd().split('/')[-2:]

dy = 0.5e-6 # default value
if mesh_lvl[1:].replace(".", "", 1).isdigit():
   dy = float(mesh_lvl[1:])*1e-6

if tcl_ext[:2] == "NO" or "TCL" not in tcl_ext:
   nmeso = 0
else:
   nmeso=int(tcl_ext[3:])

yM = dy*nmeso
idx_yM = np.argmin(abs(vy-yM))
sM = vs[idx_yM]
mM = vm[idx_yM]

plt.figure()
plt.title("mp")
plt.xlabel(r"$s$ [m]")
plt.ylabel(r"$\dot{m}$ [kg/m²/s]")
plt.plot(vs, vm,"k-", label="Simu")
plt.plot(vs, vm2,"b--", label="Simu gradT")
plt.plot(vs_ana, mp_ana,"r-", label="Analytical")
plt.plot(vs_ana, mp_anar,"r--", label=f"Ana. Ri={Ri}")
plt.plot(sM, mM,"ko", label="M (meso/macro)")
plt.legend(loc=0)
plt.savefig(os.path.join(dirTime, "mp.png"))

plt.figure()
plt.title(f"{mesh_lvl} {tcl_ext} - mp - t={time}s")
plt.xlabel(r"$s$ [m]")
plt.ylabel(r"$-\dot{m}$ [kg/m²/s]")
plt.semilogy(vs, -vm,"k-", label="Simu")
plt.semilogy(vs, -vm2,"b--", label="Simu gradT")
plt.semilogy(vs_ana, -mp_ana,"r-", label="Analytical")
plt.semilogy(vs_ana, -mp_anar,"r--", label=f"Ana. Ri={Ri}")
plt.semilogy(sM, -mM,"ko", label="M (meso/macro)")
plt.legend(loc=0)
plt.savefig(os.path.join(dirTime, "mp-logy.png"))

plt.figure()
plt.title("ai")
plt.xlabel(r"$s$ [m]")
plt.ylabel(r"$a_i$ [m²]")
plt.plot(vs, va,"k-", label="over element")
plt.plot(vsi, vai,"ko", label="at centre")
plt.legend(loc=0)
plt.savefig(os.path.join(dirTime, "ai.png"))

plt.figure()
plt.title("mp*ai")
plt.xlabel(r"$s$ [m]")
plt.ylabel(r"$\dot{m} a_i$ [kg/s]")
plt.plot(vs, vm*va,"k-", label="over element")
plt.plot(vsi, vmi*vai,"ko", label="at centre")
plt.legend(loc=0)
plt.savefig(os.path.join(dirTime, "mpai.png"))

plt.figure()
plt.title("mp*ai")
plt.xlabel(r"$s$ [m]")
plt.ylabel(r"$-\dot{m} a_i$ [kg/s]")
plt.semilogy(vs, -vm*va,"k-", label="over element")
plt.semilogy(vsi, -vmi*vai,"ko", label="at centre")
plt.legend(loc=0)
plt.savefig(os.path.join(dirTime, "mpai-logy.png"))
#plt.show()

# print("THE END")

