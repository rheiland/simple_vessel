from pyMCDS_cells import pyMCDS_cells
# from pyMCDS import pyMCDS
import numpy as np

#mcds = pyMCDS_cells('output00000001.xml','data')
#mcds = pyMCDS_cells('output00000001.xml','.') #  23123 cells
#mcds = pyMCDS_cells('output00000246.xml','.')  # 116038 cells
#mcds = pyMCDS_cells('output00000001.xml','.')  
mcds = pyMCDS_cells('output00000000.xml','.')  
# mcds = pyMCDS('output00000000.xml','.')  
tmins = mcds.get_time()
print('time (mins)=',tmins)
print('time (days)=',tmins/1440.)

print("discrete_cells.keys()= ",mcds.data['discrete_cells'].keys())
# print("ID= ",mcds.data['discrete_cells']['ID'])
ncells = len(mcds.data['discrete_cells']['ID'])
print('num cells originally = ',ncells)

#xyz = np.empty((ncells,3))
# xyz = np.zeros((ncells,3))
xvals = mcds.data['discrete_cells']['position_x']
# print("position_x=",mcds.data['discrete_cells']['position_x'])
yvals = mcds.data['discrete_cells']['position_y']
zvals = mcds.data['discrete_cells']['position_z']
print("x range: ",xvals.min(),xvals.max())
print("y range: ",yvals.min(),yvals.max())
print("z range: ",zvals.min(),zvals.max())
