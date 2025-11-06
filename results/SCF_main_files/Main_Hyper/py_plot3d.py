
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys 

nam=sys.argv[1]
# Create sample data
data=np.loadtxt(str(nam))
x=data[:,0]
y=data[:,1]
z=data[:,2]
f1=data[:,3]
#f1=data[:,4]
f2=data[:,4]

fA=data[:,3]
fB=data[:,4]

threshold=0.7

grid_size_x = np.max(x) 
grid_size_y = np.max(y)  
grid_size_z = np.max(z)   


print(len(x))
#center_x=grid_size_x/2.0
#center_y=grid_size_y/2.0
#center_z=grid_size_z/2.0

center_x =grid_size_x/2.0  # Change this to your desired X coordinate
center_y = grid_size_y/2.0  # Change this to your desired Y coordinate
center_z = grid_size_z/2.0  # Change this to your desired Z coordinate
ele=45
azi=60





#sizex=10.0
#sizey=10.0
#sizez=10.0

#sizex_half=sizex/2.0
#sizey_half=sizey/2.0
#sizez_half=sizez/2.0

#x_new=[]
#y_new=[]
#z_new=[]
#f1_new=[]

#for i in range(0,len(x)):

#    xdif=abs(center_x-x[i])
#    ydif=abs(center_y-y[i])
#    zdif=abs(center_z-z[i])
#    
#    if xdif<=sizex_half:
#        x_new.append(center_x-x[i])
#    elif abs(center_x-(x[i]+sizex))<=sizex_half:
#        x_new.append(center_x-(x[i]+sizex))
#    elif abs(center_x-(x[i]-sizex))<=sizex_half:
#        x_new.append(center_x-(x[i]-sizex))

#    if ydif<=sizey_half:
#        y_new.append(center_y-y[i])
#    elif abs(center_y-(y[i]+sizey))<=sizey_half:
#        y_new.append(center_y-(y[i]+sizey))
#    elif abs(center_x-(x[i]-sizex))<=sizey_half:
#        y_new.append(center_y-(y[i]-sizey))

#    if zdif<=sizez_half:
#        z_new.append(center_z-z[i])
#    elif abs(center_z-(z[i]+sizez))<=sizez_half:
#        x_new.append(center_x-(x[i]+sizex))
#    elif abs(center_x-(x[i]-sizex))<=sizex_half:
#        x_new.append(center_x-(x[i]-sizex))


# Calculate the shifts needed to recenter
x_shift = center_x - np.mean(x)
y_shift = center_y - np.mean(y)
z_shift = center_z - np.mean(z)



x = (x + x_shift) % grid_size_x
y = (y + y_shift) % grid_size_y
z = (z + z_shift) % grid_size_z


# Filter data above threshold
#x = x[f1 > threshold]
#y = y[f1 > threshold]
#z = z[f1 > threshold]
#f1 = f1[f1 > threshold]

#Find centre of mass in new coordiantes
mass=np.sum(f1)
x_cmass=np.sum(x*f1)/mass
y_cmass=np.sum(y*f1)/mass
z_cmass=np.sum(z*f1)/mass

print(x_cmass)
print(y_cmass)
print(z_cmass)


x_ind=np.where((abs(x-x_cmass)) == (abs(x-x_cmass)).min())[0]
y_ind=np.where((abs(y-y_cmass)) == (abs(y-y_cmass)).min())[0]
z_ind=np.where((abs(z-z_cmass)) == (abs(z-z_cmass)).min())[0]


#common=np.intersect1d(y_ind,z_ind)
#data_x=np.column_stack((x[common], f1[common]))
#plt.scatter(x[common]-x_cmass,f1[common])
#common=np.intersect1d(x_ind,z_ind)
#plt.scatter(y[common]-y_cmass,f1[common])
#common=np.intersect1d(x_ind,y_ind)
#plt.scatter(z[common]-z_cmass,f1[common])
#plt.show()


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
cmap = plt.cm.get_cmap('jet')



# Filter data above threshold
#x = x[f1 > threshold]
#y = y[f1 > threshold]
#z = z[f1 > threshold]
#f1 = f1[f1 > threshold]

f1[f1 <= threshold]=0.0

print("max vol. frac.  A " +str(np.max(fA)))
print("max vol. frac.  B " +str(np.max(fB)))

# Create plot
ax.scatter(x, y, z, c=f1, cmap=cmap, s=50*f1**2, alpha=0.5, depthshade=False)
#ax.set_aspect('equal')
ax.set_box_aspect([1,1,1])
ax.view_init(elev=ele, azim=azi)
# Set axis labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_xlim(0,grid_size_x)
ax.set_ylim(0,grid_size_y)
ax.set_zlim(0,grid_size_z)

# Add colorbar
fig.colorbar(ax.collections[0], shrink=0.5, aspect=5)

# Show plot


#plt.savefig(nam[0:-20]+"3d_plot.png",dpi=300)
#plt.savefig('3dplot.pdf')
plt.show()

