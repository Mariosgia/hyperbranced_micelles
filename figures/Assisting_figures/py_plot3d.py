import numpy as np
import pyvista as pv
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from numpy import linalg as LA

def calc_asph(fil): #calculates asphericities based on filename
    
#     cut_off=7**2
    cut_off=0.05
    dat=np.loadtxt(fil,skiprows=0)
    x=dat[:,0]
    y=dat[:,1]
    z=dat[:,2]
    phi1=dat[:,3] #Only takeing phi1
    phi2=dat[:,4] #Only takeing phi1


    grid_size_x = np.max(x) 
    grid_size_y = np.max(y)  
    grid_size_z = np.max(z)   
    center_x =grid_size_x/2.0  # Change this to your desired X coordinate
    center_y = grid_size_y/2.0  # Change this to your desired Y coordinate
    center_z = grid_size_z/2.0  # Change this to your desired Z coordinate

    x_new=x-center_x
    y_new=y-center_y
    z_new=z-center_z
    rad_squared=x_new**2+y_new**2+z_new**2
    
#     comb=[x_new[rad_squared<=cut_off],y_new[rad_squared<=cut_off],z_new[rad_squared<=cut_off]]
#     phi1=phi1[rad_squared<=cut_off]
#     phi2=phi2[rad_squared<=cut_off]
#     rad_squared=rad_squared[rad_squared<=cut_off]

    comb=[x_new[phi1>=cut_off],y_new[phi1>=cut_off],z_new[phi1>=cut_off]]
    phi2=phi2[phi1>=cut_off]
    rad_squared=rad_squared[phi1>=cut_off]
    phi1=phi1[phi1>=cut_off]

    matrx1=np.zeros((3,3))
    matrx2=np.zeros((3,3))

    for row in range(0,3): #Moment of inertia tensor
        for col in range(0,3):


            xi=comb[row]
            xj=comb[col]
            summ=-xi*xj

            if row==col:
                summ+=rad_squared

            summ1=np.sum(summ*phi1)
            summ2=np.sum(summ*phi2)

            norm1=np.sum(phi1)
            norm2=np.sum(phi2)

            matrx1[row,col]=summ1/norm1
            matrx2[row,col]=summ2/norm2


    ei1, eivec1 = LA.eig(matrx1) #eigenvalues
    ei2, eivec2 = LA.eig(matrx2)

    asph1=np.real(ei1[2]**2-0.5*(ei1[0]**2+ei1[1]**2))
    asph2=np.real(ei2[2]**2-0.5*(ei2[0]**2+ei2[1]**2))
    return asph1,asph2


#Asphericities kap : 0.0 asp = 0 , kap : 0.16 asp = 3.2275287734987046
names = ["../../results/micelles/linear/SCF_runs_3d/phi_runs/chosen_case/phibar_0.00125/stretch_runs/asp_ratio_1.0/data/volume_fractions.dat",
            "../../results/micelles/linear/SCF_runs_3d/phi_runs/chosen_case/phibar_0.00125/stretch_runs/asp_ratio_2.4/data/volume_fractions.dat"]



# Define box dimensions
# Define the fixed bounds of your simulation box
box_min, box_max = 0, 15  

# Create the box as a mesh using PyVista's Box
simulation_box = pv.Box(bounds=(box_min, box_max, box_min, box_max, box_min, box_max))

cube_size=15.0
# Create a cube with float coordinates
cube_points = np.array([
    [0.0, 0.0, 0.0], 
    [cube_size, 0.0, 0.0], 
    [cube_size, cube_size, 0.0], 
    [0.0, cube_size, 0.0],
    [0.0, 0.0, cube_size], 
    [cube_size, 0.0, cube_size], 
    [cube_size, cube_size, cube_size], 
    [0.0, cube_size, cube_size]
], dtype=float)  # Specify dtype as float

# Create a smaller inner cube (ensuring all points are inside the outer cube)
inner_cube_size = 5  # Size of the inner cube (should be less than cube_size)
inner_cube_half_size = inner_cube_size / 2  # Half size for centering
inner_cube_points = np.array([
    [cube_size / 2 - inner_cube_half_size, cube_size / 2 - inner_cube_half_size, cube_size / 2 - inner_cube_half_size],
    [cube_size / 2 + inner_cube_half_size, cube_size / 2 - inner_cube_half_size, cube_size / 2 - inner_cube_half_size],
    [cube_size / 2 + inner_cube_half_size, cube_size / 2 + inner_cube_half_size, cube_size / 2 - inner_cube_half_size],
    [cube_size / 2 - inner_cube_half_size, cube_size / 2 + inner_cube_half_size, cube_size / 2 - inner_cube_half_size],
    [cube_size / 2 - inner_cube_half_size, cube_size / 2 - inner_cube_half_size, cube_size / 2 + inner_cube_half_size],
    [cube_size / 2 + inner_cube_half_size, cube_size / 2 - inner_cube_half_size, cube_size / 2 + inner_cube_half_size],
    [cube_size / 2 + inner_cube_half_size, cube_size / 2 + inner_cube_half_size, cube_size / 2 + inner_cube_half_size],
    [cube_size / 2 - inner_cube_half_size, cube_size / 2 + inner_cube_half_size, cube_size / 2 + inner_cube_half_size]
], dtype=float)  # Specify dtype as float
# Define the edges of the cube using tuples of points

edges = [
    (inner_cube_points[0], inner_cube_points[1]),
    (inner_cube_points[1], inner_cube_points[2]),
    (inner_cube_points[2], inner_cube_points[3]),
    (inner_cube_points[3], inner_cube_points[0]),
    (inner_cube_points[4], inner_cube_points[5]),
    (inner_cube_points[5], inner_cube_points[6]),
    (inner_cube_points[6], inner_cube_points[7]),
    (inner_cube_points[7], inner_cube_points[4]),
    (inner_cube_points[0], inner_cube_points[4]),
    (inner_cube_points[1], inner_cube_points[5]),
    (inner_cube_points[2], inner_cube_points[6]),
    (inner_cube_points[3], inner_cube_points[7])
]
pv.global_theme.font.fmt = '%.2f'

# Create lines for the cube edges
cube_lines = pv.PolyData()  # Create a PolyData object to hold the lines
for start_point, end_point in edges:
    line = pv.Line(start_point, end_point)  # Create a line between each pair of points
    cube_lines += line  # Append the line to the PolyData object
for i in range(0,len(names)):
    nam=names[i]

    asph1,asph2=calc_asph(nam)
    
    out_nam='Asp_'+str(round(asph1,2))+'.jpg'
    print(out_nam)
    data = np.loadtxt(str(nam))

    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    f1 = data[:, 3]  # Use f1 for the isosurface color

    x = x.astype(float)
    y = y.astype(float)
    z = z.astype(float)
    f1 = f1.astype(float)

    grid_size_x = np.max(x) 
    grid_size_y = np.max(y)  
    grid_size_z = np.max(z)   







    center_x =grid_size_x/2.0  # Change this to your desired X coordinate
    center_y = grid_size_y/2.0  # Change this to your desired Y coordinate
    center_z = grid_size_z/2.0  # Change this to your desired Z coordinate

    threshold=0.4
    #f1[f1<=threshold]=np.nan
    f1[(x>center_x)& (y>center_y)  ]=np.nan


    x = np.unique(x)
    y = np.unique(y)
    z = np.unique(z)

    f1=f1.reshape(len(x),len(y),len(z))
    # Create a mesh grid
    x, y, z = np.meshgrid(x, y, z)

 

    # Flatten the grid to create a structured grid
    points = np.vstack((x.ravel(), y.ravel(), z.ravel())).T

    # Create a PyVista StructuredGrid
    grid = pv.StructuredGrid(x, y, z)

    # Add the function values as scalars
    grid.point_data['f1'] = f1.ravel()

    # Plotting
    plotter = pv.Plotter(off_screen=True)
    # Define specific contour values (for example: f1 = 0.5, 1.0, 1.5)
    contour_values = [0.27,0.4,0.45,0.5,0.53,0.6,0.65,0.7,0.73,0.77]  # You can change these values to your specific levels

    # Create contours only for these specified values
    contours = grid.contour(isosurfaces=contour_values)

    # Check if contours were created
    if contours.n_points == 0:
        raise ValueError("No contours created from the structured grid.")

    # Add contours to the plot
    #opacity_values = np.linspace(0.3, 1, len(contours.points))  # Example custom opacity values
    plotter.add_mesh(contours, scalars='f1', show_edges=False, opacity=1,cmap='jet',show_scalar_bar=False)
    plotter.add_scalar_bar( n_labels=5, vertical=True,position_x=0.8, position_y=0.28,label_font_size= 90)


    plotter.add_mesh(cube_lines, color='black', line_width=5)
    # Add text labels for dimensions near the bounding box (in 3D space)
#    plotter.add_text(f"Width: {box_length_x:.2f}", position="lower_left", color="black", font_size=12)
#    plotter.add_text(f"Height: {box_length_y:.2f}", position="lower_right", color="black", font_size=12)
#    plotter.add_text(f"Depth: {box_length_z:.2f}", position="upper_left", color="black", font_size=12)

    # Set camera view with specific azimuth and elevation
    azimuth = 30  # Rotate camera 45 degrees around the vertical axis
    elevation = 20 # Raise camera 30 degrees from horizontal
    plotter.view_yz()  # Set the initial view to the XY plane
    plotter.camera.azimuth = azimuth
    plotter.camera.elevation = elevation

    plotter.screenshot(out_nam, scale=10)  # Increase scale for higher quality



