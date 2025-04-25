import numpy as np
import matplotlib.pyplot as plt

# Step 1: Load the data from CSV
data = np.loadtxt("/home/user/BEACLS/beacls/sources/samples/plane4D1/all_loop.csv", delimiter=",")  # Ensure the correct delimiter is used

# Grid parameters
Nx = 11  # Number of points in x-direction


# Ensure data size matches the expected grid size
assert data.size == Nx * Nx * Nx * Nx, f"Data size {data.size} does not match expected size {Nx * Nx * Nx * Nx}."

# Reshape the data into a 3D grid
data_3d = data.reshape((Nx, Nx, Nx, Nx)).transpose(3, 2, 1, 0)

	# beacls::FloatVec mins{ (FLOAT_TYPE)-10, (FLOAT_TYPE)-5, (FLOAT_TYPE)-5,      (FLOAT_TYPE)-5 };
	# beacls::FloatVec maxs{ (FLOAT_TYPE)2,(FLOAT_TYPE)5,(FLOAT_TYPE)5,(FLOAT_TYPE)5 };
# Step 2: Map grid indices to physical coordinates
mins = [-10, -5]
maxs = [2, 5]

# Generate the physical coordinates for each dimension
x_vals = np.linspace(mins[0], maxs[0], Nx)
y_vals = np.linspace(mins[1], maxs[1], Nx)


# Step 3: Find indices with negative values in the grid
neg_indices = np.argwhere(data_3d < 0)
# print(neg_indices)
# Map indices to physical coordinates
x_coords = x_vals[neg_indices[:, 0]]
y_coords = y_vals[neg_indices[:, 1]]

# Step 4: Visualize the negative value points (3D scatter plot)
fig = plt.figure(figsize=(10, 7))
# ax = fig.add_subplot(111, projection='3d')
ax = fig.add_subplot(111)

# Scatter plot for points with negative values
ax.scatter(x_coords, y_coords, c='red', s=5, alpha=0.6)

# Labels and title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('Points with Negative Values in all_loop.csv')

# Compute ranges
x_range = x_vals.max() - x_vals.min()
y_range = y_vals.max() - y_vals.min()
xy_max_range = max(x_range, y_range)

# Compute midpoints
x_mid = (x_vals.max() + x_vals.min()) / 2
y_mid = (y_vals.max() + y_vals.min()) / 2

# Set equal limits for X and Y
ax.set_xlim(x_mid - xy_max_range / 2 - 5, x_mid + xy_max_range / 2 + 5)
ax.set_ylim(y_mid - xy_max_range / 2 - 5, y_mid + xy_max_range / 2 + 5)


# Show the plot
plt.tight_layout()
plt.show()

# # Step 5: Optionally visualize a 2D slice at a fixed z (theta) value
# fixed_z_index = Nz // 2  # You can change this index to select a specific theta slice
# slice_data = data_3d[:, :, fixed_z_index]

# # Step 6: Visualize the slice (fixed θ, 2D plot)
# plt.figure(figsize=(8, 6))
# plt.imshow(slice_data, cmap='viridis', extent=[x_vals[0], x_vals[-1], y_vals[0], y_vals[-1]], origin='lower')

# # Set equal aspect ratio for the axes in the 2D plot
# plt.gca().set_aspect('equal', adjustable='box')

# plt.colorbar(label='Value')
# plt.title(f"2D Slice at Theta = {z_vals[fixed_z_index]:.2f}")
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.show()