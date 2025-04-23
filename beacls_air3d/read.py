import numpy as np
import matplotlib.pyplot as plt

# Step 1: Load the data from CSV
data = np.loadtxt("/home/kiyong/BEACLS/beacls/sources/samples/air3D/all_loop.csv", delimiter=",")  # Ensure the correct delimiter is used

# Grid parameters
Nx = 71  # Number of points in x-direction
Ny = int(np.ceil(Nx * (10 - (-10)) / (20 - (-6))))  # Ny calculated as per your formula, should be 40
Nz = Nx - 1  # Nz = Nx - 1, should be 50

# Ensure data size matches the expected grid size
assert data.size == Nx * Ny * Nz, f"Data size {data.size} does not match expected size {Nx * Ny * Nz}."

# Reshape the data into a 3D grid
data_3d = data.reshape((Nz, Ny, Nx)).transpose(2, 1, 0)


# Step 2: Map grid indices to physical coordinates
mins = [-6, -10, 0]
maxs = [20, 10, 2 * np.pi]

# Generate the physical coordinates for each dimension
x_vals = np.linspace(mins[0], maxs[0], Nx)
y_vals = np.linspace(mins[1], maxs[1], Ny)
z_vals = np.linspace(mins[2], maxs[2], Nz)

# Step 3: Find indices with negative values in the grid
neg_indices = np.argwhere(data_3d < 0)

# Map indices to physical coordinates
x_coords = x_vals[neg_indices[:, 0]]
y_coords = y_vals[neg_indices[:, 1]]
z_coords = z_vals[neg_indices[:, 2]]

# Step 4: Visualize the negative value points (3D scatter plot)
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# Scatter plot for points with negative values
ax.scatter(x_coords, y_coords, z_coords, c='red', s=5, alpha=0.6)

# Labels and title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Theta (rad)')
ax.set_title('Points with Negative Values in all_loop.csv')

# Compute ranges
x_range = x_vals.max() - x_vals.min()
y_range = y_vals.max() - y_vals.min()
xy_max_range = max(x_range, y_range)

# Compute midpoints
x_mid = (x_vals.max() + x_vals.min()) / 2
y_mid = (y_vals.max() + y_vals.min()) / 2

# Set equal limits for X and Y
ax.set_xlim(x_mid - xy_max_range / 2, x_mid + xy_max_range / 2)
ax.set_ylim(y_mid - xy_max_range / 2, y_mid + xy_max_range / 2)

# Optional: keep original Z limits
ax.set_zlim(z_vals.min(), z_vals.max())

# Set box aspect to [1, 1, Z_ratio] to avoid squeezing Z
z_range = z_vals.max() - z_vals.min()
z_ratio = z_range / xy_max_range
ax.set_box_aspect([1, 1, 0.5])


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
