import numpy as np
import matplotlib.pyplot as plt
# x = np.linspace(0, 6000, 50)
# y = np.linspace(0, 6000, 50)
# z = np.linspace(0, 4000, 50)
dx=100
dy=100
dz=100
x = np.linspace(dx/2, 6000-dx/2, int(6000/dx))
y = np.linspace(dy/2, 6000-dy/2, int(6000/dy))
z = np.linspace(dz/2, 3500-dz/2, int(3500/dz))

X, Y, Z = np.meshgrid(x, y, z)
xf = X.flatten()
yf = Y.flatten()
zf = Z.flatten()

xc = np.array([[3000], [3000], [1250]])

angle_x = np.deg2rad(35)
Rx = np.array(
    [
        [1, 0, 0],
        [0, np.cos(angle_x), -np.sin(angle_x)],
        [0, np.sin(angle_x), np.cos(angle_x)],
    ]
)

angle_y = np.deg2rad(0)
Ry = np.array(
    [
        [np.cos(angle_y), 0, np.sin(angle_y)],
        [0, 1, 0],
        [-np.sin(angle_y), 0, np.cos(angle_y)],
    ]
)

angle_z = np.deg2rad(0)
Rz = np.array(
    [
        [np.cos(angle_z), -np.sin(angle_z), 0],
        [np.sin(angle_z), np.cos(angle_z), 0],
        [0, 0, 1],
    ]
)
R = Rz @ Ry @ Rx
a = 800
b = 1500
c = 400
A = np.array([[1 / (a**2), 0, 0], [0, 1 / (b**2), 0], [0, 0, 1 / (c**2)]])
f = []
for xi, yi, zi in zip(xf, yf, zf):
    xv = np.array([[xi], [yi], [zi]])
    fv=(xv - xc).T @ R.T @ A @ R @ (xv - xc)
    # print(fv.shape)
    f.append(fv[0,0])
f = np.array(f)
#print(f.shape)
xf = xf[f < 1]
yf = yf[f < 1]
zf = zf[f < 1]

with open("anomalies.txt","w") as f:
    for xi,yi,zi in zip(xf,yf,zf):
        f.write(f"{xi-dx/2} {xi+dx/2} {yi-dy/2} {yi+dy/2} {zi-dz/2} {zi+dz/2} 500\n")

    # f.write(f"2000 4000 1125 2875 1750 2125 300\n")
    # f.write(f"2250 3750 2250 3750 1750 1875 -400\n")
    # f.write(f"2000 4000 2000 4000 1875 2000 -300\n")
    # f.write(f"1800 4200 1800 4200 2000 2125 -200\n")
    # f.write(f"1500 4500 1500 4500 2125 2375 -100\n")


# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(xf, yf, zf)
# ax.set_xlim([min(x),max(x)])
# ax.set_ylim([min(y),max(y)])
# ax.set_zlim([min(z),max(z)])
# ax.set_xlabel("x")
# ax.set_ylabel("y")
# ax.set_zlabel("z")
# plt.show()
