# Import libraries
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

from orbit import Vector3D
from orbit import Orbit


# Creating dataset
stdGrvPrm = 3.986004418E+14
seed = 521 * 1000
pos = Vector3D(seed * 12, seed * 4, seed * 3)
vel = Vector3D(-2435, 7305, 1200)
orbit = Orbit(stdGrvPrm,  pos, vel)
preds = orbit.predict(6000)
xps = [point.x for point in preds]
yps = [point.y for point in preds]
zps = [point.z for point in preds]


# Integrating the trajectory
intgds = []
tps = 2
delTime = 1 / tps
for i in range(int(orbit.T * tps * 1.1)):
    accel = -pos * (stdGrvPrm / pos.size()**3)
    pos += vel * delTime
    vel += accel * delTime
    intgds.append(pos)
xis = [intgd.x for intgd in intgds]
yis = [intgd.y for intgd in intgds]
zis = [intgd.z for intgd in intgds]


# Creating figure
fig = plt.figure(figsize=(10,  7))
ax = plt.axes(projection="3d")


# Creating plot
ax.scatter3D(xps, yps, zps, color="green", s=0.2)
ax.scatter3D(xis, yis, zis, color="red", s=0.2)
plt.title("Orbits 3D Scatter Plot")


# show plot
plt.show()
