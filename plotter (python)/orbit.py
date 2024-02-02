import math


class Vector3D:

    def __init__(self, *args):
        if (len(args) == 3):
            self.x = args[0]
            self.y = args[1]
            self.z = args[2]
        else:
            raise IndexError

    def __neg__(self):
        return Vector3D(
            -self.x,
            -self.y,
            -self.z
        )

    def __add__(self, other):
        return Vector3D(
            self.x + other.x,
            self.y + other.y,
            self.z + other.z
        )

    def __sub__(self, other):
        return Vector3D(
            self.x - other.x,
            self.y - other.y,
            self.z - other.z
        )

    def __mul__(self,  other):
        if type(other) in [int, float]:
            return Vector3D(
                self.x * other,
                self. y * other,
                self.z * other
            )
        else:
            return self.x * other.x +\
                self.y * other.y +\
                self.z * other.z

    def __matmul__(self, other):
        return Vector3D(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x
        )

    def __repr__(self):
        return '(%s, %s, %s)' %\
            (self.x, self.y, self.z)

    def size(self):
        return (self.x**2 +
                self.y**2 +
                self.z**2)**.5

    def rotateX(self, the):
        return Vector3D(
            self.x,
            self.y * math.cos(the) - self.z * math.sin(the),
            self.z * math.cos(the) + self.y * math.sin(the)
        )

    def rotateY(self, the):
        return Vector3D(
            self.x * math.cos(the) + self.z * math.sin(the),
            self.y,
            self.z * math.cos(the) - self.x * math.sin(the)
        )

    def rotateZ(self, the):
        return Vector3D(
            self.x * math.cos(the) - self.y * math.sin(the),
            self.y * math.cos(the) + self.x * math.sin(the),
            self.z
        )


class Orbit:

    def __init__(self, mu, vec_r_0, vec_v_0):
        r_0 = vec_r_0.size()
        v_0 = vec_v_0.size()
        a = mu * r_0 / (2 * mu - r_0 * v_0**2)
        vec_H = vec_r_0 @ vec_v_0
        H = vec_H.size()
        e = (1 - H**2 / (mu * a))**.5
        T = 2 * math.pi * (a**3 / mu)**.5

        i = math.acos(vec_H.z / H)
        OM = math.atan2(vec_H.y, vec_H.x)
        vec_r_1 = vec_r_0.rotateZ(-OM).rotateY(-i)
        the_0 = math.acos((H**2 / (mu * r_0) - 1) / e)
        om = math.atan2(vec_r_1.y, vec_r_1.x)
        sgn = 1 if vec_r_0 * vec_v_0 < 0 else -1
        om += sgn * the_0

        self.a, self.e, self.T = a, e, T
        self.i, self.OM, self.om = i,  OM, om

    def predict(self, frags):
        dots = []
        for i in range(frags):
            the = 2 * math.pi * i / frags
            r = self.a * (1 - self.e**2) /\
                (1 + self.e * math. cos(the))
            dots.append(Vector3D(r * math.cos(the),
                                 r * math.sin(the),
                                 0))
            dots[-1] = dots[-1].rotateZ(self.om)\
                .rotateY(self.i).rotateZ(self.OM)
        return dots


if __name__ == "__main__":
    stdGrvPrm = 3.986004418E+14
    seed = 521 * 1000
    pos = Vector3D(seed * 12, seed * 4, seed * 3)
    vel = Vector3D(-2435, 7305, 1200)
    orbit = Orbit(stdGrvPrm,  pos, vel)
