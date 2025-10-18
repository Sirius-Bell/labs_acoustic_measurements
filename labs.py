"""Import matplotlib and math"""
from math import sqrt,pi
import matplotlib.pyplot as plt

global R, ug
R = 10**3
ug = 6

def calc_ct(f, uz, ur):
    """Calculate ct"""
    w = 2*pi*f
    return 1/(w*R*sqrt( (uz/ur)**2 - 0.25* ( ((ug**2 - uz**2)/ur**2) - 1)**2 ))

print("-----------------------")
uz_all = [4.17, 4.11, 4.03, 3.96, 3.85, 3.7, 3.49, 3.24, 3.05, 3.14, 2.91, 2.76, 3.9, 4.5, 4.53, 4.43, 4.34, 4.29, 4.22, 4.16, 4.1] # r
ur_all = [4, 4.06, 4.11, 4.17, 4.23, 4.3, 4.34, 4.29, 4.14, 3.92, 4.01, 3.6, 2.44, 2.5, 3.05, 3.37, 3.66, 3.71, 3.8, 3.86, 3.92] # r

f_all = [i for i in range(70,91)] # r
Z_all = []

for i in range(len(uz_all)):
    Z_all.append(R*uz_all[i]/ur_all[i])
    print(f"X{i+1}={Z_all[i]:.2f}")

print("-----------------------")

uz_p = 2.612 # r
ur_p = 3.82 # r

uz_a = 4.58 # r
ur_a = 3.82 # r

fp = 80.6*10**3 # r
fa = 83.5*10**3 # r
zp = R * uz_p/ur_p
za = R * uz_a/ur_a

Xp = R*sqrt( (uz_p / ur_p)**2 - 0.25* ( ((ug**2 - uz_p**2)/ur_p**2) - 1)**2 )
Xa = R*sqrt( (uz_a / ur_a)**2 - 0.25* ( ((ug**2 - uz_a**2)/ur_a**2) - 1)**2 )

gp = 1/zp
ga = 1/za

K = sqrt(1 - (fp/fa)**2)

Ct = calc_ct(70*(10**3), 4.17, 4) # r
Cd = K**2 * Ct
Cx = Ct*(1-K**2)
wp = 2*pi*fp
Ld = 1 / ( wp**2 * K**2 * Ct )

e = 1500
p = 7300
l = 20*10**(-3)

d = sqrt( (pi*e*(fa-fp)*10**(-9)) / ( 576*p*l**2 * fp**3 ) ) # r

print(f"Zp: {zp:}")
print(f"Za: {za}")
print(f"Xp: {Xp}")
print(f"Xa: {Xa}")
print(f"gp: {gp}")
print(f"ga: {ga}")
print(f"Cx: {Cx}")
print(f"Ct: {Ct}")
print(f"Cd: {Cd}")
print(f"Ld: {Ld}")
print(f"K: {K}")
print(f"d: {d}")


fig, axs = plt.subplots(3)

axs[0].plot(f_all, uz_all, linewidth=2.0)
axs[0].grid(True)
axs[0].set_ylabel('Uz, В')

axs[1].plot(f_all, ur_all, linewidth=2.0)
axs[1].grid(True)
axs[1].set_ylabel('Ur, В')

axs[2].plot(f_all, Z_all, linewidth=2.0)
axs[2].grid(True)
axs[2].set_ylabel('Z, Ом')
axs[2].set_xlabel('f, кГц')

plt.show()
