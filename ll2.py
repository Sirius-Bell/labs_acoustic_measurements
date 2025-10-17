from math import pi
from cmath import sqrt
import matplotlib.pyplot as plt

global R, ug
R = 10**3
ug = 6

def calc_ct(f, uz, ur):
    w = 2*pi*f
    return 1/(w*R*sqrt( (uz/ur)**2 - 0.25* ( ((ug**2 - uz**2)/ur**2) - 1)**2 ))

print("-----------------------")
uz_all = [1.8,1.7,1.6,1.54,1.46,1.38,1.31,1.23,1.15,1.06,0.95,0.8,0.544,0.16,3.37,1.96,1.72,1.47,1.32,1.23,1.16] # r
ur_all = [5.56,5.58,5.6,5.62,5.64,5.65,5.67,5.68,5.7,5.71,5.73, 5.74,5.77,5.72,2.6,5.21,5.53,5.6,5.65,5.67,5.68] # r

f_all = [i * 0.5 for i in range(20, 41)]# r

Z_all = []

for i in range(len(uz_all)):
    Z_all.append(R*uz_all[i]/ur_all[i])
    print(f"X{i+1}={Z_all[i]:.2f}")

print("-----------------------")

uz_p = 0.05 # r
ur_p = 5.78 # r

uz_a = 4.74 # r
ur_a = 1.41 # r

fp = 16.43*10**3 # r
fa = 17.07*10**3 # r
zp = R * uz_p/ur_p
za = R * uz_a/ur_a

Xp = R*sqrt( (uz_p / ur_p)**2 - 0.25* ( ((ug**2 - uz_p**2)/ur_p**2) - 1)**2 )
Xa = R*sqrt( (uz_a / ur_a)**2 - 0.25* ( ((ug**2 - uz_a**2)/ur_a**2) - 1)**2 )

gp = 1/zp
ga = 1/za

K = sqrt(1 - (fp/fa)**2)

Ct = calc_ct(10*(10**3), 1.8, 5.56) # r
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