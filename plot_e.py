from  matplotlib import pyplot as plt
import numpy as np
from scipy.stats import norm

# In order to run this program, please execute for all the different w in 2e.cpp
figfile = "w_all.jpeg"
dataA1 = np.loadtxt("output_e_w1_A")
dataR1 = np.loadtxt("output_e_w1_R")

dataA2 = np.loadtxt("output_e_w2_A")
dataR2 = np.loadtxt("output_e_w2_R")

dataA3 = np.loadtxt("output_e_w3_A")
dataR3 = np.loadtxt("output_e_w3_R")

dataA4 = np.loadtxt("output_e_w4_A")
dataR4 = np.loadtxt("output_e_w4_R")

N = sum(1 for line in open("output_e_w1_A"))

min_i = 1000
h = 10/(N+1);
hh = h*h;

rho = np.zeros(200)
for i in range(0, 200):
	rho[i] = (i+1)*h;
	if dataA1[i][i] < min_i:
		min_i = dataA1[i][i]
		k = i

vec = dataR1[:, k]
mag =  np.dot(vec,vec)
vec = vec/mag

plt.plot(rho, abs(vec)**2, label=r"$n=%i$"%k)

min_i = 1000
for i in range(0, 200):
	if dataA2[i][i] < min_i:
		min_i = dataA2[i][i]
		k = i

vec = dataR2[:, k]
mag =  np.dot(vec,vec)
vec = vec/mag

plt.plot(rho, abs(vec)**2, label=r"$n=%i$"%k)

min_i = 1000
for i in range(0, 200):
	if dataA3[i][i] < min_i:
		min_i = dataA3[i][i]
		k = i

vec = dataR3[:, k]
mag =  np.dot(vec,vec)
vec = vec/mag

plt.plot(rho, abs(vec)**2, label=r"$n=%i$"%k)
min_i = 1000
for i in range(0, 200):
	if dataA4[i][i] < min_i:
		min_i = dataA4[i][i]
		k = i

vec = dataR4[:, k]
mag =  np.dot(vec,vec)
vec = vec/mag

plt.plot(rho, abs(vec)**2, label=r"$n=%i$"%k)

plt.xlabel(r'$\rho$')
plt.ylabel(r'$|\psi(\rho)|^2$')
plt.title(r"Ground state for different frequencies $\omega_r$")
plt.legend([r"$\omega_1=0.01$", r"$\omega_2=0.5$" , r"$\omega_3=1.0$", r"$\omega_4=5.0$"])
plt.savefig(figfile)