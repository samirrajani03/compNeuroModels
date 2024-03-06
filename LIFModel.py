import matplotlib.pyplot as plt
import numpy as np
import math as m

EL = -0.070 #V
Rm = 5000000 #Ohm
Cm = 2e-9 #F
Vth = -0.050 #V
Vreset = -0.065 #V
tau = Rm * Cm
I0 = (1/Rm) * (Vth-EL)# + 1.746e-21 #b: Minimum Iapp to produce spikes is 4000A or 4e-9pA

#properties of time vector
t0 = 0
dt = 0.0001 # divided by 10 for 2C
tmax = 2

#intialize time vector and memebrane potential empty vector
t = np.arange(t0, tmax, dt)
Vm = np.zeros(t.size)

#set the first value of memebrane potential vector to leaky potential/resting potential
Vm[0] = EL

#intialize Applied Current Vector for 1C, Fire Rate and 1/f vectors for 1C and 1D, sigma_I for noise and the noise vector
Iapp = np.linspace(I0, I0 + 1.746e-9, 10) # +1.746e-9 gets a firing rate of 100hz 
FireRates = np.zeros(Iapp.size)
InverseF = np.zeros(Iapp.size)
sigma_I = 0.001 #value of sigma_I that noise is noticeable - causing spikes
noise = np.random.randn(t.size) * sigma_I * m.sqrt(dt)

#for loop for determining fire rates and 1/f at different applied currents + with noise included from part 2
#if statements to prevent taking the natural log of 0 or negative numbers
for j in range(0, Iapp.size):
    Vss = EL + (Iapp[j] * Rm) 
    if (Vss - Vth <= 0):
        ISI = 0
        f = 0
    else:
        ISI = -tau * m.log((Vss - Vth)/(Vss - Vreset))
        f = 1/ISI
    FireRates[j] = f
    if ((Iapp[j] * Rm) + EL - Vreset <= 0 or (Iapp[j] * Rm) + EL - Vth <= 0):
        InverseF[j] = 0
    else:
        InverseF[j] = (tau * m.log((Iapp[j] * Rm) + EL - Vreset)) - (tau * m.log((Iapp[j] * Rm) + EL - Vth))

spikes = np.zeros(FireRates.size)# initializes spikes vector for calculated rate from simulations
counts = 0 # counter for spikes per simualtion

# outer for loop iterates through Iapp values for the applied current and then adds the calculated rates to the spikes vector
# nested for loop for leaky integrate and fire model, using forward euler method, and checking and counting for spikes
for j in range(0, Iapp.size):
    for i in range(1, t.size):
        dVmdt = (((EL - Vm[i - 1])/Rm) + Iapp[j]) / Cm
        Vm[i] = Vm[i - 1] + (dt * dVmdt)# + noise[i]
        if (Vm[i] > Vth):
            Vm[i] = Vreset
            counts = counts + 1
    rate = counts/tmax
    spikes[j] = rate
    counts = 0


#code for plotting LIF model
plt.plot(t, Vm, label='Vm vs Time')
plt.xlabel('Time(s)')
plt.ylabel('Vm(V)')
plt.legend()

# code for plotting Firing Rate Curve and 1/f curve
# plt.plot(Iapp, FireRates, color='red', label='f(Iapp) - ISI Calc')
# plt.plot(Iapp, spikes, marker='+', color='red', label='f(Iapp) - simulated counts')
# # plt.plot(Iapp, InverseF, color='black', label='1/f')
# plt.xlabel('Applied Current(A)')
# plt.ylabel('Fire Rate(hz)')
# plt.legend()
plt.show()


