from pylab import *
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import wave
import struct
from scipy.signal import spectrogram
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

para = [
    ('c2', [1.23, 1/16000, 100, 160.9, 0.58, 0.25, 7.5e-5]), 
    ('c4', [0.63, 1/32000, 51, 329.6, 1.25, 1.1, 2.7e-4])]
#Parameters of different music notes stated by literature
#Function for checking the stability conditions

"""----------------------------------------------------------------------------------------------"""


def stable(X,T=1/32000,b=2.7e-4,c=329.6,k=1.25):
    sqrt_component = sqrt(16 * b**2 + 4 * (c**2 * X**2 + 4 * k**2))
    denominator = 2 * (c**2 * X**2 + 4 * k**2)
    return X**2 * ((-4 * b + sqrt_component) / denominator) - T

"""According to given conditions on stability, we get the equation for
calculating the proper spatial interval given fixed time intervals T."""

def bisection(f,a,b,tol=1e-12):
    if a > b:
        c=b
        b=a
        a=c
    if a==b:
        raise ValueError
    if f(a)*f(b) >=0.:
        raise ValueError
    while fabs(b-a)/fabs(a) > tol: 
        #print(f"a={repr(a)},   b={repr(b)},   b-a={repr(b-a)}")
        c = (a+b)/2
        sfc = sign(f(c))
        #print(f"c={repr(c)}")
        if sfc == sign(f(a)):
            a = c
        elif sfc == sign(f(b)):
            b = c
        elif sfc == 0:
            print("Here is an exact zero")
            return c
        else:
           raise ValueError("Something went horribly bad: sign(f(midpoint))={}\n".format(sfc))
        
    return c

#X1_min = bisection(stable, 0.01,0.02) #Minimum spatial step for C2

"""
By introducing bisection method,
we find X1_min approximates to 0.0117;
For simplicity, we therefore take
X = 0.0123 m for Note C2,
so that the number of grids in space are integers.
"""
"""----------------------------------------------------------------------------------------------"""

#Initial Condition Functions
def f1(x):
    if 0 < x <= 0.7:
        return (0.002/0.7)*x
    elif 0.7 < x <= 1.23:
        return (-0.002/(1.23-0.7))*(x-0.7)+0.002
    else:
        return 0
#f1 is used for simulating note C2
def f2(x):
    if 0 < x <= 0.44:
        return (0.002/0.44)*x
    elif 0.44 < x <= 0.63:
        return (-0.002/(0.63-0.44))*(x-0.44)+0.002
    else:
        return 0
#f1 is prepared for simulating note C2, while f2 is prepared for note C4.
"""----------------------------------------------------------------------------------------------"""

#Propagation Simulation for notes C2 and C4

def wave_c4(sum_of_time):
    T = 1/32000 #timestep
    fs = 1 / T #Frequency
    L = 0.63 #String length
    X = L / 51 #Space step
    #Relevant coefficients given by literature
    k = 1.25
    c = 329.6
    b1 = 1.1
    b2 = 2.7e-4
    lamda = (c * T) / X
    mu = (k * T) / X**2
    #Meshgridding
    N = int(51)
    M = int(sum_of_time / T)
    #Finite Differences
    a10 = (2 - 2*(lamda**2) - 6*(mu**2) - (4*b2*mu)/k) / (1 + b1*T)
    a11 = ((lamda**2) + 4*(mu**2) + (2*b2*mu)/k) / (1 + b1*T)
    a12 = (-(mu**2)) /(1 + b1*T)
    a20 = (-1 + ((4*b2*mu)/k) + b1*T) / (1 + b1*T)
    a21 = (-2*b2*mu/k) / (1 + b1*T)
    #Initializing a matrix
    U = zeros([N + 2, M + 1])
    #Setting up boundary conditions
    for i in arange(0, M+1):
        U[0,i] = 0
        U[N,i] = 0
    #Setting up initial conditions
    for j in arange(0, N+1):
        U[j,0] = f1(j*X)
        U[j,1] = (a10+a20)*f1(j*X) + (a11+a21)*(f1((j+1)*X) + f1((j-1)*X)) + a12*(f1((j+2)*X)+f1((j-2)*X))
    #Setting up finite-difference recurrence
    for j in arange(1, M-1):
        for i in arange(1, N):
            U[i,j+1] = a10*U[i,j] + a11*(U[i+1,j]+U[i-1,j]) + a12*(U[i+2,j]+U[i-2,j]) + a20*U[i,j-1] + a21*(U[i+1,j-1]+U[i-1,j-1])
    U = U.T
    U = U[:,:52]
    return U
"""The function returns a (32000*x) by 52 matrix U representing the displacement of the string
given the summation of time of x in seconds"""

def wave_c2(sum_of_time):
    T = 1/16000 #timestep
    fs = 1 / T
    L = 1.23 #string length
    X = L / 100 #space step
    #relevant coefficients
    k = 0.58
    c = 160.9
    b1 = 0.25
    b2 = 7.5e-5
    lamda = (c * T) / X
    mu = (k * T) / X**2
    #Meshgridding
    N = int(100)
    M = int(sum_of_time / T)
    #Finite Differences
    a10 = (2 - 2*(lamda**2) - 6*(mu**2) - (4*b2*mu)/k) / (1 + b1*T)
    a11 = ((lamda**2) + 4*(mu**2) + (2*b2*mu)/k) / (1 + b1*T)
    a12 = (-(mu**2)) /(1 + b1*T)
    a20 = (-1 + ((4*b2*mu)/k) + b1*T) / (1 + b1*T)
    a21 = (-2*b2*mu/k) / (1 + b1*T)
    #Initializing a matrix
    U = zeros([N+2, M + 1])
    #Setting up boundary conditions
    for i in arange(0, M+1):
        U[0,i] = 0
        U[N,i] = 0
    #Setting up initial conditions
    for j in arange(0, N+1):
        U[j,0] = f1(j*X)
        U[j,1] = (a10+a20)*f1(j*X) + (a11+a21)*(f1((j+1)*X) + f1((j-1)*X)) + a12*(f1((j+2)*X)+f1((j-2)*X))
    #Setting up finite-difference recurrence
    for j in arange(1, M-1):
        for i in arange(1, N):
            U[i,j+1] = a10*U[i,j] + a11*(U[i+1,j]+U[i-1,j]) + a12*(U[i+2,j]+U[i-2,j]) + a20*U[i,j-1] + a21*(U[i+1,j-1]+U[i-1,j-1])
    U = U.T
    U = U[:,:101]
    return U
"""The function returns a (16000*x) by 101 matrix U representing the displacement of the string
given the summation of time of x in seconds"""

"""----------------------------------------------------------------------------------------------------------"""

"""This 'wave_gen' function is used for generating a wav.file for reviewing the
simulation of propagation. U is the matrix for displacement returned by
'wave_prop' function. loc represents the sample location along the string.
For example, loc=59 represents the 59th grid point on string of length L.
fs represents the sample frequency.
"""
def wave_gen(U,loc=59,fs=16000):
    data = U[:,loc]
    # Normalization to 16-bit PCM range
    audio_data = np.int16(data / np.max(np.abs(data)) * 32767)
    
    #Parameters
    nchannels = 1
    sampwidth = 2
    nframes = len(audio_data)
    
    #Initialize WAV file
    current_time = datetime.datetime.now()
    formatted_time = current_time.strftime("%Y%m%d_%H%M%S")
    wav_file = wave.open('wavefile_{formatted_time}.wav','w')
    wav_file.setparams((nchannels, sampwidth, fs, nframes, 'NONE', 'not compressed'))

    #Load the data
    wav_file.writeframes(audio_data.tobytes())
    wav_file.close()

    print("WAV file created successfully.")
"""By default, sample frequency fs=16000 and sample location on the string is 59,(representing
the 59th column in matrix U of simulation on C2.
To generate that on C4, take loc=35 fs=32000 as inputs"""

"""----------------------------------------------------------------------------------------------"""
#Visualization Section

#for a real-time updating movement of note C2

figure(1) 
h = 0.0123
x = arange(101)*h
U = wave_c2(1)
framecount = 0
for u in U[0:,:]:
    framecount += 1
    clf()
    plot(x,u,'k.')
    suptitle("Wave Simulation for Note C2",fontsize=18)
    title("time={:4.3f}".format(framecount*1/16000))
    xlabel("String")
    ylim(-0.003,0.003)
    draw()
    pause(0.0001)


#for a real-time updating movement of note C4
"""
figure(2)
h = 0.01235
x = arange(52)*h
U = wave_c4(1)
framecount = 0
for u in U:
    framecount += 1
    clf()
    plot(x,u,'k.')
    suptitle("Wave Simulation for Note C4",fontsize=18)
    title("time={:4.3f}".format(framecount*1/32000))
    ylim(-0.003,0.003)
    draw()
    pause(0.001)
"""

"""---------------------------------------------------------"""
#Long-term Spectrograms for C2 simulation
"""
U = wave_c2(1)

signal = U[:8000, 59]
frequencies, times, Sxx = spectrogram(signal, 16000)
plt.figure(figsize=(10, 4))
plt.pcolormesh(times, frequencies, 10 * np.log10(Sxx), shading='gouraud')
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
plt.title('Spectrogram')
plt.colorbar(label='Intensity [dB]')
plt.ylim(0, 5000)
show()
"""
#Short-term Spectrograms for C2 simulation
"""
U = wave_c2(1)
signal2 = U[:2000, 59]
frequencies, times, Sxx = spectrogram(signal2, 16000)
plt.figure(figsize=(10, 4))
plt.pcolormesh(times, frequencies, 10 * np.log10(Sxx), shading='gouraud', vmin=-180)
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
plt.title('Spectrogram')
plt.colorbar(label='Intensity [dB]')
plt.ylim(0, 300)
show()
"""

"""---------------------------------------------------------"""

#Cumulative Spectral Decay for Note C2
"""
U = wave_c2(6)
signal = U[:, 59]
frequencies, times, Sxx = spectrogram(signal, 16000)
data = Sxx[:50,:]
Sxx_dB = 10*np.log10(Sxx)
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

time_indices = range(0, len(times), 10)  
freq_indices = range(0, len(frequencies), 10)  

X, Y = np.meshgrid(times[time_indices], frequencies[freq_indices])
Z = Sxx_dB[freq_indices, :][:, time_indices]

for i, time in enumerate(times[time_indices]):
    ys = frequencies[freq_indices]
    zs = Z[:, i]
    ax.plot(xs=np.full_like(ys, time), ys=ys, zs=zs,color='black')

ax.set_xlabel('Time [s]')
ax.set_ylabel('Frequency [Hz]')
ax.set_zlabel('Amplitude [dB]')
ax.set_title('Cumulative Spectral Decay for Note C2')
ax.view_init(elev=40, azim=-60) 

plt.show()
"""

