import numpy as np
import matplotlib.pyplot as plt

exp = np.genfromtxt('srivastava.csv', delimiter=',', skip_header=2, names=['exp_time','exp_massflow'])

cfd1 = np.genfromtxt('massflow1.csv', delimiter=',', skip_header=1, names=["massflow","Normals:0","Normals:1","Normals:2","Theta","U1:0","U1:1","U1:2","U1Mean:0","U1Mean:1","U1Mean:2","U2:0","U2:1","U2:2","U2Mean:0","U2Mean:1","U2Mean:2","alpha0","alpha1","alpha1Mean","epsilon","k","magU1Mean","p","vtkValidPointMask","Time","Point Coordinates:0","Point Coordinates:1","Point Coordinates:2"])

cfd2 = np.genfromtxt('massflow2.csv', delimiter=',', skip_header=1, names=["massflow","Normals:0","Normals:1","Normals:2","Theta","U1:0","U1:1","U1:2","U1Mean:0","U1Mean:1","U1Mean:2","U2:0","U2:1","U2:2","U2Mean:0","U2Mean:1","U2Mean:2","alpha0","alpha1","alpha1Mean","epsilon","k","magU1Mean","p","vtkValidPointMask","Time","Point Coordinates:0","Point Coordinates:1","Point Coordinates:2"])

exp_time = exp['exp_time']
exp_massflow = exp['exp_massflow']

cfd_time1 = cfd1["Time"]
cfd_massflow1 = cfd1["massflow"]

cfd_time2 = cfd2["Time"]
cfd_massflow2 = cfd2["massflow"]

plt.figure(1)

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

ax1 = plt.gca()

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

line1=ax1.plot(exp_time, exp_massflow, color='r', marker='None',markersize=7.5,linestyle='-',label=r'$Srivastava \ et \ al. \ 2003$')
line2=ax1.plot(cfd_time1, cfd_massflow1, color='r', marker='None',markersize=7.5,linestyle='--',label=r'$\Theta \ Transport \ Equation, \ Slip \ Walls$')
line3=ax1.plot(cfd_time2, cfd_massflow2, color='b', marker='None',markersize=7.5,linestyle='--',label=r'$\Theta \ Algebraic \ Equation, \ Slip \ Walls$')

lines = line1+line2+line3
labels = [l.get_label() for l in lines]
legend= ax1.legend(lines, labels, loc='upper right', fontsize='medium')

ax1.set_xlabel(r'\textit{$Time, \  t \  (s)$}')
ax1.set_ylabel(r'\textit{$Mass \ flow \ rate, \ \dot{m} \ (g/s)$}')
ax1.set_xlim([0,7])
ax1.set_ylim([0,200])
plt.savefig('massflow',bbox_inches='tight')
