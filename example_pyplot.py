import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

from palettable.colorbrewer.qualitative import Set1_9
from palettable.tableau import Tableau_10, Tableau_20
from cycler import cycler
# Set some commonly used colors
almost_black = '#262626'
light_grey = np.array([float(248) / float(255)] * 3)
brewer_set1 = Set1_9.mpl_colors
# Remove the sixth color (yellow) which is too bright
brewer_set1.pop(5)
# Swap the red and blue to let blue come first
brewer_set1[0], brewer_set1[1] = brewer_set1[1], brewer_set1[0]
# Add a decent black color to this list
brewer_set1.append(almost_black)
tableau_10 = Tableau_10.mpl_colors
# Add a decent black color
tableau_10.append(almost_black)
# Swap orange and red
tableau_10[1], tableau_10[3] = tableau_10[3], tableau_10[1]
# swap orange and purple
# now table_au has similar sequence to brewer_set1
tableau_10[3], tableau_10[4] = tableau_10[4], tableau_10[3]
# This is 20-color Tableau which contains light version of tableau_10
tableau_20 = Tableau_20.mpl_colors
default_color_cycler = cycler('color', tableau_10)
colors = default_color_cycler.by_key()['color']
plt.rcParams["axes.prop_cycle"] = default_color_cycler

plt.rcParams["figure.figsize"] = (30,20)
plt.rcParams["font.size"] = 75
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["axes.labelsize"] = 70
plt.rcParams['axes.linewidth'] = 11

# plot shell decomposition
fig, ax1 = plt.subplots()
ax1.tick_params( axis = 'y', which = 'both', bottom = False, top = False, labelbottom = False)
ax1.tick_params( axis = 'x', direction = 'out', width = 10, length = 40, which = 'both', pad = 25)
ax1.set_xticks(range(1000, 5000, 1000)) 
ax1.set_yticks([])
ax1.set_xlim(400, 4200)
#ax1.set_xlabel('Wave Number (cm$^{-1}$)') 
#ax1.set_ylabel('Solvation Shell Decompositions of rCEC Spectrum')
s = np.loadtxt("prod_ACFV0.dat"); mask = (s[:, 0] >= 600) & (s[:, 0] <= 4000)
ax1.plot(s[mask,0], s[mask,1], linewidth = 10, label = r'$\mathrm{H_3O}$ Core')
s = np.loadtxt("prod_ACFV1.dat"); mask = (s[:, 0] >= 600) & (s[:, 0] <= 4000)
ax1.plot(s[mask,0], s[mask,1] * 0.8, linewidth = 10, label = r'1$\mathrm{^{st}}$ Shell')
s = np.loadtxt("prod_ACFV2.dat"); mask = (s[:, 0] >= 600) & (s[:, 0] <= 4000)
ax1.plot(s[mask,0], s[mask,1] * 0.7, linewidth = 10, label = r'2$\mathrm{^{nd}}$ Shell')
s = np.loadtxt("prod_ACFV3.dat"); mask = (s[:, 0] >= 600) & (s[:, 0] <= 4000)
ax1.plot(s[mask,0], s[mask,1] * 0.7, linewidth = 10, label = r'3$\mathrm{^{rd}}$ Shell')
ax1.legend(fontsize = 70, ncol = 2, framealpha = 0.0, columnspacing = 1.0, handlelength = 1.2)
plt.savefig('spec_rcec_shells.png', bbox_inches = 'tight')
plt.close()

def plot_peaks(color, a = 'center'):
    pos = find_peaks(s[:,1])[0]
    ip = 0
    for p in pos:
        f = s[p, 0]
        if f >= 1100 and f <= 1900:
#            ax1.plot([f] * 50, np.linspace(-0.0001+0.0001, s[p,1], 50),
#                    c = color, linewidth = 8, linestyle = 'dashed')
            try:
                ax1.text(f, s[p, 1]+0.00001, "%d"%f, fontdict={'size':50,'color':color,'weight':'bold'}, horizontalalignment=a[ip])
            except:
                ax1.text(f, s[p, 1]+0.00001, "%d"%f, fontdict={'size':50,'color':color,'weight':'bold'}, horizontalalignment=a)
            ip += 1

# plot Eigen/Zundel decomposition (nperseg = 500)
fig, ax1 = plt.subplots()
ax1.tick_params( axis = 'y', which = 'both', bottom = False, top = False, labelbottom = False)
ax1.tick_params( axis = 'x', direction = 'out', width = 10, length = 40, which = 'both', pad = 25)
ax1.set_xticks(range(1000, 5000, 1000))
ax1.set_yticks([])
ax1.set_ylim(-0.00007, 0.00135)
#ax1.set_ylim(0, 0.00135)
ax1.set_xlim(400, 4200)
#ax1.set_xlabel('Wave Number (cm$^{-1}$)') 
#ax1.set_ylabel('Solvation Shell Decompositions of rCEC Spectrum')
s = np.loadtxt("stft_nacc1.dat"); mask = (s[:, 0] >= 600) & (s[:, 0] <= 4000)
ax1.plot(s[mask,0], s[mask,1], linewidth = 10, label = r'Zundel ($n = 1$)')
plot_peaks(colors[0])
s = np.loadtxt("stft_nacc3.dat"); mask = (s[:, 0] >= 600) & (s[:, 0] <= 4000)
ax1.plot(s[mask,0], s[mask,1], linewidth = 10, label = r'Eigen ($n = 3$)')
#plot_peaks(colors[1], a = 'left')
ax1.text(1434, 0.00047, 1434, fontdict={'size':50,'color':colors[1],'weight':'bold'}, horizontalalignment='center')
ax1.text(1834, 0.00058, 1834, fontdict={'size':50,'color':colors[1],'weight':'bold'}, horizontalalignment='left')
ax1.text(1275, 0.000415, 1275, fontdict={'size':50,'color':colors[1],'weight':'bold'}, horizontalalignment='right')
#ax1.plot([1275] * 50, np.linspace(0, 0.000415, 50), c = colors[1], linewidth = 8, linestyle = 'dashed')
s = np.loadtxt("stft_nacc2.dat"); mask = (s[:, 0] >= 600) & (s[:, 0] <= 4000)
ax1.plot(s[mask,0], s[mask,1], linewidth = 10, label = r'Unassigned ($n = 2$)')
plot_peaks(colors[2], a = 'left')
ax1.text(1216, 0.00098, 1216, fontdict={'size':50,'color':colors[2],'weight':'bold'}, horizontalalignment='right')
#ax1.plot([1216] * 50, np.linspace(0, 0.0009535, 50), c = colors[2], linewidth = 8, linestyle = 'dashed')
s = np.loadtxt("stft_total.dat"); mask = (s[:, 0] >= 600) & (s[:, 0] <= 4000)
ax1.plot(s[mask,0], s[mask,1], linewidth = 10, c = 'k', label = r'Total')
ax1.plot(np.linspace(400, 4300, 100), [0] * 100, c = 'k', linewidth = 10, linestyle = 'dashed')
plot_peaks('k')
ax1.legend(fontsize=70, framealpha = 0.0)
plt.savefig('spec_rcec_stft.png', bbox_inches = 'tight')
plt.close()

# plot Eigen/Zundel decomposition (nperseg = 350)
fig, ax1 = plt.subplots()
ax1.tick_params( axis = 'y', which = 'both', bottom = False, top = False, labelbottom = False)
ax1.tick_params( axis = 'x', direction = 'out', width = 10, length = 40, which = 'both', pad = 25)
ax1.set_xticks(range(1000, 5000, 1000))
ax1.set_yticks([])
ax1.set_ylim(-0.00007, 0.00185)
#ax1.set_ylim(0, 0.00135)
ax1.set_xlim(400, 4200)
#ax1.set_xlabel('Wave Number (cm$^{-1}$)') 
#ax1.set_ylabel('Solvation Shell Decompositions of rCEC Spectrum')
s = np.loadtxt("stft_nacc1.dat"); mask = (s[:, 0] >= 600) & (s[:, 0] <= 4000)
ax1.plot(s[mask,0], s[mask,1], linewidth = 10, label = r'Zundel ($n = 1$)')
plot_peaks(colors[0])
s = np.loadtxt("stft_nacc3.dat"); mask = (s[:, 0] >= 600) & (s[:, 0] <= 4000)
ax1.plot(s[mask,0], s[mask,1], linewidth = 10, label = r'Eigen ($n = 3$)')
plot_peaks(colors[1], a = ['center','left'])
#ax1.text(1434, 0.00047, 1434, fontdict={'size':50,'color':colors[1],'weight':'bold'}, horizontalalignment='center')
#ax1.text(1834, 0.00058, 1834, fontdict={'size':50,'color':colors[1],'weight':'bold'}, horizontalalignment='left')
s = np.loadtxt("stft_nacc2.dat"); mask = (s[:, 0] >= 600) & (s[:, 0] <= 4000)
ax1.plot(s[mask,0], s[mask,1], linewidth = 10, label = r'Unassigned ($n = 2$)')
plot_peaks(colors[2])
s = np.loadtxt("stft_total.dat"); mask = (s[:, 0] >= 600) & (s[:, 0] <= 4000)
ax1.plot(s[mask,0], s[mask,1], linewidth = 10, c = 'k', label = r'Total')
ax1.plot(np.linspace(400, 4300, 100), [0] * 100, c = 'k', linewidth = 10, linestyle = 'dashed')
plot_peaks('k')
ax1.legend(fontsize=70, framealpha = 0.0)
plt.savefig('spec_rcec_stft.png', bbox_inches = 'tight')
plt.close()

# plot ci2 VDOS
fig, ax1 = plt.subplots()
ax1.tick_params( axis = 'y', which = 'both', bottom = False, top = False, labelbottom = False)
ax1.tick_params( axis = 'x', direction = 'out', width = 10, length = 40, which = 'both', pad = 25)
ax1.set_xticks(range(1000, 5000, 1000)) 
ax1.set_yticks([])
ax1.set_xlim(400, 4200)
ax1.set_xlabel('Wavenumber (cm$^{-1}$)') 
ax1.set_ylabel('VDOS')
s = np.loadtxt("prod_ACFC1.dat"); mask = (s[:, 0] >= 600) & (s[:, 0] <= 4000)
ax1.plot(s[mask,0], s[mask,1], linewidth = 10, label = r'$c_1^2$')
s = np.loadtxt("prod_ACFC2.dat"); mask = (s[:, 0] >= 600) & (s[:, 0] <= 4000)
ax1.plot(s[mask,0], s[mask,1] * 0.8, linewidth = 10, label = r'$c_2^2$')
s = np.loadtxt("prod_ACFC3.dat"); mask = (s[:, 0] >= 600) & (s[:, 0] <= 4000)
ax1.plot(s[mask,0], s[mask,1] * 0.5, linewidth = 10, label = r'$c_3^2$')
s = np.loadtxt("prod_ACFC4.dat"); mask = (s[:, 0] >= 600) & (s[:, 0] <= 4000)
ax1.plot(s[mask,0], s[mask,1] * 0.5, linewidth = 10, label = r'$c_4^2$')
#ax1.legend(fontsize = 70, ncol = 2, framealpha = 0.0, columnspacing = 1.0, handlelength = 1.2)
ax1.legend(fontsize = 70, framealpha = 0.0)
plt.savefig('spec_rcec_ci2.png', bbox_inches = 'tight')
plt.close()

# plot Eigen/Zundel decomposition (nperseg = 350) for TOC
plt.rcParams["figure.figsize"] = (31,30)
plt.rcParams["axes.labelsize"] = 100
plt.rcParams["font.size"] = 100
fig, ax1 = plt.subplots()
ax1.set_xlabel('Wave Number (cm$^{-1}$)') 
ax1.set_ylabel('Excess Charge Spectra', labelpad = 25)
ax1.tick_params( axis = 'y', which = 'both', bottom = False, top = False, labelbottom = False, pad = 25)
ax1.tick_params( axis = 'x', direction = 'out', width = 10, length = 40, which = 'both', pad = 25)
ax1.set_xticks(range(1000, 5000, 1000))
ax1.set_yticks([])
ax1.set_ylim(-0.00007, 0.00185)
ax1.set_xlim(400, 4200)
s = np.loadtxt("stft_nacc1.dat"); mask = (s[:, 0] >= 600) & (s[:, 0] <= 4000)
ax1.plot(s[mask,0], s[mask,1], linewidth = 12.5) #, label = r'Zundel')
s = np.loadtxt("stft_nacc3.dat"); mask = (s[:, 0] >= 600) & (s[:, 0] <= 4000)
ax1.plot(s[mask,0], s[mask,1], linewidth = 12.5) #, label = r'Eigen')
s = np.loadtxt("stft_nacc2.dat"); mask = (s[:, 0] >= 600) & (s[:, 0] <= 4000)
s = np.loadtxt("stft_total.dat"); mask = (s[:, 0] >= 600) & (s[:, 0] <= 4000)
ax1.plot(s[mask,0], s[mask,1], linewidth = 10, c = 'k') #, label = r'Total')
ax1.plot(np.linspace(400, 4300, 100), [0] * 100, c = 'k', linewidth = 12.5, linestyle = 'dashed')
ax1.legend(fontsize=100, framealpha = 0.0)
plt.savefig('spec_rcec_stft_toc.png', bbox_inches = 'tight')
plt.close()
