import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d import Axes3D
import scipy
from scipy import stats
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import scipy.ndimage as ndimage
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import re

class telescopeParams:
    
    def __init__(self, nBars, d, D, L, color, name):
        """
        The constructor function of the telescopeParams class.
        
        Parameters:
        nBars (int): The number of bars in the panel.
        d (float): The pixel size in cm.
        D (float): The panel separation.
        L (float): The distance from the volcano in m.
        Color (map color): The Colormap
        """
        self.nBars = nBars
        self.d = d      
        self.D = D      
        self.L = L     
        self.color = color
        self.name = name

    def create_plot(self, data, colorbar_label):
        """
        Creates a plot with given data and colorbar label.
        
        Parameters:
        data (ndarray): 2D array data to be plotted.
        colorbar_label (str): Label for the colorbar in the plot.
        """
        
        nBars = self.nBars 
        theta = self.theta 
        Nd = self.Nd  

        fig, ax = plt.subplots(figsize=(10,7))
        extent=(-nBars*theta,nBars*theta,-nBars*theta,nBars*theta)
        im = ax.imshow(data, interpolation='nearest', extent=extent, origin='upper', cmap="jet")
        ax.set_xlabel("$\Theta_x$ [deg]", fontsize = 15)
        ax.set_ylabel("$\Theta_y$ [deg]", fontsize = 15)
        ax.set_title(self.name, fontsize = 15)

        # Color bar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)

        # Color bar
        clb = plt.colorbar(im, cax=cax)
        clb.set_label(colorbar_label, fontsize = 15)
        clb.ax.tick_params(labelsize = 12)


        labelsx = np.round(np.linspace(-nBars*theta, nBars*theta, 11),1)
        labelsy = np.round(np.linspace(-nBars*theta, nBars*theta, 11),1)
        
        output_pdf = self.name + colorbar_label + ".png"
        with PdfPages(output_pdf) as pdf:
            pdf.savefig(fig)  

        plt.show()
        
#         fig = plt.figure(figsize=(10, 7))
#         ax = fig.add_subplot(131)
#         ax.set_xlabel("$\Theta_x$ [deg]", fontsize = 15)
#         ax.set_ylabel("$\Theta_y$ [deg]", fontsize = 15)
#         plt.title(self.name, fontsize = 15)
        
        
#         ax.xaxis.set_major_locator(MaxNLocator(integer=True))
#         ax.yaxis.set_major_locator(MaxNLocator(integer=True))
   
#         extent=(-nBars*theta,nBars*theta,-nBars*theta,nBars*theta)
#         im = plt.imshow(data, interpolation='nearest', extent=extent, origin='lower', cmap=self.color)
    
#         # Color bar
     
#         clb = plt.colorbar()
#         clb.set_label(colorbar_label, fontsize = 15)
#         clb.ax.tick_params(labelsize = 12)
    
#         # Set tick labels
#         labelsx = np.round(np.linspace(-nBars*theta, nBars*theta, 11),1)
#         labelsy = np.round(np.linspace(-nBars*theta, nBars*theta, 11),1)
#         ax.set_xticklabels(labelsx.astype(int))
#         ax.set_yticklabels(labelsy.astype(int))
    
#         # Center location
#         plt.axvline(x=0, color='k', lw=1, linestyle='--')
#         plt.axhline(y=0, color='k', lw=1, linestyle='--')
    
#         # Increase tick labels
#         ax.tick_params(axis='both', which='major', labelsize=15)
        
#         plt.tight_layout()
        
#         plt.savefig(self.name + colorbar_label + '.png', bbox_inches='tight')
        
        
    def solid_angle(self):
        """
        Calculates the solid angle of the telescope.
        """
        self.A = self.d**2    #  Pixel area
        d_Omega =  self.A/(self.D**2)
        
        self.theta = int(np.rad2deg(np.arctan(self.d*self.nBars/self.D)))/float(self.nBars)
        
        self.Nd = (2*self.nBars-1)    # Number of trajectories
        C = self.nBars-1         # Shiffting index
        
        r = np.zeros((self.Nd,self.Nd))
        
        matrix_P1 = np.zeros((self.nBars, self.nBars))
        matrix_P2 = np.ones((self.nBars, self.nBars))
        
        for i in range(self.nBars):
            for j in range(self.nBars):
                for k in range(self.nBars):
                    for l in range(self.nBars):
        
                        h = k-i
                        b = l-j
                        ik = np.abs(k-i)
                        jl = np.abs(l-j)
        
                        E = np.sqrt((ik*self.d)**2 + (jl*self.d)**2)
                        r[h+C,b+C]= np.sqrt(self.D**2 + E**2)
        
        self.solid_Ang = self.A/((r/2)**2)
        self.r = r
        delta_theta = 2*np.tan(self.d/self.D)
        
     
        
        delta_x = self.L*np.tan(delta_theta)
        self.delta_theta = delta_theta
        self.delta_x = delta_x

        self.create_plot(self.solid_Ang*1000, 'Solid angle [x10$^{-3}$sr]')
        
        return self.solid_Ang
        


        
    def N_pixel(self):
        """
        Calculates the number of pixels.
        """
        
        Nd = (2*self.nBars-1)    # Number of trajectories
        C = self.nBars-1        # Shiffting index
        
        self.n_Pixel = np.zeros((Nd,Nd))
        
        matrix_P1 = np.zeros((self.nBars, self.nBars))
        matrix_P2 = np.ones((self.nBars, self.nBars))
        
        for i in range(self.nBars):
            for j in range(self.nBars):
                matrix_P1[i,j] = 1
                for k in range(self.nBars):
                    for l in range(self.nBars):
        
                        if matrix_P1[i,j]== 1:
                            iP1 = i
                            jP1 = j
                        if matrix_P2[k,l]== 1:
                            iP2 = k
                            jP2 = l
        
                        h = iP1 - iP2
                        b = jP1 - jP2
        
                        self.n_Pixel[h+C,b+C]= self.n_Pixel[h+C,b+C] + 1
                        


    def acceptance(self):
        """
        Calculates the acceptance of the telescope.
        """
        
        self.acceptance = self.n_Pixel * self.solid_Ang * self.A / 4

        self.create_plot(self.acceptance,'Acceptance [cm$^2$sr]')
        
        return self.acceptance
            
            
    def S_pixels(self):
        """
        Calculates and plots the area of pixels.
        """
        
        #self.create_plot(self.n_Pixel*self.A/1000, '$S$ [x10$^{3}$cm$^2$]')
    
    def params(self):
        
        print ("Maximum acceptance = " + str(np.max(self.acceptance)))
        print ("Maximum number of pixel = " + str(np.max(self.n_Pixel)))
        print ("Maximum number of r = " + str(np.max(self.r)))  
        print ("Angular aperture = " + str(2*self.nBars*self.theta) + " deg")
        print ("Angular aperture = " + str(np.round(2*np.arctan(self.nBars*self.d/self.D),2)) + " rad")
        print ("Angular resolution = " + str(np.round(1000*self.delta_theta,2)) + " mrad")
        print ("Spatial resolution on the volcano = %s m" % str(np.round(self.delta_x,2)))
        print ("Maximum number of solid_angle = " + str(np.max(self.solid_Ang*1000)))
        print ("Maximum area = " + str(np.max(self.A)))
        
            
    def create_subplot(self, data, colorbar_label, ax):
        """
        Creates a subplot with given data, colorbar label and axis.
        
        Parameters:
        data (ndarray): 2D array data to be plotted in the subplot.
        colorbar_label (str): Label for the colorbar in the subplot.
        ax (matplotlib.axes.Axes): Axis object to draw the plot onto.
        """
        nBars = self.nBars
        theta = self.theta
        Nd = self.Nd
        
        ax.set_xlabel("$\Theta_x$ [deg]", fontsize = 15)
        ax.set_ylabel("$\Theta_y$ [deg]", fontsize = 15)
    
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    
        extent=(-nBars*theta,nBars*theta,-nBars*theta,nBars*theta)
        im = ax.imshow(data, interpolation='nearest', extent=extent, origin='lower', cmap=self.color)
    
        # Create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
    
        # Color bar
        clb = plt.colorbar(im, cax=cax)
        clb.set_label(colorbar_label, fontsize = 15)
        clb.ax.tick_params(labelsize = 15)
    
        # Set tick labels
        labelsx = np.round(np.linspace(-nBars*theta, nBars*theta, 11),1)
        labelsy = np.round(np.linspace(-nBars*theta, nBars*theta, 11),1)
        ax.set_xticklabels(labelsx.astype(int))
        ax.set_yticklabels(labelsy.astype(int))
    
        # Center location
        ax.axvline(x=0, color='k', lw=1, linestyle='--')
        ax.axhline(y=0, color='k', lw=1, linestyle='--')
    
        # Increase tick labels
        ax.tick_params(axis='both', which='major', labelsize=20)
    
            
    
    def plot_all_params(self):
        """
        Plots all parameters by creating a grid of subplots.
        """
        
        # Create the figure with a 2x2 grid of subplots
        fig, axs = plt.subplots(2, 2, figsize=(16,16))

        # Plot each set of data in a different subplot

        self.create_subplot(self.solid_Ang*1000, 'Solid angle [x10$^{-3}$sr]',axs[0, 0])        
        self.create_subplot(self.r, '$r_{m.n}$[cm]',axs[0, 1])
        self.create_subplot(self.acceptance,'Acceptance [cm$^2$sr]',axs[1, 0])
        self.create_subplot(self.n_Pixel*self.A/1000, '$S$ [x10$^{3}$cm$^2$]',axs[1, 1])
        
        plt.tight_layout()

        output_pdf = self.name + "Params.pdf"
        with PdfPages(output_pdf) as pdf:
            pdf.savefig(fig)
            
        plt.show()


