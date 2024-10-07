
import os
import requests
import numpy as np
import math
import TopographyDownloader
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import signal
import matplotlib as mpl
from matplotlib.ticker import FuncFormatter
import timeit

# This class will make end-to-end simulation of the muography 
class FluxIntegrated(object):

    
    # Input variables are expected to be in the following units:
    # Region        = [latitude1, latitude2, longitude1, longitude2, name]
    # srtm1_data    = Module 
    # Points        = Number of resolution
    
    def __init__(self,matrix):
        """Initializing global variables 
        """
    
        self.matrix = matrix
        
    def showData(self):
        self.azimuth = self.matrix[:,1]
        self.cenith = self.matrix[:,0] * -1
        self.Np = int(np.sqrt(len(self.azimuth)))
        self.z = np.abs(self.matrix[:,2].reshape(self.Np,self.Np))/1000
        """
        print("Azimuth")
        print(self.azimuth)
        print("#########")
        print("Cenith")
        print(self.cenith)
        print("#########")
        print(self.z)
        print("#########")
        print(self.z.max())
        """
    def int_flux(self,Lm, theta_in):
        # Minimum muon energy estimation

        E = np.linspace(1,1e5,100000) # GeV
        M = 10000
        y = np.log10(E)
        l0 = 0.2549
        l1 = 0.0801
        l2 = 0.0368
        l3 = -0.0461
        l4 = 0.0154
        dE = (1e4-1e0)/M
        dEdp = -10**((l4*y**3.65) + (l3*y**3)  + (l2*y**2) + l1*y + l0)

        Op = -E*1e3/dEdp

        L = Lm*1e2 # Lenght to cm
        rho = 2.65 # Standard rock density g/cm3
        p = L*rho # Opacity g/cm2

        Eu = 105.6 # muon mass [eV/c2]

        Eminp = E[np.argmin((Op - p)**2)] - Eu/1e3 

        # Muon flux model

        cenith = theta_in
        theta = cenith*np.pi/180.0 # cenith angle / redians
        E0 = np.linspace(1e0,1e4,10000) # Muon energy / GeV

        c = 1 # Speed of light
        m0 = 0.1056 # Muon mass in GeV/c^2
        p = np.sqrt((E0**2-m0**2*c**4)/c**2) # Bugaev model

        y = np.log10(p*np.cos(theta)) # Bugaev/Reyna model
        AB = 0.00253
        a0 = 0.2455
        a1 = 1.288
        a2 = -0.2555
        a3 = 0.0209

        Phi_Bugaev_Reyna = AB*(p**(-(a3*y**3 + a2*y**2 + a1*y + a0)))*(np.cos(theta))**3

        # Integrated flux estimation

        N = len(Phi_Bugaev_Reyna)
        Int_flux = 0
        Open_Sky = 0

        for i in range(N):

            Open_Sky = Open_Sky + Phi_Bugaev_Reyna[i]  # Open sky flux
            if E0[i] >= Eminp:
                Int_flux = Int_flux + Phi_Bugaev_Reyna[i] # Traversing flux

        Open_Sky = Open_Sky*dE
        Int_flux = Int_flux*dE
                    
        return  Int_flux*86400, Open_Sky*86400
    

    def Elevation(self,Longitude1,Longitude2,Latitude1,Latitude2,Altitude1,Altitude2):
        
        cenith_matrix = self.cenith.reshape(self.Np,self.Np)
        self.cenith_matrix = cenith_matrix
        
        Ce = 40075 # Earth circunference km
        X = Ce*(Longitude1-Longitude2)*np.cos((Latitude1-Latitude2)*np.pi/180.0)/360.0 # Lenght Longitud
        Y = Ce*(Latitude1-Latitude2)/360.0 # Lenght Latitude
        H = np.sqrt(X**2+Y**2)*1000 # Horizontal distance in meters
        V = Altitude2 - Altitude1
        Elevation = np.arctan(V/H)
        print ("Elevation : %f" % (Elevation*180/np.pi))
        
        start = timeit.default_timer()

        N = self.Np
        Trav_Flux = np.zeros((N,N))
        Open_Sky_Flux = np.zeros((N,N))
        cenithal = np.zeros((N,N))

        for i in range(N):
            for j in range(N):
                cenithal[i,j] = (np.pi/2 - Elevation + self.cenith_matrix[i,j])*180/np.pi  # Cenith angle in grades

                A, B = self.int_flux(self.z[i,j]*1000, cenithal[i,j])

                Trav_Flux[i,j] = A
                Open_Sky_Flux[i,j] = B 

        stop = timeit.default_timer()
        print('Time: %f s'% (stop - start))
        
        self.Trav_Flux = Trav_Flux
        self.cenithal = cenithal
        self.Open_Sky = Open_Sky_Flux
        
        Opacity = np.zeros((N,N))
        Density = np.zeros((N,N))
        R = np.zeros((N,N))

        rho = 2.65 # Standard rock density g/cm3

        for i in range(N):
            for j in range(N):

                if self.z[i,j] > 0.1:

                    R[i,j] = Trav_Flux[i,j]/Open_Sky_Flux[i,j]
                    mu = np.log(1/R[i,j])/(self.z[i,j]*1e5)
                    kappa = mu/rho

                    Opacity[i,j] = (-1/kappa)*np.log(R[i,j])
                    Density[i,j] = Opacity[i,j]/(self.z[i,j]*1e5)
        
        self.Opacity = Opacity
        self.Density = Density
        
    def OpenSky(self):
        

        fig, ax = plt.subplots(figsize=(10,7))
        extent = (min(self.azimuth)*180/np.pi, max(self.azimuth)*180/np.pi, np.max(self.cenithal), np.min(self.cenithal))
        im = ax.imshow(self.Open_Sky, interpolation='nearest', extent=extent, origin='upper',norm=mpl.colors.LogNorm(), cmap="jet")
        ax.set_xlabel("Azimuth [degree]", fontsize = 15)
        ax.set_ylabel("Zenith [degree]", fontsize = 15)
        ax.set_title("Open Sky Flux", fontsize = 15)

        # Color bar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)

        # Color bar
        clb = plt.colorbar(im, cax=cax)
        clb.set_label('Open Sky Flux [cm$^{-2}$sr$^{-1}$day$^{-1}$]', fontsize = 15)
        clb.ax.tick_params(labelsize = 12)

        # Define contour levels
        contour_levels = np.logspace(np.log10(self.Open_Sky.min()), np.log10(self.Open_Sky.max()/2), 5)  # 4 contour lines in log scale

        # Draw contour lines
        contour_lines = ax.contour(self.Open_Sky, contour_levels, colors='k', origin='upper', extent=extent)

        # Formatter function for contour labels
        formatter = FuncFormatter(lambda x, pos: "{:.1e}".format(x))

        # Add labels to contour lines
        ax.clabel(contour_lines, inline=True, fontsize=10, colors='black', fmt=formatter)

        labelsx = np.round(np.linspace(min(self.azimuth)*180/np.pi, max(self.azimuth)*180/np.pi, 11),0)
        labelsy = np.round(np.linspace(np.max(self.cenithal), np.min(self.cenithal),  11),0)

        
        output_pdf = "Open_Sky_Flux.png"
        with PdfPages(output_pdf) as pdf:
            pdf.savefig(fig)  

        plt.show()


        
    def ShowIntegratedFlux(self):
        
        fig, ax = plt.subplots(figsize=(10,7))
        extent = (min(self.azimuth)*180/np.pi, max(self.azimuth)*180/np.pi, np.max(self.cenithal), np.min(self.cenithal))
        im = ax.imshow(self.Trav_Flux, interpolation='nearest', extent=extent, origin='upper',norm=mpl.colors.LogNorm(), cmap="jet")
        ax.set_xlabel("Azimuth [degree]", fontsize = 15)
        ax.set_ylabel("Zenith [degree]", fontsize = 15)
        ax.set_title("Integrated Flux", fontsize = 15)

        # Color bar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)

        # Color bar
        clb = plt.colorbar(im, cax=cax)
        clb.set_label('Integrated Flux [cm$^{-2}$sr$^{-1}$day$^{-1}$]', fontsize = 15)
        clb.ax.tick_params(labelsize = 12)

        # Define contour levels
        contour_levels = np.logspace(np.log10(self.Trav_Flux.min()), np.log10(self.Trav_Flux.max()/2), 5)  # 4 contour lines in log scale

        # Draw contour lines
        contour_lines = ax.contour(self.Trav_Flux, contour_levels, colors='k', origin='upper', extent=extent)

        # Formatter function for contour labels
        formatter = FuncFormatter(lambda x, pos: "{:.1e}".format(x))

        # Add labels to contour lines
        ax.clabel(contour_lines, inline=True, fontsize=10, colors='black', fmt=formatter)

        labelsx = np.round(np.linspace(min(self.azimuth)*180/np.pi, max(self.azimuth)*180/np.pi, 11),0)
        labelsy = np.round(np.linspace(np.max(self.cenithal), np.min(self.cenithal),  11),0)
        
        output_pdf = "Integrated_Flux.png"
        with PdfPages(output_pdf) as pdf:
            pdf.savefig(fig)  

        plt.show()

        return self.Trav_Flux, self.Open_Sky

    
    def Transmission(self):
        
        self.transmission = self.Trav_Flux/self.Open_Sky
        
        fig, ax = plt.subplots(figsize=(10,7))
        extent = (min(self.azimuth)*180/np.pi, max(self.azimuth)*180/np.pi, np.max(self.cenithal), np.min(self.cenithal))
        im = ax.imshow(self.transmission, interpolation='nearest', extent=extent, origin='upper',norm=mpl.colors.LogNorm(), cmap="jet")
        ax.set_xlabel("Azimuth [degree]", fontsize = 15)
        ax.set_ylabel("Zenith [degree]", fontsize = 15)
        ax.set_title("Transmission", fontsize = 15)

        # Color bar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)

        # Color bar
        clb = plt.colorbar(im, cax=cax)
        clb.set_label('Transmission', fontsize = 15)
        clb.ax.tick_params(labelsize = 12)

        labelsx = np.round(np.linspace(min(self.azimuth)*180/np.pi, max(self.azimuth)*180/np.pi, 11),0)
        labelsy = np.round(np.linspace(np.max(self.cenithal), np.min(self.cenithal),  11),0)

        output_pdf = "Transmission.png"
        
        with PdfPages(output_pdf) as pdf:
            pdf.savefig(fig)  

        plt.show()
        
        return self.transmission
        
    def ShowDensity(self):
        fig, ax = plt.subplots(figsize=(20,17))
        extent = (min(self.azimuth)*180/np.pi, max(self.azimuth)*180/np.pi, np.max(self.cenithal), np.min(self.cenithal))
        im = ax.imshow(self.Density, interpolation='nearest', extent=extent, origin='upper', norm=mpl.colors.LogNorm(), cmap="jet")
        ax.set_xlabel("Azimuth [degree]", fontsize = 25)
        ax.set_ylabel("Zenith [degree]", fontsize = 25)
        ax.set_title("Density", fontsize = 25)

        # Color bar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        # Color bar
        clb = plt.colorbar(im, cax=cax)
        clb.set_label(r'Density [g cm$^{-3}$]', fontsize = 25)
        clb.ax.tick_params(labelsize = 20)

        # Define contour levels
        contour_levels = np.logspace(self.Density.min(), self.Density.max(), 5)  # 4 contour lines in log scale

        # Draw contour lines
        contour_lines = ax.contour(self.Density, contour_levels, colors='k', origin='upper', extent=extent)

        # Formatter function for contour labels
       
        # Add labels to contour lines
        ax.clabel(contour_lines, inline=True, fontsize=10, colors='black')

        labelsx = np.round(np.linspace(min(self.azimuth)*180/np.pi, max(self.azimuth)*180/np.pi, 11),0)
        labelsy = np.round(np.linspace(np.max(self.cenithal), np.min(self.cenithal),  11),0)

        # ax.set_xticks(labelsx, fontsize = 15)
        # ax.set_yticks(labelsy, fontsize = 15)

        output_pdf = "Density.png"
        with PdfPages(output_pdf) as pdf:
            pdf.savefig(fig)  

        plt.show()

        
        
    def ShowOpacity(self):
        fig, ax = plt.subplots(figsize=(20,17))
        extent = (min(self.azimuth)*180/np.pi, max(self.azimuth)*180/np.pi, np.max(self.cenithal), np.min(self.cenithal))
        im = ax.imshow(self.Opacity/1e3, interpolation='nearest', extent=extent, origin='upper', norm=mpl.colors.LogNorm(), cmap="jet")
        ax.set_xlabel("Azimuth [degree]", fontsize = 25)
        ax.set_ylabel("Zenith [degree]", fontsize = 25)
        ax.set_title("Opacity", fontsize = 25)

        # Color bar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        # Color bar
        clb = plt.colorbar(im, cax=cax)
        clb.set_label(r'Opacity $\times 10^3$ [g cm$^{-2}$]', fontsize = 25)
        clb.ax.tick_params(labelsize = 20)

        # Define contour levels
        min_val = self.Opacity.min() if self.Opacity.min() > 0 else 1e-10  # reemplaza 1e-10 con un número adecuado para tus datos
        contour_levels = np.logspace(np.log10(min_val), np.log10(self.Opacity.max()/2), 5)


        
        # Draw contour lines
        contour_lines = ax.contour(self.Opacity, contour_levels, colors='k', origin='upper', extent=extent)

        # Formatter function for contour labels
        

        # Add labels to contour lines
        ax.clabel(contour_lines, inline=True, fontsize=10, colors='black')

        labelsx = np.round(np.linspace(min(self.azimuth)*180/np.pi, max(self.azimuth)*180/np.pi, 11),0)
        labelsy = np.round(np.linspace(np.max(self.cenithal), np.min(self.cenithal),  11),0)

        # ax.set_xticks(labelsx, fontsize = 15)
        # ax.set_yticks(labelsy, fontsize = 15)

        output_pdf = "Opacity.png"
        with PdfPages(output_pdf) as pdf:
            pdf.savefig(fig)  

        plt.show()
