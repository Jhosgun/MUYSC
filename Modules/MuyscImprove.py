#!/usr/bin/env python
# coding: utf-8
# %%

# %%

# -*- encoding: latin1 -*-

"""-------------------------------------
MUOGRAPHY SIMULATION CODE,               
  \  | |  |\ \  / __|  __|
 |\/ | |  | \  /\__ \ (   
_|  _|\__/   _| ____/\___|
The structre of MUTE    
  Jorge Jaimes, Jesus Peña  2021             
-------------------------------------
Updated: Frebrary 28, 2021
"""

"""
Full documentation at: https://github.com/Jhosgun/MUYSC/wiki
"""

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
import FluxModels


from srtm import Srtm1HeightMapCollection

srtm1_data = Srtm1HeightMapCollection()


# This class will make end-to-end simulation of the muography 
class Mute(object):

    
    # Input variables are expected to be in the following units:
    # Region        = [latitude1, latitude2, longitude1, longitude2, name]
    # srtm1_data    = Module 
    # Points        = Number of resolution
    
    def __init__(self,region,points,srtm1_data,cmap):
        """Initializing global variables 
        """
        
        self.lat1    = region[0] # First coordinates of latitude
        self.lat2    = region[1] # First coordinates of longitude
        self.lon1    = region[2] # Second coordinates of latitude
        self.lon2    = region[3] # Second coordinates of longitude
        self.name    = region[4] # Name of the study object
        
        
        self.srtm1_data = srtm1_data
        self.cmap = cmap
        self.points = points 
        self.regionPoints = region

    def inicialize_Topography(self,path):
        downloader = TopographyDownloader.TopographyData(self.regionPoints[:4])
        downloader.download_to_path(path)
        print (path)
        
    def measure(self, lat1, lon1, lat2, lon2):  
        """Converting GMS coordinates to local meters.
        """
        """
        Calculate the Haversine distance between two points on the earth's surface.

        Parameters:
        lat1, lon1, lat2, lon2 (float): Latitude and Longitude of the two points in degrees.

        Returns:
        float: Distance between the two points in meters.
        """
        EARTH_RADIUS = 6378.137  # Earth's radius in kilometers

        lat1_rad, lat2_rad = np.radians(lat1), np.radians(lat2)
        delta_lat = lat2_rad - lat1_rad
        delta_lon = np.radians(lon2) - np.radians(lon1)

        a = np.sin(delta_lat / 2) ** 2 + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(delta_lon / 2) ** 2
        c = 2 * math.atan2(np.sqrt(a), np.sqrt(1 - a))
        distance = EARTH_RADIUS * c                                                                                                                       
        return distance * 1000  # Convert distance to meters


    
    def elevation(self):
        """Calculating the elevation of the geological structure.
        """
        
        lat_points = np.linspace(self.lat1, self.lat2, num=self.points)
        lon_points = np.linspace(self.lon1, self.lon2, num=self.points)

        meters_lat = self.measure(self.lat1, self.lon1, self.lat2, self.lon1)
        meters_lon = self.measure(self.lat1, self.lon1, self.lat1, self.lon2)

        if self.lon1 and self.lon2 < 0:
            meters_lon = -meters_lon

        lat_meters_points = np.linspace(0, meters_lat, num=self.points)
        lon_meters_points = np.linspace(0, meters_lon, num=self.points)
        
       
        self.X, self.Y = np.meshgrid(lon_points, lat_points)

        self.lisZ = np.array([
            [srtm1_data.get_altitude(latitude=lat, longitude=lon) for lon in lon_points]
            for lat in lat_points
        ])

        
        
    def pointView(self, ObsPoint, RefPoint):
        """
        Set the observation point and reference point.

        Parameters:
        obs_point (tuple): Observation point as (latitude, longitude).
        ref_point (tuple or str): Reference point as (latitude, longitude) or "max" for the highest point in lis_z.
        lis_z (numpy array): Matrix of elevation data.
        x_grid (numpy array): Matrix of longitude data.
        y_grid (numpy array): Matrix of latitude data.
        measure_function (function): Function to calculate distance between two points.
        get_altitude_function (function): Function to get altitude for a given latitude and longitude.

        Returns:
        None
        """

        # Observation point
        self.obsPX = ObsPoint[1]
        self.obsPY = ObsPoint[0]
        self.obsPZ = srtm1_data.get_altitude(latitude=self.obsPY, longitude=self.obsPX)
        
        # Reference point
        if RefPoint=="max":
            self.maxpoint = np.where(self.lisZ==self.lisZ.max())
            self.RefPX = self.X[self.maxpoint[0],self.maxpoint[1]]
            self.RefPY = self.Y[self.maxpoint[0],self.maxpoint[1]]
            self.RefPZ = self.lisZ.max()
        else:
            self.RefPX=RefPoint[1]
            self.RefPY=RefPoint[0]
            self.RefPZ=srtm1_data.get_altitude(latitude=self.RefPY, longitude=self.RefPX)

        self.R = self.measure(self.obsPX, self.obsPY, self.RefPX, self.RefPY)
        self.RefZEN = np.arccos((self.RefPZ - self.obsPZ)/self.R)*180.0/np.pi
        print ("Elevation = %f" % (90.0 - self.RefZEN))
         
        
    def plot_lines(self, cenit, azimut): 
        """ Plot the geological structure with the moun's lines
        """
        
        cenitF = cenit[0]
        cenitS = cenit[1]
        cenitP = cenit[2]
        
        azimutF = azimut[0]
        azimutS = azimut[1]
        azimutP = azimut[2]
        
        # Equation creation
        t = 3 # Projection distance in km
        self.PjecX = self.obsPX+(self.RefPX-self.obsPX)*t
        self.PjecY = self.obsPY+(self.RefPY-self.obsPY)*t
        self.PjecZ = self.obsPZ+(self.RefPZ-self.obsPZ)*t
        
        
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(20,17)) 
        surf = ax.plot_surface(self.X, self.Y, self.lisZ, edgecolors='grey', alpha=0.5,  cmap=self.cmap,
                            linewidth=0, antialiased=False, label=self.name)
        
        # Plot points
        ax.scatter(self.RefPX,self.RefPY,self.RefPZ, s=20, c="r")          # Reference point
        ax.scatter(self.obsPX, self.obsPY, self.obsPZ, c="k")        # Observer point
        ax.plot([self.obsPX,self.PjecX],[self.obsPY,self.PjecY],     # Line projection                                                  
                [self.obsPZ,self.PjecZ], color="red")                      
        
        
        # Plot muon flux projections

        cenit = np.linspace(cenitS,cenitF,cenitP)
        azimut = np.linspace(azimutF,azimutS,azimutP)
        self.azimut = azimut
        self.cenit = cenit # corregir variable / creada para el funcionamiento de h flujo
        self.distances = np.zeros((cenitP,azimutP))
        t = np.linspace(0,1,1000) # Projection lines in km
        for k,i in enumerate(azimut):
            for m,j in enumerate(cenit):
                
                self.pointfX = self.obsPX + math.cos(i * np.pi / 180) * (self.PjecX - self.obsPX) - math.sin(i * np.pi / 180) * (self.PjecY - self.obsPY)
                self.pointfY = self.obsPY + math.sin(i * np.pi / 180) * (self.PjecX - self.obsPX) + math.cos(i * np.pi / 180) * (self.PjecY - self.obsPY)
                self.pointfZ = self.PjecZ + math.sin(j * np.pi / 180) * self.PjecZ 
                #=======================================================
                self.equationX = self.obsPX+(self.pointfX-self.obsPX)*t
                self.equationY = self.obsPY+(self.pointfY-self.obsPY)*t
                self.equationZ = self.obsPZ+(self.pointfZ-self.obsPZ)*t
                
                # Call function calculate distance and add to matrix distances 
                self.distances[m][k]= self.calculate_distance(self.equationX,self.equationY,self.equationZ)
                ax.plot(self.equationX,self.equationY,self.equationZ, color="black", alpha=0.05) # Line equation
             
        # Labels                                            
        ax.tick_params(axis='both', which='major',labelrotation=20, rotation=25, labelsize=20)
        #ax.xlabel('longitud', fontsize=22, labelpad=20, rotation=10)
        ax.set_xlabel("Longitude", fontsize=22, labelpad=30, rotation=10)
        ax.set_ylabel("Latitude", fontsize=22, labelpad=30, rotation=-60)
        ax.set_zlabel("Altitude", fontsize=22, labelpad=30) 
        
        ax.set_title(self.name + " topography", fontsize=25)
        
        # Save the 3D plot as a PDF file
        output_pdf = self.name + "Structure.pdf"
        with PdfPages(output_pdf) as pdf:
            pdf.savefig(fig)
            
        plt.show()
                
        Lmcenit = []
        Lmazimut = []
        for i in self.cenit:
            for j in self.azimut:
               Lmcenit.append(i)
               Lmazimut.append(j)

        matrixDatos = np.empty((self.distances.shape[0]*self.distances.shape[1],3))        
        matrixDatos[:,0] = np.radians(Lmcenit)
        matrixDatos[:,1] = np.radians(Lmazimut)
        matrixDatos[:,2] = self.distances.flatten()
        
        self.matrixDatos = matrixDatos
    
    def Flux(self):
        self.Flux = FluxModels.FluxIntegrated(self.matrixDatos)
        self.Flux.showData()
        self.Flux.Elevation(self.obsPX,self.RefPX,self.obsPY,self.RefPY,self.obsPZ,self.RefPZ)
        self.Flux.OpenSky()
        Cross_Flux, Open_Flux = self.Flux.ShowIntegratedFlux()
        # Flux.ShowOpacity()
        # Flux.ShowDensity()

        return Cross_Flux, Open_Flux
    
    def Transmission(self):
        
        return self.Flux.Transmission()
    
        
    def calculate_distance(self,equationX,equationY,equationZ):
        """ This function calculate the distance for each point
        """
        d=0
        indexlist = []
        
        self.projection=np.zeros(len(equationX))
        for i in range(len(equationX)):
            self.projection[i] =srtm1_data.get_altitude(latitude=equationY[i], longitude=equationX[i])
    
        coor = equationZ<self.projection
        
        for j,i in enumerate(coor):
            if i==True:
                if j==0:
                    if coor[j+1]==True:
                        indexlist.append(j)
                else: 
                    if coor[j+1]==False and coor[j-1]==True:
                        indexlist.append(j)
                    elif coor[j-1]==False and coor[j+1]==False:
                        dtemp = self.measure(equationY[j-1],equationX[j-1],equationY[j+1],equationX[j+1])
                        d+=dtemp/3
                    elif coor[j-1]==False:
                        indexlist.append(j)
                        
        x,y = equationX[indexlist],equationY[indexlist]
        listp = list(zip(indexlist[0::2], indexlist[1::2]))
        
        for i,j in enumerate(listp):
            dtemp = self.measure(equationY[listp[i][0]],equationX[listp[i][0]],equationY[listp[i][1]],equationX[listp[i][1]])
            d+=dtemp
        return d
    
# Plotting functions
    
    def plot_structure(self): 
        """Plot the geological structure
        """
        
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(20,17)) # Geological structure
        surf = ax.plot_surface(self.X, self.Y, self.lisZ, edgecolors='grey', alpha=0.5,  cmap=self.cmap,
                            linewidth=0, antialiased=False)
        
        
        # Labels                                            
        ax.tick_params(axis='both', which='major',labelrotation=20, rotation=25, labelsize=20)
        #ax.xlabel('longitud', fontsize=22, labelpad=20, rotation=10)
        ax.set_xlabel("Longitude", fontsize=22, labelpad=30, rotation=10)
        ax.set_ylabel("Latitude", fontsize=22, labelpad=30, rotation=-60)
        ax.set_zlabel("Altitude", fontsize=22, labelpad=30)
        ax.set_title(self.name, fontsize=25)
        
        # Save the 3D plot as a PDF file
        output_pdf = self.name + "Structure.pdf"
        with PdfPages(output_pdf) as pdf:
            pdf.savefig(fig)

        # Display the 3D plot
        plt.show()
   
        
    def plot_points(self): 
        """ Plot the geological structure with the moun's lines
        """
        
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(20,17)) 
        surf = ax.plot_surface(self.X, self.Y, self.lisZ, edgecolors='grey', alpha=0.5,  cmap=self.cmap,
                            linewidth=0, antialiased=False, label=self.name)
        
    # Plot points
        ax.scatter(self.RefPX,self.RefPY,self.RefPZ, c="black", s=100)     # Reference point
        ax.scatter(self.obsPX, self.obsPY, self.obsPZ, c="r", s=100)        # Telescope point

                            
        ax.tick_params(axis='both', which='major',labelrotation=20, rotation=25, labelsize=20)
        #ax.xlabel('longitud', fontsize=22, labelpad=20, rotation=10)
        ax.set_xlabel("Longitude", fontsize=22, labelpad=30, rotation=10)
        ax.set_ylabel("Latitude", fontsize=22, labelpad=30, rotation=-60)
        ax.set_zlabel("Altitude", fontsize=22, labelpad=30) 
        ax.set_title("Observer at "+self.name, fontsize=25)
        
        
          # Save the 3D plot as a PDF file
        output_pdf = self.name + "Structure_points.pdf"
        with PdfPages(output_pdf) as pdf:
            pdf.savefig(fig)

        plt.show()
        
        
    def section(self):
        """ This function plot the section of the geolical structure traversed
        """
        
        fig, ax = plt.subplots(figsize=(20,17))
        ax.fill_between(self.equationX, self.projection, y2=0, alpha=0.5, color="navy")
        ax.tick_params(axis='both', which='major', labelsize=20)

        ax.plot(self.equationX,self.equationZ, color="red")
         # Labels
        ax.set_ylim([2400,3000])   # Establecer el rango en el eje Y
        # Labels
  
        #ax.set_ylim([2200,3000])
        ax.set_xlabel("Latitude", fontsize=25, labelpad=20)
        ax.set_ylabel("Altitude", fontsize=25, labelpad=20)
        ax.set_title(self.name, fontsize=25)                                 # Save the 3D plot as a PDF file
        output_pdf = self.name + "section.pdf"
        with PdfPages(output_pdf) as pdf:
            pdf.savefig(fig)        
        plt.show()


        
    def show_distances(self):
        """ This function plot the distance traversed in geological structured
        """
        fig, ax = plt.subplots(figsize=(10,7))
        im = ax.imshow(self.distances/1e3, interpolation='nearest', origin='upper', cmap=self.cmap)
        ax.set_xlabel("Azimuth [degree]", fontsize = 15)
        ax.set_ylabel("Zenith [degree]", fontsize = 15)
        ax.set_title("Rock thickness : "+self.name, fontsize = 15)

        # Create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.1 inch.
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)

        # Color bar
        clb = plt.colorbar(im, cax=cax)
        clb.set_label('d [km]', fontsize = 15)
        clb.ax.tick_params(labelsize = 12)
        
        
        # Define contour levels
        contour_levels = np.linspace(self.distances.min()/1e3, self.distances.max()/1e3, 4)  # 4 contour lines

        # Draw contour lines
        contour_lines = ax.contour(self.distances/1e3, contour_levels, colors='k')

        # Add labels to contour lines
        ax.clabel(contour_lines, inline=True, fontsize=12, colors='k')

        
        labelsx_indices = np.linspace(0, self.distances.shape[1]-1, 11).astype(int)
        labelsy_indices = np.linspace(0, self.distances.shape[0]-1, 11).astype(int)

        ax.set_xticks(labelsx_indices)
        ax.set_yticks(labelsy_indices)

        # Put the actual azimuth and zenith values as labels
        ax.set_xticklabels(np.round(np.linspace(min(self.azimut), max(self.azimut), 11), 0), fontsize=12)
        ax.set_yticklabels(np.round(np.linspace(self.RefZEN + min(self.cenit),self.RefZEN + max(self.cenit), 11), 0), fontsize=12)

        output_pdf = self.name + "Rock_Thickness_Zenith.png"
        
        with PdfPages(output_pdf) as pdf:
            pdf.savefig(fig)  

        plt.show()

    def show_distances_no_extent(self):
        
        """ This function plot the distance traversed in geological structured
        """
        fig, ax = plt.subplots(figsize=(10,7))

        im = ax.imshow(self.distances/1e3, interpolation='nearest', origin='upper', cmap=self.cmap)

        ax.set_xlabel("Azimuth [degree]", fontsize = 15)
        ax.set_ylabel("Elevation [degree]", fontsize = 15)
        ax.set_title("Rock thickness : "+self.name, fontsize = 15)

        # Create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.1 inch.
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)

        # Color bar
        clb = plt.colorbar(im, cax=cax)
        clb.set_label('d [km]', fontsize = 15)
        clb.ax.tick_params(labelsize = 12)

        # Define contour levels
        contour_levels = np.linspace(self.distances.min()/1e3, self.distances.max()/1e3, 4)  # 4 contour lines

        # Draw contour lines (without specifying origin)
        contour_lines = ax.contour(self.distances/1e3, contour_levels, colors='k')

        # Add labels to contour lines
        ax.clabel(contour_lines, inline=True, fontsize=10, colors='white')

        # Adjust these to match your data's array indices
        labelsx_indices = np.linspace(0, self.distances.shape[1]-1, 11).astype(int)
        labelsy_indices = np.linspace(0, self.distances.shape[0]-1, 11).astype(int)

        ax.set_xticks(labelsx_indices)
        ax.set_yticks(labelsy_indices)

        # Put the actual azimuth and zenith values as labels
        ax.set_xticklabels(np.round(np.linspace(min(self.azimut), max(self.azimut), 11), 0), fontsize=12)
        ax.set_yticklabels(np.round(np.linspace(90.0 - self.RefZEN + max(self.cenit), 90.0 - self.RefZEN + min(self.cenit), 11), 0), fontsize=12)

        output_pdf = self.name + "Rock_Thickness_Elevation.png"
        with PdfPages(output_pdf) as pdf:
            pdf.savefig(fig)  

        plt.show()

    
    def plot_mesured_flux(self, data, plot_label, plot_units):
        """ This function plot a variable in the observation frame
        """
        fig, ax = plt.subplots(figsize=(10,7))
        im = ax.imshow(data, interpolation='nearest', origin='upper', cmap=self.cmap, norm=mpl.colors.LogNorm())
        ax.set_xlabel("Azimuth [degree]", fontsize = 15)
        ax.set_ylabel("Zenith [degree]", fontsize = 15)
        ax.set_title(plot_label, fontsize = 15)

        # Create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.1 inch.
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)

        # Color bar
        clb = plt.colorbar(im, cax=cax)
        clb.set_label(plot_label + plot_units, fontsize = 15)
        clb.ax.tick_params(labelsize = 12)
        
        
        labelsx_indices = np.linspace(0, data.shape[1]-1, 11).astype(int)
        labelsy_indices = np.linspace(0, data.shape[0]-1, 11).astype(int)

        ax.set_xticks(labelsx_indices)
        ax.set_yticks(labelsy_indices)

        # Put the actual azimuth and zenith values as labels
        ax.set_xticklabels(np.round(np.linspace(min(self.azimut), max(self.azimut), 11), 0), fontsize=12)
        ax.set_yticklabels(np.round(np.linspace(self.RefZEN + min(self.cenit),self.RefZEN + max(self.cenit), 11), 0), fontsize=12)

        output_pdf = self.name + plot_label + ".png"
        
        with PdfPages(output_pdf) as pdf:
            pdf.savefig(fig)  

        plt.show()
        
# Functions on progress
   
    def save_data(self,file_name):
        """ This function save .dat data 
        """
        
        np.savetxt(file_name+'.dat', self.matrixDatos)
        
        
        # Add first 4 lines
        
        header_lines = ["# Op " + str(self.obsPX) + " " + str(self.obsPY) + " " + str(self.obsPZ) + "\n", 
                    "# Pp " + str(self.RefPX) + " " + str(self.RefPY) + " " + str(self.RefPZ) + "\n",
                    "# A "  + str(self.lat1)  + " " + str(self.lat2)  + " " + str(self.lon1) + " " + str(self.lon2) + "\n",
                    "# cenith[rad] azimuth[rad] distance[km]\n"]
    
    
        with open(f"{file_name}.dat", "r") as file:
            file_lines = file.readlines()

        new_file_lines = header_lines + file_lines

        with open(f"{file_name}.dat", "w") as file:
            file.writelines(new_file_lines)   
            
        

    def point_show(self):
        """ This function is on progress to show multiple points
        """
        
        #a = self.lisZ==self.obsPZ
        zmin = self.obsPZ-10
        zmax = self.obsPZ+10
        indexX = []
        indexY = []
        

        # List of latitude and longitude
        lisLat = np.linspace(self.lat1, self.lat2, num=self.points)
        lisLong = np.linspace(self.lon1, self.lon2, num=self.points)
        
        # List of latitude and longitude in meteres
        meterspointLat = self.measure(self.lat1, self.lon1, self.lat2, self.lon1)
        meterspointLong = self.measure(self.lat1, self.lon1, self.lat1, self.lon2)
        
        
        self.X, self.Y = np.meshgrid(lisLong, lisLat) # Creation grid in point with meters
        
        # Calculated the elevation
        lisZ=np.zeros([self.points,self.points])
        for i in range(self.points):
            for j in range(self.points):
                lisZ[i][j]=srtm1_data.get_altitude(latitude=lisLat[i], longitude=lisLong[j])
                if self.lisZ[i][j]>zmin and self.lisZ[i][j]<zmax:
                    indexX.append(lisLat[i])
                    indexY.append(lisLong[j])
        return indexX,indexY


