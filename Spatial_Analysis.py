# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 16:29:34 2025

@author: David Hervés Pardavila

Class for spatial analysis of CFD results. 
"""

import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import scipy.stats as st
#plt.close("all")

class Spatial_Analysis:
    
    def __init__(self, path,ustar, part,time_ini,time_fin,wind_dir,field):
        
        self.path=path
        self.ustar=ustar
        self.part=part
        self.time_ini=time_ini
        self.time_fin=time_fin
        self.wind_dir=wind_dir
        self.field=field
        
        return
    
    
    def load_single_hour(self,hour,variable):
        
        df=pd.read_csv(self.path+self.ustar+self.wind_dir+str(hour)+"_"+self.field+self.part+".csv")
        matrix=np.array(df[variable])
        
        self.matrix=matrix.T
       
        X = np.array(df["X (m)"])
        self.X =X
        
        Y = np.array(df["Y (m)"])
        self.Y = Y
        
        area=np.array(df["Area: Magnitude (m^2)"])
        self.area=area
        
        return
    
    def plot_single_hour(self,hour,variable):
        
        self.load_single_hour(hour, variable)
        
        fig=plt.figure()
        ax=fig.add_subplot(111)
        im=ax.scatter(self.X,self.Y,s=self.area,c=self.matrix,marker="h",cmap="seismic",edgecolor=None)
        ax.set_aspect('equal', adjustable='box')
        ax.set_facecolor('white')
        fig.colorbar(im,location="right")
        ax.set_xlabel("X (m) ")
        ax.set_ylabel("Y (m) ")
        ax.title.set_text(variable +" "+ str(hour))
        
        #fig.savefig("imagen_scatter_numpy.png")
        plt.show()
        
        
        
    def load_data(self,variable):
        
        hours=np.arange(self.time_ini, self.time_fin + 3600,3600)
        
        matrix=[]
        
        for hour in tqdm(range(len(hours))):
            
            try:
                df=pd.read_csv(self.path+self.ustar+self.wind_dir+str(hours[hour])+"_"+self.field+self.part+".csv")
                
                matrix+=[np.array(df[variable])]
                
                
                if hour==0:
                    
                    X = np.array(df["X (m)"])
                    self.X =X
                    
                    Y = np.array(df["Y (m)"])
                    self.Y = Y
                    
                    area=np.array(df["Area: Magnitude (m^2)"])
                    self.area=area
                
            except FileNotFoundError:
                print("File Not Found for %i seconds" %hour)
                
        self.matrix=np.array(matrix).T
        
        return 
    
    def mean(self,variable):
        
        self.load_data(variable)
        
        c_mean=np.mean(self.matrix,axis=1)
        
        self.c_mean=c_mean
    
        
        return 
    
    def plot_mean(self,variable):
        
        self.mean(variable)
        
        fig=plt.figure()
        ax=fig.add_subplot(111)
        im=ax.scatter(self.X,self.Y,s=self.area,c=self.c_mean,marker="h",cmap="coolwarm",edgecolor=None)
        ax.set_aspect('equal', adjustable='box')
        ax.set_facecolor('white')
        fig.colorbar(im,location="right")
        ax.set_xlabel("X (m) ")
        ax.set_ylabel("Y (m) ")
        ax.title.set_text(variable + " Mean")
        
        #fig.savefig("imagen_scatter_numpy.png")
        plt.show()
        
        return
    
    def median(self,variable):
        
        self.load_data(variable)
        
        c_median=np.median(self.matrix,axis=1)
        
        self.c_median=c_median
    
        
        return 
    
    def plot_median(self,variable):
        
        self.median(variable)
        
        fig=plt.figure()
        ax=fig.add_subplot(111)
        im=ax.scatter(self.X,self.Y,s=self.area,c=self.c_median,marker="h",cmap="coolwarm",edgecolor=None)
        ax.set_aspect('equal', adjustable='box')
        ax.set_facecolor('white')
        fig.colorbar(im,location="right")
        ax.set_xlabel("X (m) ")
        ax.set_ylabel("Y (m) ")
        ax.title.set_text(variable + " Median")
        
        #fig.savefig("imagen_scatter_numpy.png")
        plt.show()
        
        return
    def std(self,variable):
        
        self.load_data(variable)
        
        c_std=np.std(self.matrix,axis=1)
        
        self.c_std=c_std
    
        
        return 
    
    def plot_std(self,variable):
        
        self.std(variable)
        
        fig=plt.figure()
        ax=fig.add_subplot(111)
        im=ax.scatter(self.X,self.Y,s=self.area,c=self.c_std,marker="h",cmap="coolwarm",edgecolor=None)
        ax.set_aspect('equal', adjustable='box')
        ax.set_facecolor('white')
        fig.colorbar(im,location="right")
        ax.set_xlabel("X (m) ")
        ax.set_ylabel("Y (m) ")
        ax.title.set_text(variable + " Standard Deviation")
        
        #fig.savefig("imagen_scatter_numpy.png")
        plt.show()
        
        return
        
        
    def kurtosis(self,variable):
        
        self.load_data(variable)
        
        c_kurtosis=st.kurtosis(self.matrix,axis=1)
        
        self.c_kurtosis=c_kurtosis
    
        
        return 
    
    
    def plot_kurtosis(self,variable):
        
        self.kurtosis(variable)
        
        fig=plt.figure()
        ax=fig.add_subplot(111)
        im=ax.scatter(self.X,self.Y,s=self.area,c=self.c_kurtosis,marker="h",cmap="coolwarm",edgecolor=None)
        ax.set_aspect('equal', adjustable='box')
        ax.set_facecolor('white')
        fig.colorbar(im,location="right")
        ax.set_xlabel("X (m) ")
        ax.set_ylabel("Y (m) ")
        ax.title.set_text(variable + " Kurtosis")
        
        #fig.savefig("imagen_scatter_numpy.png")
        plt.show()
        
        
        return 
        
    
if __name__ == "__main__":
    
    path_to_folder=r"C:/TFM_DAVIDH_2025/Tablas/"
    campo="PS_T_V_" #passive scalar
    part="threshold_plane_z=3m" 
    con=Spatial_Analysis(path=path_to_folder, ustar="ustar0.33_", part=part,time_ini=7200, time_fin= 176400, wind_dir="NO_", field=campo)
    
    variable = "PS1"
    
    con.mean(variable)
    con.median(variable)
    con.std(variable)
    con.kurtosis(variable)
    
    fig=plt.figure(figsize = (20,20) )

    ax1=fig.add_subplot(221)
    ax2=fig.add_subplot(222)
    ax3=fig.add_subplot(223)
    ax4=fig.add_subplot(224)
    
    im1=ax1.scatter(con.X,con.Y,s=con.area,c=con.c_mean,marker="h",cmap="seismic",edgecolor=None)
    ax1.set_aspect('equal', adjustable='box')
    ax1.set_facecolor('black')
    fig.colorbar(im1,ax=ax1, location="right")
    ax1.set_xlabel("X (m) ", fontsize=15)
    ax1.set_ylabel("Y (m) ", fontsize=15)
    ax1.set_title(variable + " Mean", fontsize=15)
    
    im2=ax2.scatter(con.X,con.Y,s=con.area,c=con.c_median,marker="h",cmap="coolwarm",edgecolor=None)
    ax2.set_aspect('equal', adjustable='box')
    ax2.set_facecolor('black')
    fig.colorbar(im2,ax=ax2, location="right")
    ax2.set_xlabel("X (m) ", fontsize=15)
    ax2.set_ylabel("Y (m) ", fontsize=15)
    ax2.set_title(variable + " Median", fontsize=15)
    
    im3=ax3.scatter(con.X,con.Y,s=con.area,c=con.c_std,marker="h",cmap="coolwarm",edgecolor=None)
    ax3.set_aspect('equal', adjustable='box')
    ax3.set_facecolor('black')
    fig.colorbar(im3,ax=ax3, location="right")
    ax3.set_xlabel("X (m) ", fontsize=15)
    ax3.set_ylabel("Y (m) ", fontsize=15)
    ax3.set_title(variable + " Standard Deviation", fontsize=15)
    
    im4=ax4.scatter(con.X,con.Y,s=con.area,c=con.c_kurtosis,marker="h",cmap="coolwarm",edgecolor=None)
    ax4.set_aspect('equal', adjustable='box')
    ax4.set_facecolor('black')
    fig.colorbar(im4,ax=ax4, location="right")
    ax4.set_xlabel("X (m) ", fontsize=15)
    ax4.set_ylabel("Y (m) ", fontsize=15)
    ax4.set_title(variable + " Kurtosis", fontsize=15)
    
    con.plot_single_hour(57600, "PS2")
    
    
    con.load_single_hour(57600, "PS2")
    las_15=con.matrix
    con.load_single_hour(79200, "Turbulent Kinetic Energy (J/kg)")
    las_21=con.matrix
    
    fig=plt.figure()
    ax=fig.add_subplot(111)
    im=ax.scatter(con.X,con.Y,s=con.area,c=las_15-las_21,marker="h",cmap="seismic",edgecolor=None)
    ax.set_aspect('equal', adjustable='box')
    ax.set_facecolor('white')
    fig.colorbar(im,location="right")
    ax.set_xlabel("X (m) ")
    ax.set_ylabel("Y (m) ")
    #ax.title.set_text(variable +" "+ str(hour))
    
    
