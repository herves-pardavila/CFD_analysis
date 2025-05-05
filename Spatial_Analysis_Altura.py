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
from scipy.interpolate import griddata
plt.close("all")

class Spatial_Analysis_Altura:
    
    def __init__(self, path,ustar, part,time_ini,time_fin,wind_dir,field,coords = ["X (m)", "Y (m)"],color_map ="seismic"):
        """
        

        Parameters
        ----------
        path : str
            Folder where the Tables are stored
        ustar : str
            Friction velocity. For example, u_star=0.33_
        part : str
            Treshold or plane
        time_ini : int
            Time of the initial table (seconds)
        time_fin : int
            Time of the last table (seconds)
        wind_dir : str
            Wind direction. For example NO_
        field : str
            An string identifying the variables the table has
        color_map : TYPE, optional
            DESCRIPTION. The default is "seismic".

        Returns
        -------
        None.

        """
        
        self.path=path
        self.ustar=ustar
        self.part=part
        self.time_ini=time_ini
        self.time_fin=time_fin
        self.wind_dir=wind_dir
        self.field=field
        self.coords = coords
        self.color_map = color_map
        
        hours=np.arange(self.time_ini, self.time_fin + 3600,3600)
        self.hours=hours
        
        
        
        return
    
    def see_variables(self,hour):
        
        df=pd.read_csv(self.path+self.ustar+self.wind_dir+str(hour)+"_"+self.field+self.part+".csv")
        
        table_variables=df.columns
        
        print(table_variables)
        
        self.table_variables=table_variables
        
        return
        
        
    
    
    def load_single_hour(self,hour,variable,norm=1):
        
        df=pd.read_csv(self.path+self.ustar+self.wind_dir+str(hour)+"_"+self.field+self.part+".csv")
        
        df.sort_values(by=self.coords[0],inplace=True,ignore_index=True)
        
        data=np.array(df[variable])
        
        self.data=data.T*norm
       
        X = np.array(df[self.coords[0]])
        self.X =X
        
        Y = np.array(df[self.coords[1]])
        self.Y = Y
        
        area=np.array(df["Area: Magnitude (m^2)"])
        self.area=area
        
        return
    
    def plot_single_hour(self,hour,variable,bounds=None,area_factor=1,norm=1,save=False,smooth=False):
        
        self.load_single_hour(hour, variable,norm=norm)
        
        if bounds != None:
            vmin,vmax=bounds[0],bounds[1]
        
        else:
            vmin,vmax=self.data.min(),self.data.max()
        
        
        
        fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(area_factor*6.4,area_factor*4.8))
        fig.subplots_adjust(top=0.95,bottom=0.1,left=0.1,right=0.95)
       
        ax1.set_aspect('equal', adjustable='box')
        ax2.set_aspect('equal', adjustable='box')
        ax1.set_title("Sur",fontsize=area_factor*10)
        ax2.set_title("Norte",fontsize=area_factor*10)
        
        ax1.set_xlim(100, 120)  # outliers only
        ax2.set_xlim(210, 230)  # most of the data
        
        im1=ax1.scatter(self.X,self.Y,s=30*area_factor*self.area**2,c=self.data,marker="h",cmap=self.color_map,vmin=vmin,vmax=vmax,edgecolor=None)
        im2=ax2.scatter(self.X,self.Y,s=30*area_factor*self.area**2,c=self.data,marker="h",cmap=self.color_map,vmin=vmin,vmax=vmax,edgecolor=None)
        ax1.set_xlabel(self.coords[0],fontsize=area_factor*10)
        ax2.set_xlabel(self.coords[0],fontsize=area_factor*10)
        ax1.set_ylabel(self.coords[1],fontsize=area_factor*10)
        ax1.set_facecolor('gray')
        ax2.set_facecolor('gray')
        
        ax1.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax1.yaxis.tick_left()
        ax1.tick_params(labelright='off')
        ax2.yaxis.tick_right()
        #ax2.set_yticklabels([])
        
        d = .015  # how big to make the diagonal lines in axes coordinates
        # arguments to pass plot, just so we don't keep repeating them
        kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
        ax1.plot((1-d, 1+d), (-d, +d), **kwargs)
        ax1.plot((1-d, 1+d), (1-d, 1+d), **kwargs)
        
        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1-d, 1+d), **kwargs)
        ax2.plot((-d, +d), (-d, +d), **kwargs)
        
       
        
        
        ax1.tick_params(axis='both', which='both', labelsize=10*area_factor)
        ax2.tick_params(axis='both', which='both', labelsize=10*area_factor)
        
        fig.colorbar(im1, ax=[ax1, ax2], orientation='horizontal', fraction=0.05, pad=0.1)
       
    
      
           
       
        plt.show()
        
        return
    def plot_single_hour_smoothed(self,hour,variable,bounds=None,area_factor=1,norm=1,save=False,smooth=False):
        
        self.load_single_hour(hour, variable,norm=norm)
        
        if bounds != None:
            vmin,vmax=bounds[0],bounds[1]
        
        else:
            vmin,vmax=self.data.min(),self.data.max()
        
        
    
        minx,maxx,dimx=min(self.X),max(self.X),len(self.X)
        miny,maxy,dimy=min(self.Y),max(self.Y),len(self.Y)

        xt=np.arange(minx,maxx,0.05)
        yt=np.arange(miny,maxy,0.05)

        X,Y=np.meshgrid(xt,yt)

        #self.data = self.data - (273.15 + 38)

        Z = griddata((self.X, self.Y), self.data, (X,Y), method='linear',rescale=True)
        
        Z[ ( (X > 115) & (X < 215) ) & (Y <15)]=np.nan
        
        fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(area_factor*6.4,area_factor*4.8))
        fig.subplots_adjust(top=0.95,bottom=0.1,left=0.1,right=0.95)
       
        ax1.set_aspect('equal', adjustable='box')
        ax2.set_aspect('equal', adjustable='box')
        ax1.set_title("Sur",fontsize=area_factor*10)
        ax2.set_title("Norte",fontsize=area_factor*10)
        
        ax1.set_xlim(100, 120)  # outliers only
        ax2.set_xlim(210, 230)  # most of the data
        
        im1=ax1.pcolormesh(X,Y,Z ,cmap=self.color_map,vmin=vmin,vmax=vmax)
        im2=ax2.pcolormesh(X,Y,Z,cmap=self.color_map,vmin=vmin,vmax=vmax)
        ax1.set_xlabel(self.coords[0],fontsize=area_factor*10)
        ax2.set_xlabel(self.coords[0],fontsize=area_factor*10)
        ax1.set_ylabel(self.coords[1],fontsize=area_factor*10)
        ax1.set_facecolor('gray')
        ax2.set_facecolor('gray')
        
        ax1.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax1.yaxis.tick_left()
        ax1.tick_params(labelright='off')
        ax2.yaxis.tick_right()
        #ax2.set_yticklabels([])
        
        d = .025  # how big to make the diagonal lines in axes coordinates
        # arguments to pass plot, just so we don't keep repeating them
        kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
        ax1.plot((1-d, 1+d), (-d, +d), **kwargs)
        ax1.plot((1-d, 1+d), (1-d, 1+d), **kwargs)
        
        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1-d, 1+d), **kwargs)
        ax2.plot((-d, +d), (-d, +d), **kwargs)
        
        
        
        
       
        fig.colorbar(im1, ax=[ax1, ax2], orientation='horizontal', fraction=0.05, pad=0.11, location = "bottom")
       
        #cbar=fig.colorbar(im1,location="right")
        ax1.tick_params(axis='both', which='both', labelsize=10*area_factor)
        ax2.tick_params(axis='both', which='both', labelsize=10*area_factor)
        #cbar.ax.tick_params(labelsize=10*area_factor) 
       
      
           
       
        plt.show()
        
        return
    
    def plot_single_hour_velocity(self,hour,bounds=None,area_factor=1,skip=2,scale=100,rotate=False,save=False):
        
        self.load_single_hour(hour,"Velocity[j] (m/s)")
        U=self.data
        self.load_single_hour(hour,"Velocity[k] (m/s)")
        V=self.data
        
        speed=np.sqrt(U*U+V*V)
        
        if bounds != None:
            vmin,vmax=bounds[0],bounds[1]
        
        else:
            vmin,vmax=speed.min(),speed.max()
        
        fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(area_factor*6.4,area_factor*4.8))
        fig.subplots_adjust(top=0.95,bottom=0.1,left=0.1,right=0.95)
        
        ax1.set_aspect('equal', adjustable='box')
        ax2.set_aspect('equal', adjustable='box')
        ax1.set_title("Sur",fontsize=area_factor*10)
        ax2.set_title("Norte",fontsize=area_factor*10)
        
        ax1.set_xlim(100, 120)  # outliers only
        ax2.set_xlim(210, 230)  # most of the data
        
        ax1.set_ylim(0, 20)  # outliers only
        ax2.set_ylim(0, 20)  # most of the data
        
        ax1.set_xlabel(self.coords[0],fontsize=area_factor*10)
        ax1.set_ylabel(self.coords[1],fontsize=area_factor*10)
        ax1.set_facecolor('gray')
        ax2.set_xlabel(self.coords[0],fontsize=area_factor*10)
        ax2.set_ylabel(self.coords[1],fontsize=area_factor*10)
        ax2.set_facecolor('gray')
       
        
        
       
        
       
        im1=ax1.quiver(self.X[::skip],self.Y[::skip],U[::skip],V[::skip],
                    speed[::skip],scale=scale,headlength=4, width=0.005,
                     cmap=self.color_map,clim=(vmin,vmax))
        
        im2=ax2.quiver(self.X[::skip],self.Y[::skip],U[::skip],V[::skip],
                    speed[::skip],scale=scale,headlength=4,width=0.005,
                     cmap=self.color_map,clim=(vmin,vmax))
        
        
        ax1.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax1.yaxis.tick_left()
        ax1.tick_params(labelright='off')
        ax2.yaxis.tick_right()
        #ax2.set_yticklabels([])
        
        d = .015  # how big to make the diagonal lines in axes coordinates
        # arguments to pass plot, just so we don't keep repeating them
        kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
        ax1.plot((1-d, 1+d), (-d, +d), **kwargs)
        ax1.plot((1-d, 1+d), (1-d, 1+d), **kwargs)
        
        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1-d, 1+d), **kwargs)
        ax2.plot((-d, +d), (-d, +d), **kwargs)
        
       
        
     
        cbar=fig.colorbar(im1, ax=[ax1, ax2], orientation='horizontal', fraction=0.05, pad=0.11, location = "bottom")
        ax1.tick_params(axis='both', which='both', labelsize=10*area_factor)
        ax2.tick_params(axis='both', which='both', labelsize=10*area_factor)
        cbar.ax.tick_params(labelsize=10*area_factor) 
    
   

        
       
       
        return
        
        
        
    def load_data(self,variable,norm=1):
        

        matrix=[]
        
        for hour in tqdm(range(len(self.hours))):
            
            try:
               
                self.load_single_hour(self.hours[hour], variable,norm=norm)
                
                matrix+=[self.data]
                
                # if hour==0:
                    
                #     X = np.array(df["X (m)"])
                #     self.X =X
                    
                #     Y = np.array(df["Y (m)"])
                #     self.Y = Y
                    
                #     area=np.array(df["Area: Magnitude (m^2)"])
                #     self.area=area
                
            except FileNotFoundError:
                print("File Not Found for %i seconds" %self.hours[hour])
                
        self.matrix=np.array(matrix).T
        
        return 
    
    
    
    def plot_any_data(self,data,title="",bounds=None,area_factor=1,smooth=False):
        
        if bounds != None:
            vmin,vmax=bounds[0],bounds[1]
        
        else:
            vmin,vmax=data.min(),data.max()
        
        
        fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(area_factor*6.4,area_factor*4.8))
        fig.subplots_adjust(top=0.95,bottom=0.1,left=0.1,right=0.95)
           
        ax1.set_aspect('equal', adjustable='box')
        ax2.set_aspect('equal', adjustable='box')
        ax1.set_title("Sur",fontsize=area_factor*10)
        ax2.set_title("Norte",fontsize=area_factor*10)
        
        ax1.set_xlim(100, 120)  # outliers only
        ax2.set_xlim(210, 230)  # most of the data
        
        im1=ax1.scatter(self.X,self.Y,s=30*area_factor*self.area**2,c=data,marker="h",cmap=self.color_map,vmin=vmin,vmax=vmax,edgecolor=None)
        im2=ax2.scatter(self.X,self.Y,s=30*area_factor*self.area**2,c=data,marker="h",cmap=self.color_map,vmin=vmin,vmax=vmax,edgecolor=None)
        ax1.set_xlabel(self.coords[0],fontsize=area_factor*10)
        ax2.set_xlabel(self.coords[0],fontsize=area_factor*10)
        ax1.set_ylabel(self.coords[1],fontsize=area_factor*10)
        ax1.set_facecolor('gray')
        ax2.set_facecolor('gray')
        
        ax1.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax1.yaxis.tick_left()
        ax1.tick_params(labelright='off')
        ax2.yaxis.tick_right()
        #ax2.set_yticklabels([])
        
        d = .025  # how big to make the diagonal lines in axes coordinates
        # arguments to pass plot, just so we don't keep repeating them
        kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
        ax1.plot((1-d, 1+d), (-d, +d), **kwargs)
        ax1.plot((1-d, 1+d), (1-d, 1+d), **kwargs)
        
        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1-d, 1+d), **kwargs)
        ax2.plot((-d, +d), (-d, +d), **kwargs)
        
           
        
        cbar=fig.colorbar(im1,location="top")
        ax1.tick_params(axis='both', which='both', labelsize=10*area_factor)
        ax2.tick_params(axis='both', which='both', labelsize=10*area_factor)
        cbar.ax.tick_params(labelsize=10*area_factor) 


     
 
        plt.show()
        
        return
        
        # if smooth==True:
            
        #     #interpolación, tarda mucho con 8GB de RAM y no aporta mucho
        #     minx,maxx,dimx=min(self.X),max(self.X),len(self.X)
        #     miny,maxy,dimy=min(self.Y),max(self.Y),len(self.Y)

        #     xt=np.linspace(minx,maxx,dimx)
        #     yt=np.linspace(miny,maxy,dimy)

        #     X,Y=np.meshgrid(xt,yt)


        #     Z = griddata((self.X, self.Y), data, (X,Y), method='linear',rescale=True)
        #     Z[ ( (X > 115) & (X < 215) ) & (Y > 115) & (Y < 215)]=np.nan

        #     fig=plt.figure(figsize=(area_factor*6.4,area_factor*4.8))
        #     im=plt.pcolormesh(X,Y,Z,cmap="seismic")
        #     plt.gca().set_aspect('equal')
        #     plt.gca().set_xlabel("X (m) ",fontsize=area_factor*10)
        #     plt.gca().set_ylabel("Y (m) ",fontsize=area_factor*10)
        #     plt.gca().set_facecolor('gray')
        #     fig.suptitle(title +self.wind_dir + " " + self.ustar,fontsize=10*area_factor)
            
        #     cbar=fig.colorbar(im,location="right")
        #     plt.gca().tick_params(axis='both', which='both', labelsize=10*area_factor)
        #     cbar.ax.tick_params(labelsize=10*area_factor) 
        #     fig.text(0.4,0.2,"South",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        #     fig.text(0.4,0.75,"North",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        #     fig.text(0.65,0.5,"East",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        #     fig.text(0.2,0.5,"West",transform=ax.transAxes,color="white",fontsize=10*area_factor)
            
        #     plt.show()
        
        return
    

    
 
        

            
            
        
    
if __name__ == "__main__":
    rho_air=1.18415 #kg/m3
    Qv=9.041e3 #kg/s
    A= 9.041e3 #m2
    
    #path_to_folder=r"C:/TFM_DAVIDH_2025/Tablas/" #PC CIEMAT
    

    path_to_folder=r"C:/Users/herve/Documents/TFM_DAVIDH_2025/Tablas/" #Asus
       
    campo="PS_T_V_" #passive scalar
    part="threshold_plane_x=165m" 
    con=Spatial_Analysis_Altura(path=path_to_folder, ustar="ustar0.33_", part=part,
                         time_ini=3600, time_fin= 176400, wind_dir="NO_", field=campo,coords = ["Y (m)", "Z (m)"],color_map="coolwarm")
    #con.see_variables(61200)
    
    con.plot_single_hour_smoothed(144000,"PS2",bounds=None,area_factor=1,norm = rho_air * 0.33 * A /Qv,save=False,smooth=False)
    PS_21h=con.data
    con.plot_single_hour_smoothed(122400,"PS2",bounds=None,area_factor=1,norm=rho_air * 0.33 * A /Qv,save=False,smooth=False)
    PS_15h=con.data
    
    con.plot_any_data(PS_15h-PS_21h,bounds=[-1,1])
    
    #con.plot_single_hour_velocity(122400,area_factor=1.5,bounds=[0,1.5],skip=2,scale=10,rotate=False)

                
    

    
