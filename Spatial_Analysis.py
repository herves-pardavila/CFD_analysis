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
#plt.close("all")

class Spatial_Analysis:
    
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
        
        
        
        fig=plt.figure(figsize=(area_factor*6.4,area_factor*4.8))
        fig.subplots_adjust(top=0.95,bottom=0.1,left=0.05,right=0.95)
        ax=fig.add_subplot(111)
        ax.set_aspect('equal', adjustable='box')
        im=ax.scatter(self.X,self.Y,s=area_factor*self.area**2,c=self.data,marker="h",cmap=self.color_map,vmin=vmin,vmax=vmax,edgecolor=None)
        ax.set_xlabel(self.coords[0],fontsize=area_factor*10)
        ax.set_ylabel(self.coords[1],fontsize=area_factor*10)
        ax.set_facecolor('gray')
        
        hour_for_title =((hour-3600)/3600)+6
    
            
        # if hour_for_title > 23:
        #     hour_for_title= str(int(hour_for_title - 24)) +"h dia2"
            
        # #ax.set_title(variable +" "+ str(hour_for_title) + "h " + self.wind_dir + " " + self.ustar,fontsize=10*area_factor)
        
        cbar=fig.colorbar(im,location="right")
        ax.tick_params(axis='both', which='both', labelsize=10*area_factor)
        cbar.ax.tick_params(labelsize=10*area_factor) 
        fig.text(0.4,0.2,"South",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        fig.text(0.4,0.75,"North",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        fig.text(0.65,0.5,"East",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        fig.text(0.2,0.5,"West",transform=ax.transAxes,color="white",fontsize=10*area_factor)
       
        # if save == True:
            
        #     fig.savefig("./savefigs/"+str(self.ustar)+"_"+str(self.wind_dir)+"_"+str(self.part)+"_"+str(hour)+".png")
        #     #fig.savefig("./savefigs/"+str(self.ustar)+"_"+str(self.wind_dir)+"_"+str(variable)+"_"+str(self.part)+"_"+str(hour)+".png")
        #     plt.close()
            
        # if smooth==True:
            
        #     #interpolación, tarda mucho con 8GB de RAM y no aporta mucho
        #     minx,maxx,dimx=min(self.X),max(self.X),len(self.X)
        #     miny,maxy,dimy=min(self.Y),max(self.Y),len(self.Y)

        #     xt=np.linspace(minx,maxx,dimx)
        #     yt=np.linspace(miny,maxy,dimy)

        #     X,Y=np.meshgrid(xt,yt)


        #     Z = griddata((self.X, self.Y), self.data, (X,Y), method='linear',rescale=True)
        #     Z[ ( (X > 115) & (X < 215) ) & (Y > 115) & (Y < 215)]=np.nan

        #     fig=plt.figure(figsize=(area_factor*6.4,area_factor*4.8))
        #     im=plt.pcolormesh(X,Y,Z,cmap="seismic")
        #     plt.gca().set_aspect('equal')
        #     plt.gca().set_xlabel("X (m) ",fontsize=area_factor*10)
        #     plt.gca().set_ylabel("Y (m) ",fontsize=area_factor*10)
        #     plt.gca().set_facecolor('gray')
        #     fig.suptitle(variable +" "+ str(hour_for_title) + "h " + self.wind_dir + " " + self.ustar,fontsize=10*area_factor)
            
        #     cbar=fig.colorbar(im,location="right")
        #     plt.gca().tick_params(axis='both', which='both', labelsize=10*area_factor)
        #     cbar.ax.tick_params(labelsize=10*area_factor) 
        #     fig.text(0.4,0.2,"Sur",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        #     fig.text(0.4,0.75,"Norte",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        #     fig.text(0.65,0.5,"Este",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        #     fig.text(0.2,0.5,"Oeste",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        #     fig.savefig("./savefigs/interpoladas/"+str(self.ustar)+"_"+str(self.wind_dir)+"_"+str(variable)+"_"+str(self.part)+"_"+str(hour)+".png")
        #     plt.close()
                
             
           
       
        plt.show()
        
        return
    
    def plot_single_hour_velocity(self,hour,bounds=None,area_factor=1,skip=2,scale=100,rotate=False,save=False):
        
        self.load_single_hour(hour,"Velocity[i] (m/s)")
        U=self.data
        self.load_single_hour(hour,"Velocity[j] (m/s)")
        V=self.data
        
        speed=np.sqrt(U*U+V*V)
        
        if bounds != None:
            vmin,vmax=bounds[0],bounds[1]
        
        else:
            vmin,vmax=speed.min(),speed.max()
        
        fig=plt.figure(figsize=(area_factor*6.4,area_factor*4.8))
        fig.subplots_adjust(top=0.95,bottom=0.1,left=0.05,right=0.95)
        #fig.suptitle("Figure 3",fontsize=area_factor*10)
        ax=fig.add_subplot(111)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel(self.coords[0],fontsize=area_factor*10)
        ax.set_ylabel(self.coords[1],fontsize=area_factor*10)
        ax.set_facecolor('gray')
       
        
        hour_for_title =((hour-3600)/3600)+6
   
           
        if hour_for_title > 23:
            hour_for_title= str(int(hour_for_title - 24)) +"h dia2"
           
        #ax.set_title("Velocity (m/s) "+ str(hour_for_title) + "h " + self.wind_dir + " " + self.ustar,fontsize=10*area_factor)
  
       
        
       
        im=ax.quiver(self.X[::skip],self.Y[::skip],U[::skip],V[::skip],
                    speed[::skip],scale=scale,headlength=4,
                     cmap=self.color_map,clim=(vmin,vmax))
        
        cbar=fig.colorbar(im,location="right")
        ax.tick_params(axis='both', which='both', labelsize=10*area_factor)
        cbar.ax.tick_params(labelsize=10*area_factor) 
        fig.text(0.4,0.2,"South",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        fig.text(0.4,0.75,"North",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        fig.text(0.65,0.5,"East",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        fig.text(0.2,0.5,"West",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        
        if save == True:
            
          
            fig.savefig("./savefigs/"+str(self.ustar)+"_"+str(self.wind_dir)+"_"+"Velocity"+"_"+str(self.part)+"_"+str(hour)+".png")
            plt.close()
       
       
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
    
    def compute_mean(self,variable):
        
        self.load_data(variable)
        
        c_mean=np.mean(self.matrix,axis=1)
        
        self.mean=c_mean
    
        
        return 
    
    def plot_mean(self,variable):
        
        self.compute_mean(variable)
        
        fig=plt.figure()
        ax=fig.add_subplot(111)
        im=ax.scatter(self.X,self.Y,s=self.area,c=self.mean,marker="h",cmap=self.color_map,edgecolor=None)
        ax.set_aspect('equal', adjustable='box')
        ax.set_facecolor('white')
        fig.colorbar(im,location="right")
        ax.set_xlabel(self.coords[0])
        ax.set_ylabel(self.coords[1])
        ax.title.set_text(variable + " Mean")
        
        #fig.savefig("imagen_scatter_numpy.png")
        plt.show()
        
        return
    
    def compute_median(self,variable):
        
        self.load_data(variable)
        
        c_median=np.median(self.matrix,axis=1)
        
        self.median=c_median
    
        
        return 
    
    def plot_median(self,variable):
        
        self.compute_median(variable)
        
        fig=plt.figure()
        ax=fig.add_subplot(111)
        im=ax.scatter(self.X,self.Y,s=self.area,c=self.median,marker="h",cmap=self.color_map,edgecolor=None)
        ax.set_aspect('equal', adjustable='box')
        ax.set_facecolor('white')
        fig.colorbar(im,location="right")
        ax.set_xlabel(self.coords[0])
        ax.set_ylabel(self.coords[1])
        ax.title.set_text(variable + " Median")
        
        #fig.savefig("imagen_scatter_numpy.png")
        plt.show()
        
        return
    def compute_std(self,variable):
        
        self.load_data(variable)
        
        c_std=np.std(self.matrix,axis=1)
        
        self.std=c_std
    
        
        return 
    
    def plot_std(self,variable):
        
        self.compute_std(variable)
        
        fig=plt.figure()
        ax=fig.add_subplot(111)
        im=ax.scatter(self.X,self.Y,s=self.area,c=self.std,marker="h",cmap=self.color_map,edgecolor=None)
        ax.set_aspect('equal', adjustable='box')
        ax.set_facecolor('white')
        fig.colorbar(im,location="right")
        ax.set_xlabel(self.coords[0])
        ax.set_ylabel(self.coords[1])
        ax.title.set_text(variable + " Standard Deviation")
        
        #fig.savefig("imagen_scatter_numpy.png")
        plt.show()
        
        return
        
        
    def compute_kurtosis(self,variable):
        
        self.load_data(variable)
        
        c_kurtosis=st.kurtosis(self.matrix,axis=1)
        
        self.kurtosis=c_kurtosis
    
        
        return 
    
    
    def plot_kurtosis(self,variable,bounds=None):
        
        self.compute_kurtosis(variable)
        
        if bounds != None:
            vmin,vmax=bounds[0],bounds[1]
        
        else:
            vmin,vmax=self.kurtosis.min(),self.kurtosis.max()
        
        fig=plt.figure()
        ax=fig.add_subplot(111)
        im=ax.scatter(self.X,self.Y,s=self.area,c=self.kurtosis,marker="h",vmin=vmin,vmax=vmax,cmap=self.color_map,edgecolor=None)
        ax.set_aspect('equal', adjustable='box')
        ax.set_facecolor('white')
        fig.colorbar(im,location="right")
        ax.set_xlabel(self.coords[0])
        ax.set_ylabel(self.coords[1])
        ax.title.set_text(variable + " Kurtosis")
        
        #fig.savefig("imagen_scatter_numpy.png")
        plt.show()
        
        
        return 
    
    def plot_any_data(self,data,title="",bounds=None,area_factor=1,smooth=False):
        
        if bounds != None:
            vmin,vmax=bounds[0],bounds[1]
        
        else:
            vmin,vmax=data.min(),data.max()
        
        
        fig=plt.figure(figsize=(area_factor*6.4,area_factor*4.8))
        ax=fig.add_subplot(111)
        ax.set_aspect('equal', adjustable='box')
        im=ax.scatter(self.X,self.Y,s=area_factor*self.area**2,c=data,marker="h",cmap="seismic",vmin=vmin,vmax=vmax,edgecolor=None)
        ax.set_xlabel(self.coords[0],fontsize=area_factor*10)
        ax.set_ylabel(self.coords[1],fontsize=area_factor*10)
        ax.set_facecolor('gray')
        #ax.set_title(title +self.wind_dir + " " + self.ustar,fontsize=10*area_factor)
        
        
        cbar=fig.colorbar(im,location="right")
        ax.tick_params(axis='both', which='both', labelsize=10*area_factor)
        cbar.ax.tick_params(labelsize=10*area_factor) 
        fig.text(0.4,0.2,"Sur",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        fig.text(0.4,0.75,"Norte",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        fig.text(0.65,0.5,"Este",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        fig.text(0.2,0.5,"Oeste",transform=ax.transAxes,color="white",fontsize=10*area_factor)
       
        plt.show()
        
        if smooth==True:
            
            #interpolación, tarda mucho con 8GB de RAM y no aporta mucho
            minx,maxx,dimx=min(self.X),max(self.X),len(self.X)
            miny,maxy,dimy=min(self.Y),max(self.Y),len(self.Y)

            xt=np.linspace(minx,maxx,dimx)
            yt=np.linspace(miny,maxy,dimy)

            X,Y=np.meshgrid(xt,yt)


            Z = griddata((self.X, self.Y), data, (X,Y), method='linear',rescale=True)
            Z[ ( (X > 115) & (X < 215) ) & (Y > 115) & (Y < 215)]=np.nan

            fig=plt.figure(figsize=(area_factor*6.4,area_factor*4.8))
            im=plt.pcolormesh(X,Y,Z,cmap="seismic")
            plt.gca().set_aspect('equal')
            plt.gca().set_xlabel("X (m) ",fontsize=area_factor*10)
            plt.gca().set_ylabel("Y (m) ",fontsize=area_factor*10)
            plt.gca().set_facecolor('gray')
            fig.suptitle(title +self.wind_dir + " " + self.ustar,fontsize=10*area_factor)
            
            cbar=fig.colorbar(im,location="right")
            plt.gca().tick_params(axis='both', which='both', labelsize=10*area_factor)
            cbar.ax.tick_params(labelsize=10*area_factor) 
            fig.text(0.4,0.2,"South",transform=ax.transAxes,color="white",fontsize=10*area_factor)
            fig.text(0.4,0.75,"North",transform=ax.transAxes,color="white",fontsize=10*area_factor)
            fig.text(0.65,0.5,"East",transform=ax.transAxes,color="white",fontsize=10*area_factor)
            fig.text(0.2,0.5,"West",transform=ax.transAxes,color="white",fontsize=10*area_factor)
            
            plt.show()
        
        return
    

    
 
        

            
            
        
    
if __name__ == "__main__":
    
    rho_air=1.18415 #kg/m3
    Qv=9.041e3 #kg/s
    A= 9.041e3 #m2
    
    #path_to_folder=r"C:/TFM_DAVIDH_2025/Tablas/" #PC CIEMAT
    
    path_to_folder=r"C:/Users/herve/Documents/TFM_DAVIDH_2025/Tablas/" #Asus
   
    campo="PS_T_V_" #passive scalar
    part="threshold_plane_z=3m" 
    con=Spatial_Analysis(path=path_to_folder, ustar="ustar0.33_", part=part,
                         time_ini=3600, time_fin= 176400, wind_dir="NO_", field=campo,color_map="coolwarm")
    con.see_variables(61200)

    #con.plot_single_hour_velocity(122400,area_factor=1.5,bounds=[1,3],skip=4,scale=60,rotate=False)
    #con.plot_single_hour_velocity(144000,area_factor=1.5,bounds=[1,3],skip=4,scale=60,rotate=False)
    #con.plot_single_hour(122400,"PS2",bounds=[0,3],area_factor=1.5,norm=rho_air * 0.33 * A /Qv,save=False,smooth=False)
    #con.plot_single_hour(144000 ,"PS2",bounds = [0,3] ,area_factor=1.5,norm=rho_air * 0.33 * A /Qv,save=False,smooth=False)
    
    con.load_single_hour(144000,"Velocity[k] (m/s)",norm = 1)
    PS_6h=con.data
    con.load_single_hour(122400,"Velocity[k] (m/s)",norm = 1)
    PS_15h=con.data
    
    con.plot_any_data(PS_15h-PS_6h,title=" ",area_factor=1.5,bounds=None)
    
   # plt.close("all")
    

    


    # for hour in con.hours:
        
    #     try:
    #         con.plot_single_hour(hour, 'Prom_w_primaT_prima_top_S (m-K/s)',bounds=[0,0.2],area_factor=1.5,norm=1, save=True,smooth=False)
    #         #con.plot_single_hour_velocity(hour,area_factor=1.5,skip=4,scale=60,rotate=False,save=True)
    #     except KeyError:
    #         continue
            
                
                
    

    
