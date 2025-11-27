import numpy as np
import matplotlib.pyplot as plt
from trustutils import run

run_dir=run.BUILD_DIRECTORY

class run_IC_Mu:

    def __init__(self,cross_sec):
       
        self.cross_sec=cross_sec
        
    def execute(self):

        self.data_dict = {
                    "variable": "U",
                    "std_dev": 0.001,#std deviation of the input
                    "cross_sec": self.cross_sec,
                    "figure": ["1","2","3","4","5"]
                }

        
        self.test_case = "Mu/Lid_driven_cavity_Poly_Chaos_"


        self.dx= np.loadtxt(run_dir+"/"+self.test_case+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[0])
        self.state_x_pcm= np.loadtxt(run_dir+"/"+self.test_case+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[1])
        self.sens_x_pcm= np.loadtxt(run_dir+"/"+self.test_case+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'_SENS.coupe', unpack=True, usecols=[1])
        self.state_y_pcm= np.loadtxt(run_dir+"/"+self.test_case+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[2])
        self.sens_y_pcm= np.loadtxt(run_dir+"/"+self.test_case+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'_SENS.coupe', unpack=True, usecols=[2])

        self.variance_x_pcm = np.abs(self.sens_x_pcm)*self.data_dict["std_dev"]
        self.state_x_pcm = self.state_x_pcm + self.sens_x_pcm*self.data_dict["std_dev"]
        self.variance_y_pcm = np.abs(self.sens_y_pcm)*self.data_dict["std_dev"]
        self.state_y_pcm = self.state_y_pcm + self.sens_y_pcm*self.data_dict["std_dev"]


        self.test_case = "Taylor/Mu/Lid_mu_"

        self.state_x= np.loadtxt(run_dir+"/"+self.test_case+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[1])
        self.sens_x= np.loadtxt(run_dir+"/"+self.test_case+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'_SENS.coupe', unpack=True, usecols=[1])
        self.state_y= np.loadtxt(run_dir+"/"+self.test_case+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[2])
        self.sens_y= np.loadtxt(run_dir+"/"+self.test_case+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'_SENS.coupe', unpack=True, usecols=[2])


        self.variance_x = np.abs(self.sens_x)*self.data_dict["std_dev"]
        self.variance_y = np.abs(self.sens_y)*self.data_dict["std_dev"]


class plot_case_IC_Mu(run_IC_Mu):

    def __init__(self,cross_sec,i):

        super().__init__(cross_sec)
         
        self.cross_sec=cross_sec
        self.i=i

    def affichage_graph(self):

        choix=self.data_dict["figure"][self.i]

        match choix:

            case "1":

                
                ref= np.loadtxt(run_dir+"/"+'Ref_sens/ux_H.txt', unpack=True, usecols=[1])
                t_ref= np.loadtxt(run_dir+"/"+'Ref_sens/ux_H.txt', unpack=True, usecols=[0])
                fig, ax = plt.subplots()
                ax.plot(t_ref, ref, self.dx, self.variance_x, self.dx, self.variance_x_pcm, 'r8', linewidth=2)
                ax.legend([ 'El-Beltagy and Wafa', 'Taylor', 'PCM'])
                plt.xlabel("x")
                plt.ylabel("Horizontal velocity component standard deviation")
                #plt.savefig('Taylor_PCM_Variance_Velocity_H_x.png')
              

            case "2":
               
                ref= np.loadtxt(run_dir+"/"+'Ref_etat/Ref5_y_velocity_horizontal.dat', unpack=True, usecols=[1])
                t_ref= np.loadtxt(run_dir+"/"+'Ref_etat/Ref5_y_velocity_horizontal.dat', unpack=True, usecols=[0])
                fig, ax = plt.subplots()
                ax.plot(t_ref, ref, self.dx, self.state_y, self.dx, self.state_y_pcm, 'r8', linewidth=2)
                ax.legend([ 'Marchi et al.', 'Taylor', 'PCM'])
                plt.xlabel("x")
                plt.ylabel("Horizontal velocity component average")
                #plt.savefig('Taylor_PCM_Mean_Velocity_H_y.png')
               

            case "3":
              
                ref= np.loadtxt(run_dir+"/"+'Ref_etat/Ref5_x_velocity_vertical.dat', unpack=True, usecols=[1])
                t_ref= np.loadtxt(run_dir+"/"+'Ref_etat/Ref5_x_velocity_vertical.dat', unpack=True, usecols=[0])

                fig, ax = plt.subplots()
                ax.plot(t_ref, ref,  self.dx, self.state_x, self.dx, self.state_x_pcm, 'r8', linewidth=2)
                ax.legend([ 'Marchi et al.', 'Taylor', 'PCM'])
                plt.xlabel("y")
                plt.ylabel("Vertical velocity component average")
                #plt.savefig('Taylor_PCM_Mean_Velocity_V_x.png')
                
            case "4":

                #plot standard deviation
                ref= np.loadtxt(run_dir+"/"+'Ref_sens/ux_V.txt', unpack=True, usecols=[1])
                t_ref= np.loadtxt(run_dir+"/"+'Ref_sens/ux_V.txt', unpack=True, usecols=[0])

                fig, ax = plt.subplots()
                ax.plot(t_ref, ref, self.dx, self.variance_x, self.dx, self.variance_x_pcm, 'r8', linewidth=2)
                ax.legend([ 'El-Beltagy and Wafa', 'Taylor', 'PCM'])
                plt.xlabel("y")
                plt.ylabel("Horizontal velocity component standard deviation")
                #plt.savefig('Taylor_PCM_Variance_Velocity_V_x.png')
             
            case "5":
            
                #plot standard deviation
                ref= np.loadtxt(run_dir+"/"+'Ref_sens/uy_V.txt', unpack=True, usecols=[1])
                t_ref= np.loadtxt(run_dir+"/"+'Ref_sens/uy_V.txt', unpack=True, usecols=[0])

                fig, ax = plt.subplots()
                ax.plot(t_ref, ref, self.dx, self.variance_y, self.dx, self.variance_y_pcm, 'r8', linewidth=2)
                ax.legend([ 'El-Beltagy and Wafa', 'Taylor', 'PCM'])
                plt.xlabel("y")
                plt.ylabel("Vertical velocity component standard deviation")
                #plt.savefig('Taylor_PCM_Variance_Velocity_V_y.png')
                
