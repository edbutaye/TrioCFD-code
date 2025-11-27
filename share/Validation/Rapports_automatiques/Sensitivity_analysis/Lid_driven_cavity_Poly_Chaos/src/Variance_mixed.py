# Sensitivity
import numpy as np
import matplotlib.pyplot as plt
from trustutils import run


run_dir=run.BUILD_DIRECTORY


class run_Variance_mixed:

    def __init__(self,cross_sec):

        self.cross_sec=cross_sec

    def execute(self):

        self.data_dict = {

            "variable"  : "U",
            "direction" : "x",            #used only for vectorial data_dict["variable"]s (i.e. velocity). x, y or z
            "cross_sec" : self.cross_sec, #"H"
            "std_dev_v" : 0.1,            #std deviation of the input
            "std_dev_mu": 0.001,          #std deviation of the input
            "figure"    :["1","2","3","4"]

        }


        self.test_case_v = "Taylor/Velocity/Lid_v_"
        self.test_case_mu = "Taylor/Mu/Lid_mu_"

        self.abscissa= np.loadtxt(run_dir+"/"+self.test_case_v+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[0])
        state_v= np.loadtxt(run_dir+"/"+self.test_case_v+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[1])
        sens_v= np.loadtxt(run_dir+"/"+self.test_case_v+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'_SENS.coupe', unpack=True, usecols=[1])
        state_mu= np.loadtxt(run_dir+"/"+self.test_case_mu+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[1])
        sens_mu= np.loadtxt(run_dir+"/"+self.test_case_mu+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'_SENS.coupe', unpack=True, usecols=[1])
        self.variance_Taylor =np.abs(sens_v*self.data_dict["std_dev_v"] + sens_mu*self.data_dict["std_dev_mu"])

########## Sensitivity section  Velocity ##########

        self.test_case_v = "Velocity/Lid_driven_cavity_Poly_Chaos_"
        self.test_case_mu = "Mu/Lid_driven_cavity_Poly_Chaos_"

        state_v_pcm= np.loadtxt(run_dir+"/"+self.test_case_v+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[1])
        sens_v_pcm= np.loadtxt(run_dir+"/"+self.test_case_v+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'_SENS.coupe', unpack=True, usecols=[1])
        state_mu_pcm= np.loadtxt(run_dir+"/"+self.test_case_mu+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[1])
        sens_mu_pcm= np.loadtxt(run_dir+"/"+self.test_case_mu+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'_SENS.coupe', unpack=True, usecols=[1])
        self.variance_pcm =np.abs(sens_v_pcm*self.data_dict["std_dev_v"] + sens_mu_pcm*self.data_dict["std_dev_mu"])

        self.abscissa= np.loadtxt(run_dir+"/"+self.test_case_v+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[0])
        state_v= np.loadtxt(run_dir+"/"+self.test_case_v+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[1])
        sens_v= np.loadtxt(run_dir+"/"+self.test_case_v+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'_SENS.coupe', unpack=True, usecols=[1])
        state_mu= np.loadtxt(run_dir+"/"+self.test_case_mu+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[1])
        sens_mu= np.loadtxt(run_dir+"/"+self.test_case_mu+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'_SENS.coupe', unpack=True, usecols=[1])

        self.variance_v =np.abs(sens_v*self.data_dict["std_dev_v"] )
        self.variance_mu =np.abs(sens_mu*self.data_dict["std_dev_mu"])
        self.variance =np.abs(sens_v*self.data_dict["std_dev_v"] + sens_mu*self.data_dict["std_dev_mu"])

        


class plot_case(run_Variance_mixed):
    def __init__(self,cross_sec,i):

        super().__init__(cross_sec)

        self.cross_sec=cross_sec
        self.i=i
    def affichage_graph(self):

        choix=self.data_dict["figure"][self.i]

        match choix:

            case "1":
               
                ref= np.loadtxt(run_dir+"/"+'Ref_sens/mixed.txt', unpack=True, usecols=[1])
                t_ref= np.loadtxt(run_dir+"/"+'Ref_sens/mixed.txt', unpack=True, usecols=[0])

                fig, ax = plt.subplots()
                ax.plot(t_ref, ref, self.abscissa, self.variance_pcm, self.abscissa, self.variance_Taylor,'r8', linewidth=2)
                ax.legend([ 'El-Beltagy and Wafa', 'Taylor', 'PCM'])
                plt.xlabel("x")
                plt.ylabel("Horizontal velocity standard deviation")
               #plt.savefig('Taylor_PCM_Variance_'+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+self.data_dict["direction"]+'mixed_v_plus_mu.png')
         
            case "2":

                ref= np.loadtxt(run_dir+"/"+'Ref_sens/mixed.txt', unpack=True, usecols=[1])
                t_ref= np.loadtxt(run_dir+"/"+'Ref_sens/mixed.txt', unpack=True, usecols=[0])
                ref_v= np.loadtxt(run_dir+"/"+'Ref_sens/velocity_H_x.txt', unpack=True, usecols=[1])
                t_ref_v= np.loadtxt(run_dir+"/"+'Ref_sens/velocity_H_x.txt', unpack=True, usecols=[0])
                ref_mu= np.loadtxt(run_dir+"/"+'Ref_sens/ux_H.txt', unpack=True, usecols=[1])
                t_ref_mu= np.loadtxt(run_dir+"/"+'Ref_sens/ux_H.txt', unpack=True, usecols=[0])

                fig, ax = plt.subplots()
                ax.plot(self.abscissa, self.variance_v, 'ro', t_ref_v, ref_v, 'r', self.abscissa, self.variance_mu,'bo', t_ref_mu, ref_mu, 'b', self.abscissa, self.variance, 'go', t_ref, ref, 'g', linewidth=2, fillstyle='none')
                ax.legend([ 'Velocity 10%', 'El-Beltagy and Wafa', 'Viscosity 10 %', 'El-Beltagy and Wafa', 'Mixed', 'El-Beltagy and Wafa'])
                plt.xlabel("x")
                plt.ylabel("Horizontal velocity standard deviation")
                #plt.savefig('Variance_'+self.data_dict["variable"]+'_'+self.data_dict["cross_sec"]+'_'+self.data_dict["direction"]+'.png')
              
            case "3":

                 vx_x =np.loadtxt(run_dir+"/"+"IC_mixed_V_x.txt") 
                 vy_x =np.loadtxt(run_dir+"/"+"IC_mixed_V_y.txt") 

                 x = vx_x[:, 0]
                 vx = vx_x[:, 1]
                 vy = vy_x[:, 1]
                 ex = np.abs(vx_x[:, 2])
                 ey = np.abs(vy_x[:, 2])
                 plt.figure(figsize=(8, 6))
                 
                 
                 plt.errorbar(x, vx, yerr=ex, fmt='-', color='red', capsize=2, markersize=1,label='_nolegend_',elinewidth=1)
                 plt.errorbar(x, vy, yerr=ey, fmt='-', color='blue',  capsize=2, markersize=1,label='_nolegend_',elinewidth=1)

                 plt.plot(x, vx, color='red', label='x-component')
                 plt.plot(x, vy, color='blue', label='y-component')
                 plt.ylim(-0.8, 1.4)
                 plt.xlim(0.0, 1)

                 plt.xlabel("x [m]")
                 plt.ylabel("Velocity [m/s]")
                 plt.legend(loc='upper center', fontsize=12, borderpad=1.5, labelspacing=1.5)
                 plt.tight_layout()
                
            case "4":

                 vx_y = np.loadtxt(run_dir+"/"+"IC_mixed_H_x.txt")
                 vy_y = np.loadtxt(run_dir+"/"+"IC_mixed_H_y.txt")

                 y = vx_y[:, 0]
                 vx = vx_y[:, 1]
                 ex = np.abs(vx_y[:, 2])
                 vy = vy_y[:, 1]
                 ey = np.abs(vy_y[:, 2])

                 plt.figure(figsize=(8,6))

                 # tracé avec barres d'erreur
                 plt.errorbar(y, vx, yerr=ex, fmt='-', color='red', capsize=2, markersize=1,label='_nolegend_',elinewidth=1)
              
                 plt.errorbar(y, vy, yerr=ey, fmt='-', color='blue',  capsize=2, markersize=1,label='_nolegend_',elinewidth=1)

                 # tracé des courbes
                 plt.plot(y, vx, color='red', label='x-component')
                 plt.plot(y, vy, color='blue', label='y-component')
                 
                 plt.ylim(-0.8, 1)
                 plt.xlim(0.0, 1)

 
                 plt.xlabel("y [m]")
                 plt.ylabel("Velocity [m/s]")
                 plt.legend(loc='upper center', fontsize=12, borderpad=1.5, labelspacing=1.5)
 
                 plt.tight_layout()
                


