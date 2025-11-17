# Sensitivity
import numpy as np
import os
import sys
import matplotlib.pyplot as plt

sys.path.append("build")
os.chdir("build")
run_dir= os.getcwd()

class run_Variance_mixed:
    def __init__(self,cross_sec,label):

        self.cross_sec=cross_sec
        self.label=label
    def execute(self):

        self.data_dict = {
                    "params": ["Tbottom", "Ttop", "Beta", "Mu", "Lambda_prime"],
                    "mean_params": [363, 313, 3e-3, 1.54e-5, 0.21847e-4],
                    "variable_1": "VITESSE",
                    "variable_2": "TEMPERATURE",
                    "variable_3": "PRESSION",
                    "confidence": 0.95,
                    "cross_sec": self.cross_sec,
                    "label":self.label,
                    "figure":["1","2","3","4","5","6"]
                    
                }

        self.var_vx  = 0
        self.var_vy = 0
        self.var_t = 0 
        self.var_p = 0
        self.pcm_var_vx=0
        self.pcm_var_vy=0
        self.pcm_var_t=0
        self.pcm_var_p=0

        sigma = self.data_dict["mean_params"]

        for i in range(len(self.data_dict["params"])):

            sigma[i] =self.data_dict["mean_params"][i]/100.

            self.test_case =run_dir+"/"+"Sensibility_Rayleigh_Bernard_Taylor/"+self.data_dict["params"][i]+"/Sensibility_Rayleigh_Bernard_"

            self.dx= np.loadtxt(self.test_case+self.data_dict["variable_1"]+'_SENSIBILITE'+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[0])
            self.sens_vx= np.loadtxt(self.test_case+self.data_dict["variable_1"]+'_SENSIBILITE'+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[1])
            self.sens_vy= np.loadtxt(self.test_case+self.data_dict["variable_1"]+'_SENSIBILITE'+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[2])
            self.sens_t= np.loadtxt(self.test_case+self.data_dict["variable_2"]+'_SENSIBILITE'+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[1])
            self.sens_p= np.loadtxt(self.test_case+self.data_dict["variable_3"]+'_SENSIBILITE'+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[1])

            self.var_vx += abs(self.sens_vx)*sigma[i]
            self.var_vy += abs(self.sens_vy)*sigma[i]
            self.var_t += abs(self.sens_t)*sigma[i]
            self.var_p += abs(self.sens_p)*sigma[i]

            self.test_case = run_dir+"/"+"Sensibility_Rayleigh_Bernard_PCM/"+self.data_dict["params"][i]+"/Sensibility_Rayleigh_Bernard_"

            self.dx= np.loadtxt(self.test_case+self.data_dict["variable_1"]+'_SENSIBILITE''_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[0])
            pcm_sens_vx= np.loadtxt(self.test_case+self.data_dict["variable_1"]+'_SENSIBILITE'+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[1])
            pcm_sens_vy= np.loadtxt(self.test_case+self.data_dict["variable_1"]+'_SENSIBILITE'+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[2])
            pcm_sens_t= np.loadtxt(self.test_case+self.data_dict["variable_2"]+'_SENSIBILITE'+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[1])
            pcm_sens_p= np.loadtxt(self.test_case+self.data_dict["variable_3"]+'_SENSIBILITE'+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[1])

            self.pcm_var_vx += abs(pcm_sens_vx)*sigma[i]
            self.pcm_var_vy += abs(pcm_sens_vy)*sigma[i]
            self.pcm_var_t += abs(pcm_sens_t)*sigma[i]
            self.pcm_var_p += abs(pcm_sens_p)*sigma[i]

        self.test_case = run_dir+"/"+"Sensibility_Rayleigh_Bernard_PCM/"+self.data_dict["params"][0]+"/Sensibility_Rayleigh_Bernard_"
        state_vx= np.loadtxt(self.test_case+self.data_dict["variable_1"]+'_ETAT'+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[1])
        upper_bound = state_vx+1/np.sqrt(1-self.data_dict["confidence"])*self.pcm_var_vx
        lower_bound = state_vx-1/np.sqrt(1-self.data_dict["confidence"])*self.pcm_var_vx
        DataOut = np.column_stack((self.dx,state_vx, upper_bound, lower_bound))
        np.savetxt('IC_mixed_'+self.data_dict["cross_sec"]+'_'+'Vx'+'.txt', DataOut)

        state_vy= np.loadtxt(self.test_case+self.data_dict["variable_1"]+'_ETAT'+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[2])
        upper_bound = state_vy+1/np.sqrt(1-self.data_dict["confidence"])*self.pcm_var_vy
        lower_bound = state_vy-1/np.sqrt(1-self.data_dict["confidence"])*self.pcm_var_vy
        DataOut = np.column_stack((self.dx,state_vy, upper_bound, lower_bound))
        np.savetxt('IC_mixed_'+self.data_dict["cross_sec"]+'_'+'Vy'+'.txt', DataOut)            
	                    

        state_p= np.loadtxt(self.test_case+self.data_dict["variable_2"]+'_ETAT'+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[1])
        upper_bound = state_p+1/np.sqrt(1-self.data_dict["confidence"])*self.pcm_var_p
        lower_bound = state_p-1/np.sqrt(1-self.data_dict["confidence"])*self.pcm_var_p
        DataOut = np.column_stack((self.dx,state_p, upper_bound, lower_bound))
        np.savetxt('IC_mixed_'+self.data_dict["cross_sec"]+'_'+'P'+'.txt', DataOut) 

        state_t= np.loadtxt(self.test_case+self.data_dict["variable_3"]+'_ETAT'+'_'+self.data_dict["cross_sec"]+'.coupe', unpack=True, usecols=[1])
        upper_bound = state_t+1/np.sqrt(1-self.data_dict["confidence"])*self.pcm_var_t
        lower_bound = state_t-1/np.sqrt(1-self.data_dict["confidence"])*self.pcm_var_t
        DataOut = np.column_stack((self.dx,state_t, upper_bound, lower_bound))
        np.savetxt('IC_mixed_'+self.data_dict["cross_sec"]+'_'+'T'+'.txt', DataOut) 


class plot_case(run_Variance_mixed):

    def __init__(self,cross_sec,label,i):

        super().__init__(cross_sec,label)

        self.i=i
        self.cross_sec=cross_sec
        self.label=label

    def affichage_graph(self):
        choix=self.data_dict["figure"][self.i]
        match choix:
            case "1":
                fig, ax = plt.subplots()
                ax.plot(self.dx, self.var_vx, self.dx, self.pcm_var_vx, 'r8', linewidth=2)
                ax.legend(['Taylor', 'PCM'])
                plt.title('Horizontal velocity component standard deviation')
                plt.xlabel(self.data_dict["label"])
                plt.ylabel("Horizontal velocity component standard deviation")
                plt.show()
                plt.close()

            case "2":
                fig, ax = plt.subplots()
                ax.plot(self.dx, self.var_vy, self.dx, self.pcm_var_vy, 'r8', linewidth=2)
                ax.legend(['Taylor', 'PCM'])
                plt.title('Vertical velocity component standard deviation')
                plt.xlabel(self.data_dict["label"])
                plt.ylabel("Vertical velocity component standard deviation")
                plt.show()
                plt.close()
            case "3":
                fig, ax = plt.subplots()
                ax.plot(self.dx, self.var_t, self.dx,self.pcm_var_t, 'r8', linewidth=2)
                ax.legend(['Taylor', 'PCM'])
                plt.title('Temperature standard deviation')
                plt.xlabel(self.data_dict["label"])
                plt.ylabel("Temperature standard deviation")
                plt.show()
                plt.close()
            case "4":

                fig, ax = plt.subplots()
                ax.plot(self.dx, self.var_p, self.dx, self.pcm_var_p, 'r8', linewidth=2)
                ax.legend(['Taylor', 'PCM'])
                plt.title('Pressure standard deviation')
                plt.xlabel(self.data_dict["label"])
                plt.ylabel("Pressure standard deviation")
         
                plt.show()
                plt.close()

            case "5":
                vx_x = np.loadtxt(run_dir+"/"+"IC_mixed_X_C_Vx.txt")
                vy_x = np.loadtxt(run_dir+"/"+"IC_mixed_X_C_Vy.txt")
                x = vx_x[:, 0]
                vx = vx_x[:, 1]
                vy = vy_x[:, 1]
                ex = np.abs(vx_x[:, 2])
                ey = np.abs(vy_x[:, 2])
                plt.figure(figsize=(7, 5))
                plt.errorbar(x, vx, yerr=ex, fmt='o', color='red', capsize=2, markersize=6,label='_nolegend_')
                plt.errorbar(x, vy, yerr=ey, fmt='o', color='blue',  capsize=2,markersize=6,label='_nolegend_')
                plt.plot(x, vx, color='red', label='x-component')
                plt.plot(x, vy, color='blue', label='y-component')

                plt.xticks(np.arange(0, np.max(x) + 0.0001, 0.005))
                plt.xlabel("x [m]")
                plt.ylabel("Velocity [m/s]")
                plt.legend(loc='upper center', fontsize=12, borderpad=1.5, labelspacing=1.5)
                plt.tight_layout()
                plt.show()
                plt.close()
            case "6":
                vx_y = np.loadtxt(run_dir+"/"+"IC_mixed_Y_C_Vx.txt")
                vy_y = np.loadtxt(run_dir+"/"+"IC_mixed_Y_C_Vy.txt")
                y = vx_y[:, 0]
                vx = vx_y[:, 1]
                ex = np.abs(vx_y[:, 2])
                vy = vy_y[:, 1]
                ey = np.abs(vy_y[:, 2])
                plt.figure(figsize=(7,5))
                # tracé avec barres d'erreur
                plt.errorbar(y, vx, yerr=ex, fmt='o', color='red', label='_nolegend_')
                plt.errorbar(y, vy, yerr=ey, fmt='o', color='blue', label='_nolegend_')
                # tracé des courbes
                plt.plot(y, vx, color='red', label='x-component')
                plt.plot(y, vy, color='blue', label='y-component')
                #plt.xticks(np.arange(0, np.max(y) + 0.0001, 0.002))
                plt.xlim(0, 0.01)
                plt.ylim(-0.08, 0.04)

                plt.xlabel("y [m]")
                plt.ylabel("Velocity [m/s]")
                plt.legend(loc='upper center', fontsize=12, borderpad=1.5, labelspacing=1.5)

                plt.tight_layout()
                plt.show()  # Affiche le graphique
                plt.savefig("IC_mixed_Y_C_V.png", dpi=300)
                plt.close()