from trustutils import run
import numpy as np

run_dir=run.BUILD_DIRECTORY


class IC_mixed_Poly_Chaos:
    data_dict = {

                    "test_case_v" : "Velocity/Lid_driven_cavity_Poly_Chaos_",
                    "test_case_mu" : "Mu/Lid_driven_cavity_Poly_Chaos_",
                    "variable" : "U",
                    "direction" : ["x","y"],
                    "dir":[1,2],
                    "cross_sec":["H","V"],
                    "confidence" : 0.95,   #level of confidence for the interval
                    "std_dev_v" : 0.1,     #std deviation of the input
                    "std_dev_mu" : 0.001   #std deviation of the input
                }
########## Sensitivity section  Velocity ##########

    for i in range(len(data_dict["dir"])):

        for j in range (len(data_dict["cross_sec"])):

            dx= np.loadtxt(run_dir+"/"+data_dict["test_case_v"]+data_dict["variable"]+'_'+data_dict["cross_sec"][i]+'.coupe', unpack=True, usecols=[0])
            state_v= np.loadtxt(run_dir+"/"+data_dict["test_case_v"]+data_dict["variable"]+'_'+data_dict["cross_sec"] [i]+'.coupe', unpack=True, usecols=[data_dict["dir"][j]])
            sens_v= np.loadtxt(run_dir+"/"+data_dict["test_case_v"]+data_dict["variable"]+'_'+data_dict["cross_sec"] [i]+'_SENS.coupe', unpack=True, usecols=[data_dict["dir"][j]])
            state_mu= np.loadtxt(run_dir+"/"+data_dict["test_case_mu"]+data_dict["variable"]+'_'+data_dict["cross_sec"] [i]+'.coupe', unpack=True, usecols=[data_dict["dir"][j]])
            sens_mu= np.loadtxt(run_dir+"/"+data_dict["test_case_mu"]+data_dict["variable"]+'_'+data_dict["cross_sec"] [i]+'_SENS.coupe', unpack=True, usecols=[data_dict["dir"][j]])

            variance =np.abs(sens_v*data_dict["std_dev_v"] + sens_mu*data_dict["std_dev_mu"])
            upper_bound = state_v+1/np.sqrt(1-data_dict["confidence"])*variance
            lower_bound = state_v-1/np.sqrt(1-data_dict["confidence"])*variance
            DataOut = np.column_stack((dx,state_v, upper_bound, lower_bound))
            np.savetxt('IC_mixed_'+data_dict["cross_sec"][i]+'_'+data_dict["direction"] [j]+'.txt', DataOut)
            np.savetxt(run_dir+"/"+'IC_mixed_'+data_dict["cross_sec"][i]+'_'+data_dict["direction"] [j]+'.txt', DataOut)