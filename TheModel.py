import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def m_inf(v):
    return 1.0/(1.0+np.exp((-24.0-v)/12.0))

def n_inf(v):
    return 1.0/(1.0+np.exp((-19.0-v)/18.0))

def h(n):
    return 1.1 - 1.0/(1.0 + np.exp(-8.0*(n-0.4)))


def model(z,t,Kb):

    #-------------------------------------------------------------
    Cm     = 1.0
    tau_n  = 0.25
    g_Cl   = 7.5  
    g_Na   = 40.0 
    g_K    = 22.0
    g_Nal  = 0.02
    g_Kl   = 0.12
    #-------------------------------------------------------------
    w_i     = 2160. 
    w_o     = 720.
    
    gamma   = 0.04
    
    rho     = 250.
    beta    = w_i/w_o 

    epsilon = 0.01
    K_bath  = Kb 

    Na_i0  = 16.0
    Na_o0  = 138.0 
    K_i0   = 140.0
    K_o0   = 4.80
    Cl_o0  = 112.0
    Cl_i0  = 5.0
    #-------------------------------------------------------------

    V    = z[0]
    n    = z[1]
    DK_i = z[2]
    Kg   = z[3]

    DNa_i = -DK_i
    DNa_o = -beta*DNa_i
    DK_o  = -beta*DK_i
    K_i  = K_i0 +DK_i                                            
    Na_i = Na_i0+DNa_i 
    Na_o = Na_o0+DNa_o 
    K_o  = K_o0 +DK_o+Kg  
    I_Na  = (g_Nal+g_Na*m_inf(V)*h(n))*(V-26.64*np.log(Na_o/Na_i))
    I_K   = (g_Kl+g_K*n)*(V-26.64*np.log(K_o/K_i))
    I_Cl  = g_Cl*(V+26.64*np.log(Cl_o0/Cl_i0))
    I_pump= rho*(1.0/(1.0+np.exp((21.0-Na_i)/2.0)))*(1.0/(1.0+np.exp((5.5-K_o))))

    dV    = (-1.0/Cm)*(I_Cl+I_Na+I_K+I_pump)  
    dn    = (n_inf(V)-n)/tau_n
    dKi   =  -(gamma/w_i)*(I_K-2.0*I_pump)
    dKg   =  epsilon*(K_bath-K_o)

    dz = [dV,dn,dKi,dKg]

    return dz




# initial condition
z0 = [-78.0, n_inf(-78.0), -0.6 ,0.8]



Np=10000
t = np.linspace(0,Np, int(Np/0.01))

KB= [4.8, 7.5, 9.5, 12.5, 17.0, 17.5, 20.0]

for kb in KB:
    print('[K]_bath =', kb)
    z = odeint(model,z0,t,args=(kb,))


    fig = plt.figure(figsize=(10, 8))
    plt.rc('xtick',labelsize=16)
    plt.rc('ytick',labelsize=16)
    plt.ylabel('Vm (mV)')
    plt.ylim(-90, 30)
    plt.xlabel('Time (ms)')
    plt.plot(t, z[:, 0])
    plt.show()

        
        





