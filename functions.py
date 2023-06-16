def E(z):
    
    O_M, O_L, O_K = 0.272, 0.728, 0.000
    
    return np.sqrt( O_M*(1+z)**3 + O_K*(1+z)**2 + O_L )


def R500(MM500, zz):
    
    rs = np.zeros(len(zz)) 
    
    # in kpc/h if M500 in M_sol/h
    
    for i in range(0, len(MM500)):
    
    	rs[i] = 557 * (MM500[i]/10**14)**(1/3) * E(zz[i])**(-2/3) * (1+zz[i])
    
    return rs


def T_X(MM500, zz):

    ts = np.zeros(len(zz))
    
    # in keV if M500 in M_sol/h
    
    for i in range(0, len(MM500)):
    
    	ts[i] = 5 * (MM500[i]/3.0/10**14)**(0.65) * E(zz[i])**(0.65) #/ np.sqrt(1+zz[i])
    
    return ts


def L_X(MM500, zz, hh):

    ls = np.zeros(len(zz))
    
    # in 10^44 erg/s if M500 in M_sol/h
    
    for i in range(0, len(MM500)):
    
    	ls[i] = 1.056 * hh**(-2) * (MM500[i]/3.9/10**14)**1.61 * E(zz[i])**(1.85) * (1+zz[i])**2
    
    return ls
