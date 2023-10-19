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
    
    	ts[i] = 5 * (MM500[i]/3.0/10**14)**(0.65) * E(zz[i])**(0.65)  / np.sqrt(1+zz[i])
    
    return ts


def L_X(MM500, zz, hh):

    ls = np.zeros(len(zz))
    
    # in 10^44 erg/s if M500 in M_sol/h
    
    for i in range(0, len(MM500)):
    
    	ls[i] = 1.056 * hh**(-2) * (MM500[i]/3.9/10**14)**1.61 * E(zz[i])**(1.85) * (1+zz[i])**2
    
    return ls
    

def L_X_from_T(temp, abund, redshift, lumin_bol):
    
    x.Xset.chatter = 0
    
    x.Model("phabs*apec", setPars={1:0.01, 2:temp, 3:abund, 4:redshift, 5:1})#, setPars={1:temp, 2:abund, 3:redshift, 4:1})
    
    x.AllModels.calcLumin(f"0.1 10.0 {redshift}")
    L_bol = x.AllModels(1).lumin[0]
    
    x.AllModels.calcLumin(f"0.5 2.0 {redshift}")
    L_05_20 = x.AllModels(1).lumin[0]
    
    x.Xset.chatter = 10
    
    return L_05_20/L_bol*lumin_bol * (1+redshift)**2
