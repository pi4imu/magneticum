# "clusters" and "binned_clusters" are external lists

# returns list of photons inside chosen radius

def extract_photons_from_cluster(current_cluster_num, r, draw=True):

    # there are several cases of SAME ihal for DIFFERENT cluster numbers
    
    current_cluster = clusters.loc[current_cluster_num]
    
    RA_c = current_cluster["x_pix"]*30-5
    DEC_c = current_cluster["y_pix"]*30-5
    R_vir = current_cluster["Rrel"]*30
    R_500 = current_cluster["R500"]*0.704  # kpc
    ztrue = current_cluster["z_true"]
    
    D_A = FlatLambdaCDM(H0=100*0.704, Om0=0.272).angular_diameter_distance(ztrue)*1000 # kpc
    R_500_rescaled = R_500/D_A.value*180/np.pi
    
    #print(R_vir, R_500_rescaled)
    
    snap_id_str = binned_clusters[current_cluster_num][1]
 #   hdul = fits.open("../data/eROSITA_30.0x30.0/Phox/phlist_"+snap_id_str+".fits")
 #   SLICE = pd.DataFrame(hdul[2].data[:])
 #   display(SLICE)
 #   hdul.close()
    
#def extract_photons(ra_cl, dec_cl, R_cl, snapID):
    
    t = Table.read("../data/eROSITA_30.0x30.0/Phox/phlist_"+snap_id_str+".fits", hdu=2)
    
    SLICE = t.to_pandas()
    
    if r == 'Rvir':
        R = R_vir

    elif r == 'R500':
        R = R_500_rescaled
    
    SLICE["check"]=np.where((SLICE["RA"]-RA_c)**2+(SLICE["DEC"]-DEC_c)**2 <= R**2, True, False)
     
    df = SLICE[SLICE['check'] == True]
    
    dddfff = df.drop("check", axis=1)
    
    if draw:
    
        #plt.figure(figsize=(6,5))
        plt.scatter(dddfff["RA"], dddfff["DEC"], c=dddfff["ENERGY"], cmap='viridis', s=1)
        plt.colorbar()
        plt.title('#'+str(current_cluster_num), fontsize=15)
        plt.tight_layout()
        #plt.show()
    
    return dddfff
    
    
# "erosita_binning" is external list
    
def create_spectrum_old(list_of_photons, create_textfile=True, create_pha=True):

    N_channels = 1024

    #dummyrsp = np.linspace(0.1, 50.0, N_channels+1)
    #dummyrsp = np.logspace(np.log10(0.1), np.log10(12.0), N_channels+1)
    dummyrsp = np.append(erosita_binning, [12.0])

    counts, energies =  np.histogram(list_of_photons["ENERGY"], bins = dummyrsp)
    
    datfile = pd.DataFrame({"COUNTS":counts})
    
    #display(datfile)
    
    if create_textfile:
        datfile.to_csv(f'../data/eROSITA_30.0x30.0/photons/photons_'+str(current_cluster_num)+'.dat', index=True, header=False, sep=' ')

    prihdr = fits.Header()

    prihdr['COMMENT'] = "FITS (Flexible Image Transport System) format is defined in 'Astronomy  and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H"

    prihdu = fits.PrimaryHDU(header=prihdr)
    
    col0 = fits.Column(name='CHANNEL',  format='I', array = np.linspace(1,1024,1024))  #  datfile["ENERGIES"])
    col3 = fits.Column(name='COUNTS',   format='I', array = datfile["COUNTS"], unit='COUNTS')

    cols_spec = fits.ColDefs([col0, col3])

    tbhdu_spec = fits.BinTableHDU.from_columns(cols_spec)
    
    thdulist1 = fits.HDUList([prihdu, tbhdu_spec])
    hdr1 = thdulist1[1].header
    
    hdr1["EXTNAME"]  = ("SPECTRUM", 'the name (i.e. type) of the extension')
    hdr1["TELESCOP"] = ("eROSITA", 'the "telescope" (i.e. mission/satellite name)')
    hdr1["INSTRUME"] = ("TM1", 'the instrument/detector')
    hdr1["FILTER"]   = ("NONE", 'the instrument filter in use (if any) ')
    hdr1["EXPOSURE"] = ("10000", 'the integration time (in seconds) for the PHA data (assumed to be corrected for deadtime, data drop-outs etc. )')
    hdr1["BACKFILE"] = ("NONE", 'the name of the corresponding background file (if any)')
    hdr1["BACKSCAL"] = ("1", 'the background scaling factor (unless included as a column)')
    hdr1["CORRFILE"] = ("NONE", 'the name of the corresponding correction file (if any)')
    hdr1["CORRSCAL"] = ("1", 'the correction scaling factor')
    hdr1["RESPFILE"] = ("NONE", 'the name of the corresponding (default) redistribution matrix file (RMF; see George et al. 1992a)')
    hdr1["ANCRFILE"] = ("NONE", 'the name of the corresponding (default) ancillary response file (ARF; see George et al. 1992a)')
    hdr1["AREASCAL"] = ("1", 'the area scaling factor (unless included as a column)')
    hdr1["HDUCLASS"] = ("OGIP", 'should contain the string "OGIP" to indicate that this is an OGIP style file')
    hdr1["HDUCLAS1"] = ("SPECTRUM", 'should contain the string "SPECTRUM" to indicate this is a spectrum')
    hdr1["HDUVERS"]  = ("1.1.0", 'the version number of the format (this document describes version 1.2.1)')
    hdr1["POISSERR"] = (True, 'whether Poissonian errors are appropriate to the data (see below)')
    hdr1["CHANTYPE"] = ("PI", 'whether the channels used in the file have been corrected in anyway (see below)')
    hdr1["DETCHANS"] = ("1024", 'the total number of detector channels available')
    
    if create_pha:
        thdulist1.writeto('../data/eROSITA_30.0x30.0/photons/spectrum'+str(current_cluster_num)+'.pha', overwrite=True)
    
    print("Spectrum succesfully created")
        
    return None
    

def create_spectrum_and_fit_it(current_cluster_num, list_of_photons, REDSHIFT, borders, Xplot=False, plot=True, also_plot_model=False):

    x.Xset.chatter = 0
    
    erosita_binning = fits.open('../erosita/erosita_pirmf_v20210719.rmf')[1].data["E_MIN"]
    #print(erosita_binning)
    
    N_channels = 1024

    # there is no big difference which binning to choose

    #dummyrsp = np.linspace(0.1, 12.0, N_channels+1)
    #dummyrsp = np.logspace(np.log10(0.1), np.log10(12.0), N_channels+1)
    dummyrsp = np.append(erosita_binning, [12.0])

    photons, energies = np.histogram(list_of_photons["ENERGY"], bins = dummyrsp)

    model_input = [a/10000/1000 for a in photons]

    # 10000 (s) is exposition time and 1000 (cm2) is nominal area
    
    if Xplot:
        x.Plot.device = "/xs"
    else:
        x.Plot.device = '/null'
        
    x.Plot.xAxis = "keV"
        
    x.AllModels.clear()

    def myModel(engs, params, flux):
        for i in range(len(engs)-1):
            if engs[i]>0.1 and engs[i]<12.0:
                val = np.interp(engs[i], dummyrsp[1:], model_input)
                #print(i, engs[i], val)
                flux[i] = val
            else:
                flux[i] = 0

    myModelParInfo = (f"par1 Number {current_cluster_num} 1 1 1 1 0.001",)

    x.AllModels.addPyMod(myModel, myModelParInfo, 'add')
    
    x.AllData.dummyrsp(lowE=0.1, highE=12.0, nBins=1024)
    #x.AllData.removeDummyrsp()
    
    mmmm = x.Model("myModel")
    
    if plot and also_plot_model:
    
        plt.subplot(121)
    
        x.Plot("model")
        xVals = x.Plot.x()[1:]
        modVals = x.Plot.model()[1:]

        plt.plot(xVals, modVals, label="My model")
        plt.xscale('log')
        plt.yscale('log')
        #plt.legend()
        
        plt.xlabel(x.Plot.labels()[0])
        plt.ylabel(x.Plot.labels()[1])
        plt.title(x.Plot.labels()[2])
               
        
    x.AllData.clear()

    fs = x.FakeitSettings(response = '../erosita/erosita_pirmf_v20210719.rmf', 
                               arf = '../erosita/tm1_arf_open_000101v02.fits', 
                        background = '', 
                          exposure = 10000, 
                        correction = '', 
                      backExposure = '', 
                          fileName = 'fakeit.pha')
    x.AllData.fakeit(nSpectra = 1, 
                     settings = fs, 
                   applyStats = True,
                   filePrefix = "",
                      noWrite = True)
                      
    if plot:
        
        if also_plot_model:
             plt.subplot(122)
        
        x.Plot("ldata")        
        
        xVals = x.Plot.x()
        xErrors = x.Plot.xErr()
        yVals = x.Plot.y()
        yErrors = x.Plot.yErr()
        modVals = x.Plot.model()
        
        every = 1

        plt.errorbar(xVals[::every], yVals[::every], yerr=yErrors[::every], xerr=xErrors[::every], linewidth=0, elinewidth=1, label = "All data")
        
        plt.xscale('log')
        plt.yscale('log')
        plt.legend(loc=1) 

        plt.xlabel(x.Plot.labels()[0])
        plt.ylabel(x.Plot.labels()[1])
        #plt.title(x.Plot.labels()[2])
       
    x.AllData.ignore(f"**-{borders[0]} {borders[1]}-**")
    #x.AllData.notice("all")

                        
    #x.AllModels.clear()

    mod = x.Model('phabs*apec')
    mod(1).values = 0.01      # n_H
    mod(1).frozen = True
    mod(3).frozen = True     # abundance
    mod(3).values = 0.3
    mod(4).values = f"{REDSHIFT}"  # of cluster, not of photon list

    mod.show()
    
    x.Fit.renorm('auto')
    x.Fit.nIterations = 100
    #x.Fit.query = 'yes'
    x.Fit.weight = 'standard'
    x.Fit.statMethod = 'cstat'
    
    x.Fit.perform()
        
    x.Xset.parallel.error = 4
    x.Fit.error('2')
    
    T_spec = mod(2).values[0]
    T_spec_left = mod(2).error[0]
    T_spec_right = mod(2).error[1]
    
    x.AllModels.calcLumin(f"0.1 10.0 {REDSHIFT}")
    luminosity = x.AllData(1).lumin
    
    # average energy:
    
    s_i = x.AllData(1).values
    ens = x.AllData(1).energies
    E_i = [(e[0]+e[1])/2 for e in ens]
    
    av_en = np.dot(E_i, s_i)/np.sum(s_i)
    
    if plot:
    
        if also_plot_model:
             plt.subplot(122)           
    
        plt.axvline(borders[0], linestyle = '--', color='black')
        plt.axvline(borders[1], linestyle = '--', color='black')
    
        x.Plot("ldata")

        xVals = x.Plot.x()
        xErrors = x.Plot.xErr()
        yVals = x.Plot.y()
        yErrors = x.Plot.yErr()
        modVals = x.Plot.model()

        plt.errorbar(xVals[::every], yVals[::every], yerr=yErrors[::every], xerr=xErrors[::every], linewidth=0, elinewidth=1, color='b', label = "Fit data")
        plt.plot(xVals, modVals, linewidth=2, color='red', label="Best-fit")
        plt.legend(loc=1, framealpha=1)
        
        plt.title(f"#{current_cluster_num}: "+"$T_{spec}="+f"{T_spec:.2f}"+f"^{{+{(T_spec-T_spec_left):.2f}}}"+f"_{{-{(T_spec_right-T_spec):.2f}}}$"+f" ($\\chi^2_{{red.}}=$ {x.Fit.statistic/x.Fit.dof:.3f})", fontsize=15)
               
        if also_plot_model:
            
            plt.subplot(121)
            
            x.Plot("model")
            xVals = x.Plot.x()
            modVals = x.Plot.model()

            plt.plot(xVals, modVals, label="Best-fit")
            plt.legend()
        
            #plt.show()  
        
    x.Xset.chatter = 10
   
    return (T_spec, T_spec_left, T_spec_right), luminosity, av_en
    

def calculate_all_and_average_it(N_usr, CLUSTERS_LIST):

    temp_usr1 = {}
    lumin_usr1 = {}
    aven_usr1 = {}
    
    for cl_num in CLUSTERS_LIST.index[:]:
    
        mean_temp = 0
        mean_lum = 0
        mean_aven = 0
        
        cl_red = CLUSTERS_LIST.loc[cl_num]["z_true"]
        cl_T500 = CLUSTERS_LIST.loc[cl_num]["T500"]
        cl_lum = CLUSTERS_LIST.loc[cl_num]["Lx500"]
        
        pho_list = extract_photons_from_cluster(cl_num, r = 'Rvir', draw=False)
        
        print(" |", cl_num,": ", end="")
        
        temps = np.zeros(N_usr)
        lumins = np.zeros(N_usr)
        avens = np.zeros(N_usr)
    
        for i in range(N_usr):
	    
            Ts = create_spectrum_and_fit_it(cl_num, pho_list, cl_red, borders=[0.4, 7.0], 
	                                    Xplot=False, plot=False, also_plot_model=False)
    
            temps[i] = Ts[0][0]
            lumins[i] = Ts[1][0]
            avens[i] = Ts[2]
	    
            print(i+1, end="")
	    
        mean_temp = np.mean(temps)
        mean_lum = np.mean(lumins)
        mean_aven = np.mean(avens)
        
        err_temp = np.std(temps)
        err_lum = np.std(lumins)
        err_aven = np.std(avens)
	    
        temp_usr1[cl_num] = [cl_T500, mean_temp, err_temp]
        lumin_usr1[cl_num] = [cl_lum, mean_lum, err_lum]
        aven_usr1[cl_num] = [mean_aven, err_aven]
        
    return temp_usr1, lumin_usr1, aven_usr1








    
