# "clusters" and "binned_clusters" are external lists

# returns list of photons inside chosen radius

def extract_photons_from_cluster(current_cluster_number, r, centroid=True, draw=True, redshifted_back=False):

    # there are several cases of SAME ihal for DIFFERENT cluster numbers
    # this is the reason for using cluster number as a counter
    
    current_cluster = clusters.loc[current_cluster_number]
    
    RA_c = current_cluster["x_pix"]*30-5
    DEC_c = current_cluster["y_pix"]*30-5
    R_vir = current_cluster["Rrel"]*30
    R_500 = current_cluster["R500"]*0.704  # kpc
    ztrue = current_cluster["z_true"]
    
    D_A = FlatLambdaCDM(H0=100*0.704, Om0=0.272).angular_diameter_distance(ztrue)*1000 # kpc
    R_500_rescaled = R_500/D_A.value*180/np.pi  # degrees
    
    snap_id_str = binned_clusters[current_cluster_number][1]   # id of photon list
        
    t = Table.read("../data/eROSITA_30.0x30.0/Phox/phlist_"+snap_id_str+".fits", hdu=2)
    
    SLICE = t.to_pandas()        # for photons extraction
    SLICE1 = t.to_pandas()       # for drawing
    SLICE2 = t.to_pandas()       # for center searching if it is not from table
    
    if r == 'Rvir':
        R = R_vir
    elif r == 'R500':
        R = R_500_rescaled
        print('OBSOLETE! Change R500 to numerical value!')
    else:
        R = r * R_500_rescaled
        
    if not centroid:
    
        SLICE["check"]=np.where((SLICE["RA"]-RA_c)**2+(SLICE["DEC"]-DEC_c)**2 <= R**2, True, False)
        df = SLICE[SLICE['check'] == True]
        dddfff = df.drop("check", axis=1)
        
        cntr = (RA_c, DEC_c)
    
    else:
    
        ang_res = 5
        half_size = 1.5*R_500_rescaled
        
        SLICE2["what"] = np.where( (np.abs(SLICE2["RA"]-RA_c) < half_size) & (np.abs(SLICE2["DEC"]-DEC_c) < half_size), True, False)
        whattodraw = SLICE2[SLICE2['what'] == True]
        whattodraw = whattodraw.drop("what", axis=1)
        nmhg, _, _ = np.histogram2d(whattodraw["RA"], whattodraw["DEC"], bins=int(2*half_size*3600/ang_res))
        
        psum = sum(nmhg.flatten())
        c_x, c_y = 0, 0
        
        for i in range(0, len(nmhg)):
            for j in range(0,  len(nmhg)):
                c_x = c_x + i*nmhg[i,j]/psum
                c_y = c_y + j*nmhg[i,j]/psum
        
        c_x_1 = RA_c - half_size + c_x*ang_res/3600
        c_y_1 = DEC_c - half_size + c_y*ang_res/3600
        cntr = (c_x_1, c_y_1)
        
        SLICE["check"]=np.where((SLICE["RA"]-c_x_1)**2+(SLICE["DEC"]-c_y_1)**2 <= R**2, True, False)
        df = SLICE[SLICE['check'] == True]
        dddfff = df.drop("check", axis=1)
           
    if draw:
    
        ang_res = 5
        half_size = 3*R_500_rescaled
            
        SLICE1["whattodraw1"] = np.where( (np.abs(SLICE1["RA"]-cntr[0]) < half_size) & (np.abs(SLICE1["DEC"]-cntr[1]) < half_size), True, False)
        whattodraw = SLICE1[SLICE1['whattodraw1'] == True]
        whattodraw = whattodraw.drop("whattodraw1", axis=1)
        nmhg, _, _, trtr = plt.hist2d(whattodraw["RA"], whattodraw["DEC"], 
                                   bins=int(2*half_size*3600/ang_res),
                                   norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1),
                                   range=np.array([(cntr[0]-half_size, cntr[0]+half_size),
                                                   (cntr[1]-half_size, cntr[1]+half_size)]))
        
  #      whereee = np.concatenate(np.where(nmhg == max(nmhg.flatten())))         
  #      reeeversed = [a*ang_res/60/60 for a in whereee]
  #      xeeec = plt.gca().get_xlim()[0] + reeeversed[0]
  #      yeeec = plt.gca().get_ylim()[0] + reeeversed[1]
        #print(xeeec[0])
        #print(yeeec[0])
            
        plt.scatter(RA_c, DEC_c, color='magenta', label = 'Catalogue')
        plt.scatter(cntr[0], cntr[1], color='orangered', label = 'Centroid')
            
        #plt.gca().add_patch(plt.Circle((RA_c, DEC_c), R_vir, color='dodgerblue', linestyle="--", lw=3, fill = False))
        plt.gca().add_patch(plt.Circle(cntr, R, color='orangered', linestyle="--", lw=3, fill = False))
        #plt.gca().add_patch(plt.Circle((xeeec, yeeec), R_500_rescaled, color='yellow', linestyle="--", lw=3, fill = False))
        
        plt.xlim(cntr[0]-half_size, cntr[0]+half_size)
        plt.ylim(cntr[1]-half_size, cntr[1]+half_size)
        plt.gca().set_aspect('equal', 'box')
        
        plt.xlabel("RA, deg")
        plt.ylabel("DEC, deg")
        plt.colorbar(trtr, label=f"Number of photons in {ang_res}''$\\times${ang_res}'' bin", fraction=0.046, pad=0.04)
        plt.title('#'+str(current_cluster_number), fontsize=15)
        #plt.tight_layout()
        
        handles, labels = plt.gca().get_legend_handles_labels()
      #  l1 = Line2D([], [], label="$R_{vir}$", color='dodgerblue', linestyle='--', linewidth=3)
        l2 = Line2D([], [], label="$R_{500}$", color='orangered', linestyle='--', linewidth=3)
        handles.extend([l2])
        plt.legend(handles=handles, loc=3)
        #plt.show()

    if redshifted_back:
        return dddfff.mul(1+ztrue)
    else:
        return dddfff
    

# returns Tspec, Lspec and Eav

def create_spectrum_and_fit_it(current_cluster_num, borders, BACKGROUND=False, inside_radius=1, Xplot=False, plot=True, draw_only=False, draw_and_save_atable_model=False):

    x.Xset.chatter = 0
    
    # binning for proper imaging of model (doesn't affect fitting)
    # erosita_binning = fits.open('../erosita/erosita_pirmf_v20210719.rmf')[1].data["E_MIN"]
  
    N_channels = 1024
    
    #binning = np.linspace(0.1, 12.0, N_channels+1)
    binning = np.logspace(np.log10(0.1), np.log10(12.0), N_channels+1)
    #binning = np.append(erosita_binning, [12.0])
    
    list_of_photons = extract_photons_from_cluster(current_cluster_num, r = inside_radius, draw=False)
    REDSHIFT = clusters.loc[current_cluster_num]["z_true"]
    D_A = FlatLambdaCDM(H0=100*0.704, Om0=0.272).angular_diameter_distance(REDSHIFT)*1000 # kpc
    R_500_rescaled = clusters.loc[current_cluster_num]["R500"]*0.704/D_A.value*180/np.pi

    # spectra from photons
    
    photons, energies_bins = np.histogram(list_of_photons["ENERGY"], bins = binning)
    
    energies = [(a+b)/2 for a, b in zip(energies_bins[:-1], energies_bins[1:])]
    dE = np.diff(energies_bins)
    
    model_input = [a/10000/1000 for a in photons]
    
    # 10000 (s) is exposition time and 1000 (cm2) is nominal area      
    
    if Xplot:
        x.Plot.device = "/xs"
    else:
        x.Plot.device = '/null'
        
    x.Plot.xAxis = "keV"
    
    # adding spectra as input model
        
    x.AllModels.clear()

    def myModel(engs, params, flux):
        for i in range(len(engs)-1):
            if engs[i]>0.1 and engs[i]<12.0:
                val = np.interp(engs[i], energies, model_input/dE)
                #print(i, engs[i], val)
                flux[i] = val*(engs[i+1]-engs[i])
            else:
                flux[i] = 0

    myModelParInfo = (f"par1 Number {current_cluster_num} 1 1 1 1 0.001",)

    x.AllModels.addPyMod(myModel, myModelParInfo, 'add')
    
    # dummy rsp for model in suitable form
    x.AllData.dummyrsp(lowE=0.1, highE=12.0, nBins=N_channels)
    
    mmmm = x.Model("myModel")
    
    if plot:
    
        model_scale = "model" # emodel, eemodel

        if draw_only!='DATA':
        
            if draw_only==False:
    	        plt.subplot(1,2,1)
            x.Plot(model_scale)
            xVals_no_bkg = x.Plot.x()[1:]
            modVals_no_bkg = x.Plot.model()[1:]
            
    # writing our model to FITS-file
    
    if draw_and_save_atable_model:
    	
        print("OBSOLETE. See database.")
        from xspec_table_models import XspecTableModelAdditive
        
        parameter = ('Number', [current_cluster_num], False, False)
        fits11 = XspecTableModelAdditive('model_atable_'+str(current_cluster_num)+'.fits', 'myAtableModel', np.array(energies), [parameter])
        # [erg/s/cm2/keV]
        atablemodel_input = [a*b*(1.6*10**(-9))/10000/1000/c for a, b, c in zip(photons, energies, np.diff(energies_bins))]
        fits11.write(0, atablemodel_input, False) 
        fits11.save()
        	
        x.Model("atable{model_atable_"+str(current_cluster_num)+".fits}")
        
        if plot:
            if draw_only==False:
    	        plt.subplot(1,2,1)
            x.Plot(model_scale)
            xVals_atable = x.Plot.x()[1:]
            modVals_atable = x.Plot.model()[1:]     
    		    
    # defining the model with background included:
    
    if BACKGROUND:
    
        df4 = pd.read_csv("utils/sky_bkg_full_arcmin_05cxb.xcm", header=None)[0]
        bkg_model_name = df4[0][6:]
        params={}
        for i in range(1,18):
            params[i+3] = df4[i]
        
        x.AllModels.clear()
        myModel_with_bkg = x.Model("myModel+const*"+bkg_model_name, setPars=params, sourceNum=1)
        myModel_with_bkg(2).values = 1           # norm for myModel
        myModel_with_bkg(2).frozen = True
        myModel_with_bkg(3).values = np.pi*(inside_radius*R_500_rescaled)**2*3600 # area of cluster = factor before background
    
    # plot initial model on the left panel

    if plot:
    
        if draw_only!='DATA':

            if draw_only==False:
                plt.subplot(121)

            plt.plot(xVals_no_bkg, modVals_no_bkg, label="Model without background", linestyle = '-', linewidth=2)
            
            if draw_and_save_atable_model:
                plt.plot(xVals_atable, modVals_atable, label="Model from atable", linestyle = '-', linewidth=2, color='g')
                
           # plt.plot(energies, [a/10000/1000/b for a,b in zip(photons, dE)], color='magenta')
            
            if BACKGROUND:
        
                x.Plot(model_scale)
                xVals_with_bkg = x.Plot.x()[1:]
                modVals_with_bkg = x.Plot.model()[1:]
                
                plt.plot(xVals_with_bkg, modVals_with_bkg, label="Model with background", alpha=0.5)

            plt.xscale('log')
            plt.yscale('log')
            #plt.legend()
        
            plt.xlabel(x.Plot.labels()[0])
            plt.ylabel(x.Plot.labels()[1])
            plt.title('Photons from circle with $R$ = '+str(inside_radius)+'$\cdot R_{500}$')
               
    # fakeit for input model (how erosita sees photons)    
        
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
    
    # plotting fakeit data on the right panel (or without left panel, if plotting all data):
                      
    if plot:
    
        if draw_only!='MODEL':
        
            if draw_only==False:
                plt.subplot(122)
        
            x.Plot("ldata")        
            
            xVals = x.Plot.x()
            xErrors = x.Plot.xErr()
            yVals = x.Plot.y()
            yErrors = x.Plot.yErr()
        
            every = 1 # how much points to draw (doesn't affect fit)

            plt.errorbar(xVals[::every], yVals[::every], yerr=yErrors[::every], xerr=xErrors[::every], linewidth=0, elinewidth=1, label = "All data")
        
            plt.xscale('log')
            plt.yscale('log')
            plt.legend(loc=1) 
            
            #plt.axvline(0.7, ls=':', color='g')
    
            plt.xlabel(x.Plot.labels()[0])
            plt.ylabel(x.Plot.labels()[1])
            #plt.title(x.Plot.labels()[2])
        
    # energy band for fitting:
          
    x.AllData.ignore(f"**-{borders[0]} {borders[1]}-**")
    #x.AllData.notice("all")
                     
    # defining the model for fitting                    
    
    x.AllModels.clear()
    
    if not BACKGROUND:
    
        mod = x.Model('phabs*apec')
        mod(1).values = 0.01      # n_H
        mod(1).frozen = True
        mod(3).frozen = True     # abundance
        mod(3).values = 0.3
        mod(4).values = f"{REDSHIFT}"  # of cluster, not of photon list

        mod.show()
        
    else:
    
        mod = x.Model('phabs*apec+const*'+bkg_model_name)

        mod(1).values = 0.01
        mod(1).frozen = True
        mod(3).frozen = True
        mod(3).values = 0.3
        mod(4).values = 0 #f"{REDSHIFT}"  #because we already redhifted energies of photons

        #mod(6).values = np.pi*(inside_radius*R_500_rescaled)**2/min2_to_deg2 # area of cluster = factor before background

        for i in range(7, 24):
            mod(i).frozen = True

        for i in range(len(params)):
            mod(i+7).values = list(params.values())[i]

        mod.show()
    
    x.Fit.renorm('auto')
    x.Fit.nIterations = 100
    x.Fit.query = 'yes'
    x.Fit.weight = 'standard'
    x.Fit.statMethod = 'cstat'
    x.Fit.perform()
    
    # extracting parameters:
        
    x.Xset.parallel.error = 4
    x.Fit.error('2')
    
    #x.Xset.parallel.steppar = 4
    #x.Fit.steppar("2 delta 0.1 5 5 delta 0.1 5")
    
    #x.Xset.parallel.goodness = 4
    #x.Fit.goodness(100)
    
    T_spec = mod(2).values[0]
    T_spec_left = mod(2).error[0]
    T_spec_right = mod(2).error[1]
    
    abund_from_fit = mod(3).values[0]
    
    if BACKGROUND:

        area_from_fit = mod(6).values[0]
        mod(6).values = 0     # to calculate luminosity just from model, excluding background
    
    x.AllModels.calcLumin(f"0.1 10.0 {REDSHIFT}")
    luminosity = x.AllData(1).lumin
    #luminosity = [0.]

    # average energy:
    
    s_i = x.AllData(1).values
    ens = x.AllData(1).energies
    E_i = [(e[0]+e[1])/2 for e in ens]
    
    av_en = np.dot(E_i, s_i)/np.sum(s_i)
    
    stats_for_header = f" ($stat/dof=$ {x.Fit.statistic/x.Fit.dof:.3f})"
    
    # plotting best-fit model at data panel:
    
    if plot:

        if draw_only!='MODEL':
            
            if draw_only==False:
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
            
            plt.title(f"#{current_cluster_num}: "+"$T_{spec}="+f"{T_spec:.2f}"+f"^{{+{(T_spec-T_spec_left):.2f}}}"+f"_{{-{(T_spec_right-T_spec):.2f}}}$", fontsize=15)
        
        # plotting best-fit model on left panel:

        if draw_only!='DATA':
        
            if draw_only==False:
                plt.subplot(121)
            
            x.Plot(model_scale)
            xVals = x.Plot.x()
            modVals = x.Plot.model()

            plt.plot(xVals, modVals, label="Best-fit model"+stats_for_header, color='red')
            plt.legend(loc=3)
            
            plt.axvline(borders[0], linestyle = '--', color='black')
            plt.axvline(borders[1], linestyle = '--', color='black')
        
            #plt.show()  
        
    x.Xset.chatter = 10

    if not BACKGROUND:
        return (T_spec, T_spec_left, T_spec_right), luminosity, av_en, abund_from_fit
    else:
        return (T_spec, T_spec_left, T_spec_right), luminosity, av_en, area_from_fit


def calculate_all_and_average_it(N_usr, bkg=False, write_to_file=False):

    temp_usr1 = {}
    lumin_usr1 = {}
    aven_usr1 = {}
    a4th_usr1 = {}   # either abundance (if no bkg) or A_from_fit (if with bkg)
    
    for cl_num in tqdm(clusters.index[:]):
    
        mean_temp = 0
        mean_lum = 0
        mean_aven = 0
        mean_a4th = 0
        
        cl_red = clusters.loc[cl_num]["z_true"]
        cl_T500 = clusters.loc[cl_num]["T500"]
        cl_lum = clusters.loc[cl_num]["Lx500"]
        
        if bkg:
            D_A = FlatLambdaCDM(H0=100*0.704, Om0=0.272).angular_diameter_distance(cl_red)*1000 # kpc
            R_500_rescaled = clusters.loc[cl_num]["R500"]*0.704/D_A.value*180/np.pi
            cl_area = np.pi*R_500_rescaled**2*3600
           
        #print(" |", cl_num,": ", end="")
        
        temps = np.zeros(N_usr)
        lumins = np.zeros(N_usr)
        avens = np.zeros(N_usr)
        a4ths = np.zeros(N_usr)
    
        for i in tqdm(range(N_usr), leave=False):
	    
            Ts = create_spectrum_and_fit_it(cl_num, borders=[0.4, 7.0], BACKGROUND=bkg, inside_radius=1,
	                                    Xplot=False, plot=False)
    
            temps[i] = Ts[0][0]
            lumins[i] = Ts[1][0]
            avens[i] = Ts[2]
            a4ths[i] = Ts[3]
	    
            #print(i+1, end="")
	    
        mean_temp = np.mean(temps)
        mean_lum = np.mean(lumins)
        mean_aven = np.mean(avens)
        mean_a4th = np.mean(a4ths)
        
        err_temp = np.std(temps)
        err_lum = np.std(lumins)
        err_aven = np.std(avens)
        err_a4th = np.std(a4ths)
	    
        temp_usr1[cl_num] = [cl_T500, mean_temp, err_temp]
        lumin_usr1[cl_num] = [cl_lum, mean_lum, err_lum]
        aven_usr1[cl_num] = [mean_aven, err_aven]
        if not bkg:
             a4th_usr1[cl_num] = [mean_a4th, err_a4th]
        else:
             a4th_usr1[cl_num] = [cl_area, mean_a4th, err_a4th]
        
    if write_to_file:

        df1 = pd.DataFrame(temp_usr1.values())
        df2 = pd.DataFrame(lumin_usr1.values())
        df3 = pd.DataFrame(aven_usr1.values())
        df4 = pd.DataFrame(a4th_usr1.values())
        
        df_all = pd.concat([df1, df2, df3, df4], axis=1)
        
        if not bkg:
             df_all.columns = ['$T_{500}$', '$T_{spec}$', '$\Delta T_{spec}$',
	                       '$L_{bol}$', '$L_{fit}$', '$\Delta L_{fit}$',
	                       '$E_{av}$', '$\Delta E_{av}$',
	                       '$Z$', '$\Delta Z$']
        else:
             df_all.columns = ['$T_{500}$', '$T_{spec}$', '$\Delta T_{spec}$',
	                       '$L_{bol}$', '$L_{fit}$', '$\Delta L_{fit}$',
	                       '$E_{av}$', '$\Delta E_{av}$',
	                       '$A_0$', '$A_{fit}$','$\Delta A_{fit}$']
                  
        df_all.index = aven_usr1.keys()
        df_all.to_csv('tables/table_'+write_to_file+'.csv', sep=' ', header=False, index=True)
        
    return temp_usr1, lumin_usr1, aven_usr1






    
