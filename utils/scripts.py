# "clusters" and "binned_clusters" are external lists

# returns list of photons inside chosen radius

def extract_photons_from_cluster(current_cluster_number, r, draw=True, draw_new=True):

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
    
    if r == 'Rvir':
        R = R_vir

    elif r == 'R500':
        R = R_500_rescaled
    
    SLICE["check"]=np.where((SLICE["RA"]-RA_c)**2+(SLICE["DEC"]-DEC_c)**2 <= R**2, True, False)
     
    df = SLICE[SLICE['check'] == True]
    
    dddfff = df.drop("check", axis=1)
    
    if draw:
    
        ang_res = 5
    
        #plt.figure(figsize=(6,5))
        if not draw_new:
            plt.scatter(dddfff["RA"], dddfff["DEC"], c=dddfff["ENERGY"], cmap='viridis', s=0.001)
        else:
            SLICE1["whattodraw1"] = np.where( (np.abs(SLICE1["RA"]-RA_c) < 1.1*R_vir) & (np.abs(SLICE1["DEC"]-DEC_c) < 1.1*R_vir), True, False)
            whattodraw = SLICE1[SLICE1['whattodraw1'] == True]
            whattodraw = whattodraw.drop("whattodraw1", axis=1)
            #plt.scatter(whattodraw["RA"], whattodraw["DEC"], c=whattodraw["ENERGY"], cmap='viridis', s=0.001) #norm=matplotlib.colors.LogNorm())
            plt.hist2d(whattodraw["RA"], whattodraw["DEC"], 
                       bins=int(2*R_500_rescaled*3600/ang_res),
                       norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1))
        
        plt.gca().add_patch(plt.Circle((RA_c, DEC_c), R_vir, color='dodgerblue', linestyle="--", lw=3, fill = False))
        plt.gca().add_patch(plt.Circle((RA_c, DEC_c), R_500_rescaled, color='orangered', linestyle="--", lw=3, fill = False))
        plt.xlim(RA_c-1.1*R_vir, RA_c+1.1*R_vir)
        plt.ylim(DEC_c-1.1*R_vir, DEC_c+1.1*R_vir)
        plt.xlabel("RA")
        plt.ylabel("DEC")
        plt.colorbar(label=f"Number of photons in {ang_res}''$\\times${ang_res}'' bin")
        plt.title('#'+str(current_cluster_number), fontsize=15)
        #plt.tight_layout()
        
        handles, labels = plt.gca().get_legend_handles_labels()
        l1 = Line2D([], [], label="$R_{vir}$", color='dodgerblue', linestyle='--', linewidth=3)
        l2 = Line2D([], [], label="$R_{500}$", color='orangered', linestyle='--', linewidth=3)
        handles.extend([l1, l2])
        plt.legend(handles=handles, loc=3)
        #plt.show()
    
    return dddfff
    
    

def create_spectrum_and_fit_it(current_cluster_num, borders, BACKGROUND=False, inside_radius="R500", Xplot=False, plot=True, draw_only=False, save_atable_model=False):

    x.Xset.chatter = 0
    
    # binning for proper imaging of model (doesn't affect fitting)
    erosita_binning = fits.open('../erosita/erosita_pirmf_v20210719.rmf')[1].data["E_MIN"]
  
    N_channels = 1024
    
    # there is no big difference which binning to choose
  
    #dummyrsp = np.linspace(0.1, 12.0, N_channels+1)
    #dummyrsp = np.logspace(np.log10(0.1), np.log10(12.0), N_channels+1)
    dummyrsp = np.append(erosita_binning, [12.0])
    
    list_of_photons = extract_photons_from_cluster(current_cluster_num, r = inside_radius, draw=False)
    REDSHIFT = clusters.loc[current_cluster_num]["z_true"]
    D_A = FlatLambdaCDM(H0=100*0.704, Om0=0.272).angular_diameter_distance(REDSHIFT)*1000 # kpc
    R_500_rescaled = clusters.loc[current_cluster_num]["R500"]*0.704/D_A.value*180/np.pi

    # spectra from photons
    photons, energies = np.histogram(list_of_photons["ENERGY"], bins = dummyrsp)

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
                val = np.interp(engs[i], dummyrsp[1:], model_input)
                #print(i, engs[i], val)
                flux[i] = val
            else:
                flux[i] = 0

    myModelParInfo = (f"par1 Number {current_cluster_num} 1 1 1 1 0.001",)

    x.AllModels.addPyMod(myModel, myModelParInfo, 'add')
    
    # dummy rsp for model in suitable form
    x.AllData.dummyrsp(lowE=0.1, highE=12.0, nBins=1024)
    
    mmmm = x.Model("myModel")
    
    if plot:
    
        model_scale = "model"

        if draw_only!='DATA':
        
            if draw_only==False:
    	        plt.subplot(1,2,1)
            x.Plot(model_scale)
            xVals_no_bkg = x.Plot.x()[1:]
            modVals_no_bkg = x.Plot.model()[1:]
            
    # writing our model to FITS-file
    
    if save_atable_model:
    	
        from xspec_table_models import XspecTableModelAdditive
        
        parameter = ('Number', [current_cluster_num], False, False)
        fits11 = XspecTableModelAdditive('atable_models/model_atable_'+str(current_cluster_num)+'.fits', 'myAtableModel', energies[1:], [parameter])
        fits11.write(0, [a*(1.6*10**(-9))/b/10000/1000 for a, b in zip(photons, np.diff(dummyrsp))], False)
        fits11.save()
	    
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
        myModel_with_bkg(3).values = np.pi*R_500_rescaled**2*3600 # area of cluster = factor before background

    
    # plot initial model on the left panel

    if plot:
    
        if draw_only!='DATA':

            if draw_only==False:
                plt.subplot(121)

            plt.plot(xVals_no_bkg, modVals_no_bkg, label="Model without background", linestyle = '--', linewidth=2)
        
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
            plt.title(x.Plot.labels()[2])
               
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
        mod(4).values = f"{REDSHIFT}"

        #mod(6).values = np.pi*R_500_rescaled**2/min2_to_deg2 # area of cluster = factor before background

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
    
    x.AllModels.calcLumin(f"0.1 10.0 {REDSHIFT}")
    luminosity = x.AllData(1).lumin
    
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
   
    return (T_spec, T_spec_left, T_spec_right), luminosity, av_en
    

def calculate_all_and_average_it(N_usr, bkg=False, write_to_file=False):

    temp_usr1 = {}
    lumin_usr1 = {}
    aven_usr1 = {}
    
    for cl_num in clusters.index[:]:
    
        mean_temp = 0
        mean_lum = 0
        mean_aven = 0
        
        cl_red = clusters.loc[cl_num]["z_true"]
        cl_T500 = clusters.loc[cl_num]["T500"]
        cl_lum = clusters.loc[cl_num]["Lx500"]
                
        print(" |", cl_num,": ", end="")
        
        temps = np.zeros(N_usr)
        lumins = np.zeros(N_usr)
        avens = np.zeros(N_usr)
    
        for i in range(N_usr):
	    
            Ts = create_spectrum_and_fit_it(cl_num, borders=[0.4, 7.0], BACKGROUND=bkg, inside_radius="R500",
	                                    Xplot=False, plot=False)
    
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
        
    if write_to_file:

        df1 = pd.DataFrame(temp_usr1.values())
        df2 = pd.DataFrame(lumin_usr1.values())
        df3 = pd.DataFrame(aven_usr1.values())
        df_all = pd.concat([df1, df2, df3], axis=1)
        df_all.columns = ['$T_{500}$', '$T_{spec}$', '$\Delta T_{spec}$',
	                  '$L_{bol}$', '$L_{fit}$', '$\Delta L_{fit}$',
	                  '$E_{av}$', '$\Delta E_{av}$']
        df_all.index = aven_usr1.keys()
        df_all.to_csv('tables/table_'+write_to_file+'.csv', sep=' ', header=False, index=True)
        
    return temp_usr1, lumin_usr1, aven_usr1


def draw_three_panels(x_array, y_array, x_label, y_label_left, y_label_right_up, y_label_right_down, clr):
    
    plt.figure(figsize=(11.5,5.5))

    plt.subplots_adjust()

    plt.subplot(121)

    xx  = [a[1] for a in x_array]#.values()]
    xxe = [a[2] for a in x_array]#.values()]
    yy  = [a[1] for a in y_array]#.values()]
    yye = [a[2] for a in y_array]#.values()]

    plt.errorbar(xx, yy, xerr=xxe, yerr=yye, linewidth=0, elinewidth=1, 
                 capsize=3, color=clr, marker='o', markersize=3)

    plt.plot([1, 9], [1, 9], color='black', linewidth=1)

    plt.xlabel(x_label, fontsize=11)
    plt.ylabel(y_label_left, fontsize=11)

    plt.xlim(1.2, 8.7)
    plt.ylim(1.2, 8.7)


    plt.subplot(222)

    plt.errorbar(xx, [YY-XX for YY, XX in zip(yy, xx)], xerr=xxe, 
                 yerr=[a+b for a, b in zip(xxe, yye)], linewidth=0, elinewidth=1, 
                 capsize=3, color=clr, marker='o', markersize=3)

    plt.axhline(0, color='black', linewidth=1)
    plt.ylabel(y_label_right_up, fontsize=11)
    
    leftb, rightb = plt.gca().get_xlim()

    plt.subplot(224)
    
    y_p = [(YY-XX)/XX for YY, XX in zip(yy, xx)]
    y_p_err = [a/b*(aa/a+bb/b) for a, aa, b, bb in zip(yy, yye, xx, xxe)]
    
    plt.errorbar(xx, y_p, xerr=xxe, yerr=y_p_err, linewidth=0, elinewidth=1, capsize=3, color=clr, marker='o', markersize=3)
    plt.scatter(xx, y_p, color=clr, marker='o', s=3)
                 
    list1, list2, list3 = zip(*sorted(zip(xx, [n-q for n, q in zip(y_p, y_p_err)], [n+q for n, q in zip(y_p, y_p_err)])))
    plt.fill_between(list1, list2, list3, interpolate=True, alpha=0.4, color=clr)

    plt.axhline(0, color='black', linewidth=1)
    plt.ylabel(y_label_right_down, fontsize=11)
    plt.xlabel(x_label, fontsize=11)
    
    plt.xlim(leftb, rightb)

    plt.show()






    
