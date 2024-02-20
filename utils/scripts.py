# "clusters" and "binned_clusters" are external lists

# returns list of photons inside chosen radius

def extract_photons_from_cluster(current_cluster_number, r, centroid=True, delete_bright_regions=True, draw=True, draw_additional=False, redshifted_back=True):

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
    SLICE3 = t.to_pandas()       # for rescaling SLICE2
    
    if r == 'Rvir':
        R = R_vir
        print('Really?')
    elif r == 'R500':
        R = R_500_rescaled
        print('OBSOLETE! Change R500 to numerical value!')
    else:
        R = r * R_500_rescaled
        
    AREA = np.pi*R**2*3600   # min2
        
    if not centroid:
    
        # taking photons from circle centered at RA_c, DEC_c
    
        SLICE["check"]=np.where((SLICE["RA"]-RA_c)**2+(SLICE["DEC"]-DEC_c)**2 <= R**2, True, False)
        df = SLICE[SLICE['check'] == True]
        dddfff = df.drop("check", axis=1)
        
        cntr = (RA_c, DEC_c) # for drawing
    
    else:
    
        #setting area and resolution for searching for center
    
        ang_res = 4
        half_size = 3*R_500_rescaled
        
        if (current_cluster_number != 13334) and (current_cluster_number != 18589):
            hs4s = half_size/3
        else:
            hs4s = half_size/1
        # making 2D histogram of size 2*half_size with center (RA_c, DEC_c) without drawing
                
        SLICE2["what"] = np.where( (np.abs(SLICE2["RA"]-RA_c) < hs4s) & (np.abs(SLICE2["DEC"]-DEC_c) < hs4s), True, False)
        whattodraw = SLICE2[SLICE2['what'] == True]
        whattodraw = whattodraw.drop("what", axis=1)
        nmhg, _, _ = np.histogram2d(whattodraw["RA"], whattodraw["DEC"], bins=int(2*hs4s*3600/ang_res))
                       
        # centroid position
        
        #nmhg = np.zeros((100,100))
        #for i in range(0, len(nmhg)):
        #    for j in range(0, len(nmhg)):
        #        if i>=70:
        #            nmhg[i,j] = 1
        #print(nmhg)
        
        psum = sum(nmhg.flatten())
        c_x, c_y = 0, 0
        
        for i in range(0, len(nmhg)):
            for j in range(0, len(nmhg)):
                c_x = c_x + i*nmhg[i,j]/psum
                c_y = c_y + j*nmhg[i,j]/psum
        
        # position of centroid in units of pixels relative to the upper left border    
        c = [int(c_x), len(nmhg)-int(c_y)]
        
        if draw_additional:
            plt.imshow(np.rot90(nmhg),
                       norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), 
                       origin='upper')
            plt.scatter(c[0], c[1], color='red')
            plt.scatter(int(len(nmhg)/2), int(len(nmhg)/2), color='magenta')
            plt.gca().set_aspect('equal', 'box')
            plt.xlim(c[0]-int(hs4s*3600/ang_res), c[0]+int(hs4s*3600/ang_res))
            plt.ylim(c[1]+int(hs4s*3600/ang_res), c[1]-int(hs4s*3600/ang_res))
            plt.gca().add_patch(plt.Circle(c, int(R*3600/ang_res), color='orangered', linestyle="--", lw=3, fill = False))
            plt.show()
                
        c_x_1 = RA_c - hs4s + c_x*ang_res/3600 
        c_y_1 = DEC_c - hs4s + c_y*ang_res/3600
        cntr = (c_x_1, c_y_1) # position in degrees

        # taking photons from circle centered at centroid
        
        SLICE["check"]=np.where((SLICE["RA"]-c_x_1)**2+(SLICE["DEC"]-c_y_1)**2 <= R**2, True, False)
        df = SLICE[SLICE['check'] == True]
        dddfff = df.drop("check", axis=1)
        
        # searching for the coordinates of maximum value
        #whereee = np.concatenate(np.where(nmhg == max(nmhg.flatten())))
        #reeeversed = [a*ang_res/60/60 for a in whereee]
        #xeeec = plt.gca().get_xlim()[0] + reeeversed[0]
        #yeeec = plt.gca().get_ylim()[0] + reeeversed[1]
        #print(xeeec)
        #print(yeeec)
        #print(np.where(nmhg == max(nmhg.flatten())))
        #print(max(nmhg.flatten()))
        #m = [int(a[0]) for a in np.where(nmhg == max(nmhg.flatten()))]
        #plt.scatter(xeeec, yeeec, color='dodgerblue', label = 'Max value')
        #plt.gca().add_patch(plt.Circle((xeeec, yeeec), R_500_rescaled, color='yellow', 
        #linestyle="--", lw=3, fill = False))                
        
        if delete_bright_regions:
        
            # recalculate nmhg relative to centroid
            SLICE3["what"] = np.where( (np.abs(SLICE3["RA"]-c_x_1) < half_size) & (np.abs(SLICE3["DEC"]-c_y_1) < half_size),
                                      True, False)
            whattodraw = SLICE3[SLICE3['what'] == True]
            whattodraw = whattodraw.drop("what", axis=1)
            nmhg, _, _ = np.histogram2d(whattodraw["RA"], whattodraw["DEC"], 
                                        bins=int(2*half_size*3600/ang_res),
                                        range=np.array([(cntr[0]-half_size, cntr[0]+half_size),
                                                        (cntr[1]-half_size, cntr[1]+half_size)]))            
            
            #fltr = Tophat2DKernel(2)
            fltr = Gaussian2DKernel(1)
            nmhg = convolve_fft(nmhg, fltr) 

            # some magic ...
            
            #shift = [int((RA_c-cntr[0]+half_size)*3600/ang_res), int((DEC_c-cntr[1]+half_size)*3600/ang_res)]
            shift = [int((half_size)*3600/ang_res), int((half_size)*3600/ang_res)]
            c1 = shift 
            
            #c2 = np.array( ) + shift #+ shift2
            
            nmhg1 = kruzhok(int(R*3600/ang_res), c1, nmhg, int(R*3600/ang_res)+1)[0]
  
            if draw_additional:
                
                plt.show()
            
                # initial plot in terms of pixels and only circle of R500
        
                plt.figure(figsize=(12,5))
                plt.subplot(121)
                plt.title("nmhg")
                plt.imshow(np.rot90(nmhg), norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), origin='upper')
                plt.scatter(c1[0], c1[1], color='orangered')
                plt.gca().add_patch(plt.Circle(c1, int(R*3600/ang_res), color='orangered', linestyle="--", lw=3, fill = False))
                plt.colorbar(fraction=0.046, pad=0.04)
                plt.subplot(122)
                plt.title("nmhg1, R = "+str(int(R*3600/ang_res))+" pixels")
                plt.imshow(np.rot90(nmhg1), norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), origin='upper')
                plt.gca().add_patch(plt.Circle((int(R*3600/ang_res), int(R*3600/ang_res)), 
                                    int(R*3600/ang_res), color='orangered', linestyle="--", lw=3, fill = False))  
                plt.colorbar(fraction=0.046, pad=0.04)
                plt.show()
            
            # searching for cutoff value in one pixel            
            
            threshold = 0.95
            
            beeens = np.geomspace(1, np.max(nmhg1.flatten()), 50)
            amount_in_bin, bin_borders = np.histogram(nmhg1.flatten(), bins = beeens)
            #bin_centers = [int(b) for b in bin_borders[:-1]]
            bin_centers = (bin_borders[:-1]+bin_borders[1:])/2
            #print(bin_centers)
            cumulative = [sum(amount_in_bin[:i]) for i in range(0, len(amount_in_bin))]
            number_cutoff = bin_centers[np.argmin(np.abs(cumulative-cumulative[-1]*threshold))]
            index_cutoff = list(bin_centers).index(number_cutoff)
            threshold_new = sum(amount_in_bin[:index_cutoff]) / cumulative[-1]
            
            #sum_amount, i = 0, 0
            #while sum_amount <= threshold * sum(amount_in_bin*bin_centers):
            #    sum_amount = sum_amount + (amount_in_bin*bin_centers)[i]
            #    #print(i, bin_centers[i], sum_amount/sum(amount_in_bin*bin_centers))
            #    i = i + 1
            #threshold = (sum_amount-(amount_in_bin*bin_centers)[i-1])/sum(amount_in_bin*bin_centers)        
            #number_cutoff = bin_centers[i-1]

            if draw_additional:
            
                plt.figure(figsize=(13,5))
                plt.subplot(121)
                plt.hist(nmhg1.flatten(), bins = beeens)
                plt.xlabel(f"Number of photons in {ang_res}''$\\times${ang_res}'' bin")
                plt.ylabel("Amount of such bins")
                plt.yscale("log")
                plt.xscale("log")
                plt.title("Flattened histogram for upper right image")
                plt.subplot(122)
                plt.scatter(bin_centers, cumulative)
                plt.xlabel(f"Number of photons in {ang_res}''$\\times${ang_res}'' bin")
                plt.ylabel("Cumulative distribution")
                plt.yscale("log")
                plt.xscale("log")         
                plt.title(f"$x$-values are added up until their sum\nis right below {threshold*100:.0f} % cutoff")    
                plt.axhline(cumulative[-1], ls='-.', color='green', label="Total sum")
                plt.axhline(cumulative[-1]*threshold, ls='-.', color='red', label=f'Total sum $\\times$ {threshold}')
                plt.axvline(number_cutoff, ls='--', color='red', 
                            label=f'Cutoff at\nnumber = {number_cutoff:.2f};\n{threshold_new*100:.2f}% reached')
                plt.legend(loc=4)
                plt.subplot(121)
                plt.axvline(number_cutoff, ls='--', color='red')
                #plt.tight_layout()
                plt.subplots_adjust()                 
                plt.show()
                        

            # making masks and applying them to images
            
            nmhg_radial = nmhg
            nmhg1_reserve = nmhg1
            
            filter_mask = nmhg <= number_cutoff
            nmhg = nmhg*filter_mask
            
            filter_mask1 = nmhg1 <= number_cutoff
            nmhg1 = nmhg1*filter_mask1
            
            #plt.imshow(np.rot90(filter_mask1))
            #plt.show()
            
            # another method of excluding bright pixels (by making rings):
            
            if False:
                
                ffltr =  nmhg1_reserve*0 #(nmhg1_reserve == 0)
                
                #plt.imshow(ffltr)
                #plt.show()
                
                r_pixels_max = int(R*3600/ang_res)
                brightness = []   #[k0[0].sum()/sum(k0[1].flatten())]
            
                ring_width = 20
                                  
                for rr in range(0, r_pixels_max+1, ring_width):
        
                    k1 = kruzhok(rr, c1, nmhg_radial, r_pixels_max+1)
                    k2 = kruzhok(rr + ring_width, c1, nmhg_radial, r_pixels_max+1)
                    ring = k2[0]-k1[0]
                    koltso = (sum(k2[1].flatten())-sum(k1[1].flatten()))
                    print('3.14 * (', (rr + ring_width)**2, '-', rr**2, ') =', koltso)
                     
                    plt.figure(figsize=(19,3))                
                    plt.subplot(151)
                    ring = np.rot90(ring)
                    plt.imshow(ring, 
                               norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), 
                               origin='upper')
                    plt.colorbar(fraction=0.046, pad=0.04)
                    
                    plt.subplot(152)
                    bbiinnss = np.geomspace(1, np.max(ring.flatten()), 50)
                    bbiinnss = sorted(bbiinnss)
                    #print(ring.flatten(), bbiinnss)
                    amount_in_bin, bin_borders = np.histogram(ring.flatten(), bins = bbiinnss)
                    #print(amount_in_bin, bin_borders)                    
                    amount_in_bin = amount_in_bin/koltso
                    bin_centers = (bin_borders[:-1]+bin_borders[1:])/2
                    #bin_centers = [int(b) for b in bin_borders[:-1]]
                    plt.plot(bin_centers, amount_in_bin)
                    
                    plt.subplot(153)
                    
                    cumulative = [sum(amount_in_bin[:i]) for i in range(0, len(amount_in_bin))]
                    number_cutoff = bin_centers[np.argmin(np.abs(cumulative-cumulative[-1]*threshold))]
                    index_cutoff = list(bin_centers).index(number_cutoff)
                    threshold_new = sum(amount_in_bin[:index_cutoff]) / cumulative[-1]
                    
                    plt.plot(bin_centers, cumulative)
                    plt.axhline(cumulative[-1], ls='-.', color='green', label="Total sum")
                    plt.axhline(cumulative[-1]*threshold, ls='-.', color='red', label=f'Total sum $\\times$ {threshold}')
                    plt.axvline(number_cutoff, ls='--', color='red', 
                                label=f'Cutoff at\nnumber = {number_cutoff:.2f};\n{threshold_new*100:.2f}% reached')
                    plt.legend(loc=4)
                    
                    plt.subplot(152)
                    plt.axvline(number_cutoff, ls='--', color='red')
                    
                    plt.subplot(154)
                    filter_mask2 = ring <= number_cutoff
                    ring = ring*filter_mask2
                    filter_mask3 = filter_mask2 == 0
                    plt.imshow(filter_mask3)
                    
                    plt.subplot(155)
                    ffltr = ffltr + ring

                    plt.imshow(ffltr, 
                               norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), 
                               origin='upper')
                    plt.colorbar(fraction=0.046, pad=0.04)
                    #plt.subplots_adjust()
                    plt.tight_layout()
                    plt.show()
            
                nmhg1 = np.rot90(ffltr, 3)
            
            
            if draw_additional:
            
                plt.figure(figsize=(11,5))
                plt.subplot(121)
                plt.title("nmhg1 (filtered)")
                plt.imshow(np.rot90(nmhg1), norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), origin='upper')
                plt.colorbar(fraction=0.046, pad=0.04)
                plt.subplot(122)
                plt.title("nmhg (filtered)")
                plt.imshow(np.rot90(nmhg), norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), origin='upper')
                plt.colorbar(fraction=0.046, pad=0.04)
                plt.show()
            
            dddfff["RA_pix"] =  (dddfff["RA"]  - cntr[0] + R)*3600/ang_res # - 1
            dddfff["DEC_pix"] = (dddfff["DEC"] - cntr[1] + R)*3600/ang_res # - 1
            
            dddfff["RA_pix"] = dddfff["RA_pix"].astype(int)
            dddfff["DEC_pix"] = dddfff["DEC_pix"].astype(int)
            
            dddfff["stay"] = filter_mask1[dddfff["RA_pix"], dddfff["DEC_pix"]]
            dddfff = dddfff[dddfff["stay"] == True]
                        
            if draw_additional:
            
                plt.figure(figsize=(11,5))
                plt.subplot(121)
                NMHG2, NMHG_X, NMHG_Y , _ = plt.hist2d(dddfff["RA"], dddfff["DEC"],
                                             bins=len(nmhg1),
                                             norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1))
                axesforsmooth = list(plt.gca().get_xlim()) + list(plt.gca().get_ylim())
                #print(axesforsmooth)
                plt.gca().set_aspect('equal', 'box')
                plt.colorbar(fraction=0.046, pad=0.04)
                plt.title('nmhg1 but in RA/DEC')
                plt.scatter(cntr[0], cntr[1], color='orangered', label = 'Centroid')
                                                
            dddfff = dddfff.drop("stay", axis=1) 
            dddfff = dddfff.drop("RA_pix", axis=1)
            dddfff = dddfff.drop("DEC_pix", axis=1)            
               
            # here we obtain brightness profile inside R_500_rescaled               
               
            if False:
            
                plt.subplot(122)
              
                r_pixels_max = int(R*3600/ang_res)
               
                #k0 = kruzhok(0, c, nmhg, r_pixels_max+1)
                brightness = []   #[k0[0].sum()/sum(k0[1].flatten())]
                brightness_filtered = []
            
                ring_width = 10
                 
                for rr in range(0, r_pixels_max+1):
        
                    k1 = kruzhok(rr, c1, nmhg, r_pixels_max+1)
                    k2 = kruzhok(rr + ring_width, c1, nmhg, r_pixels_max+1)
                    ring = k2[0]-k1[0]
                    brightness.append(ring.sum()/(sum(k2[1].flatten())-sum(k1[1].flatten())))
                    
                    k1 = kruzhok(rr, c, nmhg_radial, r_pixels_max+1)
                    k2 = kruzhok(rr + ring_width, c, nmhg_radial, r_pixels_max+1)
                    ring = k2[0]-k1[0]
                    brightness_filtered.append(ring.sum()/(sum(k2[1].flatten())-sum(k1[1].flatten())))                    
            
                r500_pix = int(R_500_rescaled*3600/ang_res)                
                plt.plot(np.linspace(0, r_pixels_max+1, r_pixels_max+1)/r500_pix, brightness)
                plt.plot(np.linspace(0, r_pixels_max+1, r_pixels_max+1)/r500_pix, brightness_filtered)
                #plt.axvline((brightness.index(max(brightness))+1)/r500_pix, ls='--', color='black')
                plt.xlabel("Radius in units of $R_{500}$")
                plt.ylabel("Brightness in relative units")
                #plt.axhline(brightness_max, ls='--', color='red', 
                #label=f'{threshold*100:.2f} % cutoff\nat brightness = {brightness_max:.2f}')
                plt.xscale("log")
                plt.yscale("log")
                #plt.legend()
                plt.subplots_adjust()                 
                plt.tight_layout()
                                
            if draw_additional and True:     
                plt.show()
                
                NMHG2 = np.rot90(NMHG2)
                plt.imshow(convolve(NMHG2, fltr), 
                           norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), 
                           origin='upper',
                           extent = axesforsmooth)
                plt.colorbar(fraction=0.046, pad=0.04)
                plt.gca().set_aspect('equal', 'box')
                plt.title('nmhg1 but in RA/DEC convolved')
                plt.show()
                   
    # this goes to final panels image centered at cntr (which is set to centroid position as default)
    
    if draw:
               
        SLICE1["todraw"] = np.where( (np.abs(SLICE1["RA"]-cntr[0]) < half_size) & (np.abs(SLICE1["DEC"]-cntr[1]) < half_size),
                                    True, False)
        SLICE1 = SLICE1[SLICE1['todraw'] == True]
        SLICE1 = SLICE1.drop("todraw", axis=1)
        
        #SLICE1 = SLICE1[SLICE1['ENERGY']>1.3]      
        #№SLICE1 = SLICE1[SLICE1['ENERGY']<2.3]
        
        if delete_bright_regions:
                
            SLICE1["RA_pix"] = ((SLICE1["RA"] - cntr[0] + half_size)*3600/ang_res).astype(int) - 1
            SLICE1["DEC_pix"] = ((SLICE1["DEC"] - cntr[1] + half_size)*3600/ang_res).astype(int) - 1
        
            #plt.hist2d(SLICE1["RA_pix"], SLICE1["DEC_pix"], 
            #           bins=len(nmhg),
            #           norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1))
            #plt.gca().set_aspect('equal', 'box')
            #plt.colorbar(fraction=0.046, pad=0.04)
            #plt.show()
        
            #print(min(SLICE1["RA_pix"]), min(SLICE1["DEC_pix"]))
        
            SLICE1['todraw'] = filter_mask[SLICE1["RA_pix"], SLICE1["DEC_pix"]]
            SLICE1 = SLICE1[SLICE1['todraw'] == True]
        
        if False:
        
            nmhg, _, _, trtr = plt.hist2d(SLICE1["RA"], SLICE1["DEC"],
                                          bins=int(2*half_size*3600/ang_res),
                                          norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1),
                                          range=np.array([(cntr[0]-half_size, cntr[0]+half_size),
                                                          (cntr[1]-half_size, cntr[1]+half_size)]))
        
        else:
        
            nmhg, nmhg_x, nmhg_y = np.histogram2d(SLICE1["RA"], SLICE1["DEC"],
                                              bins=int(2*half_size*3600/ang_res),
                                              #norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1),
                                              range=np.array([(cntr[0]-half_size, cntr[0]+half_size),
                                                              (cntr[1]-half_size, cntr[1]+half_size)]))
            
            #print([nmhg_x[0], nmhg_x[-1], nmhg_y[0], nmhg_y[-1]])                                                  
            axesforsmooth = [nmhg_x[0], nmhg_x[-1], nmhg_y[0], nmhg_y[-1]] 
           
            trtr = plt.imshow(convolve(np.rot90(nmhg), Gaussian2DKernel(1)), 
                              norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), 
                              origin='upper',
                              extent = axesforsmooth)
                        
        
        # obsolete (it was needed for estimation of correct position of circles)
        
        #m_x, m_y = RA_c - half_size + c[0]*ang_res/3600, DEC_c - half_size + c[1]*ang_res/3600
        #r_degrees = r_pixels_max*ang_res/3600
        #plt.gca().add_patch(plt.Circle((m_x, m_y), r_degrees, color='red', 
        #linestyle="-", lw=1, fill = False, label = 'Brightest'))
        #plt.gca().add_patch(plt.Rectangle((m_x-r_degrees, m_y-r_degrees), 
        #2*r_degrees, 2*r_degrees, color='red', linestyle="-", lw=1, fill = False, label = 'Brightest'))
                    
        plt.scatter(RA_c, DEC_c, color='magenta', label = 'Catalogue')
        plt.scatter(cntr[0], cntr[1], color='orangered', label = 'Centroid')
            
        #plt.gca().add_patch(plt.Circle((RA_c, DEC_c), R_vir, color='dodgerblue', linestyle="--", lw=3, fill = False))
        plt.gca().add_patch(plt.Circle(cntr, R, color='orangered', linestyle="--", lw=3, fill = False))
        
        x_s = (plt.gca().get_xlim()[1]+plt.gca().get_xlim()[0])/2
        y_s = (plt.gca().get_ylim()[1]-plt.gca().get_ylim()[0])*0.95+plt.gca().get_ylim()[0]
        y_S = (plt.gca().get_ylim()[1]-plt.gca().get_ylim()[0])*0.90+plt.gca().get_ylim()[0]
        #plt.scatter(x_s, y_s, color='red')     
        plt.plot((x_s+5/60, x_s-5/60), (y_s, y_s), color='white')
        plt.text(x_s, y_S, f'10 arcmin $\\approx$ {10/60*D_A.value*np.pi/180:.0f} kpc', 
                 color='white', ha='center', va='center')
        
        plt.xlim(cntr[0]-half_size, cntr[0]+half_size)
        plt.ylim(cntr[1]-half_size, cntr[1]+half_size)
        plt.gca().set_aspect('equal', 'box')
        
        plt.xlabel("RA, deg")
        plt.ylabel("DEC, deg")
        plt.colorbar(trtr, label=f"Number of photons in {ang_res}''$\\times${ang_res}'' bin", fraction=0.046, pad=0.04)
        plt.title(f'#{current_cluster_number}: z={ztrue:.3f}, A={AREA:.1f} min$^2$', fontsize=15)
        
        handles, labels = plt.gca().get_legend_handles_labels()
        #l1 = Line2D([], [], label="$R_{vir}$", color='dodgerblue', linestyle='--', linewidth=3)
        l2 = Line2D([], [], label=str(r)+"$\cdot R_{500}$", color='orangered', linestyle='--', linewidth=3)
        handles.extend([l2])
        plt.legend(handles=handles, loc=3)
        #plt.show()
                
              
    if redshifted_back:
        return dddfff.mul(1+ztrue)
    else:
        return dddfff

 
def kruzhok(r_pixels, mm, NMHG, d_pixels):

    kusok = np.zeros((2*d_pixels+1, 2*d_pixels+1))
        
    for i in range(mm[0]-d_pixels, mm[0]+d_pixels+1):
        for j in range(mm[1]-d_pixels, mm[1]+d_pixels+1):
            kusok[i - (mm[0]-d_pixels)][j - (mm[1]-d_pixels)] = NMHG[i][j]
        
    Y, X = np.ogrid[(mm[1]-d_pixels):(mm[1]+d_pixels+1), (mm[0]-d_pixels):(mm[0]+d_pixels+1)]
    dist_from_center = np.sqrt((X - mm[0])**2 + (Y-mm[1])**2)
    
    mask = dist_from_center <= r_pixels
        
    return mask*kusok, mask


def check_bkg():
    x.Xset.chatter = 10
    x.AllModels.show()
    x.Xset.chatter = 0
    return None    
        
# returns Tspec, Lspec and Eav

def create_spectrum_and_fit_it(current_cluster_num, borders=[0.4, 7.0], BACKGROUND=False, inside_radius=1, dbr=True, Xplot=False, plot=False, draw_only=False, draw_and_save_atable_model=False):

    x.Xset.chatter = 0
    
    # binning for proper imaging of model (doesn't affect fitting)
    # erosita_binning = fits.open('../erosita/erosita_pirmf_v20210719.rmf')[1].data["E_MIN"]
  
    N_channels = 1024
    #binning = np.linspace(0.1, 12.0, N_channels+1)
    binning = np.logspace(np.log10(0.1), np.log10(12.0), N_channels+1)
    #binning = np.append(erosita_binning, [12.0])
    
    list_of_photons = extract_photons_from_cluster(current_cluster_num, r = inside_radius, delete_bright_regions=dbr, draw=False)
    REDSHIFT = clusters.loc[current_cluster_num]["z_true"]
    D_A = FlatLambdaCDM(H0=100*0.704, Om0=0.272).angular_diameter_distance(REDSHIFT)*1000 # kpc
    R_500_rescaled = clusters.loc[current_cluster_num]["R500"]*0.704/D_A.value*180/np.pi
    AREA = np.pi*(inside_radius*R_500_rescaled)**2*3600   # min2

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
    
    #if draw_and_save_atable_model:
    	
    #    print("OBSOLETE. See database.")
    #    from xspec_table_models import XspecTableModelAdditive
        
    #    parameter = ('Number', [current_cluster_num], False, False)
    #    fits11 = XspecTableModelAdditive('model_atable_'+str(current_cluster_num)+'.fits', 
    #                                     'myAtableModel', np.array(energies), [parameter])
    #    # [erg/s/cm2/keV]
    #    atablemodel_input = [a*b*(1.6*10**(-9))/10000/1000/c for a, b, c in zip(photons, energies, np.diff(energies_bins))]
    #    fits11.write(0, atablemodel_input, False) 
    #    fits11.save()
        	
    #    x.Model("atable{model_atable_"+str(current_cluster_num)+".fits}")
    #    
    #    if plot:
    #        if draw_only==False:
    #	        plt.subplot(1,2,1)
    #        x.Plot(model_scale)
    #        xVals_atable = x.Plot.x()[1:]
    #        modVals_atable = x.Plot.model()[1:]     
    		        
    # defining the model with background included:
        
    if BACKGROUND:
    
        df4 = pd.read_csv("bkg/sky_bkg_full_arcmin_05cxb.xcm", header=None)[0]
        bkg_model_name = df4[0][6:]
        params_ph={}
        for i in range(1,18):
            params_ph[i+3] = df4[i]
        x.AllModels.clear()
        myModel_with_bkg = x.Model("myModel+const*("+bkg_model_name+")", setPars=params_ph, sourceNum=1)
        myModel_with_bkg(2).values = "1, -1.0"           # norm for myModel
        myModel_with_bkg(2).frozen = True
        myModel_with_bkg(3).values = AREA        # area of cluster = factor before background
    
        #check_bkg()
    
    # plot initial model on the left panel

    if plot:
        if draw_only!='DATA':
            if draw_only==False:
                plt.subplot(121)
            plt.plot(xVals_no_bkg, modVals_no_bkg, label="Model", linestyle = '-', linewidth=2)
            
            #if draw_and_save_atable_model:
            #    plt.plot(xVals_atable, modVals_atable, label="Model from atable", linestyle = '-', linewidth=2, color='g')    
           # plt.plot(energies, [a/10000/1000/b for a,b in zip(photons, dE)], color='magenta')
            
            if BACKGROUND:
                x.Plot(model_scale)
                xVals_with_bkg = x.Plot.x()[1:]
                modVals_with_bkg = x.Plot.model()[1:]               
                plt.plot(xVals_with_bkg, modVals_with_bkg, label="Model with background photons", alpha=0.5)
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel(x.Plot.labels()[0])
            plt.ylabel(x.Plot.labels()[1])
            plt.title('Photons from circle with $R$ = '+str(inside_radius)+'$\cdot R_{500}$')
               
    # fakeit for input model (how erosita sees photons) saved to name1
    
    rn = np.random.randint(1000)
    identificator = f'_{current_cluster_num}_{getpid()}_{rn}'
    #print(multiprocessing.current_process())
    #print(multiprocessing.Process()._identity[0])
    name1 = 'fakeit'+identificator+'.pha'
    name2 = 'pbkg'+identificator+'.pha'
    name3 = 'total'+identificator+'.pha'
        
    x.AllData.clear()
    fs = x.FakeitSettings(response = '../erosita/erosita_pirmf_v20210719.rmf', 
                               arf = '../erosita/tm1_arf_open_000101v02.fits', 
                        background = '', 
                          exposure = 10000, 
                        correction = '', 
                      backExposure = '', 
                          fileName = name1)
    x.AllData.fakeit(nSpectra = 1, 
                     settings = fs, 
                   applyStats = True,
                   filePrefix = "",
                      noWrite = False)
                      
    # particle background
    
    x.AllModels.clear()
    
    if BACKGROUND:
    
        df5 = pd.read_csv("bkg/erass_pbkg_model.xcm", header=None)[0]
        pbkg_model_name = df5[6][6:]
        params_part = {}
        for i in range(1, 17):
            params_part[i+1] = df5[i+6]
        x.AllModels.clear()
        particle_bkg_model = x.Model("const*("+pbkg_model_name+")", sourceNum=1) #setPars=params_part,
        particle_bkg_model(1).values = AREA
        for i in range(2, 18):
            particle_bkg_model(i).values = list(params_part.values())[i-2]
            particle_bkg_model(i).frozen = True              
        
        #check_bkg()

        #if plot:
        #    if draw_only!='DATA':
        #        if draw_only==False:
    	#            plt.subplot(1,2,1)
    	#            x.Plot(model_scale)
    	#            xVals_pbkg = x.Plot.x()[1:]
    	#            modVals_pbkg = np.array(x.Plot.model()[1:])
    	#            plt.plot(xVals_pbkg, modVals_pbkg, label="Model for pbkg", alpha=0.9)
        
        x.AllData.clear()
        fs = x.FakeitSettings(response = '../erosita/erosita_pirmf_v20210719.rmf', 
                                   arf = '', 
                            background = '', 
                              exposure = 10000000, 
                            correction = '', 
                          backExposure = '', 
                              fileName = name2)
        x.AllData.fakeit(nSpectra = 1, 
                         settings = fs, 
                       applyStats = True,
                       filePrefix = "",
                          noWrite = False)
        
        # clean all data and upload data for drawing only
        x.AllData.clear()	        
        x.AllData(f"1:1 {name1} 2:2 {name2}")
        s1 = x.AllData(1)
        s2 = x.AllData(2)    
    
    # plotting fakeit data on the right panel (or without left panel, if plotting only data):
    
    every = 1 # how much points to draw (doesn't affect fit)                      
    
    if plot:
        if draw_only!='MODEL':
            if draw_only==False:
                plt.subplot(122)
            x.Plot("ldata")        
            xVals = x.Plot.x()
            xErrors = x.Plot.xErr()
            yVals = x.Plot.y()
            yErrors = x.Plot.yErr()
            plt.errorbar(xVals[::every], yVals[::every], 
                         yerr=yErrors[::every], xerr=xErrors[::every], 
                         linewidth=0, elinewidth=1, label = "Photons + ph.bkg")
                
            # draw particle background separately
            
            if BACKGROUND:
            
                x.Plot.add = True
                xVals = x.Plot.x(2)
                xErrors = x.Plot.xErr(2)
                yVals = x.Plot.y(2)
                yErrors = x.Plot.yErr(2)
                plt.errorbar(xVals[::every], yVals[::every], 
                             yerr=yErrors[::every], xerr=xErrors[::every], 
                             linewidth=0, elinewidth=1, label = "Only p.bkg", color = 'darkviolet')    
        
            plt.xscale('log')
            plt.yscale('log')
            plt.legend(loc=1)   
            plt.xlabel(x.Plot.labels()[0])
            plt.ylabel(x.Plot.labels()[1])
            #plt.title(x.Plot.labels()[2])
        
    # total spectrum: photons + particles
    
    if BACKGROUND:
        
        #os.system(f"rm {name3}")
        os.system(f"mathpha {name1}+{name2} R {name3} 10000 NULL 0 chatter=0")
        #clear_output(wait=False)

        x.AllData.clear()
        #x.AllData(f"1:1 {name3} 2:2 {name2}")
        #s1 = x.AllData(1)
        #s2 = x.AllData(2)
        
        #s1.response = "../erosita/erosita_pirmf_v20210719.rmf"
        #s1.response.arf = "../erosita/tm1_arf_open_000101v02.fits"
        #s2.response = "../erosita/erosita_pirmf_v20210719.rmf"
        
        x.AllData(f"{name3}")
        s1 = x.AllData(1)
        s1.response = "../erosita/erosita_pirmf_v20210719.rmf"
        s1.response.arf = "../erosita/tm1_arf_open_000101v02.fits"
        
        s1.background = "bkg/erass_pbkg_model.pha"
        x.Xset.chatter = 10
        x.AllData.show()
        x.Xset.chatter = 0    
        s1.backScale = 1
        
        #check_bkg()   
        
        # draw total spectrum
        
        if draw_only!='MODEL':
            if draw_only==False:
                plt.subplot(122)
            x.Plot("ldata")         
            xVals = x.Plot.x()
            xErrors = x.Plot.xErr()
            yVals = x.Plot.y()
            yErrors = x.Plot.yErr()
            plt.errorbar(xVals[::every], yVals[::every], 
                         yerr=yErrors[::every], xerr=xErrors[::every], 
                         linewidth=0, elinewidth=1, label = "Total observed spectrum", color="green", alpha = 0.7)
    
    # energy band for fitting:
          
    x.AllData.ignore(f"**-{borders[0]} {borders[1]}-**")
                     
    # defining the model for fitting:                    
    
    x.AllModels.clear()
    
    if not BACKGROUND:
    
        mod = x.Model('phabs*apec')
        mod(1).values = 0.01      # n_H
        mod(1).frozen = True
        mod(3).frozen = True     # abundance
        mod(3).values = 0.3
        mod(4).values = 0 # f"{REDSHIFT}"  # of cluster, not of photon list

        mod.show()
        
    else:
        
        mod = x.Model('phabs*apec+const*'+bkg_model_name)

        mod(1).values = 0.01
        mod(1).frozen = True
        mod(3).frozen = True
        mod(3).values = 0.3
        mod(4).values = 0 #f"{REDSHIFT}"  # because we already redhifted energies of photons
        mod(6).values = AREA              # area of cluster = factor before background
        
        # photon background
                
        for i in range(7, 24):
            mod(i).frozen = True

        for i in range(len(params_ph)):
            mod(i+7).values = list(params_ph.values())[i]
        
        #mod1 = x.AllModels(1)  # data group for model
        #mod2 = x.AllModels(2)  # data group for background    
        #mod2(5).values = "0, -1.0"   # normalization in apec = 0 
        #mod2(6).values = "0, -1.0"   # const before photon bkg = 0

        #s1.multiresponse[1] = "../erosita/erosita_pirmf_v20210719.rmf"
        #s1.multiresponse[1].arf = "../erosita/tm1_arf_open_000101v02.fits"
        #s2.multiresponse[1] = "../erosita/erosita_pirmf_v20210719.rmf"
        
        #check_bkg()      
      
        #mod.show()
        
        # particle background
        #mod_pbkg = x.Model('const*('+pbkg_model_name+')', 'PBKG', 2)        
        #mod_pbkg(1).values = AREA
        #for i in range(2, 18):
        #    mod_pbkg(i).frozen = True
        #    mod_pbkg(i).values = list(params_part.values())[i-2]
        #    mod_pbkg(i).frozen = True
      
        #check_bkg()       
        
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
        
        norm_from_fit = mod(5).values[0]       
        area_from_fit = mod(6).values[0]       
        #mod(6).values = 0     # to calculate luminosity just from model, excluding background
        area_pbkg = 0 #mod_pbkg(1).values[0]
    
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
            #x.Plot.add = True #!!!!
            xVals = x.Plot.x()
            xErrors = x.Plot.xErr()
            yVals = x.Plot.y()
            yErrors = x.Plot.yErr()
            modVals = x.Plot.model()
            plt.errorbar(xVals[::every], yVals[::every], 
                         yerr=yErrors[::every], xerr=xErrors[::every], 
                         linewidth=0, elinewidth=1, color='b', label = "Data to fit", alpha=1)
            plt.plot(xVals, modVals, linewidth=2, color='red', label="Best-fit")
            
            #if BACKGROUND:
            #    xVals = x.Plot.x(2)
            #    modVals = x.Plot.model(2)
            #    plt.plot(xVals[::every], modVals[::every], linewidth=1, label = "Best-fit p.bkg", color = 'yellow', ls='-') 
                   
            plt.legend(loc=1, framealpha=1)           
            XT = [0.1, 1, 10]
            plt.gca().set_xticks(ticks=XT, labels=XT)   
            plt.title(f"#{current_cluster_num}: "+"$T_{spec}="+f"{T_spec:.2f}"+f"^{{+{(T_spec-T_spec_left):.2f}}}"+f"_{{-{(T_spec_right-T_spec):.2f}}}$", fontsize=15)
            
        # plotting best-fit model on left panel:

        if draw_only!='DATA':
            if draw_only==False:
                plt.subplot(121)
            x.Plot(model_scale)
            xVals = x.Plot.x()
            modVals = x.Plot.model()
            plt.plot(xVals, modVals, label="Best-fit model"+stats_for_header, color='red')
            plt.legend(loc=1)
            plt.gca().set_xticks(ticks=XT, labels=XT)       
            plt.axvline(borders[0], linestyle = '--', color='black')
            plt.axvline(borders[1], linestyle = '--', color='black')
        
            #plt.show()  
    
    x.Xset.chatter = 10
    
    if not BACKGROUND:
        os.system(f"rm {name1}")
        return (T_spec, T_spec_left, T_spec_right), luminosity, av_en, abund_from_fit
    else:
        os.system(f"rm {name1}")
        os.system(f"rm {name2}")
        os.system(f"rm {name3}")  
        return (T_spec, T_spec_left, T_spec_right), luminosity, av_en, area_from_fit, norm_from_fit, area_pbkg


def average_one_cluster(cl_num, N_usr=50, bkg=False):

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
        
        mean_a4th2 = 0    
        mean_a4th3 = 0
        
        a4ths2 = np.zeros(N_usr)
        a4ths3 = np.zeros(N_usr)     
                   
    #print(" |", cl_num,": ", end="")
        
    temps = np.zeros(N_usr)
    lumins = np.zeros(N_usr)
    avens = np.zeros(N_usr)
    a4ths = np.zeros(N_usr)

    for i in tqdm(range(N_usr), leave=False, desc=str(cl_num)):
	    
        Ts = create_spectrum_and_fit_it(cl_num, borders=[0.4, 7.0], BACKGROUND=bkg, inside_radius=1, dbr=True,
	                                Xplot=False, plot=False)
        temps[i] = Ts[0][0]
        lumins[i] = Ts[1][0]
        avens[i] = Ts[2]
        a4ths[i] = Ts[3]
         
        if bkg:    
            a4ths2[i] = Ts[4]
            a4ths3[i] = Ts[5]
	    
        #print(i+1, end="")
        
    print(cl_num, "done")

    mean_temp = np.mean(temps)
    mean_lum = np.mean(lumins)
    mean_aven = np.mean(avens)
    mean_a4th = np.mean(a4ths)

    err_temp = np.std(temps)
    err_lum = np.std(lumins)
    err_aven = np.std(avens)
    err_a4th = np.std(a4ths)

    if bkg:     
	
        mean_a4th2 = np.mean(a4ths2)
        mean_a4th3 = np.mean(a4ths3)
        err_a4th2 = np.std(a4ths2)
        err_a4th3 = np.std(a4ths3)    
    
    #temp_usr1[cl_num] = [cl_T500, mean_temp, err_temp]
    #lumin_usr1[cl_num] = [cl_lum, mean_lum, err_lum]
    #aven_usr1[cl_num] = [mean_aven, err_aven]
    #if not bkg:
    #     a4th_usr1[cl_num] = [mean_a4th, err_a4th]
    #else:
    #     a4th_usr1[cl_num] = [cl_area, mean_a4th, err_a4th]

    if not bkg:
        return [cl_T500, mean_temp, err_temp], [cl_lum, mean_lum, err_lum], [mean_aven, err_aven], [mean_a4th, err_a4th]
    else:
        return [cl_T500, mean_temp, err_temp], [cl_lum, mean_lum, err_lum], [mean_aven, err_aven], [cl_area, mean_a4th, err_a4th], [mean_a4th2, err_a4th2], [cl_area, mean_a4th3, err_a4th3]
        
    
def calculate_all_and_average_it(BACKGROUND, write_to_file):

    #temp_usr1 = {}
    #lumin_usr1 = {}
    #aven_usr1 = {}
    #a4th_usr1 = {}   # either abundance (if no bkg) or A_from_fit (if with bkg)
    
    df_all = pd.DataFrame()
           
    with multiprocessing.Pool(processes=6) as pool:
        output = list(tqdm(pool.map(average_one_cluster, clusters.index[:]), total=len(clusters)))
    pool.close()
        #output = average_one_cluster(cl_num, N_usr=N_USR, bkg=BACKGROUND)

    #for cl_num in tqdm(clusters.index[:]):
    
        #df1 = pd.DataFrame(temp_usr1.values())
        #df2 = pd.DataFrame(lumin_usr1.values())
        #df3 = pd.DataFrame(aven_usr1.values())
        #df4 = pd.DataFrame(a4th_usr1.values())
        
        #df_all = pd.concat([df1, df2, df3, df4], axis=1)

    for o in output:
        
        df_add = pd.DataFrame(np.concatenate(o)).T        
        df_all = pd.concat([df_all, df_add], axis=0)

    df_all.index = [clusters.index[:]]
      
    if not BACKGROUND:
        df_all.columns = ['$T_{500}$', '$T_{spec}$', '$\Delta T_{spec}$',
	                   '$L_{bol}$', '$L_{fit}$', '$\Delta L_{fit}$',
	                   '$E_{av}$', '$\Delta E_{av}$',
	                   '$Z$', '$\Delta Z$']
    else:
        df_all.columns = ['$T_{500}$', '$T_{spec}$', '$\Delta T_{spec}$',
	                  '$L_{bol}$', '$L_{fit}$', '$\Delta L_{fit}$',
	                  '$E_{av}$', '$\Delta E_{av}$',
	                  '$A_0$', '$A_{fit}$','$\Delta A_{fit}$',
	                  'NORM', '$\Delta$ NORM',
	                  '$B_0$', '$B_{fit}$','$\Delta B_{fit}$']
	                 
    display(df_all)                  
    #df_all.index = aven_usr1.keys()
    df_all.to_csv('tables/table_'+write_to_file+'.csv', sep=' ', header=False, index=True)
        
    return None # temp_usr1, lumin_usr1, aven_usr1






    
