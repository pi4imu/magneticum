### Saving original photon lists / rescaled for model / fakeit data:

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
        
    # Primary header (HDU = Header Data Unit):

    prihdr = fits.Header()

    prihdr['COMMENT'] = "FITS (Flexible Image Transport System) format is defined in 'Astronomy  and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H"

    prihdu = fits.PrimaryHDU(header=prihdr)
    
    # Table:
    
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
        
    #fits.open('../data/eROSITA_30.0x30.0/photons/spectrum'+str(current_cluster_num)+'.pha').info()
    
    #x.AllData.clear()
    #SPE = '../data/eROSITA_30.0x30.0/photons/spectrum'+str(current_cluster_num)+'.pha'
    #RMF = '../erosita/erosita_pirmf_v20210719.rmf'
    #ARF = '../erosita/tm1_arf_filter_000101v02.fits'
    #spectrum = x.Spectrum(SPE, respFile=RMF, arfFile=ARF)
    
    return None
    
    
    
def draw_additional_old()

       # additional figures to understand what's happening
        
        if delete_bright_regions and draw_additional:
        
            # some magic ...
        
            shift = [int((RA_c-cntr[0])*3600/ang_res), int((DEC_c-cntr[1])*3600/ang_res)]
            c1 = np.array(c)+shift
            nmhg1 = kruzhok(int(R*3600/ang_res), c1, nmhg, int(R*3600/ang_res))[0]
            
            plt.show()
            
            # initial plot in terms of pixels and only circle of R500
        
            plt.figure(figsize=(11,5))
            plt.subplot(121)
            plt.title("nmhg")
            plt.imshow(np.rot90(nmhg), norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), origin='upper')
            plt.colorbar(fraction=0.046, pad=0.04)
            plt.subplot(122)
            plt.title("nmhg1, R = "+str(int(R*3600/ang_res))+" pixels")
            plt.imshow(np.rot90(nmhg1), norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), origin='upper')  
            plt.colorbar(fraction=0.046, pad=0.04)
            plt.show()
            
            # searching for cutoff value in one pixel
            
            #threshold = 0.7       
        
            plt.figure(figsize=(12,5))
            plt.tight_layout()
            plt.subplot(121)
            beeens = np.geomspace(1, np.max(nmhg1.flatten()), 50)
            amount_in_bin, bin_borders, _ = plt.hist(nmhg1.flatten(), bins = beeens)
            bin_centers = [int(b) for b in bin_borders[:-1]]
            #bin_centers = (bin_borders[:-1]+bin_borders[1:])/2          
            plt.xlabel(f"Number of photons in {ang_res}''$\\times${ang_res}'' bin")
            plt.ylabel("Amount of such bins")
            plt.yscale("log")
            plt.xscale("log")
            plt.title("Flattened histogram for upper right image")
            plt.subplot(122)
            cumulative = [sum(amount_in_bin[:i]) for i in range(0, len(amount_in_bin))] #amount_in_bin*bin_centers
            plt.scatter(bin_centers, cumulative)
            plt.xlabel(f"Number of photons in {ang_res}''$\\times${ang_res}'' bin")
            #plt.ylabel("Total number of photons ($x \cdot y$ for histogram on the left)")
            plt.ylabel("Cumulative distribution")
            plt.yscale("log")
            plt.xscale("log")         
            plt.title(f"$x$-values are added up until their sum\nis right below {threshold*100:.0f} % cutoff")    
 
            plt.axhline(cumulative[-1], ls='--', color='red')
            plt.axhline(cumulative[-1]*threshold, ls='--', color='red')
            
            number_cutoff = bin_centers[np.argmin(np.abs(cumulative-cumulative[-1]*threshold))]
            
            plt.axvline(number_cutoff, ls='--', color='red')
            
            #sum_amount, i = 0, 0
            #while sum_amount <= threshold * sum(amount_in_bin*bin_centers):
            #    sum_amount = sum_amount + (amount_in_bin*bin_centers)[i]
            #    #print(i, bin_centers[i], sum_amount/sum(amount_in_bin*bin_centers))
            #    i = i + 1
            #threshold = (sum_amount-(amount_in_bin*bin_centers)[i-1])/sum(amount_in_bin*bin_centers)        
            #number_cutoff = bin_centers[i-1]
            
            #plt.axvline(number_cutoff, ls='--', color='red', label=f'{threshold*100:.2f} % cutoff\nat brightness = {number_cutoff:.2f}')
            #plt.legend()
            
            #plt.subplot(121)
            #plt.axvline(number_cutoff, ls='--', color='red', label=f'{threshold*100:.2f} % cutoff\nat brightness = {number_cutoff:.2f}')
            plt.show()
            
            # making masks and applying them to images
            
            nmhg_radial = nmhg
            
            filter_mask = nmhg <= number_cutoff
            nmhg = nmhg*filter_mask
                       
            filter_mask1 = nmhg1 <= number_cutoff
            nmhg1 = nmhg1*filter_mask1
            
            plt.figure(figsize=(11,5))
            plt.subplot(121)
            plt.title("nmhg1 (filtered), R500 is "+str(int(R*3600/ang_res))+" pixels")
            plt.imshow(np.rot90(nmhg1), norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), origin='upper')
            plt.colorbar(fraction=0.046, pad=0.04)
            plt.subplot(122)
            plt.title("nmhg (filtered)")
            plt.imshow(np.rot90(nmhg), norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1), origin='upper')
            plt.colorbar(fraction=0.046, pad=0.04)
            plt.show()
            
            #plt.imshow(np.rot90(filter_mask1))
            #plt.gca().set_aspect('equal', 'box')
            #plt.show()
            
            # obsolete section
                      
            if False:
                                 
                dddfff["RA_pix"] = (dddfff["RA"] - cntr[0] + R)*3600/ang_res 
                dddfff["DEC_pix"] = (dddfff["DEC"] - cntr[1] + R)*3600/ang_res
            
                dddfff["RA_pix"] = dddfff["RA_pix"].astype(int)
                dddfff["DEC_pix"] = dddfff["DEC_pix"].astype(int)
            
                #hj, _, _, _ = plt.hist2d(dddfff["RA_pix"], dddfff["DEC_pix"],
                #          bins=len(nmhg1),
                #          norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1)) 
                #plt.gca().set_aspect('equal', 'box')
                #plt.show()
                #plt.imshow(np.rot90(filter_mask))
                #plt.gca().set_aspect('equal', 'box')
                #plt.show()
            
                dddfff["stay"] = filter_mask1[dddfff["RA_pix"], dddfff["DEC_pix"]]
                dddfff = dddfff[dddfff["stay"] == True]
            
                plt.hist2d(dddfff["RA"], dddfff["DEC"],
                           bins=len(nmhg1),
                           norm=matplotlib.colors.SymLogNorm(linthresh=1, linscale=1))
                plt.gca().set_aspect('equal', 'box')
                plt.colorbar(fraction=0.046, pad=0.04)
                plt.show()
            
                dddfff = dddfff.drop("stay", axis=1) 
                dddfff = dddfff.drop("RA_pix", axis=1)
                dddfff = dddfff.drop("DEC_pix", axis=1)

    return None




def draw_three_panels(x_array, y_array, x_label, y_label_left, y_label_right_up, y_label_right_down, clr, NnNn, cmap=False, cmap_label="Don't forget to rename me!"):
   
    if not cmap:   
        fig = plt.figure(figsize=(11.5,5.5))
    else:
        fig = plt.figure(figsize=(11.5,6.5))
        NORM = matplotlib.colors.Normalize(vmin=min(cmap), vmax=max(cmap), clip=True)
        MAPPER = cm.ScalarMappable(norm=NORM, cmap='rainbow')
        COLOUR = np.array([(MAPPER.to_rgba(v)) for v in cmap])
    
#    plt.suptitle(f"    Mean values for {NnNn} realisations", fontsize=15)
    
    gs = GridSpec(4, 4, height_ratios=[1, 1, 0.05, 0.1], width_ratios=[5, 1.5, 4, 1], hspace = 0.3, wspace = 0.)
    ax1 = fig.add_subplot(gs[0:2, 0])
    ax2 = fig.add_subplot(gs[2])
    ax3 = fig.add_subplot(gs[6])
    ax4 = fig.add_subplot(gs[3])
    ax5 = fig.add_subplot(gs[7])
    if cmap:
        ax6 = fig.add_subplot(gs[12:16])

    #plt.subplots_adjust()
    #gs.tight_layout(figure=fig)

    #plt.subplot(121)

    xx  = [a[1] for a in x_array]
    xxe = [a[2] for a in x_array]
    yy  = [a[1] for a in y_array]
    yye = [a[2] for a in y_array]

    if not cmap:
        ax1.errorbar(xx, yy, xerr=xxe, yerr=yye, linewidth=0, elinewidth=1, 
                     capsize=3, color=clr, marker='o', markersize=3)
    else:
        for xxx, ex, yyy, ey, col in zip(xx, xxe, yy, yye, COLOUR):
            ax1.plot(xxx, yyy, 'o', color=col, markersize=3)
            ax1.errorbar(xxx, yyy, yerr=ey, xerr=ex, elinewidth=1, capsize=3, color=col)

    ax1.plot([1, 9], [1, 9], color='black', linewidth=1)

    ax1.set_xlabel(x_label, fontsize=11)
    ax1.set_ylabel(y_label_left, fontsize=11)

    if True:
    
        ax1.set_xlim(1.5, 7.2)
        ax1.set_ylim(1.5, 7.2)
    
        ax1.set_xscale("log")
        ax1.set_yscale("log")  
        ax2.set_xscale("log")
        ax3.set_xscale("log")
    
        ti = [2,3,4,5,6,7]
        ax1.set_xticks(ti, ti)
        ax1.set_yticks(ti, ti)
        ax2.set_xticks(ti, ti)
        ax3.set_xticks(ti, ti)
 
    #plt.subplot(222)
    
    y_d = [YY-XX for YY, XX in zip(yy, xx)]
    y_d_err = [a+b for a, b in zip(xxe, yye)]

    if not cmap:
        ax2.errorbar(xx, y_d, xerr=xxe, yerr=y_d_err, 
                     linewidth=0, elinewidth=1, capsize=3, color=clr, marker='o', markersize=3)
    else:
        for xxx, ex, yyy, ey, col in zip(xx, xxe, y_d, y_d_err, COLOUR):
            ax2.plot(xxx, yyy, 'o', color=col, markersize=3)
            ax2.errorbar(xxx, yyy, yerr=ey, xerr=ex, elinewidth=1, capsize=3, color=col)

    ax2.axhline(0, color='black', linewidth=1)
    ax2.set_ylabel(y_label_right_up, fontsize=11)
    #ax2.set_ylim(-3, 3)
    
    leftb, rightb = ax2.get_xlim()
    leftc, rightc = ax2.get_ylim()
        
    ba, bi, _ = ax4.hist(y_d, bins=20, histtype='stepfilled', orientation="horizontal", color=clr, density=True)
    ax4.set_ylim((leftc, rightc))
    #ax4.set_xscale("log")
    ax4.set_yticks([],[])
    ax4.set_xticks([],[])
    ax4.xaxis.set_ticks_position('none') 
    ax4.axhline(0, color='black', linewidth=1)
    
    RMS = np.sqrt( sum([(el**2) for el in y_d])/len(y_d) )  
    xxxccc = np.linspace(-3,3,1000)
    yyyccc = stats.norm.pdf(xxxccc, loc=np.mean(y_d), scale=RMS)
    ax4.plot(yyyccc, xxxccc, color='black')
    
    ax2.plot([], [], label=f"$\\mu = {np.mean([YY-XX for YY, XX in zip(yy, xx)]):.2f}$ keV")
    ax2.plot([], [], label=f"$\\sigma = {RMS:.2f}$ keV")
    ax2.legend(handlelength=0, frameon=False, fontsize=10, loc=1)

    #plt.subplot(224)
    
    y_p = [(YY-XX)/XX for YY, XX in zip(yy, xx)]
    y_p_err = [a/b*(aa/a+bb/b) for a, aa, b, bb in zip(yy, yye, xx, xxe)]

    if not cmap:
        ax3.errorbar(xx, y_p, xerr=xxe, yerr=y_p_err, linewidth=0, elinewidth=1, capsize=3, color=clr, marker='o', markersize=3)
        ax3.scatter(xx, y_p, color=clr, marker='o', s=3)
    else:
        for xxx, ex, yyy, ey, col in zip(xx, xxe, y_p, y_p_err, COLOUR):
            ax3.plot(xxx, yyy, 'o', color=col, markersize=3)
            ax3.errorbar(xxx, yyy, yerr=ey, xerr=ex, elinewidth=1, capsize=3, color=col)
                 
    #list1, list2, list3 = zip(*sorted(zip(xx, [n-q for n, q in zip(y_p, y_p_err)], [n+q for n, q in zip(y_p, y_p_err)])))
    #ax3.fill_between(list1, list2, list3, interpolate=True, alpha=0.4, color=clr)

    ax3.axhline(0, color='black', linewidth=1)
    ax3.set_ylabel(y_label_right_down, fontsize=11)
    ax3.set_xlabel(x_label, fontsize=11)
    #ax3.set_ylim(-0.8, 0.8)
    
    ax3.set_xlim(leftb, rightb)
    leftd, rightd = ax3.get_ylim()
        
    ba, bi, _ = ax5.hist(y_p, bins=20, histtype='stepfilled', orientation="horizontal", color=clr, density=True)
    ax5.set_ylim((leftd, rightd))
    #ax5.set_xscale("log")
    ax5.set_yticks([],[])
    ax5.set_xticks([],[])
    ax5.xaxis.set_ticks_position('none') 
    ax5.axhline(0, color='black', linewidth=1)
    
    RMS = np.sqrt(sum([(el**2) for el in y_p])/len(y_p))
    xxxccc = np.linspace(-1,1,1000)
    yyyccc = stats.norm.pdf(xxxccc, loc=np.mean(y_p), scale=RMS) 
    ax5.plot(yyyccc, xxxccc, color='black')
    
    ax3.plot([], [], label=f"$\\mu = {np.mean(y_p):.2f}$")
    ax3.plot([], [], label=f"$\\sigma = {RMS:.2f}$")
    ax3.legend(handlelength=0, frameon=False, fontsize=10, loc=1)

    if cmap:
        fig.colorbar(mappable=MAPPER, cax=ax6, orientation="horizontal").set_label(cmap_label, fontsize=12)
        
    #plt.show()
 
