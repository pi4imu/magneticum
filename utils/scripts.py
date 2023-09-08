# "clusters" and "binned_clusters" are external lists

# returns list of photons inside chosen radius

def extract_photons_from_cluster(ihal, r, draw=True):

    current_cluster_num = ihal
    
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
    
        plt.figure(figsize=(6,5))
        plt.scatter(dddfff["RA"], dddfff["DEC"], c=dddfff["ENERGY"], cmap='viridis', s=1)
        plt.colorbar()
        plt.show()
    
    return dddfff
    
    
# "erosita_binning" is external list
    
def create_spectrum(list_of_photons, REDSHIFT, create_textfile=True, create_pha=True):

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
    
    
def fit_spectra(filename):
 
    #fits.open(filename).info()
    x.AllData.clear()

    RMF = '../erosita/erosita_pirmf_v20210719.rmf'

    ARF = '../erosita/tm1_arf_filter_000101v02.fits'

    SPECTRUM = x.Spectrum(filename, respFile=RMF, arfFile=ARF)
    
    x.Plot.device = "/xs"
    x.Plot.xAxis = "keV"
    x.AllData.ignore("**-0.4 7.0-**")

    x.AllModels.clear()

    x.Plot("ldata")

    xVals = x.Plot.x()
    yVals = x.Plot.y()

    plt.figure(figsize=(5, 5))

    plt.plot(xVals, yVals)# np.multiply(np.multiply(xVals,xVals),yVals))
    plt.xscale('log')
    plt.yscale('log')

    plt.xlabel(x.Plot.labels()[0])
    plt.ylabel(x.Plot.labels()[1])
    #plt.title(x.Plot.labels()[2])

    #plt.ylim(10**(-1), 10**2)

    plt.show()











    
