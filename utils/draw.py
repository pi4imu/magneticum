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

    ax1.plot([1, 1.25], [1, 1.25], color='black', linewidth=1)

    ax1.set_xlabel(x_label, fontsize=11)
    ax1.set_ylabel(y_label_left, fontsize=11)

 #   ax1.set_xlim(1.5, 7.2)
 #   ax1.set_ylim(1.5, 7.2)
    
 #   ax1.set_xscale("log")
 #   ax1.set_yscale("log")  
 #   ax2.set_xscale("log")
 #   ax3.set_xscale("log")
    
 #   ti = [2,3,4,5,6,7]
 #   ax1.set_xticks(ti, ti)
 #   ax1.set_yticks(ti, ti)
 #   ax2.set_xticks(ti, ti)
 #   ax3.set_xticks(ti, ti)
 
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
    
    ax2.plot([], [], label=f"$\mu = {np.mean([YY-XX for YY, XX in zip(yy, xx)]):.2f}$ keV")
    ax2.plot([], [], label=f"$\sigma = {RMS:.2f}$ keV")
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
    
    ax3.plot([], [], label=f"$\mu = {np.mean(y_p):.2f}$")
    ax3.plot([], [], label=f"$\sigma = {RMS:.2f}$")
    ax3.legend(handlelength=0, frameon=False, fontsize=10, loc=1)

    if cmap:
        fig.colorbar(mappable=MAPPER, cax=ax6, orientation="horizontal").set_label(cmap_label, fontsize=12)
        
    #plt.show()
    
    
def draw_84_panels(mode):

    NNN = 84
    
    if mode=='IMAGE':
        size = 6
    else:
        size = 5

    plt.figure(figsize=((size)*7+6*3, 5*12+11*2.5))
    plt.tight_layout()
    
    #if mode!='IMAGE':
    
        #temp_compare = {}
        #lumin_compare = {}
        #average_ene = {}
        
    for cl_num in clusters.index[:NNN]:
        
        plt.subplot(12, 7, np.where(np.array(clusters.index[:NNN]) == cl_num)[0][0]+1)
        
        if mode=='IMAGE':
        
            pho_list = extract_photons_from_cluster(cl_num, r = 1, draw=True, delete_bright_regions=True)
        
        else:
        
            #cl_T500 = clusters.loc[cl_num]["T500"]
            #cl_lum = clusters.loc[cl_num]["Lx500"]
    
            SP = create_spectrum_and_fit_it(cl_num, borders=[0.4, 7.0], BACKGROUND=False, inside_radius="R500",
                                            dbr=False, Xplot=False, plot=True, draw_only=mode)

            #temp_compare[cl_num] = [cl_T500, SP[0][:3]]
            #lumin_compare[cl_num] = [cl_lum, SP[1][:3]]
            #average_ene[cl_num] = [SP[2]]

def func(p, x):

    a, b = p
    return a * x**b


def inv_func(y, a, b):
    return (y/a)**(1/b)
            
           
def draw_line(xs, x_es, ys, y_es, clr, l4dots, l4legend, argument, with_intervals=True, with_scatter=True):
    
    plt.errorbar(xs, ys, xerr=x_es, yerr=y_es, linewidth=0, marker='o', markersize=3, alpha=0.95,
                 elinewidth=1, capsize=2, color='r', label=l4dots)
                 
    plt.scatter(xs, ys, marker='o', s=6, color='r', alpha=0.95)

    #list1, list2, list3 = zip(*sorted(zip(xx, [n-q for n, q in zip(yy2, y2_err)], [n+q for n, q in zip(yy2, y2_err)])))
    #plt.fill_between(list1, list2, list3, interpolate=False, alpha=0.4, color=clr)
     
    #popt, pcov = curve_fit(func, xs, ys, maxfev=5000)

    power = odr.Model(func)
    #mydata = odr.Data(xs, ys, wd = 1./(np.array(x_es)+0.00001), we = 1./(np.array(y_es)+0.00001))
    mydata = odr.RealData(xs, ys, sx=(np.array(x_es)+0.00001), sy=(np.array(y_es)+0.00001))
    myodr = odr.ODR(mydata, power, beta0=[0, 0],  maxit=1000)
    output = myodr.run()
    popt = output.beta
    #print(popt)
    output.pprint()
    #print("stop reason:", output.stopreason)     
    #print("info:", output.info)
    #print("sd_beta:", output.sd_beta)
    #print("sqrt(diag(cov):", np.sqrt(np.diag(output.cov_beta)))
    
    if with_intervals:

        perr = output.sd_beta # np.sqrt(np.diagonal(pcov))
        #print(perr)  

        pp = (1. + 0.68)/2
        nstd = stats.norm.ppf(pp)

        popt_d = (popt[0]-nstd*perr[0], popt[1]-nstd*perr[1])
        popt_u = (popt[0]+nstd*perr[0], popt[1]+nstd*perr[1])

        #popt_d = popt-nstd*perr
        #popt_u = popt+nstd*perr
        
        lbl = f'${l4legend} = ({popt[0]:.2f} \pm {perr[0]:.2f}) \cdot {{{argument}}}^{{{popt[1]:.1f} \pm {perr[1]:.1f}}}$'
      
    else:
        
        lbl = f'${l4legend} = {popt[0]:.2f} \cdot {{{argument}}}^{{{popt[1]:.2f}}}$'    
    
    plt.plot(lll, [func(popt, XX) for XX in lll], color=clr, linewidth=3, linestyle='-', alpha=1, label=lbl)    
    
    if with_intervals:
        
        plt.fill_between(lll, 
                         [func(popt_u, XX) for XX in lll], 
                         [func(popt_d, XX) for XX in lll], 
                         interpolate=False, alpha=0.25, color=clr,
                         label='$1\sigma$ confidence band')      
        
    if with_scatter:
    
        ypyp = [(a-b)/b for a,b in zip(ys, [func(popt, XX) for XX in xs])]
    
        RMSp = np.sqrt( sum([(el**2) for el in ypyp])/len(ypyp))
        
        ypyp1 = [(a-b)/b for a,b in zip(xs, [inv_func(YY, *popt) for YY in ys])]
        
        RMSp1 = np.sqrt( sum([(el**2) for el in ypyp1])/len(ypyp))
        
        print(RMSp, RMSp1)
               
        plt.plot(lll, [func(popt, XX)*(1+RMSp) for XX in lll], color=clr, linewidth=3, linestyle='--', alpha=0.7, label=f'$1\sigma$ prediction band ($\pm${100*RMSp:.1f}%)')
        plt.plot(lll, [func(popt, XX)*(1-RMSp) for XX in lll], color=clr, linewidth=3, linestyle='--', alpha=0.7)
        
        if False:
        
            jj=0
            kk=0
    
            for ggg in ys:
    
                XCV = xs[ys.index(ggg)]
                       
                if ggg<=func(popt, XCV)*(1-RMSp) or (ggg>=func(popt, XCV)*(1+RMSp)):
                    jj+=1
                    plt.scatter(XCV, ggg, c='dodgerblue', marker='o')
                
                if ggg>=func(popt, XCV)*(1-RMSp) and (ggg<=func(popt, XCV)*(1+RMSp)):
                    kk+=1
                    plt.scatter(XCV, ggg, c='orangered')
            
            print(kk/84, jj/84, jj+kk)
    
    return None
    
    
def calculate_scatter(xs, ys, plot=True):
    
    #popt, pcov = curve_fit(func, xs, ys)
    power = odr.Model(func)
    mydata = odr.Data(xs, ys)#, wd = 1./(np.array(x_es)++0.00001), we = 1./(np.array(y_es)+0.00001))
    myodr = odr.ODR(mydata, power, beta0=[0, 0])
    output = myodr.run()
    popt = output.beta
            
    ypyp = [(a-b)/b for a,b in zip(ys, [func(popt, XX) for XX in xs])]
    
    RMSp = np.sqrt( sum([(el**2) for el in ypyp])/len(ypyp))
        
    if plot:
        
        plt.hist(ypyp, density=True, color='black', histtype='step', lw=2, bins=50)
        
        xxxccc = np.linspace(-0.4,0.4,100)
        yyyccc = stats.norm.pdf(xxxccc, loc=np.mean(ypyp), scale=RMSp)
        plt.plot(xxxccc, yyyccc, color='red', lw=2)
        
        plt.axvline(np.mean(ypyp), ls='--', color='r', lw=2)
        plt.axvline(RMSp, ymin=0, ymax=stats.norm.pdf(np.mean(ypyp)+RMSp, loc=np.mean(ypyp), scale=RMSp)/plt.gca().get_ylim()[1], ls='--', color='r', lw=2)
        plt.axvline(-RMSp, ymin=0, ymax=stats.norm.pdf(np.mean(ypyp)+RMSp, loc=np.mean(ypyp), scale=RMSp)/plt.gca().get_ylim()[1], ls='--', color='r', lw=2)
        plt.axhline(stats.norm.pdf(np.mean(ypyp)+RMSp, loc=np.mean(ypyp), scale=RMSp), 
                    xmin=(np.mean(ypyp)-RMSp-plt.gca().get_xlim()[0])/(plt.gca().get_xlim()[1]-plt.gca().get_xlim()[0]), 
                    xmax=(np.mean(ypyp)+RMSp-plt.gca().get_xlim()[0])/(plt.gca().get_xlim()[1]-plt.gca().get_xlim()[0]), 
                    ls='--', color='r', lw=2)

        plt.text(-0.3, 3.0, f"$1\sigma$ = {RMSp:.2f}", fontsize=12)
        #plt.text(-1.15, 1.0, f"         $1\sigma \ / \ T_{{best-fit}}$ = \n{RMS:.2f} keV / {np.mean([func(XX, *popt1) for XX in xx]):.2f} keV = {RMS/np.mean([func(XX, *popt1) for XX in xx]):.2f}", fontsize=12)
        
        plt.xlabel("$(T_{500} - T_{best-fit})/T_{best-fit}$", fontsize=12)
        plt.ylabel("Probability density", fontsize=12)
        plt.show()
    
    return RMSp
    
    
    
    
    
