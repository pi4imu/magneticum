def draw_three_panels(x_array, y_array, x_label, y_label_left, y_label_right_up, y_label_right_down, clr, NnNn=10):
    
    fig = plt.figure(figsize=(11.5,5.5))
    
        
    plt.suptitle(f"    Mean values for {NnNn} realisations", fontsize=15)
    
    gs = GridSpec(2, 4, height_ratios=[1, 1], width_ratios=[5, 1, 4, 1], hspace = 0.2, wspace = 0.)
    ax1 = fig.add_subplot(gs[:, 0])
    ax2 = fig.add_subplot(gs[2])
    ax3 = fig.add_subplot(gs[6])
    ax4 = fig.add_subplot(gs[3])
    ax5 = fig.add_subplot(gs[7])

    #plt.subplots_adjust()
    #gs.tight_layout(figure=fig)

    #plt.subplot(121)

    xx  = [a[1] for a in x_array]
    xxe = [a[2] for a in x_array]
    yy  = [a[1] for a in y_array]
    yye = [a[2] for a in y_array]

    ax1.errorbar(xx, yy, xerr=xxe, yerr=yye, linewidth=0, elinewidth=1, 
                 capsize=3, color=clr, marker='o', markersize=3)

    ax1.plot([1, 9], [1, 9], color='black', linewidth=1)

    ax1.set_xlabel(x_label, fontsize=11)
    ax1.set_ylabel(y_label_left, fontsize=11)

    ax1.set_xlim(1., 7.3)
    ax1.set_ylim(1., 7.3)
    #ax1.scatter(2.39539, 1.1103637527004389, color='red')

    #plt.subplot(222)

    ax2.errorbar(xx, [YY-XX for YY, XX in zip(yy, xx)], xerr=xxe, 
                 yerr=[a+b for a, b in zip(xxe, yye)], linewidth=0, elinewidth=1, 
                 capsize=3, color=clr, marker='o', markersize=3)

    ax2.axhline(0, color='black', linewidth=1)
    ax2.set_ylabel(y_label_right_up, fontsize=11)
#    ax2.set_ylim(-1.5, 3.5)
    
    leftb, rightb = ax2.get_xlim()
    leftc, rightc = ax2.get_ylim()
        
    ax4.hist([YY-XX for YY, XX in zip(yy, xx)], bins=20, histtype='stepfilled', orientation="horizontal", color=clr)
    ax4.set_ylim((leftc, rightc))
    #ax4.set_xscale("log")
    ax4.set_yticks([],[])
    ax4.set_xticks([],[])
    ax4.xaxis.set_ticks_position('none') 
    ax4.axhline(0, color='black', linewidth=1)
    
    RMS = np.sqrt(sum([(el**2) for el in [YY-XX for YY, XX in zip(yy, xx)]]))/len([YY-XX for YY, XX in zip(yy, xx)])    
    xxxccc = np.linspace(-3,3,100)
    ax4.plot(stats.norm.pdf(xxxccc, loc=np.mean([YY-XX for YY, XX in zip(yy, xx)]), scale=RMS), xxxccc, color='black')

    #plt.subplot(224)
    
    y_p = [(YY-XX)/XX for YY, XX in zip(yy, xx)]
    y_p_err = [a/b*(aa/a+bb/b) for a, aa, b, bb in zip(yy, yye, xx, xxe)]
    
    ax3.errorbar(xx, y_p, xerr=xxe, yerr=y_p_err, linewidth=0, elinewidth=1, capsize=3, color=clr, marker='o', markersize=3)
    ax3.scatter(xx, y_p, color=clr, marker='o', s=3)
                 
    #list1, list2, list3 = zip(*sorted(zip(xx, [n-q for n, q in zip(y_p, y_p_err)], [n+q for n, q in zip(y_p, y_p_err)])))
    #ax3.fill_between(list1, list2, list3, interpolate=True, alpha=0.4, color=clr)

    ax3.axhline(0, color='black', linewidth=1)
    ax3.set_ylabel(y_label_right_down, fontsize=11)
    ax3.set_xlabel(x_label, fontsize=11)
 #   ax2.set_ylim(-0.6, 0.8)
    
    ax3.set_xlim(leftb, rightb)
    leftd, rightd = ax3.get_ylim()
        
    ax5.hist(y_p, bins=20, histtype='stepfilled', orientation="horizontal", color=clr)
    ax5.set_ylim((leftd, rightd))
    ax5.set_xscale("log")
    ax5.set_yticks([],[])
    ax5.set_xticks([],[])
    ax5.xaxis.set_ticks_position('none') 
    ax5.axhline(0, color='black', linewidth=1)
    
  #  RMS = np.sqrt(sum([(el**2) for el in y_p]))/len(y_p)
  #  xxxccc = np.linspace(-3,3,100)
  #  ax5.plot(stats.norm.pdf(xxxccc, loc=np.mean(y_p), scale=RMS), xxxccc, color='black')

    #plt.show()
    
    
def draw_84_panels(mode):

    NNN = 84
    
    if mode=='IMAGE':
        size = 6
    else:
        size = 5

    plt.figure(figsize=((size)*7+6, 5*12+11))
    plt.tight_layout()
    
    if mode!='IMAGE':
    
        temp_compare = {}
        lumin_compare = {}
        average_ene = {}
        
    for cl_num in clusters.index[:NNN]:
        
        plt.subplot(12, 7, np.where(np.array(clusters.index[:NNN]) == cl_num)[0][0]+1)
        
        if mode=='IMAGE':
        
            pho_list = extract_photons_from_cluster(cl_num, r = 'R500', draw=True)
        
        else:
        
            cl_T500 = clusters.loc[cl_num]["T500"]
            cl_lum = clusters.loc[cl_num]["Lx500"]
    
            SP = create_spectrum_and_fit_it(cl_num, borders=[0.4, 7.0], BACKGROUND=True, inside_radius="R500",
                                            Xplot=False, plot=True, draw_only=mode)

            temp_compare[cl_num] = [cl_T500, SP[0][:3]]
            lumin_compare[cl_num] = [cl_lum, SP[1][:3]]
            average_ene[cl_num] = [SP[2]]
        
        
        
        
    
