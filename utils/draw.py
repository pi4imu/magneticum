def draw_three_panels(x_array, y_array, x_label, y_label_left, y_label_right_up, y_label_right_down, clr):
    
    plt.figure(figsize=(11.5,5.5))

    plt.subplots_adjust()
    #plt.tight_layout()

    plt.subplot(121)

    xx  = [a[1] for a in x_array]
    xxe = [a[2] for a in x_array]
    yy  = [a[1] for a in y_array]
    yye = [a[2] for a in y_array]

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
    leftc, rightc = plt.gca().get_ylim()
    
    #plt.subplot(2,5,5)
    
    #plt.hist([YY-XX for YY, XX in zip(yy, xx)], bins=30, histtype='stepfilled', orientation="horizontal", color=clr)
    #plt.ylim((leftc, rightc))

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
    leftd, rightd = plt.gca().get_ylim()
    
    #plt.subplot(2,5,10)
    #plt.hist(y_p, bins=30, histtype='stepfilled', orientation="horizontal", color=clr)
    #plt.ylim((leftd, rightd))

    plt.show()
    
    
def draw_84_panels(mode):

    NNN = 84
    
    if mode=='IMAGE':
        size=6
    else:
        size = 5

    plt.figure(figsize=((size)*7+6, size*12+11))
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
        
        
        
        
    
