def draw_three_panels(x_array, y_array, x_label, y_label_left, y_label_right_up, y_label_right_down, clr):
    
    plt.figure(figsize=(11.5,5.5))

    plt.subplots_adjust()

    plt.subplot(121)

    xx  = [a[1] for a in x_array.values()]
    xxe = [a[2] for a in x_array.values()]
    yy  = [a[1] for a in y_array.values()]
    yye = [a[2] for a in y_array.values()]

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
    
    
def plot_T_vs_avE(avE, Tsp):

	plt.figure(figsize=(7.1,7.1))

	xx = [a[0] for a in avE.values()]
	x_err = [a[1] for a in avE.values()]
	yy1 = [a[0] for a in Tsp.values()]
	yy2 = [a[1] for a in Tsp.values()]
	y2_err = [a[2] for a in Tsp.values()]

	def func(x, a, b):
	    return a * x**b

	popt1, pcov1 = curve_fit(func, xx, yy1)
	popt2, pcov2 = curve_fit(func, xx, yy2)

	#for xxx, xe, yyy, ye, col in zip(xx, x_err, yy2, y2_err, mass_colour):
	#   plt.plot(xxx, yyy, '.', color=col)
	#    plt.errorbar(xxx, yyy, xerr=xe, yerr=ye, elinewidth=1, capsize=3, color=col, label='$T_{spec}$ from fit')

	list1, list2, list3 = zip(*sorted(zip(xx, [n-q for n, q in zip(yy2, y2_err)], [n+q for n, q in zip(yy2, y2_err)])))
	plt.fill_between(list1, list2, list3, interpolate=False, alpha=0.4, color='blue')

	#plt.errorbar(xx, yy2, xerr=x_err, yerr=y2_err, linewidth=0, marker='o', markersize=4, alpha=0.1,
	#             elinewidth=1, capsize=3, color='blue', label='$T_{spec}$ from fit')

	lll = np.linspace(0.95, 1.29, 100)
	plt.plot(lll, [func(XX, *popt1) for XX in lll],  
		 color='red', linewidth=3, linestyle=':', alpha=1,
		 label=f'Best fit: $T_{{500}} = {popt1[0]:.2f} \cdot {{E_{{av}}}}^{{{popt1[1]:.2f}}}$')
	plt.plot(lll, [func(XX, *popt2) for XX in lll],  
		 color='blue', linewidth=3, linestyle='--', alpha=1,
		 label=f'Best fit: $T_{{spec}} = {popt2[0]:.2f} \cdot {{E_{{av}}}}^{{{popt2[1]:.2f}}}$')

	plt.xlabel("Average energy, keV")
	plt.ylabel("Temperature, keV")

	plt.xscale("log")
	plt.yscale("log")

	#sc = plt.scatter(xx, yy1, c='red', s=0)
	#clb = plt.colorbar(sc, label = "$M_{500}$ in units of $10^{14} M_{\odot} h^{-1}$")

	list1, list2, list3 = zip(*sorted(zip(yy1, [n-q for n, q in zip(xx, x_err)], [n+q for n, q in zip(xx, x_err)])))
	plt.gca().fill_betweenx(list1, list2, list3, interpolate=True, alpha=0.4, color='red')

	#plt.errorbar(xx, yy1, xerr=x_err, linewidth=0, elinewidth=1, capsize=3,
	#             color='red', marker='o', markersize=4, alpha=0.1, label='$T_{500}$ from simulations')

	#plt.xlim(1.4, 7.5)
	#plt.ylim(0.9, 1.3)

	plt.xticks([0.95, 1.00, 1.05, 1.10, 1.15, 1.2, 1.25, 1.3], [0.95, 1.00, 1.05, 1.10, 1.15, 1.2, 1.25, 1.3])
	plt.yticks([2,3,4,5,6,7,8,9,10], [2,3,4,5,6,7,8,9,10])

	for i in range(0, len(aven_usr)):
	    plt.plot([xx[i]+0.000454, xx[i]+0.000454], 
		     [yy1[i], yy2[i]], 
		     color='grey', alpha=0.4, marker='o', markersize=4)

	#plt.grid()
	plt.legend()

	plt.show()






    
