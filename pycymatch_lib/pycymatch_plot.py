#  @file pycymatch_plot.py
#
#  PYCYMATCH PLOTTING METHODS
#
#  Functions to plot code
#  results.
#
#  @author Samuel Farrens
#  @version 1.1
#  @date 2015
#

import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages

from pycymatch_fit import gauss, get_fits
from pycymatch_match2 import scale_matrix


##
#  Make mass-obs and P(lambda|M,z) plots.
#
#  @param[in] matrix: 3D Mass-observable matrix.
#  @param[in] hm_matrix: Mass vs z matrix.
#  @param[in] mass_hist: Halo mass histograms.
#  @param[in] ranges: Mass-observable ranges.
#  @param[in] opts: List of arguments.
#
def make_plots(matrix, hm_matrix, ranges, opts):

    mo_plt_name = opts.input_files[0] + '.pycy.mass_obs.pdf'
    nl_plt_name = opts.input_files[0] + '.pycy.n_lambda.pdf'
    mf_plt_name = opts.input_files[0] + '.pycy.mass_func.pdf'

    mass_obs_plots = PdfPages(mo_plt_name)
    n_lambda_plots = PdfPages(nl_plt_name)
    mass_func_plots = PdfPages(mf_plt_name)

    for i in range(len(ranges[0])):
        if not np.all(matrix[i] == 0.0):
            mass_obs(np.flipud(matrix[i]), ranges[0][i], opts, mass_obs_plots)
            n_lambda(matrix[i], hm_matrix, ranges, i, opts, n_lambda_plots)
            mass_func(hm_matrix[i], matrix[i], ranges[1], ranges[0][i], opts,
                      mass_func_plots)

    mass_obs_plots.close()
    n_lambda_plots.close()
    mass_func_plots.close()

    print "Plots saved to:", mo_plt_name
    print "Plots saved to:", nl_plt_name
    print "Plots saved to:", mf_plt_name


##
#  Plot mass-observable relation.
#
#  @param[in] matrix: 2D Mass-observable matrix.
#  @param[in] z_bin: Redshift bin.
#  @param[in] opts: List of arguments.
#  @param[in] page: PdfPages file.
#
def mass_obs(matrix, z_bin, opts, page):

    fig = plt.figure()

    extent = [opts.proxy_bin[0], opts.proxy_bin[1], opts.mass_bin[0],
              opts.mass_bin[1]]

    if np.any(matrix > 0.0):
        normalisation = mp.colors.LogNorm()
    else:
        normalisation = None

    plt.imshow(matrix, extent=extent, interpolation='spline16', cmap=cm.jet,
               aspect='auto', norm=normalisation)

    cbar = plt.colorbar()
    cbar.set_label(r'$n(\lambda|M,z)$', fontsize=15)

    plt.autoscale(False)
    plt.xlabel(r'$\log_{10} \lambda$', fontsize=20)
    plt.ylabel(r'$\log_{10} M$ $[M_{\odot} h^{-1}]$', fontsize=20)
    plt.title(r'z$_i$ = ' + str(z_bin), fontsize=20)

    page.savefig(fig)

    plt.close()


##
#  Plot n(L|M) distribution.
#
#  @param[in] matrix: Mass-observable matrix.
#  @param[in] opts: List of arguments.
#
def n_lambda(matrix, hm_matrix, ranges, loop_i, opts, page):

    mp.rcParams.update({'font.size': 8})

    fig = plt.figure()
    fig.set_size_inches(10.0, 15.0)

    gs = gridspec.GridSpec(len(ranges[1]) / 5, 5)
    fig.subplots_adjust(hspace=1.0)
    fig.subplots_adjust(wspace=0.0)

    fits = get_fits(matrix, ranges[2])

    for i in range(len(ranges[1])):

        ax = fig.add_subplot(gs[i])

        if i % 5:
            ax.yaxis.set_visible(False)
        else:
            ax.set_ylabel(r"$n(\lambda|M)$", fontsize=12)

        ax.set_xlabel(r"$\log_{10} \lambda$", fontsize=12)
        ax.set_title(r"$\log_{10}$ M$_i$ = " + str(ranges[1][i]) +
                     ', z$_i$ = ' + str(ranges[0][loop_i]) +
                     '\n n(M,z) = ' + str(int(hm_matrix[loop_i][i])),
                     fontsize=7)

        ax.plot(ranges[2], matrix[i], '-')
        if len(fits[i]) > 0.0 and fits[i][2] > 0.0:
            ax.plot(ranges[2], gauss(ranges[2], *fits[i]), 'r--',
                    label=r"Fit: $\sigma = $ %0.2f" % fits[i][2])
            ax.legend(prop={'size': 6})

        ax.set_xlim(opts.proxy_bin[0], opts.proxy_bin[1])
        ax.set_ylim(0.0, np.sum(matrix[i]))
        ax.xaxis.set_major_locator(MaxNLocator(2))
        ax.yaxis.set_major_locator(MaxNLocator(3))

    page.savefig(fig)

    plt.close()


##
#  Plot sigma values of P(L|M) fits.
#
#  @param[in] fits: Fits to P(L|M).
#  @param[in] opts: List of arguments.
#
def plot_sigma(fits, opts):

    fig = plt.figure(3)

    fit = fits[:, np.argsort(fits[1])]

    plt.plot(fit[1], fit[2], 'bx:')

    plt.xlabel(r"$Log_{10} \Lambda$",  fontsize=20)
    plt.ylabel(r"$\sigma$", fontsize=20)
    plt.title(str(opts.z_bin[0]) + r'$\leq z \leq$' + str(opts.z_bin[1]),
              fontsize=20)

    plt.xlim(opts.proxy_bin[0], opts.proxy_bin[1])
    plt.ylim(0.0, 0.5)

    plt_name = (opts.input_files[0] + '.' + str(opts.z_bin[0]) + 'to' +
                str(opts.z_bin[1]) + '.pycy.p_lambda_fits.pdf')
    plt.savefig(plt_name)
    print "Plot saved to:", plt_name

    plt.close()


##
#  Plot mass-function.
#
#  @param[in] hm_matrix: Histograms of mock
#  halo masses.
#  @param[in] matrix: Mass-observable matrix.
#  @param[in] x_range: X-axis data.
#  @param[in] z_bin: Redshift bin.
#  @param[in] opts: List of arguments.
#  @param[in] page: PdfPages file.
#
def mass_func(mass_hist, matrix, x_range, z_bin, opts, page):

    fig = plt.figure()

    plt.plot(x_range, mass_hist, '-', color='blue', label='Mock Haloes: ' +
             str(int(np.sum(mass_hist))))
    plt.plot(x_range, matrix.sum(axis=1), '--', color='red',
             label='Observed Clusters: ' + str(int(np.sum(matrix))))

    plt.ylabel(r'$N(M)$', fontsize=20)
    plt.xlabel(r'$\log_{10} M$ $[M_{\odot} h^{-1}]$', fontsize=20)
    plt.title(r'z$_i$ = ' + str(z_bin), fontsize=20)

    max_val = np.ceil(np.log10(np.max([np.max(mass_hist),
                      np.max(matrix.sum(axis=0))])))

    plt.ylim(10 ** 0, 10 ** max_val)
    plt.yscale('log')
    plt.legend()

    page.savefig(fig)

    plt.close()


##
#  Plot completeness.
#
#  @param[in] matrix: Completeness matrix.
#  @param[in] plot_type: Type of completeness plot.
#  @param[in] opts: List of arguments.
#
def plot_complete(matrix, plot_type, opts):

    fig_num = 5
    if plot_type == 'ngal':
        fig_num = 6

    fig = plt.figure(fig_num)

    fig.subplots_adjust(wspace=.01)
    gs = gridspec.GridSpec(1, 2, width_ratios=[12, 1])
    ax = fig.add_subplot(gs[0])
    cax = fig.add_subplot(gs[1])

    boundaries = np.arange(0.0, 1.1, 0.1)
    colormap = cm.jet

    ax.set_xlabel('$z$', fontsize=24)
    if plot_type == 'mass':
        ax.set_ylabel(r'$\log_{10} M$', fontsize=24)
        tag = '.pycy.complete_mass.pdf'
        extent = [opts.z_bin[0], opts.z_bin[1], opts.mass_bin[0],
                  opts.mass_bin[1]]
        ax.set_ylim(opts.mass_bin[0], opts.mass_bin[1])
    elif plot_type == 'ngal':
        ax.set_ylabel(r'$\log_{10}$ $N_{gal}$', fontsize=24)
        tag = '.pycy.complete_ngal.pdf'
        extent = [opts.z_bin[0], opts.z_bin[1], opts.proxy_bin[0],
                  opts.proxy_bin[1]]
        ax.set_ylim(opts.proxy_bin[0], opts.proxy_bin[1])

    im = ax.imshow(matrix, extent=extent, aspect='auto',
                   interpolation='spline16', cmap=colormap,
                   norm=mp.colors.BoundaryNorm(boundaries, colormap.N))

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(r'Completeness', fontsize=20)

    plt_name = (opts.input_files[0] + '.' + str(opts.z_bin[0]) + 'to' +
                str(opts.z_bin[1]) + tag)
    plt.savefig(plt_name)
    print 'Plot saved to:', plt_name


##
#  Plot purity.
#
#  @param[in] matrix: Purity matrix.
#  @param[in] opts: List of arguments.
#
def plot_pure(matrix, opts):

    fig = plt.figure(7)

    fig.subplots_adjust(wspace=.01)
    gs = gridspec.GridSpec(1, 2, width_ratios=[12, 1])
    ax = fig.add_subplot(gs[0])
    cax = fig.add_subplot(gs[1])

    boundaries = np.arange(0.0, 1.1, 0.1)
    colormap = cm.jet

    ax.set_xlabel('$z$', fontsize=24)
    ax.set_ylabel('$\log_{10} \lambda$', fontsize=24)

    im = ax.imshow(matrix, extent=[opts.z_bin[0], opts.z_bin[1],
                   opts.proxy_bin[0], opts.proxy_bin[1]], aspect='auto',
                   interpolation='spline16', cmap=colormap,
                   norm=mp.colors.BoundaryNorm(boundaries, colormap.N))

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label('Purity', fontsize=20)

    plt_name = (opts.input_files[0] + '.' + str(opts.z_bin[0]) + 'to' +
                str(opts.z_bin[1]) + '.pycy.purity.pdf')
    plt.savefig(plt_name)
    print 'Plot saved to:', plt_name
