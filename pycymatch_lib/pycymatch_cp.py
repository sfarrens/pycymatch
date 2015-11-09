#  @file pycymatch_cp.py
#
#  PYCYMATCH CP METHODS
#
#  Functions for computing
#  completeness and purity.
#
#  @author Samuel Farrens
#  @version 1.0
#  @date 2015
#

import numpy as np
from numpy.lib.recfunctions import append_fields
from functions import comp


##
#  Function to compute completeness of
#  the observed cluster sample.
#
#  @param[in] mock: List of mock haloes.
#  @param[in] matches: List of matches.
#  @param[in] opts: List of arguments.
#
#  @return Completeness matrix.
#
def get_completeness(mock, matches, opts):

    index, weights, mock_index = matches

    n_mass_bins = comp.num_bins(opts.mass_bin[0], opts.mass_bin[1],
                                opts.mass_bin[2])
    n_z_bins = comp.num_bins(opts.z_bin[0], opts.z_bin[1], opts.z_bin[2])
    n_ngal_bins = comp.num_bins(opts.ngal_bin[0], opts.ngal_bin[1],
                                opts.ngal_bin[2])

    mass_index = np.floor((mock.mass - opts.mass_bin[0]) /
                          opts.mass_bin[2]).astype('int')
    z_index = np.floor((mock.z - opts.z_bin[0]) / opts.z_bin[2]).astype('int')
    ngal_index = np.floor((mock.ngal - opts.ngal_bin[0]) /
                          opts.ngal_bin[2]).astype('int')

    mass_count = np.zeros((n_mass_bins, n_z_bins)).astype('float')
    match_count = np.zeros((n_mass_bins, n_z_bins)).astype('float')

    mass_count2 = np.zeros((n_ngal_bins, n_z_bins)).astype('float')
    match_count2 = np.zeros((n_ngal_bins, n_z_bins)).astype('float')

    for i in range(len(z_index)):
        mass_count[mass_index[i], z_index[i]] += 1.0
        mass_count2[ngal_index[i], z_index[i]] += 1.0
    for i in mock_index:
        match_count[mass_index[i], z_index[i]] += 1.0
        match_count2[ngal_index[i], z_index[i]] += 1.0

    return np.flipud(match_count / mass_count), np.flipud(match_count2 /
                                                          mass_count2)


##
#  Function to compute purity of
#  the observed cluster sample.
#
#  @param[in] obs: List of observed clusters.
#  @param[in] matches: List of matches.
#  @param[in] opts: List of arguments.
#
#  @return Completeness matrix.
#
def get_purity(obs, matches, opts):

    index, weights, mock_index = matches

    n_proxy_bins = comp.num_bins(opts.proxy_bin[0], opts.proxy_bin[1],
                                 opts.proxy_bin[2])
    n_z_bins = comp.num_bins(opts.z_bin[0], opts.z_bin[1], opts.z_bin[2])

    proxy_index = np.floor((obs.proxy - opts.proxy_bin[0]) /
                           opts.proxy_bin[2]).astype('int')
    z_index = np.floor((obs.z - opts.z_bin[0]) / opts.z_bin[2]).astype('int')

    proxy_count = np.zeros((n_proxy_bins, n_z_bins)).astype('float')
    match_count = np.zeros((n_proxy_bins, n_z_bins)).astype('float')

    for i in range(len(z_index)):
        proxy_count[proxy_index[i], z_index[i]] += 1.0
    for i in index:
        match_count[proxy_index[i], z_index[i]] += 1.0

    return np.flipud(match_count / proxy_count)
