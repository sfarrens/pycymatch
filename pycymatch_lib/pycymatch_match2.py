#  @file pycymatch_match.py
#
#  PYCYMATCH MATCH METHODS
#
#  Functions for finding matches
#  between observed and mock
#  clusters.
#
#  @author Samuel Farrens
#  @version 2.0
#  @date 2015
#

import sys
import numpy as np
import itertools
from functions import astro, comp


##
#  Find observed clusters that satisfy the
#  basic match criteria.
#
#  @param[in] i: Index of observed cluster.
#  @param[in] mock: List of mock haloes.
#  @param[in] obs: List of observed clusters.
#  @param[in] dz: Redshift matching treshold.
#
#  @return List of indices of observed cluster
#  that satisfy the matching criteria for a
#  given mock halo.
#
def basic_match(i, mock, obs, dz):

    return np.where((np.fabs(mock.z[i] - obs.z) <=
                    (dz * (1 + mock.z[i]))) &
                    (obs.ra >= mock.minra[i]) &
                    (obs.ra <= mock.maxra[i]) &
                    (obs.dec >= mock.mindec[i]) &
                    (obs.dec <= mock.maxdec[i]) &
                    (obs.flag == 0))[0]


#  Define weighting for matches.
#
#  @param[in] proxy: List of mass proxies of observed clusters.
#  @param[in] dist: List of distances to halo centre.
#  @param[in] r200: R_200 of mock halo.
#
#  @return Weights.
#
def match_weight(proxy, dist, r200):

    weights = (proxy ** 2) + (1.0 - (dist / r200) ** 2)

    return comp.scale(weights, 0.0, np.sum(weights))


##
#  Find observed clusters within r200 of mock
#  halo centre.
#
#  @param[in] i: Index of mock halo.
#  @param[in] index: List of indices of matches.
#  @param[in] mock: List of mock haloes.
#  @param[in] obs: List of observed clusters.
#
#  @return Indices of mathces and corresponding
#  distances of observed clusters to a given
#  mock halo.
#
def r200_match(i, index, mock, obs):

    dists = astro.ang_sep((mock.ra[i], mock.dec[i]),
                          (obs.ra[index], obs.dec[index])) * 60.0

    new_index = np.where((dists - (1.0 * mock.r200[i])) <= 0.0)[0]

    if len(new_index) > 0:
        return index[new_index], dists[new_index]
    else:
        return [False], [False]


##
#  Find cylindrical matches between mock
#  haloes and observed clusters.
#
#  @param[in] mock: List of mock haloes.
#  @param[in] obs: List of observed clusters.
#  @param[in] dz: Redshift matching treshold.
#  @param[in] opts: List of arguments.
#
#  @return Indices of mathces and corresponding
#  weights of mock haloes to observed clusters.
#
def find_matches(mock, obs, dz, opts):

    index_list = []
    weights = []
    mock_index = []

    for i in range(mock.size):
        index = basic_match(i, mock, obs, dz)

        if len(index) > 0:
            x, y = r200_match(i, index, mock, obs)
            if np.any(x):
                w = match_weight(obs.proxy[x], y, mock.r200[i])
                if opts.unique:
                    x = x[np.where(w == max(w))[0][0]]
                    w = w[np.where(w == max(w))[0][0]]
                    obs.flag[x] = 1
                index_list.append(x)
                weights.append(w)
                mock_index.extend([i])

    print 'Dected', len(np.unique(mock_index)), 'out of', len(mock), \
          'mock haloes.'

    return np.array(index_list), np.array(weights), np.array(mock_index)


##
#  Produce mass-observable matrix.
#
#  @param[in] mock: List of mock haloes.
#  @param[in] obs: List of observed clusters.
#  @param[in] matches: List of matches.
#  @param[in] opts: List of arguments.
#
#  @return Mass-observable matrix, histogram
#  of all mock halo masses, histogram of
#  matched halo masses, x-range values of
#  halo masses, x-range values of observed
#  cluster mass proxies.
#
def mo_matrix(mock, obs, matches, opts):

    index, weights, mock_index = matches

    n_z_bins = comp.num_bins(opts.z_bin[0], opts.z_bin[1], opts.z_bin[2])
    n_mass_bins = comp.num_bins(opts.mass_bin[0], opts.mass_bin[1],
                                opts.mass_bin[2])
    n_proxy_bins = comp.num_bins(opts.proxy_bin[0], opts.proxy_bin[1],
                                 opts.proxy_bin[2])

    z_index = np.floor((mock.z - opts.z_bin[0]) / opts.z_bin[2]).astype('int')
    mass_index = np.floor((mock.mass - opts.mass_bin[0]) /
                          opts.mass_bin[2]).astype('int')
    proxy_index = np.floor((obs.proxy - opts.proxy_bin[0]) /
                           opts.proxy_bin[2]).astype('int')

    z_x = comp.x_vals(n_z_bins, opts.z_bin[0], opts.z_bin[2])
    mass_x = comp.x_vals(n_mass_bins, opts.mass_bin[0], opts.mass_bin[2])
    proxy_x = comp.x_vals(n_proxy_bins, opts.proxy_bin[0], opts.proxy_bin[2])

    matrix = np.zeros((n_z_bins, n_mass_bins, n_proxy_bins))
    hm_matrix = np.zeros((n_z_bins, n_mass_bins))

    for i in range(len(mass_index)):
        hm_matrix[z_index[i], mass_index[i]] += 1

    for i in range(mock_index.size):
        if opts.unique:
            matrix[z_index[mock_index[i]], mass_index[mock_index[i]],
                   proxy_index[index[i]]] += 1
        else:
            for j in range(len(index[i])):
                matrix[z_index[mock_index[i]], mass_index[mock_index[i]],
                       proxy_index[index[i][j]]] += weights[i][j]

    return matrix, hm_matrix, np.array([z_x, mass_x, proxy_x])


##
#  Function to scale matrix lines to 1.0.
#
#  @param[in] matrix: Input matrix.
#
#  @return Scaled matrix.
#
def scale_matrix(matrix):

    return [comp.scale(matrix[i], 0.0, np.sum(matrix[i])) for i in
            range(len(matrix))]
