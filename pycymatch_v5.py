#! /Users/sfarrens/Documents/Library/anaconda/bin/python

#  @file pycymatch_v5.py
#
#  PYCYMATCH
#
#  Script for performing
#  cylindrical matching between
#  observed and mock clusters.
#
#  @author Samuel Farrens
#  @version 5.1
#  @date 2015
#

import numpy as np

from functions import errors
from pycymatch_lib import *


##
#  Code Main.
def main():

    ##
    # Read arguments
    opts = pycymatch_opts.get_opts().parse_args()

    ##
    # Check input files
    for file in opts.input_files:
        errors.file_name_error(file)

    pycymatch_extra.h_line()

    ##
    # Read mock halo catalogue
    mock = pycymatch_io.read_mock(opts)

    ##
    # Read observed catalogue
    obs = pycymatch_io.read_obs(opts)

    pycymatch_extra.h_line()

    ##
    # Find Matches
    matches = pycymatch_match2.find_matches(mock, obs, 2.0 * opts.delta_z,
                                            opts)
    pycymatch_io.print_matches(mock[matches[2]], obs, matches[0], matches[1],
                               opts)

    ##
    # Define completeness and purity of sample
    c_matrix = pycymatch_cp.get_completeness(mock, matches, opts)
    p_matrix = pycymatch_cp.get_purity(obs, matches, opts)

    ##
    # Define mass-observable matrix for matched objects
    matrix, hm_matrix, ranges = pycymatch_match2.mo_matrix(mock, obs, matches,
                                                           opts)

    pycymatch_extra.h_line()

    ##
    # Make plots
    pycymatch_plot.make_plots(matrix, hm_matrix, ranges, opts)
    pycymatch_plot.plot_complete(c_matrix[0], 'mass', opts)
    pycymatch_plot.plot_complete(c_matrix[1], 'ngal', opts)
    pycymatch_plot.plot_pure(p_matrix, opts)

    ##
    # Save matrix to file
    pycymatch_io.print_matrix(matrix, hm_matrix, ranges, opts)

    pycymatch_extra.h_line()

if __name__ == "__main__":
    main()
