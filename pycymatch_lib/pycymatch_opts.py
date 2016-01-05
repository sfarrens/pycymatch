#  @file pycymatch_opts.py
#
#  PYCYMATCH OPTION METHODS
#
#  Functions for retreiving
#  code arguments.
#
#  @author Samuel Farrens
#  @version 1.0
#  @date 2015
#

import argparse


##
#  Function to read in code arguments.
#
#  @return List of arguments.
#
def get_opts():

    m_col_help = ('Set column numbers for mock properties: ID, CENTRAL_FLAG,' +
                  'RA, DEC, Z, N_GAL, MASS, R200, MIN_RA, MAX_RA, MIN_DEC,' +
                  'MAX_DEC. [Default: 1 2 3 4 5 6 7 8 9 10 11 12]')

    o_col_help = ('Set column numbers for cluster properties: ID, RA, DEC,' +
                  'Z, NGAl, SN. [Default 1 2 3 4 5 6]')

    parser = argparse.ArgumentParser('PYCYMATCH OPTIONS:')

    parser.add_argument('-i', '--input_files', dest='input_files', nargs='+',
                        help='Input file names. (1-Observerd Clusters, ' +
                        '2-Mock Haloes)')

    parser.add_argument('-z', '--delta_z', dest='delta_z', default=0.03,
                        type=float, help='Photometric redshift error value. ' +
                        '[Default: 0.03]')

    parser.add_argument('--proxy_bin', dest='proxy_bin',
                        default=[0.0, 3.5, 0.2], nargs='+', type=float,
                        help='Mass proxy bin values: Min proxy, Max proxy, '+
                        'proxy bin size [Default: 0.0 3.5 0.2]')

    parser.add_argument('--ngal_bin', dest='ngal_bin', default=[0.0, 3.5, 0.2],
                        nargs='+', type=float, help='Halo Ngal bin values: ' +
                        'Min Ngal, Max Ngal, Ngal bin size ' +
                        '[Default: 0.0 3.5 0.2]')

    parser.add_argument('--mass_bin', dest='mass_bin',
                        default=[13.0, 15.5, 0.1], nargs='+', type=float,
                        help='Mass bin values: Min mass, Max mass, mass ' +
                        'bin size [Default: 13.0 15.5 0.1]')

    parser.add_argument('--z_bin', dest='z_bin', default=[0.0, 3.0, 0.1],
                        nargs='+', type=float, help='Redshift bin ' +
                        'values: Min z, Max z, z bin size ' +
                        '[Default: 0.0 3.0 0.1]')

    parser.add_argument('--mock_plots', action='store_true',
                        dest='make_mock_plots', help='Output plots of ' +
                        'the mock catalogue properties.')

    parser.add_argument('-u', '--unique', action='store_true', dest='unique',
                        help='Only allow unique matches.')

    parser.add_argument('--mock_cols', dest='mock_cols', default=range(1, 13),
                        nargs='+', type=int, help=m_col_help)

    parser.add_argument('--obs_cols', dest='obs_cols', default=range(1, 7),
                        nargs='+', type=int, help=o_col_help)

    check_opts(parser)

    return parser


##
#  Function to read in code arguments.
#
#  @param[in] parser: List of arguments.
#
def check_opts(parser):

    opts = parser.parse_args()

    if not opts.input_files:
        parser.error('argument --input_files: file names not provided')

    if not len(opts.input_files) == 2:
        parser.error('argument --input_files: requires 2 input file names\n' +
                     'e.g. --input_files obs_cluster_file mock_halo_file')

    if not len(opts.proxy_bin) == 3:
        parser.error('argument --proxy_bin: requires 3 input values\n' +
                     'e.g. --proxy_bin min_value max_value bin_size')

    if not len(opts.mass_bin) == 3:
        parser.error('argument --mass_bin: requires 3 input values\n' +
                     'e.g. --mass_bin min_value max_value bin_size')

    if not len(opts.z_bin) == 3:
        parser.error('argument --z_bin: requires 3 input values\n' +
                     'e.g. --z_bin min_value max_value bin_size')

    if not len(opts.mock_cols) == 12:
        parser.error('argument --mock_cols: requires 12 input values')

    if not len(opts.obs_cols) == 6:
        parser.error('argument --obs_cols: requires 6 input values')
