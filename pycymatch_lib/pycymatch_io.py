#  @file pycymatch_IO.py
#
#  PYCYMATCH IO METHODS
#
#  Functions for input
#  and output operations.
#
#  @author Samuel Farrens
#  @version 1.1
#  @date 2015
#

import numpy as np
from numpy.lib.recfunctions import append_fields


##
#  Function to read mock halo catalogue file.
#
#  @param[in] opts: List of arguments.
#
#  @return List of mock haloes.
#
def read_mock(opts):

    print 'Reading file:', opts.input_files[1],

    # read file
    mock = np.genfromtxt(opts.input_files[1], dtype='S',
                         unpack=True,
                         usecols=np.array(opts.mock_cols) - 1)

    # format array
    dtypes = [('id', 'S22'), ('cen', 'i4'), ('ra', 'f8'), ('dec', 'f8'),
              ('z', 'f8'), ('ngal', 'f8'), ('mass', 'f8'), ('r200', 'f8'),
              ('minra', 'f8'), ('maxra', 'f8'), ('mindec', 'f8'),
              ('maxdec', 'f8')]
    mock = np.core.records.fromarrays(mock, dtype=dtypes)

    # convert ngal to log_10(ngal)
    mock.ngal = np.log10(mock.ngal)

    # restrict elements to limits
    index = ((mock.z >= opts.z_bin[0]) & (mock.z <= opts.z_bin[1]) &
             (mock.mass >= opts.mass_bin[0]) & (mock.mass <= opts.mass_bin[1]))

    mock = mock[index]

    # sort by Ngal
    mock = mock[mock.argsort(order=('ngal', 'id'))[::-1]]

    print '\tComplete:', len(mock)

    return mock


##
#  Function to read observed cluster catalogue file.
#
#  @param[in] opts: List of arguments.
#
#  @return List of observed clusters.
#
def read_obs(opts):

    print 'Reading file:', opts.input_files[0],

    # read file
    obs = np.genfromtxt(opts.input_files[0], dtype='S',
                        unpack=True,
                        usecols=np.array(opts.obs_cols) - 1)

    dtypes = [('id', 'S22'), ('ra', 'f8'), ('dec', 'f8'), ('z', 'f8'),
              ('proxy', 'f8'), ('sn', 'f8')]

    obs = np.core.records.fromarrays(obs, dtype=dtypes)

    # convert proxy to log_10(proxy)
    obs.proxy = np.log10(obs.proxy)

    # restrict elements to limits
    index = ((obs.z >= opts.z_bin[0]) & (obs.z <= opts.z_bin[1]) &
             (obs.proxy >= opts.proxy_bin[0]) &
             (obs.proxy <= opts.proxy_bin[1]))

    obs = obs[index]

    # sort by mass proxy
    obs = obs[obs.argsort(order=('proxy', 'id'))[::-1]]

    ###
    # add field for flags
    obs = append_fields(obs, ['flag'], [np.zeros(obs.size)], dtypes=['int'],
                        asrecarray=True, usemask=False)
    ###

    print '\tComplete:', len(obs)

    return obs


##
#  Function to output mass-observable matrix.
#
#  @param[in] matrix: Mass-observable matrix.
#  @param[in] hm_matrix: Mass vs z matrix.
#  @param[in] ranges: Matrix ranges.
#  @param[in] opts: List of arguments.
#
def print_matrix(matrix, hm_matrix, ranges, opts):

    for i in range(len(matrix)):

        file_name = opts.input_files[0] + '.pycy.matrix.z' + \
                    str(ranges[0][i]) + '.txt'

        output = np.vstack([np.append(-1, hm_matrix[i]),
                            np.append(-1, np.sum(matrix[i], axis=1)),
                            np.transpose(np.vstack([ranges[2], matrix[i]]))])

        header = 'Lambda[1] '
        for j in range(len(ranges[1])):
            header += 'M=' + str(ranges[1][j]) + '[' + str(j + 2) + '] '

        np.savetxt(file_name, output, fmt='%.3f', header=header)

        print 'Data saved to:', file_name


##
#  Function to print matches between observed
#  clusters and mock haloes.
#
#  @param[in] mock: List of mock haloes.
#  @param[in] obs: List of observed clusters.
#  @param[in] index: List of indices of matches.
#  @param[in] weights: List of match weights.
#  @param[in] opts: List of arguments.
#
def print_matches(mock, obs, index, weights, opts):

    file_name = opts.input_files[0] + '.' + str(opts.z_bin[0]) + 'to' + \
                str(opts.z_bin[1]) + '.pycy.matches.txt'
    output = open(file_name, 'w')

    print>> output, '#MOCK_ID MOCK_Mass Mock_z [OBS_ID OBS_PROXY OBS_z' + \
                    'WEIGHT]...'

    for i in range(mock.size):
        print>> output, mock.id[i], mock.mass[i], mock.z[i],

        if isinstance(index[i], np.ndarray):
            for j in range(len(index[i])):
                    print>> output, [obs.id[index[i][j]],
                                     int(10 ** obs.proxy[index[i][j]]),
                                     obs.z[index[i][j]], weights[i][j]],
            print>> output, ''
        else:
            print>> output, [obs.id[index[i]], int(10 ** obs.proxy[index[i]]),
                             obs.z[index[i]], weights[i]]

    output.close()

    print 'Data saved to:', file_name
