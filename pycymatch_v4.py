#! /usr/bin/python

#######################
#   PYCYMATCH V.4.0   #
#######################
# Samuel Farrens 2014 #
#######################

import sys, getopt, math, optparse, warnings, pycymatch_help as pch
import os.path
from functions import astro, comp, errors
import numpy as np, matplotlib as mp, scipy.stats as ss
from numpy.lib.recfunctions import append_fields
mp.use('pdf')
import matplotlib.pyplot as plt, matplotlib.gridspec as gridspec
from matplotlib import cm as cm
from matplotlib import colors
from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter

###################
# IGNORE WARNINGS #
###################

warnings.simplefilter('ignore')

##################
# READ ARGUMENTS #
##################

def opt_callback(option, opt, value, parser):
    setattr(parser.values, option.dest, np.array(map(int, value.split(','))) - 1)

m_col_help = ("Set column numbers for mock properties: ID, CENTRAL_FLAG, RA, DEC, Z, N_GAL, MASS,"
+ "R200, MIN_RA, MAX_RA, MIN_DEC, MAX_DEC. [Default 1,2,3,4,5,6,7,8,9,10,11,12]")
c_col_help = ("Set column numbers for cluster properties: ID, RA, DEC, Z, NGAl, SN. [Default 1,2,3,4,5,6]")

parser = optparse.OptionParser()
parser.add_option("-c", "--cluster_file", dest = "cluster_file", help = "Input cluster file name.")
parser.add_option("-m", "--mock_file", dest = "mock_file", help = "Input mock file name.")
parser.add_option("-z", "--delta_z", dest = "delta_z", default = 0.03, type = "float",
                   help = "Delta_z value. [Default 0.03]")
parser.add_option("--min_rich_mock", dest = "min_rich_mock", default = 0.0, type = "float",
                  help = "Minimum mock richness value allowed. [Default 0.0]")
parser.add_option("--max_rich_mock", dest = "max_rich_mock", default = -1, type = "float",
                  help = "Maximum mock richness value allowed. [Default maximum value in list]")
parser.add_option("--min_rich_clt", dest = "min_rich_clt", default = 0.0, type = "float",
                  help = "Minimum cluster richness value allowed. [Default 0.0]")
parser.add_option("--max_rich_clt", dest = "max_rich_clt", default = -1, type = "float",
                  help = "Maximum cluster richness value allowed. [Default maximum value in list]")
parser.add_option("--rich_bin_size", dest = "rich_bin_size", default = 0.2, type = "float",
                  help = "Size of richness bins. [Default 0.2]")
parser.add_option("--min_mass_bin", dest = "min_mass_bin", default = 13.0, type = "float",
                  help = "Minimum mass value allowed. [Default 13.0]")
parser.add_option("--max_mass_bin", dest = "max_mass_bin", default = -1, type = "float",
                  help = "Maximum mass value allowed. [Default maximum value in list]")
parser.add_option("--mass_bin_size", dest = "mass_bin_size", default = 0.1, type = "float",
                   help = "Size of mass bins. [Default 0.1]")
parser.add_option("--min_z_bin", dest = "min_z_bin", default = 0.0, type = "float",
                  help = "Minimum redshift value allowed. [Default 0.0]")
parser.add_option("--max_z_bin", dest = "max_z_bin", default = -1, type = "float",
                  help = "Maximum redshift value allowed. [Default maximum value in list]")
parser.add_option("--z_bin_size", dest = "z_bin_size", default = 0.1, type = "float",
                  help = "Size of redshift bins. [Default 0.1]")
parser.add_option("--mock_plots", action = "store_true", dest = "make_mock_plots",
                  help = "Output plots of the mock catalogue properties." )
parser.add_option("--mock_cols", dest = "m_cols", type = 'string', action = 'callback',
                  callback = opt_callback, default = range(0, 12),
                  help = m_col_help)
parser.add_option("--cluster_cols", dest = "c_cols", type = 'string', action = 'callback',
                  callback = opt_callback, default = range(0, 6),
                  help = c_col_help)

(opts, args) = parser.parse_args()
    
if not opts.cluster_file:
    parser.error('Cluster file name not provided.')
if not opts.mock_file:
    parser.error('Mock file name not provided.')
    
#############
# READ DATA #
#############

#make sure files exit
errors.file_name_error(opts.mock_file)
errors.file_name_error(opts.cluster_file)

print '================================================================='

#read mock catalogue  
print 'Reading file:', opts.mock_file,        
mock = np.genfromtxt(opts.mock_file, dtype="S", unpack = True, usecols = np.array(opts.m_cols))
dtypes = [('id', 'S23'), ('cen', 'i4'), ('ra', 'f8'), ('dec', 'f8'),
          ('z', 'f8'), ('rich', 'f8'), ('mass', 'f8'), ('r200', 'f8'),
          ('minra', 'f8'), ('maxra', 'f8'), ('mindec', 'f8'), ('maxdec', 'f8')]
mock = np.core.records.fromarrays(mock, dtype = dtypes)
mock.rich = np.log10(mock.rich)
print '\tComplete:', len(mock)

#read cluster catalogue
print 'Reading file:', opts.cluster_file,        
cluster = np.genfromtxt(opts.cluster_file, dtype="S", unpack = True, usecols = np.array(opts.c_cols))

dtypes = [('id', 'S10'), ('ra', 'f8'), ('dec', 'f8'), ('z', 'f8'), ('rich', 'f8'), ('sn', 'f8')]
cluster = np.core.records.fromarrays(cluster, dtype = dtypes)
cluster.rich = np.log10(cluster.rich)
print '\tComplete:', len(cluster)

#######################
# SET BOUNDARY LIMITS #
#######################

if opts.max_rich_mock < 0: opts.max_rich_mock = max(mock.rich)
if opts.max_rich_clt < 0: opts.max_rich_clt = max(cluster.rich)
if opts.max_mass_bin < 0: opts.max_mass_bin = max(mock.mass)
if opts.max_z_bin < 0: opts.max_z_bin = max(max(cluster.z), max(mock.z))

index = ((mock.z >= opts.min_z_bin) & (mock.z  <= opts.max_z_bin) &
                 (mock.rich >= opts.min_rich_mock) & (mock.rich <= opts.max_rich_mock) &
                 (mock.mass >= opts.min_mass_bin) & (mock.mass <= opts.max_mass_bin))
mock = mock[index]
mock = append_fields(mock, ['dist'], [np.zeros(mock.size)], dtypes = ['f8'], asrecarray=True, usemask = False)
  
index = ((cluster.z >= opts.min_z_bin) & (cluster.z <= opts.max_z_bin) &
         (cluster.rich >= opts.min_rich_clt) & (cluster.rich <= opts.max_rich_clt))
cluster = cluster[index]

print '================================================================='
print 'Used Mock Haloes:', len(mock)
print 'Used Observed Clusters:', len(cluster)
print '================================================================='

####################
# SORT BY RICHNESS #
####################

mock = mock[mock.argsort(order = ('rich', 'id'))[::-1]]
cluster = cluster[cluster.argsort(order = ('rich', 'id'))[::-1]]

###################
# BIN BY RICHNESS #
###################

#set number of richness bins
n_rich_bins_mock = int(math.floor((opts.max_rich_mock - opts.min_rich_mock) / opts.rich_bin_size)) + 1
n_rich_bins_clt = int(math.floor((opts.max_rich_clt - opts.min_rich_clt) / opts.rich_bin_size)) + 1

c_rich_bin_index = np.floor((cluster.rich - opts.min_rich_clt) / opts.rich_bin_size).astype('int')
m_rich_bin_index = np.floor((mock.rich - opts.min_rich_mock) / opts.rich_bin_size).astype('int')

x_rich_vals_mock = (np.arange(n_rich_bins_mock) + 0.5) * opts.rich_bin_size + opts.min_rich_mock
x_rich_vals_clt = (np.arange(n_rich_bins_clt) + 0.5) * opts.rich_bin_size + opts.min_rich_clt

############
# BIN BY Z #
############

#set number of redshift bins
n_z_bins = int(math.floor((opts.max_z_bin - opts.min_z_bin) / opts.z_bin_size)) + 1 

c_z_bin_index = np.floor((cluster.z - opts.min_z_bin) / opts.z_bin_size).astype('int')
m_z_bin_index = np.floor((mock.z - opts.min_z_bin) / opts.z_bin_size).astype('int')

x_z_vals = (np.arange(n_z_bins) + 0.5) * opts.z_bin_size + opts.min_z_bin

################
# FIND MATCHES #
################

z_factor = 2.0 #to adjust the line-of-sight matching

c_bin_count = np.zeros((n_rich_bins_clt, n_z_bins)).astype('int')
m_bin_count = np.zeros((n_rich_bins_mock, n_z_bins)).astype('int')
c_match_bin_count = np.zeros((n_rich_bins_clt, n_z_bins)).astype('int')
m_match_bin_count = np.zeros((n_rich_bins_mock, n_z_bins)).astype('int')

m_match_flag = np.zeros(mock.size).astype('int') - 1
c_match_flag = np.zeros(cluster.size).astype('int') - 1

for i in range(cluster.size):                                           
    c_bin_count[c_rich_bin_index[i], c_z_bin_index[i]] += 1
for i in range(mock.size):    
    m_bin_count[m_rich_bin_index[i], m_z_bin_index[i]] += 1    
    if m_match_flag[i] == -1:
        #Find clusters that match primary matching conditions.
        index1 = np.where((c_match_flag == -1) &
                          (np.fabs(mock.z[i] - cluster.z) <= (z_factor * opts.delta_z * (1 + mock.z[i]))) &
                        (cluster.ra >= mock.minra[i]) & (cluster.ra <= mock.maxra[i]) &
                        (cluster.dec >= mock.mindec[i]) & (cluster.dec <= mock.maxdec[i]))[0]   
        if len(index1) > 0:       
            #Calculate the projected distance to halo centre for these clusters.
            dists = []
            for j in index1:
                dist = astro.ang_sep((mock.ra[i],  mock.dec[i]), (cluster.ra[j],cluster.dec[j])) * 60.0
                dists.extend([dist])
            dists = np.array(dists)
            #Find clusters within r200.            
            index2 = np.where(dists <= mock.r200[i])[0]
            if len(index2) > 0:
                dists = dists[index2]
                index1 = index1[index2]
                #Check if any of these have the same N_gal value.
                index3 = np.where(cluster.rich[index1] == cluster.rich[index1[0]])[0]
                if len(index3) > 0:
                    dists = dists[index3]
                    index1 = index1[index3]
                    #Choose the cluster closest to the halo.
                    index4 = np.argmin(dists)
                    index1 = index1[index4]
                    mock.dist[i] = dists[index4]
                    #Increase match count, and tag cluster and halo.     
                    c_match_bin_count[c_rich_bin_index[index1], c_z_bin_index[index1]] += 1
                    m_match_bin_count[m_rich_bin_index[i], m_z_bin_index[i]] += 1               
                    m_match_flag[i] = index1
                    c_match_flag[index1] = i
                    
pure_bin = np.flipud(np.array(c_match_bin_count).astype("float") / np.array(c_bin_count).astype("float"))
complete_bin = np.flipud(np.array(m_match_bin_count).astype("float") / np.array(m_bin_count).astype("float"))

## g_bin = 1 - np.sqrt(((1 - np.nan_to_num(pure_bin)) ** 2 + (1 - np.nan_to_num(complete_bin)) ** 2) / 2)

## g_bin_1 = np.zeros(n_z_bins).astype('float')
## for i in range(n_z_bins):
##     p = np.nansum(pure_bin[:, i]) / len(pure_bin[:, i])
##     c = np.nansum(complete_bin[:, i]) / len(complete_bin[:, i])
##     g_bin_1[i] = 1 - math.sqrt(((1 - p) ** 2 + (1 - c) ** 2) / 2)
## comp.nan2zero(g_bin_1)

## g_bin_2 = np.zeros(n_rich_bins_clt).astype('float')
## for i in range(n_rich_bins_clt):
##     p = np.nansum(pure_bin[i, :]) / len(pure_bin[i, :])
##     c = np.nansum(complete_bin[i, :]) / len(complete_bin[i, :])
##     g_bin_2[i] = 1 - math.sqrt(((1 - p) ** 2 + (1 - c) ** 2) / 2)
## comp.nan2zero(g_bin_2)
    
match_index = np.where(m_match_flag > -1)[0]

######################
# PRINT HALO MATCHES #
######################
   
match_out_file = opts.cluster_file + '_halo_matching.txt'                                          
match_out = open(match_out_file,'w')
      
print>> match_out, '#H_ID[1]               H_RA[2] H_DEC[3] H_Z[4] H_NGAL[5] H_MASS[6] H_R200[7] H_CEN[8] MATCH[9]',
print>> match_out, 'DIST[10] C_ID[11] C_RA[12] C_DEC[13] C_Z[14] C_NGAL[15] C_S/N[16]'

for i in range(mock.size):
    print>> match_out, '%-22s' % mock.id[i],'%6.3f' % mock.ra[i], '%+7.3f' % mock.dec[i], '%7.3f' % mock.z[i],
    print>> match_out, '%5i' % (10 ** mock.rich[i]), '%11.3f' % mock.mass[i], '%9.3f' % mock.r200[i], '  %2i' % mock.cen[i],
    if m_match_flag[i] > -1:
        print>> match_out, '       Y ', '%11.3f' % mock.dist[i], '   %-8s' % cluster.id[m_match_flag[i]],
        print>> match_out, '%6.3f' % cluster.ra[m_match_flag[i]], '%+8.3f' % cluster.dec[m_match_flag[i]],
        print>> match_out, '%8.3f' % cluster.z[m_match_flag[i]], '%6i' % int(round((10 ** cluster.rich[m_match_flag[i]]))),
        print>> match_out, '%14.3f' % cluster.sn[m_match_flag[i]] 
    else:
        print>> match_out, '       N'

print "Matching data saved to:", match_out_file

#########################
# PRINT CLUSTER MATCHES #
#########################

c_match_out_file = opts.cluster_file + '_cluster_matching.txt'                                          
c_match_out = open(c_match_out_file,'w')
      
print>> c_match_out, '#C_ID[1] C_RA[2] C_DEC[3] C_Z[4] C_NGAL[5] C_S/N[6] MATCH[7] DIST[8]',
print>> c_match_out, 'H_ID[9]                H_RA[10] H_DEC[11] H_Z[12] H_NGAL[13] H_MASS[14] H_R200[15] H_CEN[16]'

for i in range(cluster.size):
    print>> c_match_out, '%-8s' % cluster.id[i], '%6.3f' % cluster.ra[i], '%+7.3f' % cluster.dec[i], '%7.3f' % cluster.z[i],
    print>> c_match_out, '%5i' % int(round((10.0 ** cluster.rich[i]))), '%13.3f' % cluster.sn[i],
    if c_match_flag[i] > -1:
        print>> c_match_out, 'Y ', '%11.3f' % mock.dist[c_match_flag[i]], '  %-22s' % mock.id[c_match_flag[i]],
        print>> c_match_out, '%6.3f' % mock.ra[c_match_flag[i]], '%+8.3f' % mock.dec[c_match_flag[i]],
        print>> c_match_out, '%8.3f' % mock.z[c_match_flag[i]], '%6i  ' % int(round((10.0 ** mock.rich[c_match_flag[i]]))),
        print>> c_match_out, '%10.3f' % mock.mass[c_match_flag[i]], '%10.3f' % mock.r200[c_match_flag[i]],
        print>> c_match_out, '%5i' % mock.cen[c_match_flag[i]]
    else:
        print>> c_match_out, 'N'

print "Matching data saved to:", c_match_out_file
  
###############
# BIN BY MASS #
###############

#set number of mass bins
n_mass_bins = int(math.floor((opts.max_mass_bin - opts.min_mass_bin) / opts.mass_bin_size)) + 1

m_mass_bin_index = np.floor((mock.mass - opts.min_mass_bin) / opts.mass_bin_size).astype('int')

y_mass_vals = (np.arange(n_mass_bins) + 0.5) * opts.mass_bin_size + opts.min_mass_bin

m_bin_count2 = np.zeros((n_mass_bins, n_z_bins)).astype('int')
m_match_bin_count = np.zeros((n_mass_bins, n_z_bins)).astype('int')

for i in range(mock.size):
    m_bin_count2[m_mass_bin_index[i], m_z_bin_index[i]] += 1
    if m_match_flag[i] > -1:
        m_match_bin_count[m_mass_bin_index[i], m_z_bin_index[i]] += 1               
               
complete_bin_mass = np.flipud(np.array(m_match_bin_count).astype("float") / np.array(m_bin_count2).astype("float"))
     
################
# SPEARMAN RHO #
################
   
for i in range(n_z_bins):
    m_index_1 = mock.rich[match_index][m_z_bin_index[match_index] == i]
    m_index_2 = mock.rich[match_index][c_z_bin_index[m_match_flag[match_index]] == i]
    c_index_1 = cluster.rich[m_match_flag[match_index]][m_z_bin_index[match_index] == i]
    c_index_2 = cluster.rich[m_match_flag[match_index]][c_z_bin_index[m_match_flag[match_index]] == i]
    if m_index_1.size > 1:
        m_rho, p1 = ss.spearmanr(m_index_1, c_index_1)
        m_rho_err = 0.6325 / (len(m_index_1) - 1) ** 0.5
    else:
        m_rho = 0.0
        m_rho_err = 0.0
    if m_index_2.size > 1:
        c_rho, p2 = ss.spearmanr(m_index_2, c_index_2)
        c_rho_err = 0.6325 / (len(m_index_2) - 1) ** 0.5
    else:
        c_rho = 0.0
        c_rho_err = 0.0
    if i == 0:
        m_rhos = np.array(m_rho)
        m_rhos_err = np.array(m_rho_err)
        c_rhos = np.array(c_rho)
        c_rhos_err = np.array(c_rho_err)
    else:
        m_rhos = np.hstack((m_rhos, m_rho))
        m_rhos_err = np.hstack((m_rhos_err, m_rho_err))
        c_rhos = np.hstack((c_rhos, c_rho))
        c_rhos_err = np.hstack((c_rhos_err, c_rho_err))

#######
# S/N #
#######

euclid_bg = 30 #30 galaxies per square arcminute

m_sn_bin_sum = np.zeros((n_mass_bins, n_z_bins)).astype('float')

c_sn_bin_count = np.zeros((n_mass_bins, n_z_bins)).astype('float')
c_sn_bin_sum = np.zeros((n_mass_bins, n_z_bins)).astype('float')

mock_sn = (10.0 ** mock.rich) / np.sqrt(math.pi * mock.r200 ** 2 * euclid_bg)

for i in range(len(mock_sn)):
    m_sn_bin_sum[m_mass_bin_index[i], m_z_bin_index[i]] += mock_sn[i]

for i in range(len(cluster.sn)):
    if c_match_flag[i] > -1 and cluster.sn[i] <= 10 :
        c_sn_bin_count[m_mass_bin_index[c_match_flag[i]], m_z_bin_index[c_match_flag[i]]] += 1.0
        c_sn_bin_sum[m_mass_bin_index[c_match_flag[i]], m_z_bin_index[c_match_flag[i]]] += cluster.sn[i] - 1.0

m_sn_bin_val = np.flipud(m_sn_bin_sum / np.array(m_bin_count2).astype('float'))
c_sn_bin_val = np.flipud(c_sn_bin_sum / np.array(c_sn_bin_count).astype('float'))

if os.path.isfile("selfunct_NIPwide.dat"):
    biviano_lines = np.genfromtxt("selfunct_NIPwide.dat", unpack = True, dtype="S")
else:
    biviano_lines = 0;

##############
# MAKE PLOTS #
##############

interpolation = 'spline16'
colormap = cm.jet
boundaries = np.arange(0.0, 1.1, 0.1)
normalisation = colors.BoundaryNorm(boundaries, colormap.N)
   
#Completeness vs. z and N_true_mock
fig = plt.figure()
fig.subplots_adjust(wspace = .01)
gs = gridspec.GridSpec(1, 2, width_ratios=[12,1])
ax = fig.add_subplot(gs[0])
cax = fig.add_subplot(gs[1])
im = ax.imshow(complete_bin, extent = [opts.min_z_bin, opts.max_z_bin, opts.min_rich_mock, opts.max_rich_mock], aspect = 'auto',
                 interpolation = interpolation, cmap = colormap, norm = normalisation)
ax.set_xlabel('z', fontsize = 24)
ax.set_ylabel('log$_{10}$ N$_{true}$', fontsize = 24)
ax.set_ylim(0.0, 3.0)
cbar = fig.colorbar(im, cax=cax)
cbar.set_label('Completeness', fontsize = 20)
fig_name = opts.cluster_file + '_completeness_plot.pdf'
fig.savefig(fig_name)
print "Plots saved to:", fig_name

#Purity vs. z and N_obs
fig = plt.figure()
fig.subplots_adjust(wspace = .01)
gs = gridspec.GridSpec(1, 2, width_ratios=[12,1])
ax = fig.add_subplot(gs[0])
cax = fig.add_subplot(gs[1])
im = ax.imshow(pure_bin, extent = [opts.min_z_bin, opts.max_z_bin, opts.min_rich_clt, opts.max_rich_clt], aspect = 'auto',
                 interpolation = interpolation, cmap = colormap, norm = normalisation)
ax.set_xlabel('z', fontsize = 24)
ax.set_ylabel('log$_{10}$ N$_{obs}$', fontsize = 24)
cbar = fig.colorbar(im, cax=cax)
cbar.set_label('Purity', fontsize = 20)
fig_name = opts.cluster_file + '_purity_plot.pdf'
fig.savefig(fig_name)
print "Plots saved to:", fig_name

## #G vs. z
## fig = plt.figure()
## gs = gridspec.GridSpec(1, 1)
## ax = fig.add_subplot(gs[0])
## ax.plot(x_z_vals, g_bin_1, linewidth=1.0, c='b')
## ax.set_ylim([0.0, 1.0])
## ax.set_xlabel('z', fontsize = 24)
## ax.set_ylabel('G', fontsize = 24)
## ax.set_title('Average G vs. z')
## fig_name = opts.cluster_file + '_gvz.pdf'
## fig.savefig(fig_name)
## print "Plots saved to:", fig_name

## #G vs. rich
## fig = plt.figure()
## gs = gridspec.GridSpec(1, 1)
## ax = fig.add_subplot(gs[0])
## ax.plot(x_rich_vals_clt, g_bin_2, linewidth=1.0, c='b')
## ax.set_ylim([0.0, 1.0])
## ax.set_xlabel('log$_{10}$ N$_{gal}$', fontsize = 24)
## ax.set_ylabel('G', fontsize = 24)
## ax.set_title('Average G vs. log$_{10}$ N$_{gal}$')
## fig_name = opts.cluster_file + '_gvn.pdf'
## fig.savefig(fig_name)
## print "Plots saved to:", fig_name

## #G vs. z and N_obs
## fig = plt.figure()
## fig.subplots_adjust(wspace = .01)
## gs = gridspec.GridSpec(1, 2, width_ratios=[12,1])
## ax = fig.add_subplot(gs[0])
## cax = fig.add_subplot(gs[1])
## im = ax.imshow(g_bin, extent = [opts.min_z_bin, opts.max_z_bin, opts.min_rich_clt, opts.max_rich_clt], aspect = 'auto',
##                  interpolation = interpolation, cmap = colormap, norm = normalisation)
## ax.set_xlabel('z', fontsize = 24)
## ax.set_ylabel('log$_{10}$ N$_{gal}$', fontsize = 24)
## cbar = fig.colorbar(im, cax=cax)
## cbar.set_label('G')
## fig_name = opts.cluster_file + '_g_plot.pdf'
## fig.savefig(fig_name)
## print "Plots saved to:", fig_name
   
#z vs. spearman rho   
fig = plt.figure()
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])
ax.plot(x_z_vals, m_rhos, linewidth=1.0, c='b', label='z$_{mock}$ Bins')
ax.plot(x_z_vals, c_rhos, linewidth=1.0, linestyle='dashed', c='r', label='z$_{obs}$ Bins')
ax.errorbar(x_z_vals, m_rhos, yerr=m_rhos_err, fmt='bx')
ax.errorbar(x_z_vals, c_rhos, yerr=c_rhos_err, fmt='rx')
ax.set_ylim(0.0, 1.0)
ax.set_xlabel('z', fontsize = 24)
ax.set_ylabel(r'$\rho$', fontsize = 24)
ax.set_title('Spearman Rank Order Coeffcient')
ax.legend()
fig_name = opts.cluster_file + '_spearman.pdf'
fig.savefig(fig_name)
print "Plots saved to:", fig_name
   
#Completeness vs. z and Mass
fig = plt.figure()
fig.subplots_adjust(wspace = .01)
gs = gridspec.GridSpec(1, 2, width_ratios=[12,1])
ax = fig.add_subplot(gs[0])
cax = fig.add_subplot(gs[1])
im = ax.imshow(complete_bin_mass, extent = [opts.min_z_bin, opts.max_z_bin, opts.min_mass_bin, opts.max_mass_bin], aspect = 'auto',
                 interpolation=interpolation, cmap = colormap, norm = normalisation)
ax.set_xlabel('z', fontsize = 24)
ax.set_ylabel('log$_{10}$ M', fontsize = 24)
cbar = fig.colorbar(im, cax=cax)
cbar.set_label('Completeness', fontsize = 20)
fig_name = opts.cluster_file + '_completeness_mass_plot.pdf'
fig.savefig(fig_name)
print "Plots saved to:", fig_name

#Signal-to-noise vs. z and N_obs
sn_bounds = np.arange(0.0, 10.0, 1.0)
sn_norm = colors.BoundaryNorm(sn_bounds, colormap.N)
fig = plt.figure()
fig.subplots_adjust(wspace = .01)
gs = gridspec.GridSpec(1, 2, width_ratios=[12,1])
ax = fig.add_subplot(gs[0])
cax = fig.add_subplot(gs[1])
little_h_factor = math.log10(0.7)
mass_min_limit = opts.min_mass_bin - little_h_factor
mass_max_limit = opts.max_mass_bin - little_h_factor
im = ax.imshow(c_sn_bin_val, extent = [opts.min_z_bin, opts.max_z_bin, mass_min_limit, mass_max_limit], aspect = 'auto',
               interpolation = interpolation, cmap = colormap, norm = sn_norm)
if biviano_lines != 0:
    ax.plot(biviano_lines[0], biviano_lines[1], linewidth=1.0, c = 'k', ls = '--', label = r'3 $\sigma$')
    ax.plot(biviano_lines[0], biviano_lines[2], linewidth=1.0, c = 'k', ls = '-', label = r'5 $\sigma$')
ax.set_xlabel('z', fontsize = 24)
ax.set_ylabel('log$_{10}$ M', fontsize = 24)
ax.legend()
cbar = fig.colorbar(im, cax=cax)
cbar.set_label('S/N', fontsize = 20)
fig_name = opts.cluster_file + '_sn_plot.pdf'
fig.savefig(fig_name)
print "Plots saved to:", fig_name

#MOCK PLOTS

if opts.make_mock_plots:

    #Mock signal-to-noise vs. z and Mass
    fig = plt.figure()
    fig.subplots_adjust(wspace = .01)
    gs = gridspec.GridSpec(1, 2, width_ratios=[12,1])
    ax = fig.add_subplot(gs[0])
    cax = fig.add_subplot(gs[1])
    im = ax.imshow(m_sn_bin_val, extent = [opts.min_z_bin, opts.max_z_bin, opts.min_mass_bin, opts.max_mass_bin],
                   aspect = 'auto', interpolation = interpolation, cmap = colormap)
    ax.set_xlabel('z', fontsize = 24)
    ax.set_ylabel('log$_{10}$ M', fontsize = 24)
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label('S/N', fontsize = 20)
    fig_name = opts.cluster_file + '_mock_sn_plot.pdf'
    fig.savefig(fig_name)
    print "Plots saved to:", fig_name

    #Mass vs. N_true_mock
    heat, yedges, xedges = np.histogram2d(mock.rich, mock.mass, bins=(n_rich_bins_mock, n_mass_bins))
    heat = np.log10(heat)   
    fig3 = plt.figure()
    fig3.subplots_adjust(wspace = .01)
    gs3 = gridspec.GridSpec(1, 2, width_ratios=[12,1])
    ax31 = fig3.add_subplot(gs3[0])
    cax31 = fig3.add_subplot(gs3[1])
    im = ax31.imshow(np.flipud(heat), extent = [opts.min_mass_bin, opts.max_mass_bin, opts.min_rich_mock, opts.max_rich_mock], aspect = 'auto',
                    interpolation = interpolation, cmap = colormap)
    ax31.set_xlabel('log$_{10}$ M', fontsize = 24)
    ax31.set_ylabel('log$_{10}$ N$_{true}$', fontsize = 24)
    cbar = fig3.colorbar(im, cax=cax31)
    cbar.set_label('log$_{10}$ N', fontsize = 20)
    fig3_name = opts.cluster_file + '_mock_plot1.pdf'
    fig3.savefig(fig3_name)
    print "Plots saved to:", fig3_name

    #z vs. N_true_mock
    heat, yedges, xedges = np.histogram2d(mock.rich, mock.z, bins=(n_rich_bins_mock, n_z_bins))
    heat = np.log10(heat)   
    fig4 = plt.figure()
    fig4.subplots_adjust(wspace = .01)
    gs4 = gridspec.GridSpec(1, 2, width_ratios=[12,1])
    ax41 = fig4.add_subplot(gs4[0])
    cax41 = fig4.add_subplot(gs4[1])
    im = ax41.imshow(np.flipud(heat), extent = [opts.min_z_bin, opts.max_z_bin, opts.min_rich_mock, opts.max_rich_mock], aspect = 'auto',
                    interpolation = interpolation, cmap = colormap)
    ax41.set_xlabel('z', fontsize = 24)
    ax41.set_ylabel('log$_{10}$ N$_{true}$', fontsize = 24)
    cbar = fig4.colorbar(im, cax=cax41)
    cbar.set_label('log$_{10}$ N', fontsize = 20)
    fig4_name = opts.cluster_file + '_mock_plot2.pdf'
    fig4.savefig(fig4_name)
    print "Plots saved to:", fig4_name

#########################################################################
