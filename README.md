# PYCYMATCH

> Author: **Samuel Farrens**

> Year: **2016**

> Version: **4.0, 5.0**

> Email: **[samuel.farrens@gmail.com](mailto:samuel.farrens@gmail.com)**

## Contents
1. [Introduction](#intro_anchor)
1. [Method](#method_anchor)
1. [Dependenies](#depend_anchor)
1. [Input](#in_anchor)
1. [Execution](#exe_anchor)
1. [Output](#out_anchor)
1. [Troubleshooting](#trb_anchor)


<a name="intro_anchor"></a>
## Introduction

Pycymatch (versions 4 & 5) is a cylindrical matching code for identifying
matches between a catalogue of simulated dark matter haloes populated with
galaxies and the results of a cluster detection algorithm run on said
catalogue. The code defines the "cylinder" as the R200 distance from the halo
centre and `2 x dz(1 + z)` (where dz is the expected photometric redshift error).
The code additionally requires that any potential matches be within the
boundaries defined by the members of the mock halo (*i.e.* the farthest extent in
RA and Dec of the members).

Version 4 of the code rank orders the mock haloes and the detections by Ngal
(*i.e.* the number of galaxy members) and enforces unique matches. This ensures
that a given halo is always matched to the largest available detection that
satisfies the matching criteria. (*Note:* This method was implemented for the
Euclid Cluster Finder Challenges (CFC) 1 and 2)

Version 5 of the code does not rank order the objects. Instead this code finds
all possible detections (*i.e.* those that satisfy the matching criteria). A
weighting is then assigned to each match to either find the best possible
unique match (providing identical results to V4) or to find the probability
distribution of all of the matches.

Full Doxygen documentation is available for both versions of the code in the
docs folder. Additionally, version 5 and all its dependencies adhere to PEP8
style guidelines.

<a name="method_anchor"></a>
## Method

This section describes the cylindrical matching method in both versions (4 & 5) of Pycymatch.

### Version 4

This version of the code uses the following procedure for determining matches between mock haloes and detected clusters:

* The mock halo catalogue is rank ordered by *N_gal* (the number of galaxy members contained within the cluster/halo).
* The detected cluster catalogue is rank ordered by *lambda_obs* (the observed mass proxy, which can also be *N_gal*).
* For each mock halo a matching threshold region is defined around the halo centre using the following criteria:
 * The radial limit set by the *R_200* of the mock halo.
 * The maximum extent (max RA/Dec) set by the galaxy member farthest from the halo centre.
 * The line-of-sight limit set by *2 x sigma_z(1 + z)*, where *sigma_z* is the photometric redshift error of the mock halo catalogue.
* For each halo (starting with the highest ranked) a search is performed for the highest ranked detection within the matching threshold region.
* Only unique matches are permitted. In the case of multiple detections with equal rank the object closest to the halo centre is chosen as the match.
* All unmatched detections count as impurities.
* Completeness is measured as *N_matches/N_haloes*.
* Purity is measured as *N_matches/N_detections*.

### Version 5

This version of the code implements some minor modifications with respect to the previous version.

* The input catalogues are no longer rank ordered.
* The same matching threshold region is used.
* For each halo a list of all detections within the matching threshold region is found.
* Each match is assigned a weight *w = 1 - (d_proj/R_200)*, where *d_proj* is the projected distance between the halo and the detection. (Note that more sophisticated weights could be implemented)
* For unique matching only the match with the highest weight is kept. In this case the results are identical to those of version 4.
* For non-unique matching the weights can be used to make a probability distribution of *lambda_obs* for the mass of the halo.
* Completeness and purity are measured in the same way.

<a name="depend_anchor"></a>
## Dependencies

The code requires the following Python packages:

* <a href="http://www.numpy.org/" target="_blank">Numpy</a> (version 1.7 or greater)

* <a href="http://www.scipy.org/" target="_blank">SciPy</a>

* <a href="http://matplotlib.org/" target="_blank">Matplotlib</a>

<a name="in_anchor"></a>
## Input

### Input Files

The expected inputs for both versions of the code are the following two files (in ASCII):

1. The mock halo catalogue.
2. A catalogue of cluster detections.

The column values required for each of these files are described in the following subsections. Note, however, that the order of the properties provided are the code defaults. If needed the order can be altered using the code options described further down.

### Mock Halo Catalogue Format

For the mock halo catalogue the codes expect to find the following properties
in the following order:

1. **ID** : `(A unique string of numbers and/or characters to identify the halo)`
2. **Central Galaxy Flag** : `(A boolean, i.e. 0 or 1)`
3. **RA** : `(The halo Right Ascension in degrees)`
4. **Dec** : `(The halo Declination in degrees)`
5. **z** : `(The halo redshift)`
6. **Ngal** : `(The number of galaxy members in the halo)`
7. **Mass** : `(Log10 of the halo mass)`
8. **R200** : `(R200 of the halo in arcminutes)`
9. **Min RA** : `(Minimum halo member RA)`
10. **Max RA** : `(Maximum halo member RA)`
11. **Min Dec** : `(Minimum halo member Dec)`
12. **Max Dec** : `(Maximum halo member Dec)`

### Cluster Detection Catalogue Format

For the detection catalogue the codes expect to find the following properties
in the following order:

1. **ID** : `(A unique string of numbers and/or characters to identify the cluster)`
2. **RA** : `(The cluster Right Ascension in degrees)`
3. **Dec** : `(The cluster Declination in degrees)`
4. **z** : `(The cluster redshift)`
5. **Ngal** : `(The number of galaxy members in the cluster)`
6. **SNR** : `(The signal-to-noise ratio of the cluster)`

<a name="exe_anchor"></a>
## Execution

### Running the Codes

The codes can be run as executables by changing the file permissions
and specifying the path to Python (this is first line of the
pycymatch_v4.py and pycymatch_v5.py files and the default is /usr/bin/python)
e.g.:

> \>\> chmod +x pycymatch_v5.py

> \>\> pycymatch_v5.py

Otherwise to run the code symply run with Python e.g:

> \>\> python pycymatch_v5.py

Help and a list of arguments are provided with the `--help` option e.g:

> \>\> python pycymatch_v5.py --help

### Examples

To find matches with dz = 0.05 with version 4:

> \>\> python pycymatch_v4.py -c DETECTIONS_FILE -m MOCK_FILE -z 0.05

To find unique matches with version 5:

> \>\> python pycymatch_v5.py -i DETECTIONS_FILE MOCK_FILE -z 0.05 -u

To find non-unique matches with version 5 within 0.1<=z<2.6 and bin size 0.3:

> \>\> python pycymatch_v5.py -i DETECTIONS_FILE MOCK_FILE -z 0.05 --z_bin 0.1 2.6 0.3

### Code Options

The following screenshots show the available options for both versions of the
code.

<img src=figures/v4_options.jpg width=500>

<img src=figures/v5_options.jpg width=500>

*NOTE: A full description of code options will be added in the future.*

<a name="out_anchor"></a>
# Output

### Basic Outputs

Both versions of the code output the following plots:

* **Completeness vs. mass and redshift:** A 2D plot of the completeness `(N_Matches/N_Haloes)` in bins of halo mass and redshift.

*e.g.*

<img src=figures/completeness_mass_plot.jpg width=500>

* **Completeness vs. Ngal (halo) and redshift:** A 2D plot of the completeness `(N_Matches/N_Haloes)` in bins of halo Ngal (number of galaxy members in the halo) and redshift.

*e.g.*

<img src=figures/completeness_ngal_plot.jpg width=500>

* **Purity vs. Ngal (detection) and redshift:** A 2D plot of the purity `(N_Matches/N_Detections)` in bins of detection Ngal (number of galaxy members in the detection) and redshift.

*e.g.*

<img src=figures/purity_plot.jpg width=500>

### Additional Version 4 Outputs

Plots not included in version 5:

* **SNR vs. mass and redshift:** A 2D plot of the signal-to-noise ratio of the detections in bins of halo mass and redshift.

*e.g.*

<img src=figures/sn_plot.jpg width=500>

* **Spearman's rank order correlation coefficient vs. redshift:** A nonparametric measure of the correlation between the halo mass and the detection observable mass proxy.

*e.g.*

<img src=figures/spear_plot.jpg width=500>

Version 4 also outputs the following text files with the matches between mock haloes and detections:

* A list of all mock haloes and corresponding detection matches.
* A list of all detections and corresponding halo matches.

### Additional Version 5 Outputs

New plots provided in version 5:

* **Cluster mass observable (lambda) vs. halo mass in bins of redshift:** A set of 2D histograms (one for each redshift bin) in bins of lambda and halo mass.

*e.g.*

<img src=figures/mass_obs_plot.jpg width=500>

* **Histograms of halo mass in bins of redshift:** A set of histograms (one for each redshift bin) of the halo mass (blue solid line) and the corresponding matched detections (red dashed line).

*e.g.*

<img src=figures/mass_func_plot.jpg width=500>

* **Histograms of lambda in bins of mass and redshift.** A set of histograms (one for each mass and redshift bin) of the mass observable lambda (blue solid line) and a gaussian fit (red dashed line).

*e.g.*

<img src=figures/nlambda_plot.jpg width=500>

Version 5 outputs the following text files that differ from version 4:

* A list of matched haloes and corresponding detections.
* Mass-observable matrices in bins of redshift (*i.e.* the data required to reproduce all of the plots.)

<a name="trb_anchor"></a>
# Troubleshooting

### General

* Make sure that both of the input files contain all of the required columns with the correct units.
* Make sure that the input file columns are in the default order or that the appropriate options have been used to adjust the input order.

### Version 5 Specific

* Make sure that the observed detections file and then the mock halo file are provided following the `--i` option.
* Make sure that the default bin ranges (*e.g.* mass bins, redshift bins or Ngal bins) cover the full extent of the input data. For example, the default redshift range is `0<=z<=3.0` and if the input contains an object with `z=3.2` it will cause problems.
* Make sure to use the `--u` option to obtain unique matches.
