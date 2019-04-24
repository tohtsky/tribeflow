# TribeFlow C++ re-implementation

Contains a C++ Implementation of TribeFlow ([origina implementation](https://github.com/flaviovdf/tribeflow)) ``plearn`` functionality.
Instead of using MPI, parallel collapsed Gibbs sampling are performed in a single process with multiple threads. 
Thanks to pybind11, the installation procedure has been made simpler(I believe), if you have sufficiently new version of GCC or Clang.

# Prerequisites

To compile the module, you need recent (C++11 compatible) version of Clang or GCC and Eigen3.3.
Download Eigen3.3 from the official site.

The python dependencies are:

* pybind11
* numpy


# How to install 

```
EIGEN3_INCLUDE_DIR=/path/to/eigen python setup.py install
```

How to use
----------

Either use `python setup.py install` to install the packager or just use it from
the package folder using the `run_script.sh` command.

*How to parse datasets:* Use the `scripts/trace_converter.py` script. It has a help.

For command line help:

```bash
$ python scripts/trace_converter.py -h
$ python main.py -h
```

Running with mpi

```bash
$ mpiexec -np 4 python main.py [OPTIONS]
```

Running TribeFlow from other python code:

```
Check the api_singlecore_example.py file
```

Example
-------

**Converting the Trace**

Let's assume we have a trace like the Last.FM trace from [Oscar
Celma](http://www.dtic.upf.edu/~ocelma/MusicRecommendationDataset/lastfm-1K.html).
In this example, each line is of the form:

```bash
userid \t timestamp \t musicbrainz-artist-id \t artist-name \t
musicbrainz-track-id \t track-name
```

For instance:

```bash
user_000001 2009-05-01T09:17:36Z    c74ee320-1daa-43e6-89ee-f71070ee9e8f
Impossible Beings   952f360d-d678-40b2-8a64-18b4fa4c5f8Dois Pólos
```

First, we want to convert this file to our input format. We do this with the
`scripts/trace_converter.py` script. Let's have a look at the options from
this script:

```bash
$ python scripts/trace_converter.py -h
usage: trace_converter.py [-h] [-d DELIMITER] [-l LOOPS] [-r SORT] [-f FMT]
                          [-s SCALE] [-k SKIP_HEADER] [-m MEM_SIZE]
                          original_trace tstamp_column hypernode_column
                          obj_node_column

positional arguments:
  original_trace        The name of the original trace
  tstamp_column         The column of the time stamp
  hypernode_column      The column of the time hypernode
  obj_node_column       The column of the object node

optional arguments:
  -h, --help            show this help message and exit
  -d DELIMITER, --delimiter DELIMITER
                        The delimiter
  -l LOOPS, --loops LOOPS
                        Consider loops
  -r SORT, --sort SORT  Sort the trace
  -f FMT, --fmt FMT     The format of the date in the trace
  -s SCALE, --scale SCALE
                        Scale the time by this value
  -k SKIP_HEADER, --skip_header SKIP_HEADER
                        Skip these first k lines
  -m MEM_SIZE, --mem_size MEM_SIZE
                        Memory Size (the markov order is m - 1)
```

The positional (obrigatory) arguments are:

   * *original_trace* is the input file
   * *hypernode_column* represents the users (called hypernodes since it can 
     be playlists as well)
   * *tstamp_column* the column of the time stamp
   * *obj_node_column* the objects of interest

We can convert the file with the following line:

```bash
python scripts/trace_converter.py scripts/test_parser.dat 1 0 2 -d$'\t' \
        -f'%Y-%m-%dT%H:%M:%SZ' > trace.dat
```

Here, we are saying that column 1 are the timestamps, 0 is the user, and 2 are the
objects (artist ids). The delimiter *-d* is a tab. The time stamp format is
`'%Y-%m-%dT%H:%M:%SZ'`.

***Adding memory***

Use the -m argument to increase the burst (B parameter in the paper) size. 

```bash
python scripts/trace_converter.py scripts/test_parser.dat 1 0 2 -d$'\t' \
        -f'%Y-%m-%dT%H:%M:%SZ' -m 3 > trace.dat
```

**Learning the Model**

The example below is the same code used for every result in the paper. It runs
TribeFlow with the options used in every result in the paper. Explaining the
parameters:

   * *-np 20* Number of cores for execution.
   * *100* topics.
   * *output.h5* model file.
   * *--kernel eccdf* The kernel heuristic for inter-event time estimation. ECCDF
     based as per described on the paper. We also have a t-student kernel.
   * *--residency_priors 1 99* The priors for the inter-event time estimation.
   * *--leaveout 0.3* Number of transitions to leaveout.
   * *--num_iter 2000* Number of iterations.
   * *--dynamic True* Try and find the number of latent spaces from the data
   * *--num_batches 20* Number of split/merge moves.

*The example below uses 20 cores*
```bash
$ mpiexec -np 20 python main.py trace.dat 100 output.h5 \
    --kernel eccdf --residency_priors 1 99 --dynamic True \
    --leaveout 0.3 --num_iter 2000 --num_batches 20
```

**What if I don't want to explore inter-event times? Like Tribeflow-nt in the paper?**

Use the noop kernel

```bash
mpiexec -np 20 python main.py trace.dat 100 output.h5 \
    --kernel noop \
    --leaveout 0.3 --num_iter 2000
```

Notice that this line does not use *--dynamic* nor *--num_batches*. They are useless
and may lead to strange results in this case.

**Predictions**

The mean reciprocal rank script will generate predictions and save them to the 
given files. Just run:

```bash
$ PYTHONPATH=. python scripts/mrr.py output.h5 rss.dat predictions.dat
```

`output.h5` is the model trained.

**Other useful scripts**

Similar to the script above, you can use the scripts:

1. `view_topics.py` to print a summary of the topics with most likely objects
2. `printmat.py` to print either an O by O matrix or a Z by Z matrix 
3. `plotmat-toyplot.py` to generate the Z by Z matrix in the ISMIR jazz paper
4. `fancyplot.py` to generate the Miles Davis plot in the ISMIR jazz paper

<a name="data"></a>
Datasets
========

Below we have the list of datasets explored on the paper. We also curated links
to various other timestamp datasets that can be exploited by TribeFlow and 
future efforts.

Datasets used on the paper:

1. [LastFM-1k](http://www.dtic.upf.edu/~ocelma/MusicRecommendationDataset/lastfm-1K.html)
2. [LastFM-Our](https://www.dropbox.com/s/49utmj05rmv55wz/lastfm_groups.dat.gz?dl=0) 
3. [FourSQ](https://archive.org/details/201309_foursquare_dataset_umn)
    This dataset was removed from the original website. Still available on
    archive. Other, more recent, FourSQ datasets are available. See below.
4. [Brightkite](https://snap.stanford.edu/data/loc-brightkite.html)
5. [Yes](http://www.cs.cornell.edu/people/tj/playlists/index.html)

List of other, some more recent, datasets that can be explored by TribeFlow.

1. [Newer FourSQ](https://sites.google.com/site/yangdingqi/home/foursquare-dataset)
2. [Million Music Tweet](http://www.cp.jku.at/datasets/MMTD/)
3. [Movie Ratings](https://github.com/sidooms/MovieTweetings)
4. [Twitter](https://snap.stanford.edu/data/twitter7.html)
5. [Gowalla](https://snap.stanford.edu/data/loc-gowalla.html)
6. [Yelp](https://www.yelp.com/dataset_challenge)
7. [Best Buy](https://www.kaggle.com/c/acm-sf-chapter-hackathon-big/data)

Basically, anything with users (playlists, actors, etc also work), objects and 
timestamps.

On the `example` folder we have some sub-sampled datasets that can be used to
better understand the method.

<a name="reproductibility"></a>
Reproducibility
===============

The current version of the code may not be the exact version used in any of the
papers that employ Tribeflow. However, and most importantly, I am tagging the
commits closest to each paper. Please check the tags if you want to run an exact
version of tribeflow used in a given paper.

<a name="competition"></a>
Competing Methods
=================

* [HMMs](https://github.com/guyz/HMM) - or any other implementation
* [PRLME](http://github.com/flaviovdf/plme)
* [FPMC](http://github.com/flaviovdf/fpmc)
* [LME](http://www.cs.cornell.edu/people/tj/playlists/index.html)
* [Gravity Model](https://github.com/flaviovdf/tribeflow/blob/master/scripts/gravity_model.py)
* [TMLDA](https://github.com/flaviovdf/tribeflow/blob/master/scripts/tmlda.py)
* [StagesModel](http://infolab.stanford.edu/~crucis/code/stages-package.zip)
