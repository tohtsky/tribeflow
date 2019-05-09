# TribeFlow C++ re-implementation

Contains a C++ Implementation of TribeFlow ([original implementation](https://github.com/flaviovdf/tribeflow)). Currently ``plearn`` functionality has been re-written.  

The most notable difference is that instead of using MPI, parallel collapsed Gibbs sampling are performed in a single process with multiple ``std::thread``.

Also, thanks to richer standard library of C++, several handmade functionality present in the original implementation are simply replaced.

The port to Python is done by using pybind11, which makes the installation procedure somewhat simpler(I believe), if you have a C++11 compatible version of GCC or Clang.

# Prerequisites

To compile the module, you need recent (C++11 compatible) version of Clang or GCC and Eigen3.3.
Download Eigen3.3 from the official site.

Other notable dependency in the C++ side is [TartanLlama/optional](https://github.com/TartanLlama/optional), but it has been included in the source and no action is required.

# How to install 
To compile and install the module, run
```
EIGEN3_INCLUDE_DIR=/path/to/eigen python setup.py install
```
which also installs pybind11 and numpy (if not present).

# To Dos

* implement eccdf kernel
* implement single thread mode
* clean unnecessary python codes
* dynamic topic size


How to use
----------

*How to parse datasets:* Use the `scripts/trace_converter.py` script. It has a help.

For command line help:

```bash
$ python scripts/trace_converter.py -h
```

Example
-------

### **Converting the Trace**

See the original implementation.

### **Learning the Model**

TBD...


