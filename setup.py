#!/usr/bin/env python
# *************************************************************************
# * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
# * Produced at the Lawrence Livermore National Laboratory
# * Written by M. Black amd S. C. Wilks, LLNL
# * LLNL-CODE-739358
# * All rights reserved.
# *
# * This file is part of prorad.   For details, see https://github/LLNL/prorad.
# * Please also read this link:  https://github/LLNL/prorad/AdditionalBSDNotice.
# *
# * Redistribution and use in source and binary forms, with or without 
# * modification, are permitted provided that the following conditions are met:
# *
# *   *  Redistributions of source code must retain the above copyright notice, 
# *      this list of conditions and the disclaimer below.
# *   *  Redistributions in binary form must reproduce the above copyright 
# *      notice, this list of conditions and the disclaimer (as noted below) 
# *      in the documentation and/or other materials provided with the 
# *      distribution.
# *   *  Neither the name of the LLNS/LLNL nor the names of its contributorsa
# *      may be used to endorse or promote products derived from this softwarea
# *      without specific prior written permission.
# *
# * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# * ARE DISCLAIMED.  IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, 
# * LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
# * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
# * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
# * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
# * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
# ***************************************************************************/

# Use this script to compile the pmover module, or to rebuild it after changing anything in pmover.c or _pmover.c
# To compile, run:
#   
#       python setup.py build
# 
# If your system has OpenMP, you can compile to use it by instead typing:
#
#       python setup.py build -omp
#
# In either case you will probably get a warning about deprecated Numpy API, that is fine.
# TODO: Figure out how to suppress warning

from distutils.core import setup, Extension
import numpy.distutils.misc_util
import sys
import os
import shutil

# Get the directory of this script and remove the previous compiled module and build folder
path = os.path.dirname(os.path.realpath(__file__))
if '_pmover.so' in os.listdir(path):
    os.remove(os.path.join(path,"_pmover.so"))
if 'build' in os.listdir(path):
    shutil.rmtree(os.path.join(path,"build"))

# Determine compiler arguments
link_args = []
compile_args = ["-std=c99", "-O2"]
defines = []
if '-omp' in sys.argv:
    sys.argv.remove('-omp')
    link_args.append("-lgomp")
    compile_args.append("-fopenmp")
    defines.append(('USE_OMP','1'))
    print("Compiling to use OpenMP.")
else:
    print("Not compiling to use OpenMP.")

# Build Extension object
module = Extension("_pmover", 
                   sources=["_pmover.c", "pmover.c"],
                   extra_compile_args=compile_args, extra_link_args=link_args,
                   define_macros=defines)

# Compile module
setup(
    ext_modules=[module],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(), 
)
