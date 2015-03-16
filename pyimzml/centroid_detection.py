# -*- coding: UTF-8 -*-

# Copyright 2015 Andrew Palmer
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import numpy as np

"""
Performs gradient-based centroid detection on the specified spectrum.

input
mzs:
    m/z array
intensities:
    intensity array
max_output (Default: -1):
    limits the amount of centroids to be returned. If not set, all centroids are returned
weighted_bins (Default: 1):

min_intensity (Default: 0):
    threshold beneath which intensity values are ignored. If not set, all intensity values will be considered.

output
mzs_list:
    those m/z values from the input array that correspond to a centroid
intensities_list:
    those intensity values from the input array that correspond to a centroid
indices_list:
    position of the centroids within the original array
"""
def gradient(mzs,intensities,**opt_args):
    function_args = {'max_output': -1, 'weighted_bins':1,'min_intensity':0}
    for key, val in opt_args.iteritems():
        if key in function_args.keys():
            function_args[key] = val
        else:
            print ('possible arguments:')
            for i in function_args.keys():
                   print (i)
            raise NameError('gradient does not take argument: %s'%(key))
    mzMaxNum = function_args['max_output']
    weighted_bins = function_args['weighted_bins']
    min_intensity = function_args['min_intensity']
    #calc first differential
    MZgrad = np.gradient(intensities)
    #calc second differential
    MZgrad2 = np.gradient(MZgrad)
    # detect crossing points
    indices = MZgrad[0:-1]*MZgrad[1:]<=0
    indices = np.concatenate((indices,[False]))
    # bool->list of indices
    indices_list = np.transpose(np.where(indices==True))
    indices_list = indices_list[MZgrad2[indices_list]<0]#pick maxima



    #Remove any 'peaks' that aren't real
    indices_list = indices_list[intensities[indices_list] > min_intensity]

    #Select the peaks
    intensities_list = intensities[indices_list]
    mzs_list = mzs[indices_list]

    # Tidy up if required
    if mzMaxNum > 0:
        if len(mzs_list) > mzMaxNum:
            sort_idx = np.argsort(intensities_list)
            intensities_list = intensities_list[sort_idx[-mzMaxNum:]]
            mzs_list = mzs_list[sort_idx[-mzMaxNum:]]
            indices_list = indices_list[sort_idx[-mzMaxNum:]]
        elif len(mzs) < mzMaxNum:
            lengthDiff = mzMaxNum - len(indices_list)
            mzs_list = np.concatenate((mzs_list, np.zeros((lengthDiff,1))))
            intensities_list  = np.concatenate((intensities_list, np.zeros((lengthDiff,1))))
            indices_list = np.concatenate((indices_list, np.zeros((lengthDiff,1))))

    if weighted_bins > 0:
        # check no peaks within bin width of spectrum edge
        good_idx = (indices_list>weighted_bins) & (indices_list<(len(mzs)-weighted_bins))
        mzs_list = mzs_list[good_idx]
        intensities_list = intensities_list[good_idx]
        indices_list = indices_list[good_idx]
        bin_shift = range(-weighted_bins,weighted_bins+1)
        for ii in range(0,len(mzs_list)):
            bin_idx = bin_shift+indices_list[ii]
            mzs_list[ii] = np.average(mzs[bin_idx],weights=intensities[bin_idx])
            indices_list[ii] = indices_list[ii] + np.argmax(intensities[bin_idx]) - weighted_bins
            intensities_list[ii] = intensities[indices_list[ii]]
            #mzs_list[ii] = np.mean(mzs[bin_idx])
    return  (mzs_list,intensities_list,indices_list)
