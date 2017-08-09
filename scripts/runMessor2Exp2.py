#
#	author: Michal Nowicki
#
from subprocess import call
import sys
import os
import fileinput


# Path to save main results - create if needed
if not os.path.exists("results"):
    os.makedirs("results");
else:
    call('rm results/*', shell=True);


sequences = [
    'outdoor1', \
    'outdoor2', \
    'outdoor3', \
    'testWalk2', \
    'tripodFreeVel0p35', \
    'tripodFreeVel0p7', \
    'tripodRampVel0p2', \
    'tripodRampVel0p25', \
    'tripodRampVel0p25_1', \
    'tripodSquareRotVel0p1', \
    'tripodSquareRotVel0p35', \
    'tripodSquareRotVel0p7', \
    'tripodSquareRotVel1p0', \
    'tripodSquareVel0p1', \
    'tripodSquareVel0p35', \
    'tripodSquareVel0p7', \
    'tripodSquareVel1p0', \
    'tripodXYZvel0p1', \
    'tripodXYZvel0p35', \
    'tripodXYZvel0p7', \
    'tripodXYZvel1p0' \
];

runsPerSequence = 1;


mainDatasetPath = '/mnt/data/MessorExp2';

print 'mainDatasetPath: ' + mainDatasetPath


# For all selected sequences
for seq in sequences:
    print("Current sequence: " + seq);

    if not os.path.exists("results/" + str(seq)):
        os.makedirs("results/" + str(seq));
    else:
        call('rm results/' + str(seq) + '/*', shell=True);

    # We call this command
    print('./Examples/RGB-D/rgbd_tum Vocabulary/ORBvoc.txt Examples/RGB-D/Messor2Exp2.yaml ' + str(mainDatasetPath) + '/' + str(seq) +' ' + str(mainDatasetPath) +'/orbslamAssociate.txt');

    # Run code
    call('./Examples/RGB-D/rgbd_tum Vocabulary/ORBvoc.txt Examples/RGB-D/Messor2Exp2.yaml ' + str(mainDatasetPath) + '/' + str(seq) +' ' + str(mainDatasetPath) +'/orbslamAssociate.txt', shell=True);

    # Copy results
    #CameraTrajectory.txt
    #call('mv KeyFrameTrajectory.txt results/' + detector+ "_harrisK_" + str(harrisK) + '/sequence_' + str(seq) + '_' + str(runId) + '.txt', shell=True);
