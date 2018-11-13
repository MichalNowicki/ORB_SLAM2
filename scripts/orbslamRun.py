#
#	author: Michal Nowicki
#
from subprocess import call
import sys
import os
import fileinput

def setYamlFile(filePathWithName, parameter, value):
    print "Setting " + parameter + " = " + str(value);

    outLines = []
    with open(filePathWithName, 'r') as f_in:
        for line in f_in:
            if parameter in line:
                outLines.append(parameter + str(value) +  "\n");
            else:
                outLines.append("%s" % line);
    with open(filePathWithName, 'w') as f_out:
        for line in outLines:
            f_out.write(line);



# Path to save main results - create if needed
if not os.path.exists("results"):
	os.makedirs("results");	
else:
	call('rm results/*', shell=True);


sequences = [
'00', \
'01', \
'02', \
'03', \
'04', \
'05', \
'06', \
'07', \
'08', \
'09', \
'10', \
'11', \
'12', \
    ];

runsPerSequence = 1;


if len(sys.argv) != 2:
    print('Please provide the name of the host machine, e.g. python orbslamRun.py lrm2')
else:

    mainDatasetPath = '';

    # # PPCM
    # if "ppcm" in sys.argv[1]:
    #     mainDatasetPath = '/home/michal/dsoDataset';

    # LRM2
    if "lrm2" in sys.argv[1]:
        mainDatasetPath = '/mnt/data/Datasets/kitti/sequences';

    # Toshiba Portege Z30
    if "local" in sys.argv[1]:
        mainDatasetPath = '/home/mnowicki/Projects/Datasets/KITTI/sequences';

    print 'mainDatasetPath: ' + mainDatasetPath

    # -------------------------------------------
    # Parameters:
    #   - tracking.kltTrack: 1
    #   - tracking.kltMaxIterations: 30
    #   - tracking.kltEPS: 0.01
    #   - tracking.kltZNCCThreshold: 0.900000
    #   - tracking.kltPatchSize: 11
    #   - tracking.kltMaxMovement: 5
    #
    kltTrackings = [1, 1, 1, 1];
    kltZNCCThresholds = [0.85, 0.9, 0.92, 0.95];
    kltPatchSizes = [11, 11, 11, 11];
    kltMaxMovements = [5, 5, 5, 5];

    # For chosen detector
    for (kltTracking,kltZNCCThreshold,kltPatchSize, kltMaxMovement) in zip(kltTrackings, kltZNCCThresholds,kltPatchSizes, kltMaxMovements):

        # For all selected sequences
        for seq in sequences:

            yamlName = "";
            # Depending on the KITTI yaml
            if seq in ['00', '01', '02']:
                yamlName = "KITTI00-02.yaml";
            elif seq in ['03']:
                yamlName = "KITTI03.yaml";
            else:
                yamlName = "KITTI04-12.yaml";

            setYamlFile("Examples/Stereo/" + yamlName, "tracking.kltTrack: ", kltTracking);
            setYamlFile("Examples/Stereo/" + yamlName, "tracking.kltZNCCThreshold: ", kltZNCCThreshold);
            setYamlFile("Examples/Stereo/" + yamlName, "tracking.kltPatchSize: ", kltPatchSize);
            setYamlFile("Examples/Stereo/" + yamlName, "tracking.kltMaxMovement: ", kltMaxMovement);

            # Patch depending on the parameters
            dir = "klt_" + str(kltTracking) + "_znccThr_" + str(kltZNCCThreshold) + "_patchSize_" + str(kltPatchSize) + "_kltMaxMovement_" + str(kltMaxMovement);

            # Create dir for chosen detector
            if not os.path.exists("results/" + dir + "/data"):
                os.makedirs("results/" + dir + "/data");
            else:
                call('rm results/' + dir + '/data/*', shell=True);

            if not os.path.exists("results/" + dir + "/inliers/sequence_" + str(seq)):
                os.makedirs("results/" + dir + "/inliers/sequence_" + str(seq));
            else:
                call('rm results/' + dir + 'inliers/sequence_' + str(seq) + '/*', shell=True);


            # We call this command
            print('./Examples/Stereo/stereo_kitti Vocabulary/ORBvoc.txt Examples/Stereo/' + yamlName + ' ' + str(mainDatasetPath) +'/' + seq +'/');

            # Run code
            call('./Examples/Stereo/stereo_kitti Vocabulary/ORBvoc.txt Examples/Stereo/' + yamlName + ' '  + str(mainDatasetPath) +'/' + seq +'/', shell=True);

            # Copy results
            call('mv CameraTrajectory.txt results/' + dir + '/data/' + str(seq) + '.txt', shell=True);
            call('mv logs/*.txt results/' + dir +  '/inliers/sequence_' + str(seq) + '/', shell=True);
            call('rm logs/*.txt', shell=True);