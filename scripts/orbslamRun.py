#
#	author: Michal Nowicki
#
from subprocess import call
import sys
import os
import fileinput

def setYamlFile(detector, harrisK):
    print "Setting detector = " + detector + " harrisK = " + str(harrisK);

    outLines = []
    with open("Examples/Monocular/TUM_MONO_DSO.yaml", 'r') as f_in:
        for line in f_in:
            if "ORBextractor.detectorType" in line:
                outLines.append("ORBextractor.detectorType: %s\n" % detector);
            elif "ORBextractor.HarrisK" in line:
                outLines.append("ORBextractor.HarrisK: %f\n" % harrisK);
            else:
                outLines.append("%s" % line);
    with open("Examples/Monocular/TUM_MONO_DSO.yaml", 'w') as f_out:
        for line in outLines:
            f_out.write(line);



# Path to save main results - create if needed
if not os.path.exists("results"):
	os.makedirs("results");	
else:
	call('rm results/*', shell=True);


sequences = [
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
'13', \
'14', \
'15', \
'16', \
'17', \
'18', \
'19', \
'20', \
'21', \
'22', \
'23', \
'24', \
'25', \
'26', \
'27', \
'28', \
'29', \
'30', \
'31', \
'32', \
'33', \
'34', \
'35', \
'36', \
'37', \
'38', \
'39', \
'40', \
'41', \
'42', \
'43', \
'44', \
'45', \
'46', \
'47', \
'48', \
'49', \
'50'];

runsPerSequence = 1;


if len(sys.argv) != 2:
    print('Please provide the name of the host machine, e.g. python orbslamRun.py lrm2')
else:

    mainDatasetPath = '';

    # PPCM
    if "ppcm" in sys.argv[1]:
        mainDatasetPath = '/home/michal/dsoDataset';

    # LRM2
    if "lrm2" in sys.argv[1]:
        mainDatasetPath = '/mnt/data/Datasets/dso';

    # Toshiba Portege Z30
    if "local" in sys.argv[1]:
        mainDatasetPath = '/media/michalnowicki/MNowicki-Private/DSO/sequences';

    print 'mainDatasetPath: ' + mainDatasetPath

    # -------------------------------------------
    # Parameters
    detectorTypes = ["HARRIS", "HARRIS"];#, "ShiTomasi"];
    #harrisK = [0.002, 0.005, 0.01, 0.02, 0.04, 0.08]; Those were used in DSO
    #harrisKs = [0.002, 0.005]
    harrisKs = [0.02, 0.04]

    # For chosen detector
    for (detector, harrisK) in zip(detectorTypes, harrisKs):

        setYamlFile(detector, harrisK);

        # Create dir for chosen detector
        if not os.path.exists("results/" + detector + "_harrisK_" + str(harrisK)):
            os.makedirs("results/" + detector+ "_harrisK_" + str(harrisK));
        else:
            call('rm results/' + detector+ "_harrisK_" + str(harrisK) + '/*', shell=True);


        # For all selected sequences
        for seq in sequences:

            # For number of runs
            for runId in range(0, runsPerSequence):
                print("Current sequence: " + seq);

                # We call this command
                print('./Examples/Monocular/mono_tum_dso Vocabulary/ORBvoc.txt Examples/Monocular/TUM_MONO_DSO.yaml ' + str(mainDatasetPath) +'/sequence_' + seq +'/');

                # Run code
                call('./Examples/Monocular/mono_tum_dso Vocabulary/ORBvoc.txt Examples/Monocular/TUM_MONO_DSO.yaml ' + str(mainDatasetPath) +'/sequence_' + seq +'/', shell=True);

                # Copy results
                call('mv KeyFrameTrajectory.txt results/' + detector + '/sequence_' + str(seq) + '_' + str(runId) + '.txt', shell=True);



