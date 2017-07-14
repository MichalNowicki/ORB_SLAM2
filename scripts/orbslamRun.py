#
#	author: Michal Nowicki
#
from subprocess import call
import sys
import os
import fileinput

def setYamlFile(detector, harrisK, lambdaThreshold, sigma):
    print "Setting detector = " + detector + " harrisK = " + str(harrisK);

    outLines = []
    with open("Examples/Monocular/TUM_MONO_DSO.yaml", 'r') as f_in:
        for line in f_in:
            if "ORBextractor.detectorType" in line:
                outLines.append("ORBextractor.detectorType: %s\n" % detector);
            elif "ORBextractor.HarrisK" in line:
                outLines.append("ORBextractor.HarrisK: %f\n" % harrisK);
            elif "ORBextractor.lambdaThreshold" in line:
                outLines.append("ORBextractor.lambdaThreshold: %f\n" % lambdaThreshold);
            elif "Optimization.sigma" in line:
                outLines.append("Optimization.sigma: %f\n" % sigma);
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
    detectorTypes = ["FAST"]#["FAST", "FAST", "FAST", "FAST"]#, "Harris", "HarrisCE"];#, "ShiTomasi"];
    harrisKs = [0]#[0, 0, 0 ,0]#, 0.01, 0.01];
    lambdaThresholds = [0]#[0, 0, 0 ,0]#, 0.001, 0.001];
    sigmas = [1.0]#[0.35, 0.5, 0.75, 1.0]
    #harrisKs = [0.002, 0.005, 0.01, 0.02, 0.04]; # Those were used in DSO

    # For chosen detector
    for (detector, harrisK, lambdaThreshold, sigma) in zip(detectorTypes, harrisKs, lambdaThresholds, sigmas):

        setYamlFile(detector, harrisK, lambdaThreshold, sigma);

        # Create dir for chosen detector
        if "HarrisCE" in detector:
            if not os.path.exists("results/" + detector + "_lambdaThreshold_" + str(lambdaThreshold)):
                os.makedirs("results/" + detector+ "_lambdaThreshold_" + str(lambdaThreshold));
            else:
                call('rm results/' + detector+ "_lambdaThreshold_" + str(lambdaThreshold) + '/*', shell=True);
        elif "Harris" in detector:
            if not os.path.exists("results/" + detector + "_harrisK_" + str(harrisK)):
                os.makedirs("results/" + detector+ "_harrisK_" + str(harrisK));
            else:
                call('rm results/' + detector+ "_harrisK_" + str(harrisK) + '/*', shell=True);
        else:
            if not os.path.exists("results/" + detector +"_sigma_" + str(sigma) ):
                os.makedirs("results/" + detector +"_sigma_" + str(sigma));
            else:
                call('rm results/' + detector +"_sigma_" + str(sigma), shell=True);


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
                if "HarrisCE" in detector:
                    call('mv KeyFrameTrajectory.txt results/' + detector+ "_harrisK_" + str(harrisK) + '/sequence_' + str(seq) + '_' + str(runId) + '.txt', shell=True);
                    call('mv GBA_chi2.txt results/' + detector+ "_harrisK_" + str(harrisK) + '/sequence_' + str(seq) + '_GBA_chi2.txt', shell=True);
                    call('mv GBA_reproj.txt results/' + detector+ "_harrisK_" + str(harrisK) + '/sequence_' + str(seq) + '_GBA_reproj.txt', shell=True);
                elif "Harris" in detector:
                    call('mv KeyFrameTrajectory.txt results/' + detector+ "_lambdaThreshold_" + str(lambdaThreshold) + '/sequence_' + str(seq) + '_' + str(runId) + '.txt', shell=True);
                    call('mv GBA_chi2.txt results/' + detector+ "_lambdaThreshold_" + str(lambdaThreshold) + '/sequence_' + str(seq) + '_GBA_chi2.txt', shell=True);
                    call('mv GBA_reproj.txt results/' + detector+ "_lambdaThreshold_" + str(lambdaThreshold) + '/sequence_' + str(seq) + '_GBA_reproj.txt', shell=True);
                else:
                    call('mv KeyFrameTrajectory.txt results/' + detector +"_sigma_" + str(sigma)+ '/sequence_' + str(seq) + '_' + str(runId) + '.txt', shell=True);
                    call('mv GBA_chi2.txt results/' + detector +"_sigma_" + str(sigma)+ '/sequence_' + str(seq) + '_GBA_chi2.txt', shell=True);
                    call('mv GBA_reproj.txt results/' + detector +"_sigma_" + str(sigma)+ '/sequence_' + str(seq) + '_GBA_reproj.txt', shell=True);

