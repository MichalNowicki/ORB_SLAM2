#
#	author: Michal Nowicki
#
from subprocess import call
import subprocess
import sys
import glob
import os

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


mainDatasetPath = '/mnt/data/Datasets/dso';
#mainDatasetPath = '/media/michalnowicki/MNowicki-Private/DSO/sequences/';

for seq in sequences:


	for runId in range(0, runsPerSequence):
		print("Current sequence: " + seq);

		print('./Examples/Monocular/mono_tum_dso Vocabulary/ORBvoc.txt Examples/Monocular/TUM_MONO_DSO.yaml ' + str(mainDatasetPath) +'/sequence_' + seq +'/');

		# Copy to currently used settings
		call('./Examples/Monocular/mono_tum_dso Vocabulary/ORBvoc.txt Examples/Monocular/TUM_MONO_DSO.yaml ' + str(mainDatasetPath) +'/sequence_' + seq +'/', shell=True);
		
		# Run software
		call('mv KeyFrameTrajectory.txt results/sequence_' + str(seq) + '_' + str(runId) + '.txt', shell=True);

