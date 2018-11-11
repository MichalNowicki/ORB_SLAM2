#
#	author: Michal Nowicki
#
from subprocess import call
import sys
import os
import fileinput

import matplotlib.pyplot as plt
import numpy as np



dirName = sys.argv[1];

call('./evaluate_odometry /home/mnowicki/Projects/Tardos/ORB_SLAM2/results/' + str(dirName), shell=True);

t = np.arange(0.0, 2.0, 0.01)
s = 1 + np.sin(2*np.pi*t)
plt.plot(t, s)

plt.xlabel('time (s)')
plt.ylabel('voltage (mV)')
plt.title('About as simple as it gets, folks')
plt.grid(True)
# plt.savefig("test.png")
plt.show()