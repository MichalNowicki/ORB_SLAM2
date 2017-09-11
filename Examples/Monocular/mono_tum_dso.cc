
#include<iostream>
#include<algorithm>
#include<fstream>
#include<chrono>

#include<opencv2/core/core.hpp>

#include<System.h>

#include<string>

// Code from DSO
#include <opencv2/highgui/highgui.hpp>
#include "Modifications/DatasetReader.h"


using namespace std;

void LoadImages(const string &strFile, vector<string> &vstrImageFilenames,
                vector<double> &vTimestamps);

int main(int argc, char **argv)
{
    if(argc != 4)
    {
        cerr << endl << "Usage: ./mono_tum_dso path_to_vocabulary path_to_settings path_to_sequence" << endl;
        return 1;
    }


    std::string basePath = std::string(argv[3]);
    std::string source = basePath + "images.zip";
    std::string calib = basePath + "camera.txt";
    std::string gammaCalib = basePath + "pcalib.txt";
    std::string vignette = basePath + "vignette.png";

    ImageFolderReader* reader = new ImageFolderReader(source,calib, gammaCalib, vignette);
    reader->setGlobalCalibration();

    std::vector<int> idsToPlay;
    std::vector<double> timesToPlayAt, timestamps;
    int lstart=0;
    int lend=100000;
    int linc = 1;
    for(int i=lstart;i>= 0 && i< reader->getNumImages() && linc*i < linc*lend;i+=linc)
    {
        idsToPlay.push_back(i);
        double tsThis = reader->getTimestamp(idsToPlay[idsToPlay.size()-1]);
        timestamps.push_back(tsThis);

        if(timesToPlayAt.size() == 0)
        {
            timesToPlayAt.push_back((double)0);
        }
        else
        {
            double tsThis = reader->getTimestamp(idsToPlay[idsToPlay.size()-1]);
            double tsPrev = reader->getTimestamp(idsToPlay[idsToPlay.size()-2]);
            timesToPlayAt.push_back(timesToPlayAt.back() +  fabs(tsThis-tsPrev));
            //printf("tsThis = %f\n", tsThis);
        }
    }


    // Create SLAM system. It initializes all system threads and gets ready to process frames.
    ORB_SLAM2::System SLAM(argv[1],argv[2],ORB_SLAM2::System::MONOCULAR,true);

    Eigen::Matrix3f K;
    int w;
    int h;
    reader->getCalibMono(K, w, h);

    SLAM.ChangeCalibration(K(0,0), K(1,1), K(0,2), K(1,2), 0, 0, 0, 0, 0);
    

    // TODO: Modified so if we skip every 2nd frame it is still 2.5 times slower than real-time
    const float playbackSpeed = 0.02;

    struct timeval tv_start;
    gettimeofday(&tv_start, NULL);
    double sInitializerOffset=timesToPlayAt[0];

    for(int ii=0;ii<(int)idsToPlay.size(); ii++) {

        int i = idsToPlay[ii];

        ImageAndExposure *img;
        cv::Mat imgOpenCV;

        img = reader->getImage(i, imgOpenCV);
        delete img;

        printf ("Processing image : %d \tTimestamp %f\r\n", ii, timestamps[ii]);

//        // tODO: Let's skip every 2nd frame
//        if ( ii%2 == 1 )
//            continue;

        // Pass the image to the SLAM system
//        SLAM.TrackMonocular(imgOpenCV,timesToPlayAt[ii]);
        SLAM.TrackMonocular(imgOpenCV,timestamps[ii]);


        // TODO: We will stop sequence after 1500 frames
//        if ( ii == 1500 )
//            break;

        bool skipFrame = false;
        struct timeval tv_now;
        gettimeofday(&tv_now, NULL);
        double sSinceStart = sInitializerOffset + ((tv_now.tv_sec - tv_start.tv_sec) +
                                                   (tv_now.tv_usec - tv_start.tv_usec) / (1000.0f * 1000.0f));

        if (sSinceStart*playbackSpeed < timesToPlayAt[ii])
            usleep((int) ((timesToPlayAt[ii] - sSinceStart*playbackSpeed) * 1000 * 1000));
    }



    // Stop all threads
    SLAM.Shutdown();

    // Save camera trajectory
    SLAM.SaveKeyFrameTrajectoryTUM("KeyFrameTrajectory.txt");

    return 0;
}

void LoadImages(const string &strFile, vector<string> &vstrImageFilenames, vector<double> &vTimestamps)
{
    ifstream f;
    f.open(strFile.c_str());

    // skip first three lines
    string s0;
    getline(f,s0);
    getline(f,s0);
    getline(f,s0);

    while(!f.eof())
    {
        string s;
        getline(f,s);
        if(!s.empty())
        {
            stringstream ss;
            ss << s;
            double t;
            string sRGB;
            ss >> t;
            vTimestamps.push_back(t);
            ss >> sRGB;
            vstrImageFilenames.push_back(sRGB);
        }
    }
}
