clear;
clc;
close all;

% dir = '/home/mnowicki/Projects/Tardos/KITTI/orbslam2_klt/klt0.85_patchSize9/inliers/sequence_10/';
% dir = '/home/mnowicki/Projects/Tardos/KITTI/orbslam2_klt/klt0.9_run1/inliers/sequence_10/';
% dir = '/home/mnowicki/Projects/Tardos/KITTI/orbslam2_klt/klt0.95_run1/inliers/sequence_01/';
% dir = '/home/mnowicki/Projects/Tardos/KITTI/orbslam2_klt/orig_run1/inliers/sequence_01/';
%dir = '/home/mnowicki/Projects/Tardos/ORB_SLAM2/logs/';

% dir = '/home/mnowicki/Projects/Tardos/KITTI/orbslam2_klt/klt0.85_patchSize_11_kltMaxMovement_5/inliers/sequence_01/';
% dir = '/home/mnowicki/Projects/Tardos/KITTI/orbslam2_klt/klt0.9_patchSize_11_kltMaxMovement_5/inliers/sequence_01/';
% dir = '/home/mnowicki/Projects/Tardos/KITTI/orbslam2_klt/klt0.92_patchSize_11_kltMaxMovement_5/inliers/sequence_01/';


% dirVersion = '/home/mnowicki/Projects/Tardos/KITTI/orbslam2_woLC/klt0.8_patchSize_11_kltMaxMovement_5/';
%dirVersion = '/home/mnowicki/Projects/Tardos/KITTI/orbslam2_woLC/orig_woLC_run1/';
dirVersion = '/home/mnowicki/Projects/Tardos/KITTI/orbslam2_woLC/test/';

seq = 'sequence_01';

dir = strcat(dirVersion, 'inliers/', seq);


%
% candidates (1), matches (2), tracks (3), extra tracks over matches (4), inliers (5),
% inliersTracking (6), zncc (7), avgTravelDist (8)
%
voInliers = importdata(strcat(dir, '/voInlierCount.txt'));

%
% matches, tracks, candidates, inliers, total
%
mapInliers = importdata(strcat(dir, '/mapInlierCount.txt'));


clc
fileID = fopen(strcat(dirVersion, seq, '_inlierStats.txt'),'w');

fprintf(['VO\n\tCand: %.2f \n\tMatches: %.2f \n\tTracks: %.2f \n\tMatches + tracks: %.2f (%.2f) ',...
    '\n\tInliers: %.2f (t:%.2f, m:%.2f)\n\tInlier perc. - t:%.2f m:%.2f\n\tAvg znnc of inliners: %.2f\n'], ...
    mean(voInliers(:,1)), mean(voInliers(:,2)), mean(voInliers(:,3)), mean(voInliers(:,2) + voInliers(:,4)), ...p[
    mean(voInliers(:,4)), mean(voInliers(:,5)), mean(voInliers(:,6)), mean(voInliers(:,5) - voInliers(:,6)), ...
    mean(voInliers(:,6)) / mean(voInliers(:,4)), ...
    (mean(voInliers(:,5)) - mean(voInliers(:,6))) / mean(voInliers(:,2)), ...
    mean(voInliers(:,7)));

fprintf(['Map\n\tAdd cand: %f \n\tAdd matches: %f \n\tAdd tracks: %f \n\tInliers: %f \n\tTotal: %f', ...
    '\n\tInlier perc. - m:%.2f\n\n'], ...
    mean(mapInliers(:,3)), mean(mapInliers(:,1)), mean(mapInliers(:,2)), ...
    mean(mapInliers(:,4)), mean(mapInliers(:,5)), mean(mapInliers(:,4))/mean(mapInliers(:,5)));

fprintf(fileID,['VO\n\tCand: %.2f \n\tMatches: %.2f \n\tTracks: %.2f \n\tMatches + tracks: %.2f (%.2f) ',...
    '\n\tInliers: %.2f (t:%.2f, m:%.2f)\n\tInlier perc. - t:%.2f m:%.2f\n\tAvg znnc of inliners: %.2f\n'], ...
    mean(voInliers(:,1)), mean(voInliers(:,2)), mean(voInliers(:,3)), mean(voInliers(:,2) + voInliers(:,4)), ...p[
    mean(voInliers(:,4)), mean(voInliers(:,5)), mean(voInliers(:,6)), mean(voInliers(:,5) - voInliers(:,6)), ...
    mean(voInliers(:,6)) / mean(voInliers(:,4)), ...
    (mean(voInliers(:,5)) - mean(voInliers(:,6))) / mean(voInliers(:,2)), ...
    mean(voInliers(:,7)));

fprintf(fileID,['Map\n\tAdd cand: %f \n\tAdd matches: %f \n\tAdd tracks: %f \n\tInliers: %f \n\tTotal: %f', ...
    '\n\tInlier perc. - m:%.2f\n\n'], ...
    mean(mapInliers(:,3)), mean(mapInliers(:,1)), mean(mapInliers(:,2)), ...
    mean(mapInliers(:,4)), mean(mapInliers(:,5)), mean(mapInliers(:,4))/mean(mapInliers(:,5)));

fclose(fileID);

%% Figures

% VO candidates
figure;
plot(voInliers(:,1), 'r', 'LineWidth',1);
hold on;
plot(voInliers(:,2), 'g', 'LineWidth',1);
plot(voInliers(:,3), 'b', 'LineWidth',1);
plot(voInliers(:,2) + voInliers(:,4), 'c', 'LineWidth',1);
legend('Last inliers (candidates)', 'VO Matches', 'VO Trackings', 'VO Matches + extra tracks');
title('VO - frame2frame');
saveas(gcf,strcat(dirVersion, seq, '_VO_candidates.png'));

% VO inliers
figure;
plot(voInliers(:,5), 'r', 'LineWidth',1);
hold on;
plot(voInliers(:,5) - voInliers(:,6), 'g', 'LineWidth',1);
plot(voInliers(:,6), 'b', 'LineWidth',1);
legend('Total inliers', 'Inliers matches', 'Inliers trackings');
title('VO - inlier analysis');
saveas(gcf,strcat(dirVersion, seq, '_VO_inliers.png'));

% VO inliers vs ZNCC
figure;
plot(voInliers(:,6), 'r', 'LineWidth',1);
hold on;
yyaxis right;
plot(voInliers(:,7), 'g', 'LineWidth',1);
legend('Inliers tracking', 'ZNCC');
title('Tracking inliers vs avg zncc');
saveas(gcf,strcat(dirVersion, seq, '_VO_inliersVsZNCC.png'));

% Map frame2localFrames
figure;
plot(mapInliers(:,1), 'r', 'LineWidth',1);
hold on;
plot(mapInliers(:,2), 'g', 'LineWidth',1);
plot(mapInliers(:,3), 'b', 'LineWidth',1);
plot(mapInliers(:,4), 'c', 'LineWidth',1);
plot(mapInliers(:,5), 'm', 'LineWidth',1);
legend('Map matches', 'Map tracks', 'Map add cand', 'Map inliers', 'Map total cand');
title('Map - frame2localFrames');
saveas(gcf,strcat(dirVersion, seq, '_Map_frame2localFrames.png'));
